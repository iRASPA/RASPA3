module;

#define CL_TARGET_OPENCL_VERSION 120
#define CL_SILENCE_DEPRECATION
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

module energy_opencl_molecular_energy_grid;

import std;

import opencl;
import int3;
import uint3;
import double3;
import double3x3;
import unit_cell;
import crystal;
import pair_interactions;
import units;
import energy_shared_linear_probe;
import energy_shared_probe_energy_grid;
import energy_shared_electrostatic_potential_grid;
import energy_shared_molecular_energy_grid;


MolecularEnergyGrid MolecularEnergyGridOpenCL::compute(const PairInteractions &interactions, const Crystal &framework,
                                                 const LinearProbe &probe, uint3 gridSize,
                                                 std::size_t numberOfOrientations, double temperature,
                                                 const ElectrostaticPotentialGrid *potential)
{
  if (!OpenCL::clContext.has_value() || !OpenCL::clDeviceId.has_value() || !OpenCL::clCommandQueue.has_value())
  {
    throw std::runtime_error("MolecularEnergyGridOpenCL: no OpenCL device found, the grid route needs a GPU\n");
  }
  if (gridSize.x == 0 || gridSize.y == 0 || gridSize.z == 0)
  {
    throw std::runtime_error("MolecularEnergyGridOpenCL: the grid must have at least one point along each axis\n");
  }
  if (probe.sites.empty())
  {
    throw std::runtime_error("MolecularEnergyGridOpenCL: the probe has no sites\n");
  }
  if (numberOfOrientations == 0)
  {
    throw std::runtime_error("MolecularEnergyGridOpenCL: at least one orientation is needed\n");
  }
  if (temperature <= 0.0)
  {
    throw std::runtime_error("MolecularEnergyGridOpenCL: the free energy needs a temperature above zero\n");
  }

  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  MolecularEnergyGrid grid;
  grid.gridSize = gridSize;
  grid.backend = "gpu";
  grid.unitCell = framework.unitCell;
  grid.probe = probe;
  grid.overHemisphere = probe.headTailSymmetric;
  grid.temperature = temperature;
  grid.cutOff = interactions.cutOffVDW;
  grid.ceiling = probeEnergyCeilingInKelvin * Units::KelvinToEnergy;

  // Electrostatics are acted on only when the molecule has charges to act on and a potential was supplied to
  // act with. A charged molecule without one is the case worth being loud about, so it is recorded rather
  // than quietly treated as neutral.
  grid.chargesIncluded = probe.isCharged() && potential != nullptr && potential->numberOfVoxels() > 0;
  grid.chargesIgnored = probe.isCharged() && !grid.chargesIncluded;
  if (grid.chargesIncluded)
  {
    grid.ewaldAlpha = potential->alpha;
    grid.numberOfWaveVectors = potential->numberOfWaveVectors;
    if (potential->gridSize.x != gridSize.x || potential->gridSize.y != gridSize.y ||
        potential->gridSize.z != gridSize.z)
    {
      throw std::runtime_error("MolecularEnergyGridOpenCL: the potential is on a different grid than this one\n");
    }
  }

  // A molecule that is the same end for end gives the same energy for a direction and its opposite, so half
  // the sphere holds every distinct orientation and sampling the whole of it would only duplicate.
  std::vector<double3> directions = orientationSet(numberOfOrientations, probe.headTailSymmetric);
  grid.numberOfOrientations = directions.size();

  std::size_t numberOfVoxels = gridSize.x * gridSize.y * gridSize.z;
  grid.freeEnergy.assign(numberOfVoxels, 0.0f);
  grid.minimumEnergy.assign(numberOfVoxels, 0.0f);

  std::vector<double3> fractionalPositions = framework.fractionalPositions;
  if (fractionalPositions.empty()) return grid;

  const std::size_t numberOfAtoms = fractionalPositions.size();
  const std::size_t numberOfSites = probe.sites.size();

  double3x3 cell = framework.unitCell.cell;
  double3x3 inverseCell = framework.unitCell.inverseCell;

  // Where each site sits relative to the centre of mass, as a fractional displacement. It depends on the
  // orientation and the site but not on the grid point, so it is built once here and the kernel needs no
  // inverse cell of its own.
  std::vector<cl_float4> siteOffsets(directions.size() * numberOfSites);
  for (std::size_t o = 0; o < directions.size(); ++o)
  {
    for (std::size_t s = 0; s < numberOfSites; ++s)
    {
      double3 displacement = inverseCell * (probe.sites[s].offset * directions[o]);
      siteOffsets[o * numberOfSites + s] = {{cl_float(displacement.x), cl_float(displacement.y),
                                             cl_float(displacement.z), 0.0f}};
    }
  }

  // The mixed size and strength for every site against every framework atom, taken from the force field so
  // that its own rule applies rather than one repeated here.
  std::vector<cl_float4> positions(numberOfAtoms);
  std::vector<cl_float> epsilonTimesFour(numberOfSites * numberOfAtoms);
  std::vector<cl_float> sigmaSquared(numberOfSites * numberOfAtoms);
  std::vector<cl_float> chargeProduct(numberOfSites * numberOfAtoms, 0.0f);
  std::vector<cl_float> shiftValue(numberOfSites * numberOfAtoms, 0.0f);
  std::vector<cl_float> siteCharge(numberOfSites);
  for (std::size_t s = 0; s < numberOfSites; ++s) siteCharge[s] = cl_float(probe.sites[s].charge);
  for (std::size_t i = 0; i < numberOfAtoms; ++i)
  {
    positions[i] = {{cl_float(fractionalPositions[i].x), cl_float(fractionalPositions[i].y),
                     cl_float(fractionalPositions[i].z), 0.0f}};

    std::size_t atomType = framework.atoms[i].type;
    for (std::size_t s = 0; s < numberOfSites; ++s)
    {
      double sigma = interactions(probe.sites[s].type, atomType).sizeParameter;
      double epsilon = interactions(probe.sites[s].type, atomType).strengthParameter;
      epsilonTimesFour[s * numberOfAtoms + i] = cl_float(4.0 * epsilon);
      sigmaSquared[s * numberOfAtoms + i] = cl_float(sigma * sigma);

      // The near half of the lattice sum is done pair by pair alongside the dispersion, since the walk over
      // neighbours is already being made and an erfc is all that is added to it.
      chargeProduct[s * numberOfAtoms + i] =
          cl_float(Units::CoulombicConversionFactor * probe.sites[s].charge * framework.atoms[i].charge);

      // The force field's own truncation convention, applied here as everywhere else in this route.
      shiftValue[s * numberOfAtoms + i] = cl_float(interactions(probe.sites[s].type, atomType).shift);
    }
  }

  // The molecule reaches beyond its own centre, so an atom that is out of range of the centre may still be
  // within range of a site. The reach is widened by the half-length of the molecule before the image range
  // is settled, which keeps the same inequality valid for every site at once.
  double coulombCutOff = grid.chargesIncluded ? potential->cutOff : 0.0;
  double longestCutOff = std::max(grid.cutOff, coulombCutOff);

  double3 widths = framework.unitCell.perpendicularWidths();
  double reach = longestCutOff + 0.5 * probe.length();
  grid.numberOfImageShells = int3(static_cast<std::int32_t>(std::floor(reach / widths.x + 0.5)),
                                  static_cast<std::int32_t>(std::floor(reach / widths.y + 0.5)),
                                  static_cast<std::int32_t>(std::floor(reach / widths.z + 0.5)));

  cl_int err;

  const char *source = MolecularEnergyGridOpenCL::molecularKernelSource;
  cl_program program = clCreateProgramWithSource(OpenCL::clContext.value(), 1, &source, nullptr, &err);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(
        std::format("MolecularEnergyGridOpenCL: OpenCL clCreateProgramWithSource failed at {}\n", __LINE__));
  }

  err = clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    std::size_t length;
    char buffer[2048];
    clGetProgramBuildInfo(program, OpenCL::clDeviceId.value(), CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &length);
    throw std::runtime_error(std::format("MolecularEnergyGridOpenCL: OpenCL failed to build program at {} (error: {})\n",
                                         __LINE__, std::string(buffer)));
  }

  cl_kernel kernel = clCreateKernel(program, "MolecularEnergyField", &err);
  if (err != CL_SUCCESS)
  {
    clReleaseProgram(program);
    throw std::runtime_error(std::format("MolecularEnergyGridOpenCL: OpenCL clCreateKernel failed at {}\n", __LINE__));
  }

  cl_int positionError, epsilonError, sigmaError, offsetError, freeError, minimumError;
  cl_int chargeError, siteChargeError, potentialError, shiftError;
  cl_mem positionBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY,
                                         sizeof(cl_float4) * positions.size(), nullptr, &positionError);
  cl_mem epsilonBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY,
                                        sizeof(cl_float) * epsilonTimesFour.size(), nullptr, &epsilonError);
  cl_mem sigmaBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY,
                                      sizeof(cl_float) * sigmaSquared.size(), nullptr, &sigmaError);
  cl_mem offsetBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY,
                                       sizeof(cl_float4) * siteOffsets.size(), nullptr, &offsetError);
  cl_mem freeBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_WRITE_ONLY, sizeof(cl_float) * numberOfVoxels,
                                     nullptr, &freeError);
  cl_mem minimumBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_WRITE_ONLY,
                                        sizeof(cl_float) * numberOfVoxels, nullptr, &minimumError);
  cl_mem chargeBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY,
                                       sizeof(cl_float) * chargeProduct.size(), nullptr, &chargeError);
  cl_mem shiftBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY,
                                      sizeof(cl_float) * shiftValue.size(), nullptr, &shiftError);
  cl_mem siteChargeBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY,
                                           sizeof(cl_float) * siteCharge.size(), nullptr, &siteChargeError);
  cl_mem potentialBuffer =
      clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY,
                     sizeof(cl_float) * (grid.chargesIncluded ? numberOfVoxels : 1), nullptr, &potentialError);
  if (positionError != CL_SUCCESS || epsilonError != CL_SUCCESS || sigmaError != CL_SUCCESS ||
      offsetError != CL_SUCCESS || freeError != CL_SUCCESS || minimumError != CL_SUCCESS ||
      chargeError != CL_SUCCESS || shiftError != CL_SUCCESS || siteChargeError != CL_SUCCESS || potentialError != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("MolecularEnergyGridOpenCL: OpenCL clCreateBuffer failed at {}\n", __LINE__));
  }

  err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), positionBuffer, CL_TRUE, 0,
                             sizeof(cl_float4) * positions.size(), positions.data(), 0, nullptr, nullptr);
  err |= clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), epsilonBuffer, CL_TRUE, 0,
                              sizeof(cl_float) * epsilonTimesFour.size(), epsilonTimesFour.data(), 0, nullptr, nullptr);
  err |= clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), sigmaBuffer, CL_TRUE, 0,
                              sizeof(cl_float) * sigmaSquared.size(), sigmaSquared.data(), 0, nullptr, nullptr);
  err |= clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), offsetBuffer, CL_TRUE, 0,
                              sizeof(cl_float4) * siteOffsets.size(), siteOffsets.data(), 0, nullptr, nullptr);
  err |= clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), chargeBuffer, CL_TRUE, 0,
                              sizeof(cl_float) * chargeProduct.size(), chargeProduct.data(), 0, nullptr, nullptr);
  err |= clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), shiftBuffer, CL_TRUE, 0,
                              sizeof(cl_float) * shiftValue.size(), shiftValue.data(), 0, nullptr, nullptr);
  err |= clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), siteChargeBuffer, CL_TRUE, 0,
                              sizeof(cl_float) * siteCharge.size(), siteCharge.data(), 0, nullptr, nullptr);
  if (grid.chargesIncluded)
  {
    err |= clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), potentialBuffer, CL_TRUE, 0,
                                sizeof(cl_float) * numberOfVoxels, potential->smoothPotential.data(), 0, nullptr,
                                nullptr);
  }
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("MolecularEnergyGridOpenCL: OpenCL clEnqueueWriteBuffer failed at {}\n", __LINE__));
  }

  cl_float4 cellRow0 = {{cl_float(cell[0][0]), cl_float(cell[1][0]), cl_float(cell[2][0]), 0.0f}};
  cl_float4 cellRow1 = {{cl_float(cell[0][1]), cl_float(cell[1][1]), cl_float(cell[2][1]), 0.0f}};
  cl_float4 cellRow2 = {{cl_float(cell[0][2]), cl_float(cell[1][2]), cl_float(cell[2][2]), 0.0f}};
  cl_float4 clWidths = {{cl_float(widths.x), cl_float(widths.y), cl_float(widths.z), 0.0f}};
  cl_int3 clShells = {{cl_int(grid.numberOfImageShells.x), cl_int(grid.numberOfImageShells.y),
                       cl_int(grid.numberOfImageShells.z), 0}};
  cl_int3 clGridSize = {{cl_int(gridSize.x), cl_int(gridSize.y), cl_int(gridSize.z), 0}};
  cl_int clNumberOfAtoms = static_cast<cl_int>(numberOfAtoms);
  cl_int clNumberOfSites = static_cast<cl_int>(numberOfSites);
  cl_int clNumberOfOrientations = static_cast<cl_int>(directions.size());
  cl_float clCutOff = cl_float(longestCutOff);
  cl_float clCutOffSquared = cl_float(grid.cutOff * grid.cutOff);
  cl_float clCoulombCutOffSquared = cl_float(coulombCutOff * coulombCutOff);
  cl_float clLongestCutOffSquared = cl_float(longestCutOff * longestCutOff);
  cl_float clAlpha = cl_float(grid.chargesIncluded ? potential->alpha : 0.0);
  cl_int clUseCharges = grid.chargesIncluded ? 1 : 0;
  cl_float clCeiling = cl_float(grid.ceiling);
  cl_float clBeta = cl_float(1.0 / (Units::KB * temperature));
  cl_float clKT = cl_float(Units::KB * temperature);

  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &positionBuffer);
  err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &epsilonBuffer);
  err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &sigmaBuffer);
  err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &offsetBuffer);
  err |= clSetKernelArg(kernel, 4, sizeof(cl_mem), &freeBuffer);
  err |= clSetKernelArg(kernel, 5, sizeof(cl_mem), &minimumBuffer);
  err |= clSetKernelArg(kernel, 6, sizeof(cl_float4), &cellRow0);
  err |= clSetKernelArg(kernel, 7, sizeof(cl_float4), &cellRow1);
  err |= clSetKernelArg(kernel, 8, sizeof(cl_float4), &cellRow2);
  err |= clSetKernelArg(kernel, 9, sizeof(cl_float4), &clWidths);
  err |= clSetKernelArg(kernel, 10, sizeof(cl_int), &clNumberOfAtoms);
  err |= clSetKernelArg(kernel, 11, sizeof(cl_int), &clNumberOfSites);
  err |= clSetKernelArg(kernel, 12, sizeof(cl_int), &clNumberOfOrientations);
  err |= clSetKernelArg(kernel, 13, sizeof(cl_int3), &clShells);
  err |= clSetKernelArg(kernel, 14, sizeof(cl_int3), &clGridSize);
  err |= clSetKernelArg(kernel, 15, sizeof(cl_float), &clCutOff);
  err |= clSetKernelArg(kernel, 16, sizeof(cl_float), &clCutOffSquared);
  err |= clSetKernelArg(kernel, 17, sizeof(cl_float), &clCeiling);
  err |= clSetKernelArg(kernel, 18, sizeof(cl_float), &clBeta);
  err |= clSetKernelArg(kernel, 19, sizeof(cl_float), &clKT);
  err |= clSetKernelArg(kernel, 20, sizeof(cl_mem), &chargeBuffer);
  err |= clSetKernelArg(kernel, 21, sizeof(cl_mem), &siteChargeBuffer);
  err |= clSetKernelArg(kernel, 22, sizeof(cl_mem), &potentialBuffer);
  err |= clSetKernelArg(kernel, 23, sizeof(cl_float), &clCoulombCutOffSquared);
  err |= clSetKernelArg(kernel, 24, sizeof(cl_float), &clLongestCutOffSquared);
  err |= clSetKernelArg(kernel, 25, sizeof(cl_float), &clAlpha);
  err |= clSetKernelArg(kernel, 26, sizeof(cl_int), &clUseCharges);
  err |= clSetKernelArg(kernel, 27, sizeof(cl_mem), &shiftBuffer);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("MolecularEnergyGridOpenCL: OpenCL clSetKernelArg failed at {}\n", __LINE__));
  }

  std::size_t globalWorkSize[3] = {static_cast<std::size_t>(gridSize.x), static_cast<std::size_t>(gridSize.y),
                                   static_cast<std::size_t>(gridSize.z)};
  err = clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), kernel, 3, nullptr, globalWorkSize, nullptr, 0, nullptr,
                               nullptr);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(
        std::format("MolecularEnergyGridOpenCL: OpenCL clEnqueueNDRangeKernel failed at {}\n", __LINE__));
  }

  err = clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), freeBuffer, CL_TRUE, 0, sizeof(cl_float) * numberOfVoxels,
                            grid.freeEnergy.data(), 0, nullptr, nullptr);
  err |= clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), minimumBuffer, CL_TRUE, 0,
                             sizeof(cl_float) * numberOfVoxels, grid.minimumEnergy.data(), 0, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("MolecularEnergyGridOpenCL: OpenCL clEnqueueReadBuffer failed at {}\n", __LINE__));
  }

  clFinish(OpenCL::clCommandQueue.value());

  clReleaseMemObject(positionBuffer);
  clReleaseMemObject(epsilonBuffer);
  clReleaseMemObject(sigmaBuffer);
  clReleaseMemObject(offsetBuffer);
  clReleaseMemObject(freeBuffer);
  clReleaseMemObject(minimumBuffer);
  clReleaseMemObject(chargeBuffer);
  clReleaseMemObject(shiftBuffer);
  clReleaseMemObject(siteChargeBuffer);
  clReleaseMemObject(potentialBuffer);
  clReleaseKernel(kernel);
  clReleaseProgram(program);

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  grid.seconds = elapsed.count();

  return grid;
}


std::vector<std::int32_t> MolecularEnergyGridOpenCL::strongestAtoms(const PairInteractions &interactions,
                                                                    const Crystal &framework,
                                                                    const LinearProbe &probe, uint3 gridSize,
                                                                    std::size_t numberOfOrientations,
                                                                    double temperature,
                                                                    const ElectrostaticPotentialGrid *potential)
{
  if (!OpenCL::clContext.has_value() || !OpenCL::clDeviceId.has_value() || !OpenCL::clCommandQueue.has_value())
  {
    throw std::runtime_error("MolecularEnergyGridOpenCL: no OpenCL device found, the grid route needs a GPU\n");
  }
  if (gridSize.x == 0 || gridSize.y == 0 || gridSize.z == 0)
  {
    throw std::runtime_error("MolecularEnergyGridOpenCL: the grid must have at least one point along each axis\n");
  }
  if (probe.sites.empty()) throw std::runtime_error("MolecularEnergyGridOpenCL: the probe has no sites\n");
  if (numberOfOrientations == 0)
  {
    throw std::runtime_error("MolecularEnergyGridOpenCL: at least one orientation is needed\n");
  }
  if (temperature <= 0.0)
  {
    throw std::runtime_error("MolecularEnergyGridOpenCL: the weighting over orientations needs a temperature\n");
  }

  const std::size_t numberOfVoxels = gridSize.x * gridSize.y * gridSize.z;
  std::vector<std::int32_t> labels(numberOfVoxels, -1);

  std::vector<double3> fractionalPositions = framework.fractionalPositions;
  if (fractionalPositions.empty()) return labels;

  const bool useCharges = probe.isCharged() && potential != nullptr && potential->numberOfVoxels() > 0;
  if (useCharges && (potential->gridSize.x != gridSize.x || potential->gridSize.y != gridSize.y ||
                     potential->gridSize.z != gridSize.z))
  {
    throw std::runtime_error("MolecularEnergyGridOpenCL: the potential is on a different grid than this one\n");
  }

  std::vector<double3> directions = orientationSet(numberOfOrientations, probe.headTailSymmetric);

  const std::size_t numberOfAtoms = fractionalPositions.size();
  const std::size_t numberOfSites = probe.sites.size();
  const std::size_t numberOfDirections = directions.size();

  double3x3 cell = framework.unitCell.cell;
  double3x3 inverseCell = framework.unitCell.inverseCell;

  std::vector<cl_float4> siteOffsets(numberOfDirections * numberOfSites);
  for (std::size_t o = 0; o < numberOfDirections; ++o)
  {
    for (std::size_t s = 0; s < numberOfSites; ++s)
    {
      double3 displacement = inverseCell * (probe.sites[s].offset * directions[o]);
      siteOffsets[o * numberOfSites + s] = {
          {cl_float(displacement.x), cl_float(displacement.y), cl_float(displacement.z), 0.0f}};
    }
  }

  std::vector<cl_float4> positions(numberOfAtoms);
  std::vector<cl_float> epsilonTimesFour(numberOfSites * numberOfAtoms);
  std::vector<cl_float> sigmaSquared(numberOfSites * numberOfAtoms);
  std::vector<cl_float> chargeProduct(numberOfSites * numberOfAtoms, 0.0f);
  std::vector<cl_float> shiftValue(numberOfSites * numberOfAtoms, 0.0f);
  std::vector<cl_float> siteCharge(numberOfSites);
  for (std::size_t s = 0; s < numberOfSites; ++s) siteCharge[s] = cl_float(probe.sites[s].charge);
  for (std::size_t i = 0; i < numberOfAtoms; ++i)
  {
    positions[i] = {{cl_float(fractionalPositions[i].x), cl_float(fractionalPositions[i].y),
                     cl_float(fractionalPositions[i].z), 0.0f}};

    std::size_t atomType = framework.atoms[i].type;
    for (std::size_t s = 0; s < numberOfSites; ++s)
    {
      double sigma = interactions(probe.sites[s].type, atomType).sizeParameter;
      double epsilon = interactions(probe.sites[s].type, atomType).strengthParameter;
      epsilonTimesFour[s * numberOfAtoms + i] = cl_float(4.0 * epsilon);
      sigmaSquared[s * numberOfAtoms + i] = cl_float(sigma * sigma);
      chargeProduct[s * numberOfAtoms + i] =
          cl_float(Units::CoulombicConversionFactor * probe.sites[s].charge * framework.atoms[i].charge);
      shiftValue[s * numberOfAtoms + i] = cl_float(interactions(probe.sites[s].type, atomType).shift);
    }
  }

  const double coulombCutOff = useCharges ? potential->cutOff : 0.0;
  const double vdwCutOff = interactions.cutOffVDW;
  const double longestCutOff = std::max(vdwCutOff, coulombCutOff);

  double3 widths = framework.unitCell.perpendicularWidths();
  double reach = longestCutOff + 0.5 * probe.length();
  const int3 shells(static_cast<std::int32_t>(std::floor(reach / widths.x + 0.5)),
                    static_cast<std::int32_t>(std::floor(reach / widths.y + 0.5)),
                    static_cast<std::int32_t>(std::floor(reach / widths.z + 0.5)));

  // The energies of every orientation have to be in hand before any of them can be weighed, and holding them
  // for the whole cell at once would be a gigabyte at a middling grid. A slab of z-planes at a time is held
  // instead, sized to a few hundred megabytes, and the pair of kernels is run over each slab in turn.
  const std::size_t planeVoxels = static_cast<std::size_t>(gridSize.x) * gridSize.y;
  const std::size_t perPlane = planeVoxels * numberOfDirections * sizeof(cl_float);
  const std::size_t budget = std::size_t{192} * 1024 * 1024;
  std::size_t planesPerSlab = std::max(std::size_t{1}, budget / std::max(perPlane, std::size_t{1}));
  planesPerSlab = std::min(planesPerSlab, static_cast<std::size_t>(gridSize.z));

  cl_int err;

  const char *source = MolecularEnergyGridOpenCL::tessellationKernelSource;
  cl_program program = clCreateProgramWithSource(OpenCL::clContext.value(), 1, &source, nullptr, &err);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(
        std::format("MolecularEnergyGridOpenCL: OpenCL clCreateProgramWithSource failed at {}\n", __LINE__));
  }

  err = clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    std::size_t length;
    char buffer[4096];
    clGetProgramBuildInfo(program, OpenCL::clDeviceId.value(), CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &length);
    clReleaseProgram(program);
    throw std::runtime_error(std::format(
        "MolecularEnergyGridOpenCL: OpenCL failed to build program at {} (error: {})\n", __LINE__,
        std::string(buffer)));
  }

  cl_kernel energyKernel = clCreateKernel(program, "MolecularOrientationEnergies", &err);
  cl_int labelError;
  cl_kernel labelKernel = clCreateKernel(program, "MolecularStrongestAtom", &labelError);
  if (err != CL_SUCCESS || labelError != CL_SUCCESS)
  {
    clReleaseProgram(program);
    throw std::runtime_error(std::format("MolecularEnergyGridOpenCL: OpenCL clCreateKernel failed at {}\n", __LINE__));
  }

  auto buffer = [&](cl_mem_flags flags, std::size_t bytes)
  {
    cl_int bufferError;
    cl_mem memory = clCreateBuffer(OpenCL::clContext.value(), flags, std::max(bytes, std::size_t{1}), nullptr,
                                   &bufferError);
    if (bufferError != CL_SUCCESS)
    {
      throw std::runtime_error(
          std::format("MolecularEnergyGridOpenCL: OpenCL clCreateBuffer failed at {}\n", __LINE__));
    }
    return memory;
  };

  cl_mem positionBuffer = buffer(CL_MEM_READ_ONLY, sizeof(cl_float4) * positions.size());
  cl_mem epsilonBuffer = buffer(CL_MEM_READ_ONLY, sizeof(cl_float) * epsilonTimesFour.size());
  cl_mem sigmaBuffer = buffer(CL_MEM_READ_ONLY, sizeof(cl_float) * sigmaSquared.size());
  cl_mem offsetBuffer = buffer(CL_MEM_READ_ONLY, sizeof(cl_float4) * siteOffsets.size());
  cl_mem chargeBuffer = buffer(CL_MEM_READ_ONLY, sizeof(cl_float) * chargeProduct.size());
  cl_mem shiftBuffer = buffer(CL_MEM_READ_ONLY, sizeof(cl_float) * shiftValue.size());
  cl_mem siteChargeBuffer = buffer(CL_MEM_READ_ONLY, sizeof(cl_float) * siteCharge.size());
  cl_mem potentialBuffer = buffer(CL_MEM_READ_ONLY, sizeof(cl_float) * (useCharges ? numberOfVoxels : 1));
  cl_mem energyBuffer = buffer(CL_MEM_READ_WRITE, sizeof(cl_float) * planeVoxels * planesPerSlab * numberOfDirections);
  cl_mem labelBuffer = buffer(CL_MEM_WRITE_ONLY, sizeof(cl_int) * numberOfVoxels);

  err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), positionBuffer, CL_TRUE, 0,
                             sizeof(cl_float4) * positions.size(), positions.data(), 0, nullptr, nullptr);
  err |= clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), epsilonBuffer, CL_TRUE, 0,
                              sizeof(cl_float) * epsilonTimesFour.size(), epsilonTimesFour.data(), 0, nullptr,
                              nullptr);
  err |= clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), sigmaBuffer, CL_TRUE, 0,
                              sizeof(cl_float) * sigmaSquared.size(), sigmaSquared.data(), 0, nullptr, nullptr);
  err |= clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), offsetBuffer, CL_TRUE, 0,
                              sizeof(cl_float4) * siteOffsets.size(), siteOffsets.data(), 0, nullptr, nullptr);
  err |= clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), chargeBuffer, CL_TRUE, 0,
                              sizeof(cl_float) * chargeProduct.size(), chargeProduct.data(), 0, nullptr, nullptr);
  err |= clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), shiftBuffer, CL_TRUE, 0,
                              sizeof(cl_float) * shiftValue.size(), shiftValue.data(), 0, nullptr, nullptr);
  err |= clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), siteChargeBuffer, CL_TRUE, 0,
                              sizeof(cl_float) * siteCharge.size(), siteCharge.data(), 0, nullptr, nullptr);
  if (useCharges)
  {
    err |= clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), potentialBuffer, CL_TRUE, 0,
                                sizeof(cl_float) * numberOfVoxels, potential->smoothPotential.data(), 0, nullptr,
                                nullptr);
  }
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(
        std::format("MolecularEnergyGridOpenCL: OpenCL clEnqueueWriteBuffer failed at {}\n", __LINE__));
  }

  cl_float4 cellRow0 = {{cl_float(cell[0][0]), cl_float(cell[1][0]), cl_float(cell[2][0]), 0.0f}};
  cl_float4 cellRow1 = {{cl_float(cell[0][1]), cl_float(cell[1][1]), cl_float(cell[2][1]), 0.0f}};
  cl_float4 cellRow2 = {{cl_float(cell[0][2]), cl_float(cell[1][2]), cl_float(cell[2][2]), 0.0f}};
  cl_float4 clWidths = {{cl_float(widths.x), cl_float(widths.y), cl_float(widths.z), 0.0f}};
  cl_int3 clShells = {{cl_int(shells.x), cl_int(shells.y), cl_int(shells.z), 0}};
  cl_int3 clGridSize = {{cl_int(gridSize.x), cl_int(gridSize.y), cl_int(gridSize.z), 0}};
  cl_int clNumberOfAtoms = static_cast<cl_int>(numberOfAtoms);
  cl_int clNumberOfSites = static_cast<cl_int>(numberOfSites);
  cl_int clNumberOfOrientations = static_cast<cl_int>(numberOfDirections);
  cl_float clCutOff = cl_float(longestCutOff);
  cl_float clCutOffSquared = cl_float(vdwCutOff * vdwCutOff);
  cl_float clCoulombCutOffSquared = cl_float(coulombCutOff * coulombCutOff);
  cl_float clLongestCutOffSquared = cl_float(longestCutOff * longestCutOff);
  cl_float clAlpha = cl_float(useCharges ? potential->alpha : 0.0);
  cl_int clUseCharges = useCharges ? 1 : 0;
  cl_float clCeiling = cl_float(probeEnergyCeilingInKelvin * Units::KelvinToEnergy);
  cl_float clBeta = cl_float(1.0 / (Units::KB * temperature));

  err = clSetKernelArg(energyKernel, 0, sizeof(cl_mem), &positionBuffer);
  err |= clSetKernelArg(energyKernel, 1, sizeof(cl_mem), &epsilonBuffer);
  err |= clSetKernelArg(energyKernel, 2, sizeof(cl_mem), &sigmaBuffer);
  err |= clSetKernelArg(energyKernel, 3, sizeof(cl_mem), &offsetBuffer);
  err |= clSetKernelArg(energyKernel, 4, sizeof(cl_mem), &energyBuffer);
  err |= clSetKernelArg(energyKernel, 5, sizeof(cl_float4), &cellRow0);
  err |= clSetKernelArg(energyKernel, 6, sizeof(cl_float4), &cellRow1);
  err |= clSetKernelArg(energyKernel, 7, sizeof(cl_float4), &cellRow2);
  err |= clSetKernelArg(energyKernel, 8, sizeof(cl_float4), &clWidths);
  err |= clSetKernelArg(energyKernel, 9, sizeof(cl_int), &clNumberOfAtoms);
  err |= clSetKernelArg(energyKernel, 10, sizeof(cl_int), &clNumberOfSites);
  err |= clSetKernelArg(energyKernel, 11, sizeof(cl_int), &clNumberOfOrientations);
  err |= clSetKernelArg(energyKernel, 12, sizeof(cl_int3), &clShells);
  err |= clSetKernelArg(energyKernel, 13, sizeof(cl_int3), &clGridSize);
  err |= clSetKernelArg(energyKernel, 14, sizeof(cl_float), &clCutOff);
  err |= clSetKernelArg(energyKernel, 15, sizeof(cl_float), &clCutOffSquared);
  err |= clSetKernelArg(energyKernel, 16, sizeof(cl_float), &clCeiling);
  err |= clSetKernelArg(energyKernel, 17, sizeof(cl_mem), &chargeBuffer);
  err |= clSetKernelArg(energyKernel, 18, sizeof(cl_mem), &siteChargeBuffer);
  err |= clSetKernelArg(energyKernel, 19, sizeof(cl_mem), &potentialBuffer);
  err |= clSetKernelArg(energyKernel, 20, sizeof(cl_float), &clCoulombCutOffSquared);
  err |= clSetKernelArg(energyKernel, 21, sizeof(cl_float), &clLongestCutOffSquared);
  err |= clSetKernelArg(energyKernel, 22, sizeof(cl_float), &clAlpha);
  err |= clSetKernelArg(energyKernel, 23, sizeof(cl_int), &clUseCharges);
  err |= clSetKernelArg(energyKernel, 24, sizeof(cl_mem), &shiftBuffer);

  err |= clSetKernelArg(labelKernel, 0, sizeof(cl_mem), &positionBuffer);
  err |= clSetKernelArg(labelKernel, 1, sizeof(cl_mem), &epsilonBuffer);
  err |= clSetKernelArg(labelKernel, 2, sizeof(cl_mem), &sigmaBuffer);
  err |= clSetKernelArg(labelKernel, 3, sizeof(cl_mem), &offsetBuffer);
  err |= clSetKernelArg(labelKernel, 4, sizeof(cl_mem), &energyBuffer);
  err |= clSetKernelArg(labelKernel, 5, sizeof(cl_mem), &labelBuffer);
  err |= clSetKernelArg(labelKernel, 6, sizeof(cl_float4), &cellRow0);
  err |= clSetKernelArg(labelKernel, 7, sizeof(cl_float4), &cellRow1);
  err |= clSetKernelArg(labelKernel, 8, sizeof(cl_float4), &cellRow2);
  err |= clSetKernelArg(labelKernel, 9, sizeof(cl_float4), &clWidths);
  err |= clSetKernelArg(labelKernel, 10, sizeof(cl_int), &clNumberOfAtoms);
  err |= clSetKernelArg(labelKernel, 11, sizeof(cl_int), &clNumberOfSites);
  err |= clSetKernelArg(labelKernel, 12, sizeof(cl_int), &clNumberOfOrientations);
  err |= clSetKernelArg(labelKernel, 13, sizeof(cl_int3), &clShells);
  err |= clSetKernelArg(labelKernel, 14, sizeof(cl_int3), &clGridSize);
  err |= clSetKernelArg(labelKernel, 15, sizeof(cl_float), &clCutOff);
  err |= clSetKernelArg(labelKernel, 16, sizeof(cl_float), &clCutOffSquared);
  err |= clSetKernelArg(labelKernel, 17, sizeof(cl_float), &clCeiling);
  err |= clSetKernelArg(labelKernel, 18, sizeof(cl_float), &clBeta);
  err |= clSetKernelArg(labelKernel, 19, sizeof(cl_mem), &chargeBuffer);
  err |= clSetKernelArg(labelKernel, 20, sizeof(cl_float), &clCoulombCutOffSquared);
  err |= clSetKernelArg(labelKernel, 21, sizeof(cl_float), &clLongestCutOffSquared);
  err |= clSetKernelArg(labelKernel, 22, sizeof(cl_float), &clAlpha);
  err |= clSetKernelArg(labelKernel, 23, sizeof(cl_int), &clUseCharges);
  err |= clSetKernelArg(labelKernel, 24, sizeof(cl_mem), &shiftBuffer);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("MolecularEnergyGridOpenCL: OpenCL clSetKernelArg failed at {}\n", __LINE__));
  }

  for (std::size_t first = 0; first < gridSize.z; first += planesPerSlab)
  {
    const std::size_t planes = std::min(planesPerSlab, static_cast<std::size_t>(gridSize.z) - first);
    cl_int clFirstPlane = static_cast<cl_int>(first);
    cl_int clNumberOfPlanes = static_cast<cl_int>(planes);

    err = clSetKernelArg(energyKernel, 25, sizeof(cl_int), &clFirstPlane);
    err |= clSetKernelArg(energyKernel, 26, sizeof(cl_int), &clNumberOfPlanes);
    err |= clSetKernelArg(labelKernel, 25, sizeof(cl_int), &clFirstPlane);
    err |= clSetKernelArg(labelKernel, 26, sizeof(cl_int), &clNumberOfPlanes);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(
          std::format("MolecularEnergyGridOpenCL: OpenCL clSetKernelArg failed at {}\n", __LINE__));
    }

    std::size_t globalWorkSize[3] = {static_cast<std::size_t>(gridSize.x), static_cast<std::size_t>(gridSize.y),
                                     planes};
    err = clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), energyKernel, 3, nullptr, globalWorkSize, nullptr, 0,
                                 nullptr, nullptr);
    err |= clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), labelKernel, 3, nullptr, globalWorkSize, nullptr, 0,
                                  nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(
          std::format("MolecularEnergyGridOpenCL: OpenCL clEnqueueNDRangeKernel failed at {}\n", __LINE__));
    }

    clFinish(OpenCL::clCommandQueue.value());
  }

  err = clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), labelBuffer, CL_TRUE, 0, sizeof(cl_int) * numberOfVoxels,
                            labels.data(), 0, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(
        std::format("MolecularEnergyGridOpenCL: OpenCL clEnqueueReadBuffer failed at {}\n", __LINE__));
  }

  clFinish(OpenCL::clCommandQueue.value());

  clReleaseMemObject(positionBuffer);
  clReleaseMemObject(epsilonBuffer);
  clReleaseMemObject(sigmaBuffer);
  clReleaseMemObject(offsetBuffer);
  clReleaseMemObject(chargeBuffer);
  clReleaseMemObject(shiftBuffer);
  clReleaseMemObject(siteChargeBuffer);
  clReleaseMemObject(potentialBuffer);
  clReleaseMemObject(energyBuffer);
  clReleaseMemObject(labelBuffer);
  clReleaseKernel(energyKernel);
  clReleaseKernel(labelKernel);
  clReleaseProgram(program);

  return labels;
}


const char *MolecularEnergyGridOpenCL::molecularKernelSource = R"foo(
// The far half of the electrostatic sum, read off between grid points. It is built from waves no shorter
// than a few Ångström, so it bends very little over one spacing and reading it off straight lines costs
// almost nothing. The near half, which has the 1/r in it and bends a great deal, is never read off but
// summed exactly below.
static float smoothPotentialAt(__global const float *smoothPotential, float3 s, int3 gridSize)
{
  s -= floor(s);

  float3 scaled = s * (float3)((float)(gridSize.x), (float)(gridSize.y), (float)(gridSize.z));
  int3 low = convert_int3(floor(scaled));
  float3 frac = scaled - convert_float3(low);

  low.x = low.x % gridSize.x;
  low.y = low.y % gridSize.y;
  low.z = low.z % gridSize.z;

  int hx = (low.x + 1) % gridSize.x;
  int hy = (low.y + 1) % gridSize.y;
  int hz = (low.z + 1) % gridSize.z;

  float c00 = mix(smoothPotential[(low.z * gridSize.y + low.y) * gridSize.x + low.x],
                  smoothPotential[(low.z * gridSize.y + low.y) * gridSize.x + hx], frac.x);
  float c10 = mix(smoothPotential[(low.z * gridSize.y + hy) * gridSize.x + low.x],
                  smoothPotential[(low.z * gridSize.y + hy) * gridSize.x + hx], frac.x);
  float c01 = mix(smoothPotential[(hz * gridSize.y + low.y) * gridSize.x + low.x],
                  smoothPotential[(hz * gridSize.y + low.y) * gridSize.x + hx], frac.x);
  float c11 = mix(smoothPotential[(hz * gridSize.y + hy) * gridSize.x + low.x],
                  smoothPotential[(hz * gridSize.y + hy) * gridSize.x + hx], frac.x);

  return mix(mix(c00, c10, frac.y), mix(c01, c11, frac.y), frac.z);
}

__kernel void MolecularEnergyField(__global const float4 *position,
                                   __global const float *epsilonTimesFour,
                                   __global const float *sigmaSquared,
                                   __global const float4 *siteOffset,
                                   __global float *freeEnergy,
                                   __global float *minimumEnergy,
                                   const float4 cellRow0,
                                   const float4 cellRow1,
                                   const float4 cellRow2,
                                   const float4 widths,
                                   const int numberOfAtoms,
                                   const int numberOfSites,
                                   const int numberOfOrientations,
                                   const int3 shells,
                                   const int3 gridSize,
                                   const float cutOff,
                                   const float cutOffSquared,
                                   const float ceiling,
                                   const float beta,
                                   const float kT,
                                   __global const float *chargeProduct,
                                   __global const float *siteCharge,
                                   __global const float *smoothPotential,
                                   const float coulombCutOffSquared,
                                   const float longestCutOffSquared,
                                   const float alpha,
                                   const int useCharges,
                                   __global const float *shiftValue)
{
  int ix = get_global_id(0);
  int iy = get_global_id(1);
  int iz = get_global_id(2);

  if(ix >= gridSize.x || iy >= gridSize.y || iz >= gridSize.z) return;

  // Endpoint-exclusive sampling, as in the other fields of this route, so that what is read off one may be
  // read off another. The point is the molecule's centre of mass; its sites hang off it.
  float4 centre = (float4)((float)(ix) / (float)(gridSize.x),
                           (float)(iy) / (float)(gridSize.y),
                           (float)(iz) / (float)(gridSize.z),
                           0.0f);

  // The average of exp(-U/kT) over orientations is accumulated as it goes, rescaled whenever a deeper
  // orientation turns up. Taken plainly the exponential would overflow at a strong site long before the
  // physics did; carried this way every term is at most one and the least is harmless.
  float least = MAXFLOAT;
  float sum = 0.0f;

  for(int o = 0; o < numberOfOrientations; o++)
  {
    float total = 0.0f;

    for(int site = 0; site < numberOfSites; site++)
    {
      float4 s = centre + siteOffset[o * numberOfSites + site];
      s.w = 0.0f;

      __global const float *epsilonForSite = epsilonTimesFour + site * numberOfAtoms;
      __global const float *sigmaForSite = sigmaSquared + site * numberOfAtoms;
      __global const float *chargeForSite = chargeProduct + site * numberOfAtoms;
      __global const float *shiftForSite = shiftValue + site * numberOfAtoms;

      if(useCharges != 0) total += siteCharge[site] * smoothPotentialAt(smoothPotential, s.xyz, gridSize);

      for(int iatom = 0; iatom < numberOfAtoms; iatom++)
      {
        float4 ds = s - position[iatom];
        ds -= rint(ds);
        ds.w = 0.0f;

        float epsilon4 = epsilonForSite[iatom];
        float sigma2 = sigmaForSite[iatom];
        float charges = chargeForSite[iatom];

        for(int a = -shells.x; a <= shells.x; a++)
        {
          for(int b = -shells.y; b <= shells.y; b++)
          {
            for(int c = -shells.z; c <= shells.z; c++)
            {
              float4 t = ds + (float4)((float)(a), (float)(b), (float)(c), 0.0f);

              // A fractional difference is at least |t_k| w_k long whatever the shape of the cell, so an
              // image beyond the cutoff is thrown out here for three multiplies rather than three dot
              // products and a square root's worth of work afterwards.
              float far = fmax(fmax(fabs(t.x) * widths.x, fabs(t.y) * widths.y), fabs(t.z) * widths.z);
              if(far > cutOff) continue;


              float4 dr;
              dr.x = dot(cellRow0, t);
              dr.y = dot(cellRow1, t);
              dr.z = dot(cellRow2, t);
              dr.w = 0.0f;

              float rr = dot(dr, dr);
              if(rr >= longestCutOffSquared) continue;

              // A site can land on an atom's centre, where the pair energy has no value. Holding the
              // separation off zero keeps the sum a number; such an orientation is buried far above the
              // ceiling and contributes nothing to the average either way.
              rr = fmax(rr, 1.0e-6f);

              if(rr < cutOffSquared)
              {
                float ratio = sigma2 / rr;
                float ratio3 = ratio * ratio * ratio;
                total += min(epsilon4 * ratio3 * (ratio3 - 1.0f) - shiftForSite[iatom], ceiling);
              }

              // The near half of the electrostatic sum. The dispersion has already walked to this neighbour
              // and worked out how far away it is, so all this costs is the erfc.
              if(useCharges != 0 && rr < coulombCutOffSquared && charges != 0.0f)
              {
                float r = sqrt(rr);
                total += min(charges * erfc(alpha * r) / r, ceiling);
              }
            }
          }
        }
      }
    }

    total = min(total, ceiling);

    if(total < least)
    {
      // A deeper orientation resets what the sum is measured from, so everything gathered so far is scaled
      // down to the new floor before this one is added at its full weight.
      sum = sum * exp(-beta * (least - total)) + 1.0f;
      least = total;
    }
    else
    {
      sum += exp(-beta * (total - least));
    }
  }

  int index = (iz * gridSize.y + iy) * gridSize.x + ix;
  minimumEnergy[index] = least;

  // -kT ln <exp(-U/kT)>, written from the floor upwards. The average can only be at most one when measured
  // from the least, so the logarithm is never positive and the free energy never falls below the minimum:
  // the gap between them is what the molecule pays for having to be turned a particular way.
  freeEnergy[index] = min(least - kT * log(sum / (float)(numberOfOrientations)), ceiling);
}
)foo";


const char *MolecularEnergyGridOpenCL::tessellationKernelSource = R"foo(
static float smoothPotentialAt(__global const float *smoothPotential, float3 s, int3 gridSize)
{
  s -= floor(s);

  float3 scaled = s * (float3)((float)(gridSize.x), (float)(gridSize.y), (float)(gridSize.z));
  int3 low = convert_int3(floor(scaled));
  float3 frac = scaled - convert_float3(low);

  low.x = low.x % gridSize.x;
  low.y = low.y % gridSize.y;
  low.z = low.z % gridSize.z;

  int hx = (low.x + 1) % gridSize.x;
  int hy = (low.y + 1) % gridSize.y;
  int hz = (low.z + 1) % gridSize.z;

  float c00 = mix(smoothPotential[(low.z * gridSize.y + low.y) * gridSize.x + low.x],
                  smoothPotential[(low.z * gridSize.y + low.y) * gridSize.x + hx], frac.x);
  float c10 = mix(smoothPotential[(low.z * gridSize.y + hy) * gridSize.x + low.x],
                  smoothPotential[(low.z * gridSize.y + hy) * gridSize.x + hx], frac.x);
  float c01 = mix(smoothPotential[(hz * gridSize.y + low.y) * gridSize.x + low.x],
                  smoothPotential[(hz * gridSize.y + low.y) * gridSize.x + hx], frac.x);
  float c11 = mix(smoothPotential[(hz * gridSize.y + hy) * gridSize.x + low.x],
                  smoothPotential[(hz * gridSize.y + hy) * gridSize.x + hx], frac.x);

  return mix(mix(c00, c10, frac.y), mix(c01, c11, frac.y), frac.z);
}

// What one framework atom, over all of its images in range, does to the molecule held one way. It is the
// part of the energy that can be laid at a single atom's door, and it is what the label is decided on. The
// nearest of those images is handed back beside it, for the wall, where no atom pulls at all and the energy
// has nothing to say.
static float contributionOfAtom(__global const float4 *position,
                                __global const float *epsilonTimesFour,
                                __global const float *sigmaSquared,
                                __global const float *chargeProduct,
                                __global const float *shiftValue,
                                __global const float4 *siteOffset,
                                float4 centre,
                                int orientation,
                                int iatom,
                                const float4 cellRow0,
                                const float4 cellRow1,
                                const float4 cellRow2,
                                const float4 widths,
                                const int numberOfAtoms,
                                const int numberOfSites,
                                const int3 shells,
                                const float cutOff,
                                const float cutOffSquared,
                                const float coulombCutOffSquared,
                                const float longestCutOffSquared,
                                const float alpha,
                                const int useCharges,
                                const float ceiling,
                                float *closestSquared)
{
  float contribution = 0.0f;
  float closest = MAXFLOAT;

  for(int site = 0; site < numberOfSites; site++)
  {
    float4 s = centre + siteOffset[orientation * numberOfSites + site];
    s.w = 0.0f;

    float4 ds = s - position[iatom];
    ds -= rint(ds);
    ds.w = 0.0f;

    float epsilon4 = epsilonTimesFour[site * numberOfAtoms + iatom];
    float sigma2 = sigmaSquared[site * numberOfAtoms + iatom];
    float charges = chargeProduct[site * numberOfAtoms + iatom];
    float shift = shiftValue[site * numberOfAtoms + iatom];

    for(int a = -shells.x; a <= shells.x; a++)
    {
      for(int b = -shells.y; b <= shells.y; b++)
      {
        for(int c = -shells.z; c <= shells.z; c++)
        {
          float4 t = ds + (float4)((float)(a), (float)(b), (float)(c), 0.0f);

          float far = fmax(fmax(fabs(t.x) * widths.x, fabs(t.y) * widths.y), fabs(t.z) * widths.z);
          if(far > cutOff) continue;

          float4 dr;
          dr.x = dot(cellRow0, t);
          dr.y = dot(cellRow1, t);
          dr.z = dot(cellRow2, t);
          dr.w = 0.0f;

          float rr = dot(dr, dr);
          if(rr >= longestCutOffSquared) continue;

          rr = fmax(rr, 1.0e-6f);
          closest = fmin(closest, rr);

          if(rr < cutOffSquared)
          {
            float ratio = sigma2 / rr;
            float ratio3 = ratio * ratio * ratio;
            contribution += min(epsilon4 * ratio3 * (ratio3 - 1.0f) - shift, ceiling);
          }

          if(useCharges != 0 && rr < coulombCutOffSquared && charges != 0.0f)
          {
            float r = sqrt(rr);
            contribution += min(charges * erfc(alpha * r) / r, ceiling);
          }
        }
      }
    }
  }

  *closestSquared = closest;
  return contribution;
}

// The energy of every orientation at every point of one slab, written out rather than reduced. It is the
// same arithmetic as MolecularEnergyField and is kept apart from it because what is wanted here is the whole
// set and not the two reductions of it.
__kernel void MolecularOrientationEnergies(__global const float4 *position,
                                           __global const float *epsilonTimesFour,
                                           __global const float *sigmaSquared,
                                           __global const float4 *siteOffset,
                                           __global float *orientationEnergy,
                                           const float4 cellRow0,
                                           const float4 cellRow1,
                                           const float4 cellRow2,
                                           const float4 widths,
                                           const int numberOfAtoms,
                                           const int numberOfSites,
                                           const int numberOfOrientations,
                                           const int3 shells,
                                           const int3 gridSize,
                                           const float cutOff,
                                           const float cutOffSquared,
                                           const float ceiling,
                                           __global const float *chargeProduct,
                                           __global const float *siteCharge,
                                           __global const float *smoothPotential,
                                           const float coulombCutOffSquared,
                                           const float longestCutOffSquared,
                                           const float alpha,
                                           const int useCharges,
                                           __global const float *shiftValue,
                                           const int firstPlane,
                                           const int numberOfPlanes)
{
  int ix = get_global_id(0);
  int iy = get_global_id(1);
  int slab = get_global_id(2);

  if(ix >= gridSize.x || iy >= gridSize.y || slab >= numberOfPlanes) return;

  int iz = firstPlane + slab;
  if(iz >= gridSize.z) return;

  float4 centre = (float4)((float)(ix) / (float)(gridSize.x),
                           (float)(iy) / (float)(gridSize.y),
                           (float)(iz) / (float)(gridSize.z),
                           0.0f);

  int slabVoxel = (slab * gridSize.y + iy) * gridSize.x + ix;

  for(int o = 0; o < numberOfOrientations; o++)
  {
    float total = 0.0f;
    for(int iatom = 0; iatom < numberOfAtoms; iatom++)
    {
      float closest;
      total += contributionOfAtom(position, epsilonTimesFour, sigmaSquared, chargeProduct, shiftValue, siteOffset,
                                  centre, o, iatom, cellRow0, cellRow1, cellRow2, widths, numberOfAtoms,
                                  numberOfSites, shells, cutOff, cutOffSquared, coulombCutOffSquared,
                                  longestCutOffSquared, alpha, useCharges, ceiling, &closest);
    }

    // The reciprocal half of the electrostatic sum belongs to no one atom, so it takes no part in the label.
    // It does belong here, since this is what decides how the molecule sits.
    if(useCharges != 0)
    {
      for(int site = 0; site < numberOfSites; site++)
      {
        float4 s = centre + siteOffset[o * numberOfSites + site];
        total += siteCharge[site] * smoothPotentialAt(smoothPotential, s.xyz, gridSize);
      }
    }

    orientationEnergy[slabVoxel * numberOfOrientations + o] = min(total, ceiling);
  }
}

// The atom that pulls hardest, averaged over the orientations the molecule actually takes up.
__kernel void MolecularStrongestAtom(__global const float4 *position,
                                     __global const float *epsilonTimesFour,
                                     __global const float *sigmaSquared,
                                     __global const float4 *siteOffset,
                                     __global const float *orientationEnergy,
                                     __global int *label,
                                     const float4 cellRow0,
                                     const float4 cellRow1,
                                     const float4 cellRow2,
                                     const float4 widths,
                                     const int numberOfAtoms,
                                     const int numberOfSites,
                                     const int numberOfOrientations,
                                     const int3 shells,
                                     const int3 gridSize,
                                     const float cutOff,
                                     const float cutOffSquared,
                                     const float ceiling,
                                     const float beta,
                                     __global const float *chargeProduct,
                                     const float coulombCutOffSquared,
                                     const float longestCutOffSquared,
                                     const float alpha,
                                     const int useCharges,
                                     __global const float *shiftValue,
                                     const int firstPlane,
                                     const int numberOfPlanes)
{
  int ix = get_global_id(0);
  int iy = get_global_id(1);
  int slab = get_global_id(2);

  if(ix >= gridSize.x || iy >= gridSize.y || slab >= numberOfPlanes) return;

  int iz = firstPlane + slab;
  if(iz >= gridSize.z) return;

  float4 centre = (float4)((float)(ix) / (float)(gridSize.x),
                           (float)(iy) / (float)(gridSize.y),
                           (float)(iz) / (float)(gridSize.z),
                           0.0f);

  int slabVoxel = (slab * gridSize.y + iy) * gridSize.x + ix;
  __global const float *energies = orientationEnergy + (size_t)(slabVoxel) * (size_t)(numberOfOrientations);

  // Measured from the deepest orientation, so the largest weight is one whatever the site is worth and
  // nothing overflows. Inside a wall every orientation is at the ceiling and the weights come out equal,
  // which is the right answer there: no orientation is preferred because none is possible.
  float least = MAXFLOAT;
  for(int o = 0; o < numberOfOrientations; o++) least = fmin(least, energies[o]);

  float normalisation = 0.0f;
  for(int o = 0; o < numberOfOrientations; o++) normalisation += exp(-beta * (energies[o] - least));

  float strongestPull = 0.0f;
  float nearestDistanceSquared = MAXFLOAT;
  int pullingAtom = -1;
  int nearestAtom = -1;

  for(int iatom = 0; iatom < numberOfAtoms; iatom++)
  {
    float pull = 0.0f;
    float closestOverOrientations = MAXFLOAT;

    for(int o = 0; o < numberOfOrientations; o++)
    {
      float closest;
      float contribution =
          contributionOfAtom(position, epsilonTimesFour, sigmaSquared, chargeProduct, shiftValue, siteOffset,
                             centre, o, iatom, cellRow0, cellRow1, cellRow2, widths, numberOfAtoms, numberOfSites,
                             shells, cutOff, cutOffSquared, coulombCutOffSquared, longestCutOffSquared, alpha,
                             useCharges, ceiling, &closest);

      pull += (exp(-beta * (energies[o] - least)) / normalisation) * contribution;
      closestOverOrientations = fmin(closestOverOrientations, closest);
    }

    if(pull < strongestPull)
    {
      strongestPull = pull;
      pullingAtom = iatom;
    }
    if(closestOverOrientations < nearestDistanceSquared)
    {
      nearestDistanceSquared = closestOverOrientations;
      nearestAtom = iatom;
    }
  }

  int index = (iz * gridSize.y + iy) * gridSize.x + ix;
  label[index] = (pullingAtom >= 0) ? pullingAtom : nearestAtom;
}
)foo";
