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

module energy_opencl_well_field;

import std;

import opencl;
import int3;
import uint3;
import double3;
import double3x3;
import unit_cell;
import crystal;
import pair_interactions;
import blocking_spheres;
import energy_shared_linear_probe;
import energy_shared_electrostatic_potential_grid;
import energy_shared_well_field;

namespace
{
constexpr std::size_t maxOrientations = 128;

void check(cl_int err, const char *what, int line)
{
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("WellFieldOpenCL: OpenCL {} failed at {} (err {})\n", what, line, err));
  }
}

cl_program buildProgram()
{
  if (!OpenCL::clContext.has_value() || !OpenCL::clDeviceId.has_value() || !OpenCL::clCommandQueue.has_value())
  {
    throw std::runtime_error("WellFieldOpenCL: no OpenCL device found, the well field needs a GPU\n");
  }

  std::string source = std::string(WellFieldOpenCL::doubleFloatSource) + WellFieldOpenCL::kernelSource;
  const char *ptr = source.c_str();
  cl_int err = CL_SUCCESS;
  cl_program program = clCreateProgramWithSource(OpenCL::clContext.value(), 1, &ptr, nullptr, &err);
  check(err, "clCreateProgramWithSource", __LINE__);

  err = clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    std::size_t length = 0;
    clGetProgramBuildInfo(program, OpenCL::clDeviceId.value(), CL_PROGRAM_BUILD_LOG, 0, nullptr, &length);
    std::string log(length, '\0');
    clGetProgramBuildInfo(program, OpenCL::clDeviceId.value(), CL_PROGRAM_BUILD_LOG, length, log.data(), nullptr);
    clReleaseProgram(program);
    throw std::runtime_error(std::format("WellFieldOpenCL: OpenCL failed to build program (error: {})\n", log));
  }
  return program;
}

struct GpuPack
{
  std::vector<cl_float4> positions;
  std::vector<cl_float> epsilonTimesFour;
  std::vector<cl_float> sigmaSquared;
  std::vector<cl_float> contact;
  std::vector<cl_float> shiftValue;
  std::vector<cl_float> chargeProduct;
  std::vector<cl_float4> siteOffset;
  std::vector<cl_float4> axes;
  std::vector<cl_float> siteCharge;
  std::vector<cl_int> siteDispersion;
  std::vector<cl_float> siteLength;
  std::vector<cl_float> screenedTable;
  std::vector<cl_float4> spheres;
  std::span<const float> smoothPotential{};

  cl_float4 cellRow0{};
  cl_float4 cellRow1{};
  cl_float4 cellRow2{};
  cl_float4 invRow0{};
  cl_float4 invRow1{};
  cl_float4 invRow2{};
  cl_float4 widths{};
  cl_int3 shells{};
  cl_int numberOfAtoms{0};
  cl_int numberOfSites{0};
  cl_int numberOfOrientations{0};
  cl_int numberOfSpheres{0};
  cl_int useCharges{0};
  cl_float cutOffSquared{0.0f};
  cl_float coulombCutOffSquared{0.0f};
  cl_float longestCutOff{0.0f};
  cl_float longestCutOffSquared{0.0f};
  cl_float ceiling{0.0f};
  cl_float beta{0.0f};
  cl_float kT{0.0f};
  cl_float screenedScale{0.0f};
  cl_int screenedBins{static_cast<cl_int>(ScreenedCoulomb::bins)};
  cl_float alpha{0.0f};
  cl_float softminTau{static_cast<cl_float>(wellSoftminTau)};
  cl_float blockedEnergyPerAngstrom{0.0f};
  cl_float extraReach{0.0f};
  bool chargesIncluded{false};
  bool chargesIgnored{false};
  double ewaldAlpha{0.0};
  std::size_t numberOfWaveVectors{0};
  double cutOff{0.0};
  uint3 potentialGridSize{0, 0, 0};
};

GpuPack packNeighbourhood(const PairInteractions &interactions, const Crystal &framework, const LinearProbe &probe,
                          std::size_t numberOfOrientations, double thermalEnergy,
                          std::span<const BlockingSphere> blockingSpheres, double blockedEnergyPerAngstrom,
                          double ceiling, const ElectrostaticPotentialGrid *potential, double coulombFactor,
                          double extraReach)
{
  GpuPack pack;
  pack.ceiling = cl_float(ceiling);
  pack.blockedEnergyPerAngstrom = cl_float(blockedEnergyPerAngstrom);
  pack.extraReach = cl_float(extraReach);
  pack.cutOff = interactions.cutOffVDW;
  pack.cutOffSquared = cl_float(pack.cutOff * pack.cutOff);
  pack.kT = cl_float(std::max(thermalEnergy, 0.0));
  pack.beta = (thermalEnergy > 0.0) ? cl_float(1.0 / thermalEnergy) : 0.0f;

  pack.chargesIncluded = probe.isCharged() && potential != nullptr && potential->numberOfVoxels() > 0;
  pack.chargesIgnored = probe.isCharged() && !pack.chargesIncluded;
  pack.useCharges = pack.chargesIncluded ? 1 : 0;

  double coulombCutOff = 0.0;
  ScreenedCoulomb screened;
  if (pack.chargesIncluded)
  {
    if (coulombFactor == 0.0)
    {
      throw std::runtime_error(
          "Well field: the near half of the electrostatic sum needs the same charge-to-energy conversion "
          "the far half was built with\n");
    }
    coulombCutOff = potential->cutOff;
    pack.coulombCutOffSquared = cl_float(coulombCutOff * coulombCutOff);
    screened.build(potential->alpha, coulombCutOff * coulombCutOff);
    pack.screenedScale = cl_float(screened.scale);
    pack.alpha = cl_float(potential->alpha);
    pack.ewaldAlpha = potential->alpha;
    pack.numberOfWaveVectors = potential->numberOfWaveVectors;
    pack.smoothPotential = std::span<const float>(potential->smoothPotential);
    pack.potentialGridSize = potential->gridSize;
    pack.screenedTable.resize(screened.table.size());
    for (std::size_t i = 0; i < screened.table.size(); ++i) pack.screenedTable[i] = cl_float(screened.table[i]);
  }
  else
  {
    pack.screenedTable = {0.0f, 0.0f};
  }

  double longest = std::max(pack.cutOff, coulombCutOff);
  pack.longestCutOff = cl_float(longest);
  pack.longestCutOffSquared = cl_float(longest * longest);

  const std::size_t numberOfAtoms = framework.atoms.size();
  pack.numberOfAtoms = static_cast<cl_int>(numberOfAtoms);
  pack.positions.resize(std::max<std::size_t>(numberOfAtoms, 1));
  for (std::size_t i = 0; i < numberOfAtoms; ++i)
  {
    double3 fractional = framework.fractionalPositions.empty()
                             ? framework.unitCell.inverseCell * framework.atoms[i].position
                             : framework.fractionalPositions[i];
    pack.positions[i] = {{cl_float(fractional.x), cl_float(fractional.y), cl_float(fractional.z), 0.0f}};
  }

  std::vector<std::size_t> kept;
  std::vector<ProbeSite> sites;
  for (std::size_t s = 0; s < probe.sites.size(); ++s)
  {
    const LinearProbe::Site &site = probe.sites[s];
    double charge = pack.chargesIncluded ? site.charge : 0.0;
    bool carriesDispersion = false;
    for (std::size_t i = 0; i < numberOfAtoms; ++i)
    {
      carriesDispersion =
          carriesDispersion || interactions(site.type, framework.atoms[i].type).strengthParameter != 0.0;
    }
    if (!carriesDispersion && charge == 0.0) continue;
    kept.push_back(s);
    sites.push_back(ProbeSite{site.offset, charge, carriesDispersion});
  }
  if (!std::ranges::any_of(sites, &ProbeSite::dispersion))
  {
    throw std::runtime_error(
        std::format("Well field: probe '{}' has no site with any dispersion against this framework\n", probe.name));
  }

  const std::size_t numberOfSites = sites.size();
  pack.numberOfSites = static_cast<cl_int>(numberOfSites);
  pack.siteCharge.resize(numberOfSites);
  pack.siteDispersion.resize(numberOfSites);
  pack.siteLength.resize(numberOfSites);
  double halfSpan = 0.0;
  for (std::size_t s = 0; s < numberOfSites; ++s)
  {
    pack.siteCharge[s] = cl_float(sites[s].charge);
    pack.siteDispersion[s] = sites[s].dispersion ? 1 : 0;
    pack.siteLength[s] = cl_float(sites[s].offset);
    halfSpan = std::max(halfSpan, std::abs(sites[s].offset));
  }

  pack.epsilonTimesFour.resize(numberOfSites * std::max<std::size_t>(numberOfAtoms, 1));
  pack.sigmaSquared.resize(pack.epsilonTimesFour.size());
  pack.contact.resize(pack.epsilonTimesFour.size());
  pack.shiftValue.resize(pack.epsilonTimesFour.size());
  pack.chargeProduct.resize(pack.epsilonTimesFour.size());
  for (std::size_t s = 0; s < numberOfSites; ++s)
  {
    for (std::size_t i = 0; i < numberOfAtoms; ++i)
    {
      const PairParameters &pair = interactions(probe.sites[kept[s]].type, framework.atoms[i].type);
      const std::size_t index = s * numberOfAtoms + i;
      pack.epsilonTimesFour[index] = cl_float(4.0 * pair.strengthParameter);
      pack.sigmaSquared[index] = cl_float(pair.sizeParameter * pair.sizeParameter);
      pack.contact[index] = cl_float(wellContactPrefactor * pair.sizeParameter);
      pack.shiftValue[index] = cl_float(pair.shift);
      pack.chargeProduct[index] = cl_float(coulombFactor * sites[s].charge * framework.atoms[i].charge);
    }
  }
  const float dummyCore2 = cl_float(wellDummyCoreRadius * wellDummyCoreRadius);
  for (std::size_t i = 0; i < numberOfAtoms; ++i)
  {
    for (std::size_t s = 0; s < numberOfSites; ++s)
    {
      if (!sites[s].dispersion) pack.sigmaSquared[s * numberOfAtoms + i] = dummyCore2;
    }
  }

  std::vector<double3> directions = (halfSpan > 0.0 && numberOfOrientations > 1)
                                        ? orientationSet(numberOfOrientations, probe.headTailSymmetric)
                                        : std::vector<double3>{double3(0.0, 0.0, 1.0)};
  if (directions.size() > maxOrientations)
  {
    throw std::runtime_error(std::format("WellFieldOpenCL: {} orientations exceeds the GPU cap of {}\n",
                                         directions.size(), maxOrientations));
  }
  pack.numberOfOrientations = static_cast<cl_int>(directions.size());
  pack.axes.resize(directions.size());
  pack.siteOffset.resize(directions.size() * numberOfSites);
  const double3x3 inverseCell = framework.unitCell.inverseCell;
  for (std::size_t o = 0; o < directions.size(); ++o)
  {
    pack.axes[o] = {{cl_float(directions[o].x), cl_float(directions[o].y), cl_float(directions[o].z), 0.0f}};
    for (std::size_t s = 0; s < numberOfSites; ++s)
    {
      double3 displacement = inverseCell * (sites[s].offset * directions[o]);
      pack.siteOffset[o * numberOfSites + s] = {
          {cl_float(displacement.x), cl_float(displacement.y), cl_float(displacement.z), 0.0f}};
    }
  }

  pack.spheres.resize(std::max<std::size_t>(blockingSpheres.size(), 1));
  pack.numberOfSpheres = static_cast<cl_int>(blockingSpheres.size());
  for (std::size_t i = 0; i < blockingSpheres.size(); ++i)
  {
    pack.spheres[i] = {{cl_float(blockingSpheres[i].centerFractional.x),
                        cl_float(blockingSpheres[i].centerFractional.y),
                        cl_float(blockingSpheres[i].centerFractional.z), cl_float(blockingSpheres[i].radius)}};
  }

  double3 width = framework.unitCell.perpendicularWidths();
  pack.widths = {{cl_float(width.x), cl_float(width.y), cl_float(width.z), 0.0f}};
  double reach = longest + halfSpan + extraReach;
  pack.shells = {{cl_int(std::floor(reach / width.x + 0.5)), cl_int(std::floor(reach / width.y + 0.5)),
                  cl_int(std::floor(reach / width.z + 0.5)), 0}};

  double3x3 cell = framework.unitCell.cell;
  pack.cellRow0 = {{cl_float(cell[0][0]), cl_float(cell[1][0]), cl_float(cell[2][0]), 0.0f}};
  pack.cellRow1 = {{cl_float(cell[0][1]), cl_float(cell[1][1]), cl_float(cell[2][1]), 0.0f}};
  pack.cellRow2 = {{cl_float(cell[0][2]), cl_float(cell[1][2]), cl_float(cell[2][2]), 0.0f}};
  double3x3 inverse = framework.unitCell.inverseCell;
  pack.invRow0 = {{cl_float(inverse[0][0]), cl_float(inverse[1][0]), cl_float(inverse[2][0]), 0.0f}};
  pack.invRow1 = {{cl_float(inverse[0][1]), cl_float(inverse[1][1]), cl_float(inverse[2][1]), 0.0f}};
  pack.invRow2 = {{cl_float(inverse[0][2]), cl_float(inverse[1][2]), cl_float(inverse[2][2]), 0.0f}};
  return pack;
}

struct DeviceBuffers
{
  cl_mem position{};
  cl_mem epsilon{};
  cl_mem sigma{};
  cl_mem contact{};
  cl_mem shift{};
  cl_mem charge{};
  cl_mem siteOffset{};
  cl_mem axes{};
  cl_mem siteCharge{};
  cl_mem siteDispersion{};
  cl_mem siteLength{};
  cl_mem screened{};
  cl_mem potential{};
  cl_mem spheres{};

  void release() const
  {
    if (position) clReleaseMemObject(position);
    if (epsilon) clReleaseMemObject(epsilon);
    if (sigma) clReleaseMemObject(sigma);
    if (contact) clReleaseMemObject(contact);
    if (shift) clReleaseMemObject(shift);
    if (charge) clReleaseMemObject(charge);
    if (siteOffset) clReleaseMemObject(siteOffset);
    if (axes) clReleaseMemObject(axes);
    if (siteCharge) clReleaseMemObject(siteCharge);
    if (siteDispersion) clReleaseMemObject(siteDispersion);
    if (siteLength) clReleaseMemObject(siteLength);
    if (screened) clReleaseMemObject(screened);
    if (potential) clReleaseMemObject(potential);
    if (spheres) clReleaseMemObject(spheres);
  }
};

cl_mem bufferOf(cl_mem_flags flags, std::size_t bytes, const void *host)
{
  cl_mem mem = OpenCL::createBuffer(flags, bytes);
  if (host != nullptr && bytes > 0)
  {
    OpenCL::writeBuffer(mem, bytes, host);
  }
  return mem;
}

DeviceBuffers uploadPack(const GpuPack &pack, std::size_t potentialVoxels)
{
  DeviceBuffers buffers;
  buffers.position = bufferOf(CL_MEM_READ_ONLY, sizeof(cl_float4) * pack.positions.size(), pack.positions.data());
  buffers.epsilon =
      bufferOf(CL_MEM_READ_ONLY, sizeof(cl_float) * pack.epsilonTimesFour.size(), pack.epsilonTimesFour.data());
  buffers.sigma = bufferOf(CL_MEM_READ_ONLY, sizeof(cl_float) * pack.sigmaSquared.size(), pack.sigmaSquared.data());
  buffers.contact = bufferOf(CL_MEM_READ_ONLY, sizeof(cl_float) * pack.contact.size(), pack.contact.data());
  buffers.shift = bufferOf(CL_MEM_READ_ONLY, sizeof(cl_float) * pack.shiftValue.size(), pack.shiftValue.data());
  buffers.charge = bufferOf(CL_MEM_READ_ONLY, sizeof(cl_float) * pack.chargeProduct.size(), pack.chargeProduct.data());
  buffers.siteOffset = bufferOf(CL_MEM_READ_ONLY, sizeof(cl_float4) * pack.siteOffset.size(), pack.siteOffset.data());
  buffers.axes = bufferOf(CL_MEM_READ_ONLY, sizeof(cl_float4) * pack.axes.size(), pack.axes.data());
  buffers.siteCharge = bufferOf(CL_MEM_READ_ONLY, sizeof(cl_float) * pack.siteCharge.size(), pack.siteCharge.data());
  buffers.siteDispersion =
      bufferOf(CL_MEM_READ_ONLY, sizeof(cl_int) * pack.siteDispersion.size(), pack.siteDispersion.data());
  buffers.siteLength = bufferOf(CL_MEM_READ_ONLY, sizeof(cl_float) * pack.siteLength.size(), pack.siteLength.data());
  buffers.screened = bufferOf(CL_MEM_READ_ONLY, sizeof(cl_float) * pack.screenedTable.size(), pack.screenedTable.data());
  buffers.spheres = bufferOf(CL_MEM_READ_ONLY, sizeof(cl_float4) * pack.spheres.size(), pack.spheres.data());
  if (pack.chargesIncluded && potentialVoxels > 0)
  {
    buffers.potential =
        bufferOf(CL_MEM_READ_ONLY, sizeof(cl_float) * potentialVoxels, pack.smoothPotential.data());
  }
  else
  {
    float zero = 0.0f;
    buffers.potential = bufferOf(CL_MEM_READ_ONLY, sizeof(cl_float), &zero);
  }
  return buffers;
}

int setCommonArgs(cl_kernel kernel, int start, const DeviceBuffers &buffers, const GpuPack &pack)
{
  cl_int err = CL_SUCCESS;
  int a = start;
  err |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &buffers.position);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &buffers.epsilon);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &buffers.sigma);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &buffers.contact);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &buffers.shift);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &buffers.charge);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &buffers.siteOffset);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &buffers.axes);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &buffers.siteCharge);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &buffers.siteDispersion);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &buffers.siteLength);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &buffers.screened);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &buffers.potential);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &buffers.spheres);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_float4), &pack.cellRow0);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_float4), &pack.cellRow1);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_float4), &pack.cellRow2);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_float4), &pack.invRow0);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_float4), &pack.invRow1);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_float4), &pack.invRow2);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_float4), &pack.widths);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_int), &pack.numberOfAtoms);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_int), &pack.numberOfSites);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_int), &pack.numberOfOrientations);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_int3), &pack.shells);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_int), &pack.numberOfSpheres);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_float), &pack.cutOffSquared);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_float), &pack.coulombCutOffSquared);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_float), &pack.longestCutOff);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_float), &pack.longestCutOffSquared);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_float), &pack.ceiling);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_float), &pack.beta);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_float), &pack.kT);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_float), &pack.screenedScale);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_int), &pack.screenedBins);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_float), &pack.alpha);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_int), &pack.useCharges);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_float), &pack.softminTau);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_float), &pack.blockedEnergyPerAngstrom);
  err |= clSetKernelArg(kernel, a++, sizeof(cl_float), &pack.extraReach);
  check(err, "clSetKernelArg", __LINE__);
  return a;
}
}  // namespace


WellField WellFieldOpenCL::compute(const PairInteractions &interactions, const Crystal &framework,
                                   const LinearProbe &probe, uint3 gridSize, std::size_t numberOfOrientations,
                                   double thermalEnergy, std::span<const BlockingSphere> blockingSpheres,
                                   double blockedEnergyPerAngstrom, double ceiling,
                                   const ElectrostaticPotentialGrid *potential, double coulombFactor)
{
  if (gridSize.x == 0 || gridSize.y == 0 || gridSize.z == 0)
  {
    throw std::runtime_error("Well field: the grid must have at least one point along each axis\n");
  }
  if (potential != nullptr && potential->numberOfVoxels() > 0 &&
      (potential->gridSize.x != gridSize.x || potential->gridSize.y != gridSize.y ||
       potential->gridSize.z != gridSize.z))
  {
    throw std::runtime_error("Well field: the electrostatic potential is on a different grid than this one\n");
  }

  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  GpuPack pack = packNeighbourhood(interactions, framework, probe, numberOfOrientations, thermalEnergy,
                                   blockingSpheres, blockedEnergyPerAngstrom, ceiling, potential, coulombFactor, 0.0);

  WellField field;
  field.gridSize = gridSize;
  field.unitCell = framework.unitCell;
  field.probe = probe;
  field.probeName = probe.name;
  field.numberOfOrientations = static_cast<std::size_t>(pack.numberOfOrientations);
  field.thermalEnergy = thermalEnergy;
  field.cutOff = pack.cutOff;
  field.ceiling = ceiling;
  field.chargesIncluded = pack.chargesIncluded;
  field.chargesIgnored = pack.chargesIgnored;
  field.ewaldAlpha = pack.ewaldAlpha;
  field.numberOfWaveVectors = pack.numberOfWaveVectors;

  const std::size_t numberOfVoxels = static_cast<std::size_t>(gridSize.x) * static_cast<std::size_t>(gridSize.y) *
                                    static_cast<std::size_t>(gridSize.z);
  field.energy.assign(numberOfVoxels, 0.0f);
  field.distance.assign(numberOfVoxels, 1.0e10f);
  field.reliability.assign(numberOfVoxels, 1.0f);
  if (framework.atoms.empty())
  {
    field.seconds = 0.0;
    return field;
  }

  cl_program program = buildProgram();
  cl_int err = CL_SUCCESS;
  cl_kernel kernel = clCreateKernel(program, "WellField", &err);
  check(err, "clCreateKernel WellField", __LINE__);

  DeviceBuffers buffers = uploadPack(pack, numberOfVoxels);
  cl_mem energyBuffer = bufferOf(CL_MEM_WRITE_ONLY, sizeof(cl_float) * numberOfVoxels, nullptr);
  cl_mem distanceBuffer = bufferOf(CL_MEM_WRITE_ONLY, sizeof(cl_float) * numberOfVoxels, nullptr);
  cl_mem reliabilityBuffer = bufferOf(CL_MEM_WRITE_ONLY, sizeof(cl_float) * numberOfVoxels, nullptr);

  cl_int3 clGridSize = {{cl_int(gridSize.x), cl_int(gridSize.y), cl_int(gridSize.z), 0}};
  int arg = 0;
  err = clSetKernelArg(kernel, arg++, sizeof(cl_mem), &energyBuffer);
  err |= clSetKernelArg(kernel, arg++, sizeof(cl_mem), &distanceBuffer);
  err |= clSetKernelArg(kernel, arg++, sizeof(cl_mem), &reliabilityBuffer);
  err |= clSetKernelArg(kernel, arg++, sizeof(cl_int3), &clGridSize);
  check(err, "clSetKernelArg outputs", __LINE__);
  setCommonArgs(kernel, arg, buffers, pack);

  std::size_t globalWorkSize[3] = {static_cast<std::size_t>(gridSize.x), static_cast<std::size_t>(gridSize.y),
                                   static_cast<std::size_t>(gridSize.z)};
  err = clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), kernel, 3, nullptr, globalWorkSize, nullptr, 0, nullptr,
                               nullptr);
  check(err, "clEnqueueNDRangeKernel WellField", __LINE__);

  OpenCL::readBuffer(energyBuffer, sizeof(cl_float) * numberOfVoxels, field.energy.data());
  OpenCL::readBuffer(distanceBuffer, sizeof(cl_float) * numberOfVoxels, field.distance.data());
  OpenCL::readBuffer(reliabilityBuffer, sizeof(cl_float) * numberOfVoxels, field.reliability.data());
  clFinish(OpenCL::clCommandQueue.value());

  clReleaseMemObject(energyBuffer);
  clReleaseMemObject(distanceBuffer);
  clReleaseMemObject(reliabilityBuffer);
  buffers.release();
  clReleaseKernel(kernel);
  clReleaseProgram(program);

  field.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - time_begin).count();
  return field;
}


void WellFieldOpenCL::refineVertices(std::vector<double3> &corners, std::vector<double> &energies,
                                    const PairInteractions &interactions, const Crystal &framework,
                                    const NeighbourhoodParameters &parameters,
                                    std::span<const BlockingSphere> blockingSpheres, double iso)
{
  energies.assign(corners.size(), 0.0);
  if (corners.empty()) return;

  GpuPack pack =
      packNeighbourhood(interactions, framework, parameters.probe, parameters.numberOfOrientations,
                        parameters.thermalEnergy, blockingSpheres, parameters.blockedEnergyPerAngstrom,
                        parameters.ceiling, parameters.potential, parameters.coulombFactor,
                        parameters.extraReach > 0.0 ? parameters.extraReach : wellRefinementReach);

  cl_program program = buildProgram();
  cl_int err = CL_SUCCESS;
  cl_kernel kernel = clCreateKernel(program, "RefineWellVertices", &err);
  check(err, "clCreateKernel RefineWellVertices", __LINE__);

  std::vector<cl_float4> cornerData(corners.size());
  std::vector<cl_float> energyData(corners.size(), 0.0f);
  for (std::size_t i = 0; i < corners.size(); ++i)
  {
    cornerData[i] = {{cl_float(corners[i].x), cl_float(corners[i].y), cl_float(corners[i].z), 0.0f}};
  }

  DeviceBuffers buffers = uploadPack(pack, pack.chargesIncluded ? pack.smoothPotential.size() : 0);
  cl_mem cornerBuffer =
      bufferOf(CL_MEM_READ_WRITE, sizeof(cl_float4) * cornerData.size(), cornerData.data());
  cl_mem energyBuffer = bufferOf(CL_MEM_WRITE_ONLY, sizeof(cl_float) * energyData.size(), nullptr);

  cl_int numberOfVertices = static_cast<cl_int>(corners.size());
  cl_float clIso = cl_float(iso);
  cl_int3 clPotentialSize = {{cl_int(pack.potentialGridSize.x), cl_int(pack.potentialGridSize.y),
                              cl_int(pack.potentialGridSize.z), 0}};

  int arg = 0;
  err = clSetKernelArg(kernel, arg++, sizeof(cl_mem), &cornerBuffer);
  err |= clSetKernelArg(kernel, arg++, sizeof(cl_mem), &energyBuffer);
  err |= clSetKernelArg(kernel, arg++, sizeof(cl_int), &numberOfVertices);
  err |= clSetKernelArg(kernel, arg++, sizeof(cl_float), &clIso);
  err |= clSetKernelArg(kernel, arg++, sizeof(cl_int3), &clPotentialSize);
  check(err, "clSetKernelArg refine outputs", __LINE__);
  arg = setCommonArgs(kernel, arg, buffers, pack);
  std::size_t globalWork = corners.size();
  err = clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), kernel, 1, nullptr, &globalWork, nullptr, 0, nullptr,
                               nullptr);
  check(err, "clEnqueueNDRangeKernel RefineWellVertices", __LINE__);

  err = clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), cornerBuffer, CL_TRUE, 0,
                            sizeof(cl_float4) * cornerData.size(), cornerData.data(), 0, nullptr, nullptr);
  err |= clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), energyBuffer, CL_TRUE, 0,
                             sizeof(cl_float) * energyData.size(), energyData.data(), 0, nullptr, nullptr);
  check(err, "clEnqueueReadBuffer refine", __LINE__);
  clFinish(OpenCL::clCommandQueue.value());

  for (std::size_t i = 0; i < corners.size(); ++i)
  {
    corners[i] = double3(static_cast<double>(cornerData[i].s[0]), static_cast<double>(cornerData[i].s[1]),
                         static_cast<double>(cornerData[i].s[2]));
    energies[i] = static_cast<double>(energyData[i]);
  }

  clReleaseMemObject(cornerBuffer);
  clReleaseMemObject(energyBuffer);
  buffers.release();
  clReleaseKernel(kernel);
  clReleaseProgram(program);
}


const char *WellFieldOpenCL::doubleFloatSource = R"foo(
// Unevaluated hi/lo pair, about 48 mantissa bits. Only the accumulator is compensated; each pair term
// stays a single float.
typedef float2 df;

static df twoSum(float a, float b)
{
  float s = a + b;
  float v = s - a;
  float e = (a - (s - v)) + (b - v);
  return (df)(s, e);
}

static df quickTwoSum(float a, float b)
{
  float s = a + b;
  float e = b - (s - a);
  return (df)(s, e);
}

static df df_add(df a, df b)
{
  df s = twoSum(a.x, b.x);
  df t = twoSum(a.y, b.y);
  s.y += t.x;
  s = quickTwoSum(s.x, s.y);
  s.y += t.y;
  return quickTwoSum(s.x, s.y);
}

static df df_add_f(df a, float b)
{
  df s = twoSum(a.x, b);
  s.y += a.y;
  return quickTwoSum(s.x, s.y);
}

static float df_hi(df a) { return a.x + a.y; }

// Ordered comparison of the unevaluated pair. A float compare of df_hi is only
// good to a ulp of the energy (~2e-4 K at a few thousand K); this one keeps the
// extra bits so a twenty-step golden section can still tell the two interior
// samples apart.
static int df_lt(df a, df b)
{
  return (a.x < b.x) || (a.x == b.x && a.y < b.y);
}
)foo";


const char *WellFieldOpenCL::kernelSource = R"foo(
#define MAX_ORIENT 128
#define NEIGHBOUR_CAP 1024
#define COARSE_STEPS 14
#define GOLDEN_STEPS 20
#define GOLDEN_INVPHI 0.6180339887f

static float smoothPotentialAt(__global const float *smoothPotential, float3 s, int3 gridSize)
{
  if(gridSize.x <= 0 || gridSize.y <= 0 || gridSize.z <= 0) return 0.0f;
  s -= floor(s);
  float3 scaled = s * (float3)((float)(gridSize.x), (float)(gridSize.y), (float)(gridSize.z));
  int3 low = convert_int3(floor(scaled));
  float3 frac = scaled - convert_float3(low);
  low.x = ((low.x % gridSize.x) + gridSize.x) % gridSize.x;
  low.y = ((low.y % gridSize.y) + gridSize.y) % gridSize.y;
  low.z = ((low.z % gridSize.z) + gridSize.z) % gridSize.z;
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

static float screenedAt(float rr, __global const float *table, float scale, int bins, float alpha)
{
  if(rr < 0.25f)
  {
    float r = sqrt(fmax(rr, 1.0e-8f));
    return erfc(alpha * r) / r;
  }
  float x = (rr - 0.25f) * scale;
  int i = (int)(x);
  int last = bins + 1;
  if(i >= last) return table[last];
  if(i < 0) i = 0;
  float f = x - (float)(i);
  return table[i] + f * (table[i + 1] - table[i]);
}

static float4 cartesianOf(float4 t, float4 cellRow0, float4 cellRow1, float4 cellRow2)
{
  float4 dr;
  dr.x = dot(cellRow0, t);
  dr.y = dot(cellRow1, t);
  dr.z = dot(cellRow2, t);
  dr.w = 0.0f;
  return dr;
}

static float blockingPocket(float4 centre, __global const float4 *spheres, int numberOfSpheres,
                            float4 cellRow0, float4 cellRow1, float4 cellRow2)
{
  float nearest = 1.0e10f;
  for(int i = 0; i < numberOfSpheres; i++)
  {
    float4 ds = centre - spheres[i];
    ds -= rint(ds);
    ds.w = 0.0f;
    float4 dr = cartesianOf(ds, cellRow0, cellRow1, cellRow2);
    nearest = fmin(nearest, sqrt(dot(dr, dr)) - spheres[i].w);
  }
  return nearest;
}

static void addPair(float rr, float epsilon4, float sigma2, float contact, float shift, float charges,
                    int dispersion, int wantClearance, int withLJ, int useCharges, float cutOffSquared,
                    float coulombCutOffSquared, float ceiling, __global const float *screenedTable,
                    float screenedScale, int screenedBins, float alpha, df *energy, float *nearest,
                    int *skipSmooth, int *overlapHit)
{
  if(withLJ != 0 && dispersion != 0 && rr < cutOffSquared)
  {
    float clamped = fmax(rr, 1.0e-8f);
    float ratio = sigma2 / clamped;
    float ratio3 = ratio * ratio * ratio;
    *energy = df_add_f(*energy, fmin(epsilon4 * ratio3 * (ratio3 - 1.0f) - shift, ceiling));
    if(overlapHit != 0 && rr < sigma2) *overlapHit = 1;
    if(wantClearance != 0) *nearest = fmin(*nearest, sqrt(clamped) - contact);
  }
  int inCore = (dispersion == 0 && sigma2 > 0.0f && rr < sigma2) ? 1 : 0;
  if(inCore != 0 && skipSmooth != 0) *skipSmooth = 1;
  if(useCharges != 0 && inCore == 0 && rr < coulombCutOffSquared && charges != 0.0f)
  {
    *energy = df_add_f(*energy, fmin(charges * screenedAt(rr, screenedTable, screenedScale, screenedBins, alpha),
                                     ceiling));
  }
}

static void visitPairs(float4 centre,
                       __global const float4 *position,
                       __global const float *epsilonTimesFour,
                       __global const float *sigmaSquared,
                       __global const float *contact,
                       __global const float *shiftValue,
                       __global const float *chargeProduct,
                       __global const float4 *siteOffset,
                       __global const float *siteCharge,
                       __global const int *siteDispersion,
                       __global const float *screenedTable,
                       __global const float *smoothPotential,
                       int3 potentialGridSize,
                       int numberOfAtoms,
                       int numberOfSites,
                       int3 shells,
                       float4 widths,
                       float4 cellRow0,
                       float4 cellRow1,
                       float4 cellRow2,
                       float cutOffSquared,
                       float coulombCutOffSquared,
                       float longestCutOff,
                       float longestCutOffSquared,
                       float extraReach,
                       float ceiling,
                       float screenedScale,
                       int screenedBins,
                       float alpha,
                       int o,
                       int withLJ,
                       int withCharges,
                       int wantClearance,
                       df *energy,
                       float *nearest,
                       int *dummyHit,
                       int *overlapHit)
{
  float reach = longestCutOff + extraReach;
  for(int site = 0; site < numberOfSites; site++)
  {
    float4 s = centre + siteOffset[o * numberOfSites + site];
    s.w = 0.0f;
    int skipSmooth = 0;
    __global const float *epsilonForSite = epsilonTimesFour + site * numberOfAtoms;
    __global const float *sigmaForSite = sigmaSquared + site * numberOfAtoms;
    __global const float *contactForSite = contact + site * numberOfAtoms;
    __global const float *shiftForSite = shiftValue + site * numberOfAtoms;
    __global const float *chargeForSite = chargeProduct + site * numberOfAtoms;
    int dispersion = siteDispersion[site];
    for(int iatom = 0; iatom < numberOfAtoms; iatom++)
    {
      float4 ds = s - position[iatom];
      ds -= rint(ds);
      ds.w = 0.0f;
      for(int a = -shells.x; a <= shells.x; a++)
      {
        for(int b = -shells.y; b <= shells.y; b++)
        {
          for(int c = -shells.z; c <= shells.z; c++)
          {
            float4 t = ds + (float4)((float)(a), (float)(b), (float)(c), 0.0f);
            float far = fmax(fmax(fabs(t.x) * widths.x, fabs(t.y) * widths.y), fabs(t.z) * widths.z);
            if(far > reach) continue;
            float4 dr = cartesianOf(t, cellRow0, cellRow1, cellRow2);
            float rr = dot(dr, dr);
            if(rr >= longestCutOffSquared && !(dispersion != 0)) continue;
            addPair(rr, epsilonForSite[iatom], sigmaForSite[iatom], contactForSite[iatom], shiftForSite[iatom],
                    chargeForSite[iatom], dispersion, wantClearance, withLJ, withCharges, cutOffSquared,
                    coulombCutOffSquared, ceiling, screenedTable, screenedScale, screenedBins, alpha, energy,
                    nearest, &skipSmooth, overlapHit);
          }
        }
      }
    }
    if(dispersion == 0 && skipSmooth != 0) *dummyHit = 1;
    if(withCharges != 0 && siteCharge[site] != 0.0f && skipSmooth == 0)
    {
      *energy = df_add_f(*energy, siteCharge[site] * smoothPotentialAt(smoothPotential, s.xyz, potentialGridSize));
    }
  }
}

static df reduceOrientations(df *energies, float *clearances, int *dummyHit, int numberOfOrientations, float beta,
                             float kT, float ceiling, int *bestOrientation, float *bestClearance)
{
  df leastDf = (df)(MAXFLOAT, 0.0f);
  float least = MAXFLOAT;
  float sum = 0.0f;
  float best = -1.0e10f;
  int bestO = 0;
  for(int o = 0; o < numberOfOrientations; o++)
  {
    df energyDf = energies[o];
    if(dummyHit != 0 && dummyHit[o] != 0) energyDf = (df)(ceiling, 0.0f);
    else if(df_hi(energyDf) > ceiling) energyDf = (df)(ceiling, 0.0f);
    float energy = df_hi(energyDf);
    if(o == 0 || clearances[o] > best)
    {
      best = clearances[o];
      bestO = o;
    }
    if(beta > 0.0f)
    {
      if(energy < least)
      {
        sum = sum * exp(-beta * (least - energy)) + 1.0f;
        least = energy;
        leastDf = energyDf;
      }
      else
      {
        sum += exp(-beta * (energy - least));
      }
    }
    else if(df_lt(energyDf, leastDf))
    {
      leastDf = energyDf;
      least = energy;
    }
  }
  *bestOrientation = bestO;
  *bestClearance = best;
  if(beta > 0.0f)
  {
    return df_add_f(leastDf, -kT * log(sum / (float)(numberOfOrientations)));
  }
  return leastDf;
}

static float3 wrapUnit(float3 s)
{
  return s - floor(s);
}

static float3 quantizeFractional(float3 s)
{
  float3 r = rint(s * 1048576.0f);
  if(r.x == 1048576.0f) r.x = 0.0f;
  if(r.y == 1048576.0f) r.y = 0.0f;
  if(r.z == 1048576.0f) r.z = 0.0f;
  return r / 1048576.0f;
}

static float smoothstepf(float edge0, float edge1, float x)
{
  if(edge1 <= edge0) return x >= edge1 ? 1.0f : 0.0f;
  float t = clamp((x - edge0) / (edge1 - edge0), 0.0f, 1.0f);
  return t * t * (3.0f - 2.0f * t);
}

__kernel void WellField(__global float *energyOut,
                        __global float *distanceOut,
                        __global float *reliabilityOut,
                        const int3 gridSize,
                        __global const float4 *position,
                        __global const float *epsilonTimesFour,
                        __global const float *sigmaSquared,
                        __global const float *contact,
                        __global const float *shiftValue,
                        __global const float *chargeProduct,
                        __global const float4 *siteOffset,
                        __global const float4 *axes,
                        __global const float *siteCharge,
                        __global const int *siteDispersion,
                        __global const float *siteLength,
                        __global const float *screenedTable,
                        __global const float *smoothPotential,
                        __global const float4 *spheres,
                        const float4 cellRow0,
                        const float4 cellRow1,
                        const float4 cellRow2,
                        const float4 invRow0,
                        const float4 invRow1,
                        const float4 invRow2,
                        const float4 widths,
                        const int numberOfAtoms,
                        const int numberOfSites,
                        const int numberOfOrientations,
                        const int3 shells,
                        const int numberOfSpheres,
                        const float cutOffSquared,
                        const float coulombCutOffSquared,
                        const float longestCutOff,
                        const float longestCutOffSquared,
                        const float ceiling,
                        const float beta,
                        const float kT,
                        const float screenedScale,
                        const int screenedBins,
                        const float alpha,
                        const int useCharges,
                        const float softminTau,
                        const float blockedEnergyPerAngstrom,
                        const float extraReach)
{
  int ix = get_global_id(0);
  int iy = get_global_id(1);
  int iz = get_global_id(2);
  if(ix >= gridSize.x || iy >= gridSize.y || iz >= gridSize.z) return;

  float4 centre = (float4)((float)(ix) / (float)(gridSize.x),
                           (float)(iy) / (float)(gridSize.y),
                           (float)(iz) / (float)(gridSize.z),
                           0.0f);

  df energies[MAX_ORIENT];
  float clearances[MAX_ORIENT];
  int dummyHit[MAX_ORIENT];
  int overlapHit[MAX_ORIENT];
  for(int o = 0; o < numberOfOrientations; o++)
  {
    energies[o] = (df)(0.0f, 0.0f);
    clearances[o] = 1.0e10f;
    dummyHit[o] = 0;
    overlapHit[o] = 0;
  }

  for(int o = 0; o < numberOfOrientations; o++)
  {
    visitPairs(centre, position, epsilonTimesFour, sigmaSquared, contact, shiftValue, chargeProduct, siteOffset,
               siteCharge, siteDispersion, screenedTable, smoothPotential, gridSize, numberOfAtoms, numberOfSites,
               shells, widths, cellRow0, cellRow1, cellRow2, cutOffSquared, coulombCutOffSquared, longestCutOff,
               longestCutOffSquared, extraReach, ceiling, screenedScale, screenedBins, alpha, o, 1, 0, 1,
               &energies[o], &clearances[o], &dummyHit[o], &overlapHit[o]);
  }
  if(useCharges != 0)
  {
    for(int o = 0; o < numberOfOrientations; o++)
    {
      if(dummyHit[o] != 0 || overlapHit[o] != 0) continue;
      visitPairs(centre, position, epsilonTimesFour, sigmaSquared, contact, shiftValue, chargeProduct, siteOffset,
                 siteCharge, siteDispersion, screenedTable, smoothPotential, gridSize, numberOfAtoms, numberOfSites,
                 shells, widths, cellRow0, cellRow1, cellRow2, cutOffSquared, coulombCutOffSquared, longestCutOff,
                 longestCutOffSquared, extraReach, ceiling, screenedScale, screenedBins, alpha, o, 0, 1, 0,
                 &energies[o], &clearances[o], &dummyHit[o], 0);
    }
  }

  int bestOrientation = 0;
  float bestClearance = 1.0e10f;
  float value = df_hi(reduceOrientations(energies, clearances, dummyHit, numberOfOrientations, beta, kT, ceiling,
                                         &bestOrientation, &bestClearance));

  float reliability = 1.0f;
  if(bestClearance < 1.0e9f)
  {
    float4 axis = axes[bestOrientation];
    float4 directionSum = (float4)(0.0f);
    float weightSum = 0.0f;
    for(int site = 0; site < numberOfSites; site++)
    {
      if(siteDispersion[site] == 0) continue;
      float4 s = centre + siteOffset[bestOrientation * numberOfSites + site];
      s.w = 0.0f;
      __global const float *contactForSite = contact + site * numberOfAtoms;
      for(int iatom = 0; iatom < numberOfAtoms; iatom++)
      {
        float4 ds = s - position[iatom];
        ds -= rint(ds);
        ds.w = 0.0f;
        for(int a = -shells.x; a <= shells.x; a++)
        {
          for(int b = -shells.y; b <= shells.y; b++)
          {
            for(int c = -shells.z; c <= shells.z; c++)
            {
              float4 t = ds + (float4)((float)(a), (float)(b), (float)(c), 0.0f);
              float4 dr = cartesianOf(t, cellRow0, cellRow1, cellRow2);
              float r = sqrt(fmax(dot(dr, dr), 1.0e-16f));
              if(!(r > 1.0e-6f)) continue;
              float weighted = r - contactForSite[iatom];
              if(weighted - bestClearance > 6.0f * softminTau) continue;
              float w = exp(-(weighted - bestClearance) / softminTau);
              directionSum += (-dr / r) * w;
              weightSum += w;
            }
          }
        }
      }
      (void)axis;
    }
    if(weightSum > 0.0f)
    {
      reliability = sqrt(dot(directionSum, directionSum)) / weightSum;
    }
    else
    {
      reliability = 0.0f;
    }
  }

  float pocket = blockingPocket(centre, spheres, numberOfSpheres, cellRow0, cellRow1, cellRow2);
  int index = (iz * gridSize.y + iy) * gridSize.x + ix;
  float energy = pocket < 0.0f ? fmin(-pocket * blockedEnergyPerAngstrom, ceiling) : value;
  energyOut[index] = energy;
  distanceOut[index] = fmin(bestClearance, pocket);
  reliabilityOut[index] = reliability;
  (void)invRow0; (void)invRow1; (void)invRow2;
}

static df energyAtShift(__global const float4 *position,
                           __global const float *epsilonTimesFour,
                           __global const float *sigmaSquared,
                           __global const float *contact,
                           __global const float *shiftValue,
                           __global const float *chargeProduct,
                           __global const float4 *siteOffset,
                           __global const float *siteCharge,
                           __global const int *siteDispersion,
                           __global const float *screenedTable,
                           __global const float *smoothPotential,
                           float4 centre,
                           float4 shiftCartesian,
                           float4 shiftFractional,
                           int3 potentialGridSize,
                           int numberOfAtoms,
                           int numberOfSites,
                           int numberOfOrientations,
                           int3 shells,
                           float4 widths,
                           float4 cellRow0,
                           float4 cellRow1,
                           float4 cellRow2,
                           float cutOffSquared,
                           float coulombCutOffSquared,
                           float longestCutOff,
                           float longestCutOffSquared,
                           float ceiling,
                           float beta,
                           float kT,
                           float screenedScale,
                           int screenedBins,
                           float alpha,
                           int useCharges,
                           float extraReach)
{
  df energies[MAX_ORIENT];
  float clearances[MAX_ORIENT];
  int dummyHit[MAX_ORIENT];
  int overlapHit[MAX_ORIENT];
  for(int o = 0; o < numberOfOrientations; o++)
  {
    energies[o] = (df)(0.0f, 0.0f);
    clearances[o] = 1.0e10f;
    dummyHit[o] = 0;
    overlapHit[o] = 0;
  }

  float4 query = centre + shiftFractional;
  query.w = 0.0f;

  for(int o = 0; o < numberOfOrientations; o++)
  {
    visitPairs(query, position, epsilonTimesFour, sigmaSquared, contact, shiftValue, chargeProduct, siteOffset,
               siteCharge, siteDispersion, screenedTable, smoothPotential, potentialGridSize, numberOfAtoms,
               numberOfSites, shells, widths, cellRow0, cellRow1, cellRow2, cutOffSquared, coulombCutOffSquared,
               longestCutOff, longestCutOffSquared, extraReach, ceiling, screenedScale, screenedBins, alpha, o, 1, 0,
               0, &energies[o], &clearances[o], &dummyHit[o], &overlapHit[o]);
  }
  if(useCharges != 0)
  {
    for(int o = 0; o < numberOfOrientations; o++)
    {
      if(dummyHit[o] != 0 || overlapHit[o] != 0) continue;
      visitPairs(query, position, epsilonTimesFour, sigmaSquared, contact, shiftValue, chargeProduct, siteOffset,
                 siteCharge, siteDispersion, screenedTable, smoothPotential, potentialGridSize, numberOfAtoms,
                 numberOfSites, shells, widths, cellRow0, cellRow1, cellRow2, cutOffSquared, coulombCutOffSquared,
                 longestCutOff, longestCutOffSquared, extraReach, ceiling, screenedScale, screenedBins, alpha, o, 0,
                 1, 0, &energies[o], &clearances[o], &dummyHit[o], 0);
    }
  }

  int bestOrientation = 0;
  float bestClearance = 0.0f;
  return reduceOrientations(energies, clearances, dummyHit, numberOfOrientations, beta, kT, ceiling, &bestOrientation,
                            &bestClearance);
}

__kernel void RefineWellVertices(__global float4 *corners,
                                 __global float *energies,
                                 const int numberOfVertices,
                                 const float iso,
                                 const int3 potentialGridSize,
                                 __global const float4 *position,
                                 __global const float *epsilonTimesFour,
                                 __global const float *sigmaSquared,
                                 __global const float *contact,
                                 __global const float *shiftValue,
                                 __global const float *chargeProduct,
                                 __global const float4 *siteOffset,
                                 __global const float4 *axes,
                                 __global const float *siteCharge,
                                 __global const int *siteDispersion,
                                 __global const float *siteLength,
                                 __global const float *screenedTable,
                                 __global const float *smoothPotential,
                                 __global const float4 *spheres,
                                 const float4 cellRow0,
                                 const float4 cellRow1,
                                 const float4 cellRow2,
                                 const float4 invRow0,
                                 const float4 invRow1,
                                 const float4 invRow2,
                                 const float4 widths,
                                 const int numberOfAtoms,
                                 const int numberOfSites,
                                 const int numberOfOrientations,
                                 const int3 shells,
                                 const int numberOfSpheres,
                                 const float cutOffSquared,
                                 const float coulombCutOffSquared,
                                 const float longestCutOff,
                                 const float longestCutOffSquared,
                                 const float ceiling,
                                 const float beta,
                                 const float kT,
                                 const float screenedScale,
                                 const int screenedBins,
                                 const float alpha,
                                 const int useCharges,
                                 const float softminTau,
                                 const float blockedEnergyPerAngstrom,
                                 const float extraReach)
{
  int vertex = (int)get_global_id(0);
  if(vertex >= numberOfVertices) return;

  float4 raw = corners[vertex];
  float3 q = quantizeFractional(wrapUnit(raw.xyz));
  float4 point = (float4)(q.x, q.y, q.z, 0.0f);

  df energyHereDf = energyAtShift(position, epsilonTimesFour, sigmaSquared, contact, shiftValue, chargeProduct,
                                  siteOffset, siteCharge, siteDispersion, screenedTable, smoothPotential, point,
                                  (float4)(0.0f), (float4)(0.0f), potentialGridSize, numberOfAtoms, numberOfSites,
                                  numberOfOrientations, shells, widths, cellRow0, cellRow1, cellRow2, cutOffSquared,
                                  coulombCutOffSquared, longestCutOff, longestCutOffSquared, ceiling, beta, kT,
                                  screenedScale, screenedBins, alpha, useCharges, extraReach);
  float energyHere = df_hi(energyHereDf);
  energies[vertex] = energyHere;
  if(energyHere > iso - 0.02f * fabs(iso)) return;

  float pocket = blockingPocket(point, spheres, numberOfSpheres, cellRow0, cellRow1, cellRow2);

  df energiesO[MAX_ORIENT];
  float clearancesO[MAX_ORIENT];
  int dummyHitO[MAX_ORIENT];
  int overlapHitO[MAX_ORIENT];
  for(int o = 0; o < numberOfOrientations; o++)
  {
    energiesO[o] = (df)(0.0f, 0.0f);
    clearancesO[o] = 1.0e10f;
    dummyHitO[o] = 0;
    overlapHitO[o] = 0;
  }
  for(int o = 0; o < numberOfOrientations; o++)
  {
    visitPairs(point, position, epsilonTimesFour, sigmaSquared, contact, shiftValue, chargeProduct, siteOffset,
               siteCharge, siteDispersion, screenedTable, smoothPotential, potentialGridSize, numberOfAtoms,
               numberOfSites, shells, widths, cellRow0, cellRow1, cellRow2, cutOffSquared, coulombCutOffSquared,
               longestCutOff, longestCutOffSquared, extraReach, ceiling, screenedScale, screenedBins, alpha, o, 1, 0,
               1, &energiesO[o], &clearancesO[o], &dummyHitO[o], &overlapHitO[o]);
  }
  if(useCharges != 0)
  {
    for(int o = 0; o < numberOfOrientations; o++)
    {
      if(dummyHitO[o] != 0 || overlapHitO[o] != 0) continue;
      visitPairs(point, position, epsilonTimesFour, sigmaSquared, contact, shiftValue, chargeProduct, siteOffset,
                 siteCharge, siteDispersion, screenedTable, smoothPotential, potentialGridSize, numberOfAtoms,
                 numberOfSites, shells, widths, cellRow0, cellRow1, cellRow2, cutOffSquared, coulombCutOffSquared,
                 longestCutOff, longestCutOffSquared, extraReach, ceiling, screenedScale, screenedBins, alpha, o, 0,
                 1, 0, &energiesO[o], &clearancesO[o], &dummyHitO[o], 0);
    }
  }
  int bestOrientation = 0;
  float nearest = 1.0e10f;
  reduceOrientations(energiesO, clearancesO, dummyHitO, numberOfOrientations, beta, kT, ceiling, &bestOrientation,
                     &nearest);
  if(pocket < nearest) return;

  float4 directionSum = (float4)(0.0f);
  float weightSum = 0.0f;
  for(int site = 0; site < numberOfSites; site++)
  {
    if(siteDispersion[site] == 0) continue;
    float4 s = point + siteOffset[bestOrientation * numberOfSites + site];
    s.w = 0.0f;
    __global const float *contactForSite = contact + site * numberOfAtoms;
    for(int iatom = 0; iatom < numberOfAtoms; iatom++)
    {
      float4 ds = s - position[iatom];
      ds -= rint(ds);
      ds.w = 0.0f;
      for(int a = -shells.x; a <= shells.x; a++)
      {
        for(int b = -shells.y; b <= shells.y; b++)
        {
          for(int c = -shells.z; c <= shells.z; c++)
          {
            float4 t = ds + (float4)((float)(a), (float)(b), (float)(c), 0.0f);
            float4 dr = cartesianOf(t, cellRow0, cellRow1, cellRow2);
            float r = sqrt(fmax(dot(dr, dr), 1.0e-16f));
            if(!(r > 1.0e-6f)) continue;
            float weighted = r - contactForSite[iatom];
            if(weighted - nearest > 6.0f * softminTau) continue;
            float w = exp(-(weighted - nearest) / softminTau);
            directionSum += (-dr / r) * w;
            weightSum += w;
          }
        }
      }
    }
  }
  float reliability = (weightSum > 0.0f) ? sqrt(dot(directionSum, directionSum)) / weightSum : 0.0f;
  float span = 0.7f * smoothstepf(0.25f, 0.6f, reliability);
  float dir2 = dot(directionSum, directionSum);
  if(span < 0.05f || dir2 < 1.0e-16f) return;
  float4 direction = directionSum * (1.0f / sqrt(dir2));
  direction.w = 0.0f;

  float4 rayFractional;
  rayFractional.x = dot(invRow0, direction);
  rayFractional.y = dot(invRow1, direction);
  rayFractional.z = dot(invRow2, direction);
  rayFractional.w = 0.0f;

  float sBest = 0.0f;
  df uBest = energyHereDf;
  for(int i = -COARSE_STEPS; i <= COARSE_STEPS; i++)
  {
    float s = span * (float)(i) / (float)(COARSE_STEPS);
    df u = energyAtShift(position, epsilonTimesFour, sigmaSquared, contact, shiftValue, chargeProduct, siteOffset,
                         siteCharge, siteDispersion, screenedTable, smoothPotential, point, direction * s,
                         rayFractional * s, potentialGridSize, numberOfAtoms, numberOfSites, numberOfOrientations,
                         shells, widths, cellRow0, cellRow1, cellRow2, cutOffSquared, coulombCutOffSquared,
                         longestCutOff, longestCutOffSquared, ceiling, beta, kT, screenedScale, screenedBins, alpha,
                         useCharges, extraReach);
    if(df_lt(u, uBest))
    {
      uBest = u;
      sBest = s;
    }
  }
  if(fabs(sBest) >= span - 0.5f * span / (float)(COARSE_STEPS)) return;

  // Same 20-step golden section as the processor. The two interior samples are
  // compared as double-floats, so the last iterations still see the slope.
  float a = sBest - span / (float)(COARSE_STEPS);
  float b = sBest + span / (float)(COARSE_STEPS);
  float x1 = b - GOLDEN_INVPHI * (b - a);
  float x2 = a + GOLDEN_INVPHI * (b - a);
  df f1 = energyAtShift(position, epsilonTimesFour, sigmaSquared, contact, shiftValue, chargeProduct, siteOffset,
                        siteCharge, siteDispersion, screenedTable, smoothPotential, point, direction * x1,
                        rayFractional * x1, potentialGridSize, numberOfAtoms, numberOfSites, numberOfOrientations,
                        shells, widths, cellRow0, cellRow1, cellRow2, cutOffSquared, coulombCutOffSquared,
                        longestCutOff, longestCutOffSquared, ceiling, beta, kT, screenedScale, screenedBins, alpha,
                        useCharges, extraReach);
  df f2 = energyAtShift(position, epsilonTimesFour, sigmaSquared, contact, shiftValue, chargeProduct, siteOffset,
                        siteCharge, siteDispersion, screenedTable, smoothPotential, point, direction * x2,
                        rayFractional * x2, potentialGridSize, numberOfAtoms, numberOfSites, numberOfOrientations,
                        shells, widths, cellRow0, cellRow1, cellRow2, cutOffSquared, coulombCutOffSquared,
                        longestCutOff, longestCutOffSquared, ceiling, beta, kT, screenedScale, screenedBins, alpha,
                        useCharges, extraReach);
  for(int iteration = 0; iteration < GOLDEN_STEPS; iteration++)
  {
    if(df_lt(f1, f2))
    {
      b = x2;
      x2 = x1;
      f2 = f1;
      x1 = b - GOLDEN_INVPHI * (b - a);
      f1 = energyAtShift(position, epsilonTimesFour, sigmaSquared, contact, shiftValue, chargeProduct, siteOffset,
                         siteCharge, siteDispersion, screenedTable, smoothPotential, point, direction * x1,
                         rayFractional * x1, potentialGridSize, numberOfAtoms, numberOfSites, numberOfOrientations,
                         shells, widths, cellRow0, cellRow1, cellRow2, cutOffSquared, coulombCutOffSquared,
                         longestCutOff, longestCutOffSquared, ceiling, beta, kT, screenedScale, screenedBins, alpha,
                         useCharges, extraReach);
    }
    else
    {
      a = x1;
      x1 = x2;
      f1 = f2;
      x2 = a + GOLDEN_INVPHI * (b - a);
      f2 = energyAtShift(position, epsilonTimesFour, sigmaSquared, contact, shiftValue, chargeProduct, siteOffset,
                         siteCharge, siteDispersion, screenedTable, smoothPotential, point, direction * x2,
                         rayFractional * x2, potentialGridSize, numberOfAtoms, numberOfSites, numberOfOrientations,
                         shells, widths, cellRow0, cellRow1, cellRow2, cutOffSquared, coulombCutOffSquared,
                         longestCutOff, longestCutOffSquared, ceiling, beta, kT, screenedScale, screenedBins, alpha,
                         useCharges, extraReach);
    }
  }
  float s = 0.5f * (a + b);
  float4 moved = point + rayFractional * s;
  float3 wrapped = wrapUnit(moved.xyz);
  corners[vertex] = (float4)(wrapped.x, wrapped.y, wrapped.z, 0.0f);
  energies[vertex] = df_hi(energyAtShift(position, epsilonTimesFour, sigmaSquared, contact, shiftValue, chargeProduct,
                                         siteOffset, siteCharge, siteDispersion, screenedTable, smoothPotential, point,
                                         direction * s, rayFractional * s, potentialGridSize, numberOfAtoms,
                                         numberOfSites, numberOfOrientations, shells, widths, cellRow0, cellRow1,
                                         cellRow2, cutOffSquared, coulombCutOffSquared, longestCutOff,
                                         longestCutOffSquared, ceiling, beta, kT, screenedScale, screenedBins, alpha,
                                         useCharges, extraReach));
  (void)siteLength;
  (void)axes;
  (void)blockedEnergyPerAngstrom;
}
)foo";
