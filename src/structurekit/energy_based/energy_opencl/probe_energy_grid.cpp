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

module energy_opencl_probe_energy_grid;

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

import energy_shared_probe_energy_grid;


ProbeEnergyGrid ProbeEnergyGridOpenCL::compute(const PairInteractions &interactions, const Crystal &framework,
                                               std::string probePseudoAtom, uint3 gridSize)
{
  if (!OpenCL::clContext.has_value() || !OpenCL::clDeviceId.has_value() || !OpenCL::clCommandQueue.has_value())
  {
    throw std::runtime_error("ProbeEnergyGridOpenCL: no OpenCL device found, the grid route needs a GPU\n");
  }
  if (gridSize.x == 0 || gridSize.y == 0 || gridSize.z == 0)
  {
    throw std::runtime_error("ProbeEnergyGridOpenCL: the grid must have at least one point along each axis\n");
  }

  std::optional<std::size_t> probeType = interactions.findType(probePseudoAtom);
  if (!probeType.has_value())
  {
    throw std::runtime_error(
        std::format("ProbeEnergyGridOpenCL: probe atom '{}' not found in the force field\n", probePseudoAtom));
  }

  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  ProbeEnergyGrid grid;
  grid.gridSize = gridSize;
  grid.unitCell = framework.unitCell;
  grid.backend = "gpu";
  grid.probeName = probePseudoAtom;
  grid.probeEpsilon = interactions[probeType.value()].strengthParameter;
  grid.probeSigma = interactions[probeType.value()].sizeParameter;
  grid.cutOff = interactions.cutOffVDW;

  grid.ceiling = probeEnergyCeilingInKelvin * Units::KelvinToEnergy;

  std::size_t numberOfVoxels = gridSize.x * gridSize.y * gridSize.z;
  grid.energy.assign(numberOfVoxels, 0.0f);
  grid.strongestAtom.assign(numberOfVoxels, -1);

  std::vector<double3> fractionalPositions = framework.fractionalPositions;
  if (fractionalPositions.empty()) return grid;

  // The mixing is the force field's own, taken for the pair rather than assembled here, so that the field
  // agrees with what a simulation with this probe would feel.
  std::vector<cl_float4> positions(fractionalPositions.size());
  std::vector<cl_float> epsilonTimesFour(fractionalPositions.size());
  std::vector<cl_float> sigmaSquared(fractionalPositions.size());
  std::vector<cl_float> shiftValue(fractionalPositions.size());
  for (std::size_t i = 0; i < fractionalPositions.size(); ++i)
  {
    std::size_t atomType = framework.atoms[i].type;
    double sigma = interactions(probeType.value(), atomType).sizeParameter;
    double epsilon = interactions(probeType.value(), atomType).strengthParameter;

    positions[i] = {{cl_float(fractionalPositions[i].x), cl_float(fractionalPositions[i].y),
                     cl_float(fractionalPositions[i].z), 0.0f}};
    epsilonTimesFour[i] = cl_float(4.0 * epsilon);
    sigmaSquared[i] = cl_float(sigma * sigma);

    // A truncated potential is usually shifted so that it reaches the cutoff at zero rather than stepping
    // there, and the force field says whether it is. Leaving the step in puts the energy out by however many
    // neighbours are in range, which varies from place to place and so bends the landscape rather than
    // raising it.
    shiftValue[i] = cl_float(interactions(probeType.value(), atomType).shift);
  }

  double3x3 cell = framework.unitCell.cell;

  // A fractional difference t is at least |t_k| w_k long in Cartesian terms, w_k being the width of the cell
  // perpendicular to the plane of the other two axes. Reducing into [-1/2, 1/2) leaves |t_k| at least
  // |n_k| - 1/2 for an image n_k cells out, so an image can only come within the cutoff while
  // |n_k| <= cutoff / w_k + 1/2, and that settles how far to look without any reference to the cell's shape.
  double3 widths = framework.unitCell.perpendicularWidths();
  grid.numberOfImageShells =
      int3(static_cast<std::int32_t>(std::floor(grid.cutOff / widths.x + 0.5)),
           static_cast<std::int32_t>(std::floor(grid.cutOff / widths.y + 0.5)),
           static_cast<std::int32_t>(std::floor(grid.cutOff / widths.z + 0.5)));

  cl_int err;

  const char *source = ProbeEnergyGridOpenCL::energyKernelSource;
  cl_program program = clCreateProgramWithSource(OpenCL::clContext.value(), 1, &source, nullptr, &err);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("ProbeEnergyGridOpenCL: OpenCL clCreateProgramWithSource failed at {}\n", __LINE__));
  }

  err = clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    std::size_t length;
    char buffer[2048];
    clGetProgramBuildInfo(program, OpenCL::clDeviceId.value(), CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &length);
    throw std::runtime_error(std::format("ProbeEnergyGridOpenCL: OpenCL failed to build program at {} (error: {})\n",
                                         __LINE__, std::string(buffer)));
  }

  cl_kernel kernel = clCreateKernel(program, "ProbeEnergyField", &err);
  if (err != CL_SUCCESS)
  {
    clReleaseProgram(program);
    throw std::runtime_error(std::format("ProbeEnergyGridOpenCL: OpenCL clCreateKernel failed at {}\n", __LINE__));
  }

  cl_int shiftError;
  cl_int positionError, epsilonError, sigmaError, energyError;
  cl_mem positionBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY,
                                         sizeof(cl_float4) * positions.size(), nullptr, &positionError);
  cl_mem epsilonBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY,
                                        sizeof(cl_float) * epsilonTimesFour.size(), nullptr, &epsilonError);
  cl_mem sigmaBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY,
                                      sizeof(cl_float) * sigmaSquared.size(), nullptr, &sigmaError);
  cl_mem shiftBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY,
                                      sizeof(cl_float) * shiftValue.size(), nullptr, &shiftError);
  cl_mem energyBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_WRITE_ONLY, sizeof(cl_float) * numberOfVoxels,
                                       nullptr, &energyError);
  cl_int strongestError;
  cl_mem strongestBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_WRITE_ONLY,
                                          sizeof(cl_int) * numberOfVoxels, nullptr, &strongestError);
  if (positionError != CL_SUCCESS || epsilonError != CL_SUCCESS || sigmaError != CL_SUCCESS ||
      energyError != CL_SUCCESS || strongestError != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("ProbeEnergyGridOpenCL: OpenCL clCreateBuffer failed at {}\n", __LINE__));
  }

  err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), positionBuffer, CL_TRUE, 0,
                             sizeof(cl_float4) * positions.size(), positions.data(), 0, nullptr, nullptr);
  err |= clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), epsilonBuffer, CL_TRUE, 0,
                              sizeof(cl_float) * epsilonTimesFour.size(), epsilonTimesFour.data(), 0, nullptr, nullptr);
  err |= clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), sigmaBuffer, CL_TRUE, 0,
                              sizeof(cl_float) * sigmaSquared.size(), sigmaSquared.data(), 0, nullptr, nullptr);
  err |= clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), shiftBuffer, CL_TRUE, 0,
                              sizeof(cl_float) * shiftValue.size(), shiftValue.data(), 0, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("ProbeEnergyGridOpenCL: OpenCL clEnqueueWriteBuffer failed at {}\n", __LINE__));
  }

  cl_float4 cellRow0 = {{cl_float(cell[0][0]), cl_float(cell[1][0]), cl_float(cell[2][0]), 0.0f}};
  cl_float4 cellRow1 = {{cl_float(cell[0][1]), cl_float(cell[1][1]), cl_float(cell[2][1]), 0.0f}};
  cl_float4 cellRow2 = {{cl_float(cell[0][2]), cl_float(cell[1][2]), cl_float(cell[2][2]), 0.0f}};
  cl_float4 clWidths = {{cl_float(widths.x), cl_float(widths.y), cl_float(widths.z), 0.0f}};
  cl_int3 clShells = {{cl_int(grid.numberOfImageShells.x), cl_int(grid.numberOfImageShells.y),
                       cl_int(grid.numberOfImageShells.z), 0}};
  cl_int3 clGridSize = {{cl_int(gridSize.x), cl_int(gridSize.y), cl_int(gridSize.z), 0}};
  cl_int numberOfAtoms = static_cast<cl_int>(positions.size());
  cl_float clCutOffSquared = cl_float(grid.cutOff * grid.cutOff);
  cl_float clCutOff = cl_float(grid.cutOff);
  cl_float clCeiling = cl_float(grid.ceiling);

  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &positionBuffer);
  err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &epsilonBuffer);
  err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &sigmaBuffer);
  err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &shiftBuffer);
  err |= clSetKernelArg(kernel, 4, sizeof(cl_mem), &energyBuffer);
  err |= clSetKernelArg(kernel, 5, sizeof(cl_mem), &strongestBuffer);
  err |= clSetKernelArg(kernel, 6, sizeof(cl_float4), &cellRow0);
  err |= clSetKernelArg(kernel, 7, sizeof(cl_float4), &cellRow1);
  err |= clSetKernelArg(kernel, 8, sizeof(cl_float4), &cellRow2);
  err |= clSetKernelArg(kernel, 9, sizeof(cl_float4), &clWidths);
  err |= clSetKernelArg(kernel, 10, sizeof(cl_int), &numberOfAtoms);
  err |= clSetKernelArg(kernel, 11, sizeof(cl_int3), &clShells);
  err |= clSetKernelArg(kernel, 12, sizeof(cl_int3), &clGridSize);
  err |= clSetKernelArg(kernel, 13, sizeof(cl_float), &clCutOff);
  err |= clSetKernelArg(kernel, 14, sizeof(cl_float), &clCutOffSquared);
  err |= clSetKernelArg(kernel, 15, sizeof(cl_float), &clCeiling);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("ProbeEnergyGridOpenCL: OpenCL clSetKernelArg failed at {}\n", __LINE__));
  }

  std::size_t globalWorkSize[3] = {static_cast<std::size_t>(gridSize.x), static_cast<std::size_t>(gridSize.y),
                                   static_cast<std::size_t>(gridSize.z)};
  err = clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), kernel, 3, nullptr, globalWorkSize, nullptr, 0, nullptr,
                               nullptr);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("ProbeEnergyGridOpenCL: OpenCL clEnqueueNDRangeKernel failed at {}\n", __LINE__));
  }

  err = clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), energyBuffer, CL_TRUE, 0, sizeof(cl_float) * numberOfVoxels,
                            grid.energy.data(), 0, nullptr, nullptr);
  err |= clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), strongestBuffer, CL_TRUE, 0,
                             sizeof(cl_int) * numberOfVoxels, grid.strongestAtom.data(), 0, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("ProbeEnergyGridOpenCL: OpenCL clEnqueueReadBuffer failed at {}\n", __LINE__));
  }

  clFinish(OpenCL::clCommandQueue.value());

  clReleaseMemObject(positionBuffer);
  clReleaseMemObject(epsilonBuffer);
  clReleaseMemObject(sigmaBuffer);
  clReleaseMemObject(shiftBuffer);
  clReleaseMemObject(energyBuffer);
  clReleaseMemObject(strongestBuffer);
  clReleaseKernel(kernel);
  clReleaseProgram(program);

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  grid.seconds = elapsed.count();

  return grid;
}


const char *ProbeEnergyGridOpenCL::energyKernelSource = R"foo(
__kernel void ProbeEnergyField(__global const float4 *position,
                               __global const float *epsilonTimesFour,
                               __global const float *sigmaSquared,
                               __global const float *shiftValue,
                               __global float *energy,
                               __global int *strongestAtom,
                               const float4 cellRow0,
                               const float4 cellRow1,
                               const float4 cellRow2,
                               const float4 widths,
                               const int numberOfAtoms,
                               const int3 shells,
                               const int3 gridSize,
                               const float cutOff,
                               const float cutOffSquared,
                               const float ceiling)
{
  int ix = get_global_id(0);
  int iy = get_global_id(1);
  int iz = get_global_id(2);

  if(ix >= gridSize.x || iy >= gridSize.y || iz >= gridSize.z) return;

  // Endpoint-exclusive sampling, as in the clearance field: fractional 0 and 1 are the same periodic point,
  // so dividing by the grid size rather than one less keeps the spacing uniform and every sample distinct.
  // The two fields have to be sampled alike for anything to be read off both.
  float4 s = (float4)((float)(ix) / (float)(gridSize.x),
                      (float)(iy) / (float)(gridSize.y),
                      (float)(iz) / (float)(gridSize.z),
                      0.0f);

  float total = 0.0f;

  // The atom pulling hardest, and the nearest one in case none of them pulls. An atom's claim is its own
  // contribution to the sum, over all of its images within the cutoff, so that an atom is not credited twice
  // for being near in a small cell.
  //
  // Inside a wall every term is held at the ceiling, which makes the atoms there exactly equal and the
  // strongest of them a matter of which was looked at first. Distance separates them and energy does not, so
  // that region falls back on the nearest atom, which is also what the geometric route would say about it.
  float strongestPull = 0.0f;
  float nearestDistanceSquared = MAXFLOAT;
  int pullingAtom = -1;
  int nearestAtom = -1;

  for(int iatom = 0; iatom < numberOfAtoms; iatom++)
  {
    float contribution = 0.0f;
    float closest = MAXFLOAT;
    float4 ds = s - position[iatom];
    ds -= rint(ds);
    ds.w = 0.0f;

    float epsilon4 = epsilonTimesFour[iatom];
    float sigma2 = sigmaSquared[iatom];
    float shift = shiftValue[iatom];

    for(int a = -shells.x; a <= shells.x; a++)
    {
      for(int b = -shells.y; b <= shells.y; b++)
      {
        for(int c = -shells.z; c <= shells.z; c++)
        {
          float4 t = ds + (float4)((float)(a), (float)(b), (float)(c), 0.0f);

          // The same inequality that sizes the shells above, used again to throw out an image before it is
          // built. Most of the shell lies outside the cutoff on all but the smallest cells, and this settles
          // that with three multiplies rather than three dot products and a square.
          float reach = fmax(fmax(fabs(t.x) * widths.x, fabs(t.y) * widths.y), fabs(t.z) * widths.z);
          if(reach > cutOff) continue;

          float4 dr;
          dr.x = dot(cellRow0, t);
          dr.y = dot(cellRow1, t);
          dr.z = dot(cellRow2, t);
          dr.w = 0.0f;

          float rr = dot(dr, dr);
          if(rr >= cutOffSquared) continue;

          // A grid point can land on an atom's centre, where the pair energy has no value. The separation is
          // held off zero so that the sum stays a number; the point is buried far above the ceiling either
          // way, so nothing that is read off the field can turn on the figure chosen here.
          rr = fmax(rr, 1.0e-6f);
          closest = fmin(closest, rr);

          float ratio = sigma2 / rr;
          float ratio3 = ratio * ratio * ratio;

          // Each term is held down before it is added rather than the sum afterwards, so that no term ever
          // overflows to something a later addition would turn into a value that is not a number.
          contribution += min(epsilon4 * ratio3 * (ratio3 - 1.0f) - shift, ceiling);
        }
      }
    }

    total += contribution;

    if(contribution < strongestPull)
    {
      strongestPull = contribution;
      pullingAtom = iatom;
    }
    if(closest < nearestDistanceSquared)
    {
      nearestDistanceSquared = closest;
      nearestAtom = iatom;
    }
  }

  int voxel = (iz * gridSize.y + iy) * gridSize.x + ix;
  energy[voxel] = min(total, ceiling);
  strongestAtom[voxel] = (pullingAtom >= 0) ? pullingAtom : nearestAtom;
}
)foo";
