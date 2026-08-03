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

module energy_opencl_electrostatic_potential_grid;

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

import energy_shared_electrostatic_potential_grid;
import energy_shared_ewald;



ElectrostaticPotentialGrid ElectrostaticPotentialGridOpenCL::compute(const PairInteractions &interactions, const Crystal &framework,
                                                               uint3 gridSize, double relativePrecision,
                                                               std::optional<double> alphaOverride)
{
  if (!OpenCL::clContext.has_value() || !OpenCL::clDeviceId.has_value() || !OpenCL::clCommandQueue.has_value())
  {
    throw std::runtime_error("ElectrostaticPotentialGridOpenCL: no OpenCL device found, the grid route needs a GPU\n");
  }

  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();

  ElectrostaticPotentialGrid grid;
  grid.gridSize = gridSize;
  grid.backend = "gpu";
  grid.unitCell = framework.unitCell;
  grid.cutOff = interactions.cutOffCoulomb;
  grid.relativePrecision = relativePrecision;

  EwaldSplit split = ewaldSplit(grid.cutOff, relativePrecision, alphaOverride);
  grid.alpha = split.alpha;
  grid.largestWaveVector = split.largestWaveVector;

  std::size_t numberOfVoxels = gridSize.x * gridSize.y * gridSize.z;
  grid.smoothPotential.assign(numberOfVoxels, 0.0f);

  std::vector<double3> fractionalPositions = framework.fractionalPositions;
  if (fractionalPositions.empty()) return grid;

  for (const CrystalAtom &atom : framework.atoms)
  {
    grid.netCharge += atom.charge;
    grid.largestFrameworkCharge = std::max(grid.largestFrameworkCharge, std::abs(atom.charge));
  }

  std::vector<WaveVector> vectors = waveVectors(framework, fractionalPositions, grid.alpha, grid.largestWaveVector);
  grid.numberOfWaveVectors = vectors.size();

  double background = neutralisingBackground(grid.netCharge, grid.alpha, framework.unitCell.volume);

  std::vector<cl_float4> waves(vectors.size());
  std::vector<cl_float> imaginary(vectors.size());
  for (std::size_t i = 0; i < vectors.size(); ++i)
  {
    waves[i] = {{cl_float(vectors[i].h), cl_float(vectors[i].k), cl_float(vectors[i].l),
                 cl_float(Units::CoulombicConversionFactor * vectors[i].weightedReal)}};
    imaginary[i] = cl_float(Units::CoulombicConversionFactor * vectors[i].weightedImaginary);
  }

  cl_int err;

  const char *source = ElectrostaticPotentialGridOpenCL::reciprocalKernelSource;
  cl_program program = clCreateProgramWithSource(OpenCL::clContext.value(), 1, &source, nullptr, &err);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(
        std::format("ElectrostaticPotentialGridOpenCL: OpenCL clCreateProgramWithSource failed at {}\n", __LINE__));
  }

  err = clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    std::size_t length;
    char buffer[2048];
    clGetProgramBuildInfo(program, OpenCL::clDeviceId.value(), CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &length);
    throw std::runtime_error(std::format(
        "ElectrostaticPotentialGridOpenCL: OpenCL failed to build program at {} (error: {})\n", __LINE__,
        std::string(buffer)));
  }

  cl_kernel kernel = clCreateKernel(program, "ReciprocalPotentialField", &err);
  if (err != CL_SUCCESS)
  {
    clReleaseProgram(program);
    throw std::runtime_error(
        std::format("ElectrostaticPotentialGridOpenCL: OpenCL clCreateKernel failed at {}\n", __LINE__));
  }

  cl_int waveError, imaginaryError, potentialError;
  cl_mem waveBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY,
                                     sizeof(cl_float4) * std::max<std::size_t>(1, waves.size()), nullptr, &waveError);
  cl_mem imaginaryBuffer =
      clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY,
                     sizeof(cl_float) * std::max<std::size_t>(1, imaginary.size()), nullptr, &imaginaryError);
  cl_mem potentialBuffer = clCreateBuffer(OpenCL::clContext.value(), CL_MEM_WRITE_ONLY,
                                          sizeof(cl_float) * numberOfVoxels, nullptr, &potentialError);
  if (waveError != CL_SUCCESS || imaginaryError != CL_SUCCESS || potentialError != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("ElectrostaticPotentialGridOpenCL: OpenCL clCreateBuffer failed at {}\n", __LINE__));
  }

  if (!waves.empty())
  {
    err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), waveBuffer, CL_TRUE, 0, sizeof(cl_float4) * waves.size(),
                               waves.data(), 0, nullptr, nullptr);
    err |= clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), imaginaryBuffer, CL_TRUE, 0,
                                sizeof(cl_float) * imaginary.size(), imaginary.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(
          std::format("ElectrostaticPotentialGridOpenCL: OpenCL clEnqueueWriteBuffer failed at {}\n", __LINE__));
    }
  }

  cl_int3 clGridSize = {{cl_int(gridSize.x), cl_int(gridSize.y), cl_int(gridSize.z), 0}};
  cl_int clNumberOfWaves = static_cast<cl_int>(waves.size());
  cl_float clBackground = cl_float(background);

  err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &waveBuffer);
  err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &imaginaryBuffer);
  err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &potentialBuffer);
  err |= clSetKernelArg(kernel, 3, sizeof(cl_int), &clNumberOfWaves);
  err |= clSetKernelArg(kernel, 4, sizeof(cl_int3), &clGridSize);
  err |= clSetKernelArg(kernel, 5, sizeof(cl_float), &clBackground);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("ElectrostaticPotentialGridOpenCL: OpenCL clSetKernelArg failed at {}\n", __LINE__));
  }

  std::size_t globalWorkSize[3] = {static_cast<std::size_t>(gridSize.x), static_cast<std::size_t>(gridSize.y),
                                   static_cast<std::size_t>(gridSize.z)};
  err = clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), kernel, 3, nullptr, globalWorkSize, nullptr, 0, nullptr,
                               nullptr);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(
        std::format("ElectrostaticPotentialGridOpenCL: OpenCL clEnqueueNDRangeKernel failed at {}\n", __LINE__));
  }

  err = clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), potentialBuffer, CL_TRUE, 0,
                            sizeof(cl_float) * numberOfVoxels, grid.smoothPotential.data(), 0, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(
        std::format("ElectrostaticPotentialGridOpenCL: OpenCL clEnqueueReadBuffer failed at {}\n", __LINE__));
  }

  clFinish(OpenCL::clCommandQueue.value());

  clReleaseMemObject(waveBuffer);
  clReleaseMemObject(imaginaryBuffer);
  clReleaseMemObject(potentialBuffer);
  clReleaseKernel(kernel);
  clReleaseProgram(program);

  std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - time_begin;
  grid.seconds = elapsed.count();

  return grid;
}


const char *ElectrostaticPotentialGridOpenCL::reciprocalKernelSource = R"foo(
__kernel void ReciprocalPotentialField(__global const float4 *wave,
                                       __global const float *weightedImaginary,
                                       __global float *potential,
                                       const int numberOfWaves,
                                       const int3 gridSize,
                                       const float background)
{
  int ix = get_global_id(0);
  int iy = get_global_id(1);
  int iz = get_global_id(2);

  if(ix >= gridSize.x || iy >= gridSize.y || iz >= gridSize.z) return;

  float3 s = (float3)((float)(ix) / (float)(gridSize.x),
                      (float)(iy) / (float)(gridSize.y),
                      (float)(iz) / (float)(gridSize.z));

  // Holding the wave vectors as whole numbers of reciprocal cells means the phase is just their dot product
  // with the fractional position, so nothing here needs to know the shape of the cell.
  float total = background;
  for(int i = 0; i < numberOfWaves; i++)
  {
    float4 w = wave[i];
    float phase = 2.0f * M_PI_F * (w.x * s.x + w.y * s.y + w.z * s.z);
    total += w.w * cos(phase) - weightedImaginary[i] * sin(phase);
  }

  potential[(iz * gridSize.y + iy) * gridSize.x + ix] = total;
}
)foo";
