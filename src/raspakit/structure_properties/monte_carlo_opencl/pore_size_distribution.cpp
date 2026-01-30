module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <print>
#include <string>
#include <vector>
#include <chrono>
#include <thread>
#pragma push_macro("__SSE3__")
#undef __SSE3__
#include <random>
#pragma pop_macro("__SSE3__")
#endif

#define CL_TARGET_OPENCL_VERSION 120
#define CL_SILENCE_DEPRECATION
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

module mc_opencl_pore_size_distribution;

#ifdef USE_STD_IMPORT
import std;
#endif

import opencl;
import double2;
import double3;
import double4;
import float4;
import double3x3;
import randomnumbers;
import framework;
import forcefield;
import atom;


MC_OpenCL_PoreSizeDistribution::MC_OpenCL_PoreSizeDistribution(std::size_t numberOfBins) : numberOfBins(numberOfBins),
             histogram(numberOfBins), histogram_cummulative(numberOfBins)
{
  if (OpenCL::clContext.has_value() && OpenCL::clDeviceId.has_value())
  {
    cl_int err;

    const char* poreSizeDistributionShaderSourceCode = MC_OpenCL_PoreSizeDistribution::poreSizeDistributionKernelSource;
    poreSizeDistributionProgram =
        clCreateProgramWithSource(OpenCL::clContext.value(), 1, &poreSizeDistributionShaderSourceCode, nullptr, &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("MC_OpenCL_PoreSizeDistribution: OpenCL clCreateProgramWithSource failed at {}\n", __LINE__));
    }

    err = clBuildProgram(poreSizeDistributionProgram, 0, nullptr, nullptr, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      size_t len;
      char buffer[2048];
      clGetProgramBuildInfo(poreSizeDistributionProgram, OpenCL::clDeviceId.value(), CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer,
                            &len);
      std::string message =
          std::format("MC_OpenCL_PoreSizeDistribution: OpenCL Failed to build program at {} (line {} error: {})\n", __FILE__,
                      __LINE__, std::string(buffer));
      throw std::runtime_error(message);
    }

    poreSizeDistributionKernel = clCreateKernel(poreSizeDistributionProgram, "ComputePoreSizeDistribution", &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("MC_OpenCL_PoreSizeDistribution: OpenCL clCreateKernel failed {} : {}\n", __FILE__, __LINE__));
    }
    err = clGetKernelWorkGroupInfo(poreSizeDistributionKernel, OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE,
                                   sizeof(size_t), &poreSizeDistributionWorkGroupSize, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("MC_OpenCL_PoreSizeDistribution: OpenCL clGetKernelWorkGroupInfo failed at {} : {}\n", __FILE__, __LINE__));
    }
  }
};

MC_OpenCL_PoreSizeDistribution::~MC_OpenCL_PoreSizeDistribution()
{
  if (OpenCL::clContext.has_value() && OpenCL::clDeviceId.has_value())
  {
    clReleaseKernel(poreSizeDistributionKernel);
    clReleaseProgram(poreSizeDistributionProgram);
  }
}


void MC_OpenCL_PoreSizeDistribution::run(const ForceField &forceField,
                                         const Framework &framework,
                                         double wellDepthFactor,
                                         std::size_t numberOfIterations,
                                         std::optional<std::size_t> numberOfInnerSteps)
{
  RandomNumber random{std::nullopt};
  cl_int err;
  std::chrono::system_clock::time_point time_begin, time_end;
  std::size_t number_of_inner_steps = numberOfInnerSteps.value_or(10000);
  std::size_t histogram_size{ 1000 };

  size_t global_work_size = (number_of_inner_steps + poreSizeDistributionWorkGroupSize - 1) & ~(poreSizeDistributionWorkGroupSize - 1);

  double3x3 unit_cell = framework.simulationBox.cell;
  double3x3 inverse_unit_cell = framework.simulationBox.inverseCell;
  std::vector<double3> positions = framework.cartesianAtomPositionsUnitCell();

  std::size_t numberOfAtoms = positions.size();

  std::vector<cl_float4> framework_positions(numberOfAtoms);
  std::vector<cl_float> sigma(numberOfAtoms);


  time_begin = std::chrono::system_clock::now();

  for (size_t i = 0; i < numberOfAtoms; i++)
  {
    std::size_t atomType = static_cast<std::size_t>(framework.unitCellAtoms[i].type);
    double size_parameter = forceField[atomType].sizeParameter();
    double equilibrium_distance = wellDepthFactor * size_parameter;

    // fill in the Cartesian position
    double3 position = positions[i];
    framework_positions[i] = {{cl_float(position.x), cl_float(position.y), cl_float(position.z), 0.0f}};

    sigma[i] = cl_float(equilibrium_distance);
  }


  // upload framework position array
  cl_mem framework_positions_mem =
        clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY, sizeof(float4) * framework_positions.size(), nullptr, &err);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));
  }

  err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), framework_positions_mem, CL_TRUE, 0, sizeof(float4) * framework_positions.size(),
                               framework_positions.data(), 0, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("OpenCL clCommandQueue failed {} : {}\n", __FILE__, __LINE__));
  }


  // upload equilibrium size array
  cl_mem inputSigma =
      clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY, sizeof(cl_float) * sigma.size(), nullptr, &err);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));
  }

  err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), inputSigma, CL_TRUE, 0, sizeof(cl_float) * sigma.size(),
                             sigma.data(), 0, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    throw std::runtime_error(std::format("OpenCL clCommandQueue failed {} : {}\n", __FILE__, __LINE__));
  }



  // setup unit cell matrix
  cl_float4 clCella = {{cl_float(unit_cell[0][0]), cl_float(unit_cell[1][0]), cl_float(unit_cell[2][0]), cl_float(0.0)}};
  cl_float4 clCellb = {{cl_float(unit_cell[0][1]), cl_float(unit_cell[1][1]), cl_float(unit_cell[2][1]), cl_float(0.0)}};
  cl_float4 clCellc = {{cl_float(unit_cell[0][2]), cl_float(unit_cell[1][2]), cl_float(unit_cell[2][2]), cl_float(0.0)}};

  cl_float4 inverse_clCella = {{cl_float(inverse_unit_cell[0][0]), cl_float(inverse_unit_cell[1][0]), cl_float(inverse_unit_cell[2][0]), cl_float(0.0)}};
  cl_float4 inverse_clCellb = {{cl_float(inverse_unit_cell[0][1]), cl_float(inverse_unit_cell[1][1]), cl_float(inverse_unit_cell[2][1]), cl_float(0.0)}};
  cl_float4 inverse_clCellc = {{cl_float(inverse_unit_cell[0][2]), cl_float(inverse_unit_cell[1][2]), cl_float(inverse_unit_cell[2][2]), cl_float(0.0)}};

  cl_int number_of_atoms = cl_int(numberOfAtoms);
  cl_int number_of_inner_steps_parameter = cl_int(number_of_inner_steps);

  double delta_r = 10.0 / static_cast<double>(numberOfBins);

  for(std::size_t i = 0; i < numberOfIterations; ++i)
  {
    double3 fractional_position = {random.uniform(), random.uniform(), random.uniform()};
    double3 cartesian_position =  unit_cell * fractional_position;
    cl_float4 posA = {{cl_float(cartesian_position.x), cl_float(cartesian_position.y), cl_float(cartesian_position.z), 0.0f}};

    // upload random Cartesian positions
    std::size_t number_of_random_position{global_work_size};
    std::vector<cl_float4> random_cartesian_positions(number_of_random_position);
    for (size_t j = 0; j < number_of_random_position; j++)
    {
      fractional_position = {random.uniform(), random.uniform(), random.uniform()};
      cartesian_position =  unit_cell * fractional_position;
      random_cartesian_positions[j] = {{cl_float(cartesian_position.x), cl_float(cartesian_position.y), cl_float(cartesian_position.z), 0.0f}};
    }
    cl_mem random_cartesian_positions_mem =
          clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_ONLY, sizeof(float4) *  random_cartesian_positions.size(), nullptr, &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));
    }

    err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), random_cartesian_positions_mem, CL_TRUE, 0, sizeof(float4) * random_cartesian_positions.size(),
                                 random_cartesian_positions.data(), 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCommandQueue failed {} : {}\n", __FILE__, __LINE__));
    }

    // setup output array
    cl_mem largest_radius_mem =
          clCreateBuffer(OpenCL::clContext.value(), CL_MEM_READ_WRITE, sizeof(cl_float), nullptr, &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clCreateBuffer failed {} : {}\n", __FILE__, __LINE__));
    }

    float largest_radius[1] = {-1.0};
    err = clEnqueueWriteBuffer(OpenCL::clCommandQueue.value(), largest_radius_mem, CL_TRUE, 0,
                               sizeof(cl_float), largest_radius, 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clEnqueueWriteBuffer failed {} : {}\n", __FILE__, __LINE__));
    }

    err = clSetKernelArg(poreSizeDistributionKernel,   0, sizeof(cl_mem), &framework_positions_mem);
    err |= clSetKernelArg(poreSizeDistributionKernel,  1, sizeof(cl_mem), &random_cartesian_positions_mem);
    err |= clSetKernelArg(poreSizeDistributionKernel,  2, sizeof(cl_mem), &inputSigma);
    err |= clSetKernelArg(poreSizeDistributionKernel,  3, sizeof(cl_mem), &largest_radius_mem);
    err |= clSetKernelArg(poreSizeDistributionKernel,  4, sizeof(cl_float4), &posA);
    err |= clSetKernelArg(poreSizeDistributionKernel,  5, sizeof(cl_int), &number_of_atoms);
    err |= clSetKernelArg(poreSizeDistributionKernel,  6, sizeof(cl_int), &number_of_inner_steps_parameter);
    err |= clSetKernelArg(poreSizeDistributionKernel,  7, sizeof(cl_float4), &clCella);
    err |= clSetKernelArg(poreSizeDistributionKernel,  8, sizeof(cl_float4), &clCellb);
    err |= clSetKernelArg(poreSizeDistributionKernel,  9, sizeof(cl_float4), &clCellc);
    err |= clSetKernelArg(poreSizeDistributionKernel, 10, sizeof(cl_float4), &inverse_clCella);
    err |= clSetKernelArg(poreSizeDistributionKernel, 11, sizeof(cl_float4), &inverse_clCellb);
    err |= clSetKernelArg(poreSizeDistributionKernel, 12, sizeof(cl_float4), &inverse_clCellc);

    err |= clEnqueueNDRangeKernel(OpenCL::clCommandQueue.value(), poreSizeDistributionKernel, 1, nullptr, &global_work_size,
                                  &poreSizeDistributionWorkGroupSize, 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("OpenCL clEnqueueNDRangeKernel failed (err: {}) {} : {}\n", err, __FILE__, __LINE__));
    }

    // Read the buffer back to the array
    err = clEnqueueReadBuffer(OpenCL::clCommandQueue.value(), largest_radius_mem, CL_TRUE, 0, sizeof(float),
                              &largest_radius, 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error("MC_OpenCL_SurfaceArea: error in clEnqueueReadBuffer");
    }

    if(largest_radius[0] >= 0.0f)
    {
      std::size_t index = static_cast<std::size_t>(static_cast<double>(largest_radius[0]) / delta_r);
      if(index >= 0 && index < numberOfBins)
      {
        histogram[index]++;
        for(std::size_t k = 0; k <= index; ++k)
        {
           histogram_cummulative[k]++;
        }
      }
    }

    clReleaseMemObject(random_cartesian_positions_mem);
    clReleaseMemObject(largest_radius_mem);
  }
  clFinish(OpenCL::clCommandQueue.value());

  clReleaseMemObject(inputSigma);
  clReleaseMemObject(framework_positions_mem);

  time_end = std::chrono::system_clock::now();

  std::chrono::duration<double> timing = time_end - time_begin;

  std::ofstream myfile;
  myfile.open(framework.name + ".mc.psd.gpu.txt");
  std::print(myfile, "# Pore-size distribution using Mont Carlo-based method\n");
  std::print(myfile, "# GPU Timing: {} [s]\n", timing.count());

  myfile << "# column 1: diameter d [A]\n";
  myfile << "# column 2: pore size distribution\n";
  myfile << "# column 3: cummulative pore volume\n";
  myfile << "# value at d=0 is related to the void-fraction\n";

  double normalization = 1.0 / static_cast<double>(numberOfIterations);
  for(std::size_t index = 0; index < histogram_size; ++index)
  {
    std::print(myfile, "{} {} {}\n", 2.0 * delta_r * (static_cast<double>(index) + 0.5), 
               histogram[index] * normalization,
               histogram_cummulative[index] * normalization);
  }

  myfile.close();
}

const char* MC_OpenCL_PoreSizeDistribution::poreSizeDistributionKernelSource = R"foo(

// Function to perform the atomic max
// https://stackoverflow.com/questions/18950732/atomic-max-for-floats-in-opencl
//
inline void AtomicMax(volatile __global float *source, const float operand)
{
  union 
  {
    unsigned int intVal;
    float floatVal;
  } newVal;
  union
  {
    unsigned int intVal;
    float floatVal;
  } prevVal;
  do 
  {
     prevVal.floatVal = *source;
     newVal.floatVal = max(prevVal.floatVal,operand);
  } while (atomic_cmpxchg((volatile __global unsigned int *)source, prevVal.intVal, newVal.intVal) != prevVal.intVal);
}

// check overlap with other atoms
int CheckSphereOverLap(float4 posA, float *smallest_radius, int numberOfFrameworkAtoms, __global float4 *frameworkPositions, 
                       float4 inverse_cella, float4 inverse_cellb, float4 inverse_cellc,
                       float4 cella, float4 cellb, float4 cellc, __global float *sigma)
{
  float4 ds;
  *smallest_radius = 1000.0;
  for(int jatom = 0; jatom < numberOfFrameworkAtoms; ++jatom)
  {
    float4 dr = posA - frameworkPositions[jatom];

    ds.x = dot(inverse_cella, dr);
    ds.y = dot(inverse_cellb, dr);
    ds.z = dot(inverse_cellc, dr);
    ds.w = 0.0f;

    float4 t = ds - rint(ds);

    dr.x = dot(cella, t);
    dr.y = dot(cellb, t);
    dr.z = dot(cellc, t);
    dr.w = 0.0f;

    float rr = dot(dr, dr);

    if(rr < (0.5 * sigma[jatom]) * (0.5 * sigma[jatom]))
    {
      return true;
    }

    float radius = sqrt(rr) - 0.5 * sigma[jatom];
    if(radius < *smallest_radius)
    {
      *smallest_radius = radius;
    }
  }
  return false;
}

__kernel void ComputePoreSizeDistribution(__global float4 *frameworkPositions,
                                          __global float4 *randomCartesianPositions,
                                          __global float *sigma,
                                          __global float *largest_radius,
                                          const float4 posA,
                                          const int numberOfFrameworkAtoms,
                                          const int numberOfInnerSteps,
                                          const float4 cella,
                                          const float4 cellb,
                                          const float4 cellc,
                                          const float4 inverse_cella,
                                          const float4 inverse_cellb,
                                          const float4 inverse_cellc)
{
  int global_id = get_global_id(0);
  float4 ds;

  // check overlap with other atoms
  float radius = 0.0;
  if(CheckSphereOverLap(posA, &radius, numberOfFrameworkAtoms, frameworkPositions, 
                     inverse_cella, inverse_cellb, inverse_cellc,
                     cella, cellb, cellc, sigma)) return;

  float4 posB = randomCartesianPositions[global_id];

  if(CheckSphereOverLap(posB, &radius, numberOfFrameworkAtoms, frameworkPositions, 
                     inverse_cella, inverse_cellb, inverse_cellc,
                     cella, cellb, cellc, sigma)) return;

  float4 dr = posA - posB;

  ds.x = dot(inverse_cella, dr);
  ds.y = dot(inverse_cellb, dr);
  ds.z = dot(inverse_cellc, dr);
  ds.w = 0.0f;

  float4 t = ds - rint(ds);

  dr.x = dot(cella, t);
  dr.y = dot(cellb, t);
  dr.z = dot(cellc, t);
  dr.w = 0.0f;

  float rr = dot(dr, dr);

  if(rr > radius * radius) return;

  AtomicMax(largest_radius, radius);
}
)foo";


