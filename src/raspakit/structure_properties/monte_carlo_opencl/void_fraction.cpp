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
#include <tuple>
#include <vector>
#include <chrono>
#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif
#endif

module mc_opencl_void_fraction;

#ifdef USE_STD_IMPORT
import std;
#endif

import opencl;
import double2;
import double3;
import float4;
import double3x3;
import randomnumbers;
import atom;
import framework;
import forcefield;
import component;
import system;
import mc_moves_widom;

MC_OpenCL_VoidFraction::MC_OpenCL_VoidFraction()
{
  if (OpenCL::clContext.has_value() && OpenCL::clDeviceId.has_value())
  {
    cl_int err;

    const char* voidFractionShaderSourceCode = MC_OpenCL_VoidFraction::voidFractionKernelSource;
    voidFractionProgram =
        clCreateProgramWithSource(OpenCL::clContext.value(), 1, &voidFractionShaderSourceCode, nullptr, &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("MC_OpenCL_VoidFraction: OpenCL clCreateProgramWithSource failed at {}\n", __LINE__));
    }

    err = clBuildProgram(voidFractionProgram, 0, nullptr, nullptr, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      size_t len;
      char buffer[2048];
      clGetProgramBuildInfo(voidFractionProgram, OpenCL::clDeviceId.value(), CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer,
                            &len);
      std::string message =
          std::format("MC_OpenCL_VoidFraction: OpenCL Failed to build program at {} (line {} error: {})\n", __FILE__,
                      __LINE__, std::string(buffer));
      throw std::runtime_error(message);
    }

    voidFractionKernel = clCreateKernel(voidFractionProgram, "ComputeVoidFraction", &err);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("MC_OpenCL_VoidFraction: OpenCL clCreateKernel failed {} : {}\n", __FILE__, __LINE__));
    }
    err = clGetKernelWorkGroupInfo(voidFractionKernel, OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE,
                                   sizeof(size_t), &voidFractionWorkGroupSize, nullptr);
    if (err != CL_SUCCESS)
    {
      throw std::runtime_error(std::format("MC_OpenCL_VoidFraction: OpenCL clGetKernelWorkGroupInfo failed at {} : {}\n", __FILE__, __LINE__));
    }
  }

}

void MC_OpenCL_VoidFraction::run([[maybe_unused]] const ForceField &forceField,
                                 [[maybe_unused]] const Framework &framework,
                                 [[maybe_unused]] std::size_t number_of_iterations)
{
  RandomNumber random{std::nullopt};
}

const char* MC_OpenCL_VoidFraction::voidFractionKernelSource = R"foo(
__kernel void ComputeVoidFraction(__global float4 *position,
                                  __global float *sigma,
                                  __global float *output,
                                  const int numberOfAtoms,
                                  const int numberOfSlices,
                                  const float4 cella,
                                  const float4 cellb,
                                  const float4 cellc,
                                  const float4 inverse_cella,
                                  const float4 inverse_cellb,
                                  const float4 inverse_cellc)
{
}
)foo";

