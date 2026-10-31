module;

#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

module mc_opencl_void_fraction;

import std;

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

namespace
{
// The build log of a program for a device, sized by a first query (the log of a failed build regularly exceeds a
// fixed buffer, and a query with a too-small buffer fails without writing anything).
std::string programBuildLog(cl_program program, cl_device_id device)
{
  size_t length = 0;
  if (clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, nullptr, &length) != CL_SUCCESS || length == 0)
  {
    return std::string("(no build log)");
  }
  std::string log(length, '\0');
  if (clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, length, log.data(), nullptr) != CL_SUCCESS)
  {
    return std::string("(no build log)");
  }
  while (!log.empty() && log.back() == '\0') log.pop_back();
  return log;
}
}  // namespace

MC_OpenCL_VoidFraction::MC_OpenCL_VoidFraction()
{
  // idempotent: selects the device and creates the shared context and queue when not done yet
  OpenCL::initialize();

  if (OpenCL::clContext.has_value() && OpenCL::clDeviceId.has_value())
  {
    cl_int err;

    const char* voidFractionShaderSourceCode = MC_OpenCL_VoidFraction::voidFractionKernelSource;
    voidFractionProgram =
        clCreateProgramWithSource(OpenCL::clContext.value(), 1, &voidFractionShaderSourceCode, nullptr, &err);
    if (err != CL_SUCCESS)
    {
      voidFractionProgram = nullptr;
      throw std::runtime_error(
          std::format("MC_OpenCL_VoidFraction: OpenCL clCreateProgramWithSource failed at {}\n", __LINE__));
    }

    err = clBuildProgram(voidFractionProgram, 0, nullptr, nullptr, nullptr, nullptr);
    if (err != CL_SUCCESS)
    {
      const std::string log = programBuildLog(voidFractionProgram, OpenCL::clDeviceId.value());
      clReleaseProgram(voidFractionProgram);
      voidFractionProgram = nullptr;
      throw std::runtime_error(
          std::format("MC_OpenCL_VoidFraction: OpenCL Failed to build program at {} (line {} error: {})\n", __FILE__,
                      __LINE__, log));
    }

    voidFractionKernel = clCreateKernel(voidFractionProgram, "ComputeVoidFraction", &err);
    if (err != CL_SUCCESS)
    {
      voidFractionKernel = nullptr;
      clReleaseProgram(voidFractionProgram);
      voidFractionProgram = nullptr;
      throw std::runtime_error(
          std::format("MC_OpenCL_VoidFraction: OpenCL clCreateKernel failed {} : {}\n", __FILE__, __LINE__));
    }
    err = clGetKernelWorkGroupInfo(voidFractionKernel, OpenCL::clDeviceId.value(), CL_KERNEL_WORK_GROUP_SIZE,
                                   sizeof(size_t), &voidFractionWorkGroupSize, nullptr);
    if (err != CL_SUCCESS)
    {
      clReleaseKernel(voidFractionKernel);
      voidFractionKernel = nullptr;
      clReleaseProgram(voidFractionProgram);
      voidFractionProgram = nullptr;
      throw std::runtime_error(std::format(
          "MC_OpenCL_VoidFraction: OpenCL clGetKernelWorkGroupInfo failed at {} : {}\n", __FILE__, __LINE__));
    }
  }
}

MC_OpenCL_VoidFraction::~MC_OpenCL_VoidFraction()
{
  if (voidFractionKernel != nullptr) clReleaseKernel(voidFractionKernel);
  if (voidFractionProgram != nullptr) clReleaseProgram(voidFractionProgram);
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

