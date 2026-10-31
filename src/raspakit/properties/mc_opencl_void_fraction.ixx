module;

#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

export module mc_opencl_void_fraction;

import std;

import framework;
import forcefield;

export struct MC_OpenCL_VoidFraction
{
  cl_program voidFractionProgram{nullptr};
  cl_kernel voidFractionKernel{nullptr};
  static const char* voidFractionKernelSource;
  size_t voidFractionWorkGroupSize{0};

  std::vector<double> data;

  MC_OpenCL_VoidFraction();
  ~MC_OpenCL_VoidFraction();

  // the kernel and program handles are owned (released once in the destructor)
  MC_OpenCL_VoidFraction(const MC_OpenCL_VoidFraction &) = delete;
  MC_OpenCL_VoidFraction &operator=(const MC_OpenCL_VoidFraction &) = delete;

  void run(const ForceField &forceField, const Framework &framework, std::size_t number_of_iterations);
};
