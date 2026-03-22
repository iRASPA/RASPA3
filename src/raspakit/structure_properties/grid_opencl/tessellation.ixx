module;

#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

export module tessellation;

import std;

import int3;
import uint3;
import framework;
import forcefield;

export struct Tessellation
{
  uint3 grid_size;

  cl_program tessellationProgram;
  cl_kernel tessellationKernel;
  static const char *tessellationKernelSource;
  std::size_t tessellationWorkGroupSize;

  Tessellation(uint3 grid_size);
  ~Tessellation();

  void run(const ForceField &forceField, const Framework &framework);
};
