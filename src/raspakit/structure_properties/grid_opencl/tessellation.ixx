module;

#ifdef USE_LEGACY_HEADERS
#include <stdint.h>

#include <cstddef>
#include <print>
#include <string>
#include <vector>
#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif
#endif

#ifndef USE_LEGACY_HEADERS
#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif
#endif

export module tessellation;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import int3;
import framework;
import forcefield;

export struct Tessellation
{
  int3 grid_size;

  cl_program tessellationProgram;
  cl_kernel tessellationKernel;
  static const char *tessellationKernelSource;
  std::size_t tessellationWorkGroupSize;

  Tessellation(int3 grid_size);
  ~Tessellation();

  void run(const ForceField &forceField, const Framework &framework);
};
