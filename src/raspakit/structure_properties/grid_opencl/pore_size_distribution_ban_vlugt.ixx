module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <stdint.h>
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


export module pore_size_distribution_ban_vlugt;

#ifndef USE_LEGACY_HEADERS
import std;
#endif


import int3;
import framework;
import forcefield;

export struct BanVlugtPoreSizeDistribution
{
  int3 grid_size;

  cl_program PSDProgram;
  cl_kernel PSDKernel;
  static const char* PSDKernelSource;
  size_t PSDWorkGroupSize;

  BanVlugtPoreSizeDistribution(int3 grid_size);
  ~BanVlugtPoreSizeDistribution();

  void run(const ForceField &forceField, const Framework &framework);

};

