module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <stdint.h>
#include <print>
#include <string>
#include <vector>
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

#ifdef USE_STD_IMPORT
#define CL_TARGET_OPENCL_VERSION 120
#define CL_SILENCE_DEPRECATION
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif
#endif


export module pore_size_distribution_ban_vlugt;

#ifdef USE_STD_IMPORT
import std;
#endif


import int3;
import uint3;
import framework;
import forcefield;

export struct BanVlugtPoreSizeDistribution
{
  uint3 grid_size;

  cl_program PSDProgram;
  cl_kernel PSDKernel;
  static const char* PSDKernelSource;
  size_t PSDWorkGroupSize;

  BanVlugtPoreSizeDistribution(uint3 grid_size);
  ~BanVlugtPoreSizeDistribution();

  void run(const ForceField &forceField, const Framework &framework);

};

