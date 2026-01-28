module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <print>
#include <string>
#include <vector>
#endif

#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

export module mc_opencl_pore_size_distribution;

#ifdef USE_STD_IMPORT
import std;
#endif

import framework;
import forcefield;

export struct MC_OpenCL_PoreSizeDistribution
{
  cl_program poreSizeDistributionProgram;
  cl_kernel poreSizeDistributionKernel;
  static const char* poreSizeDistributionKernelSource;
  size_t poreSizeDistributionWorkGroupSize;

  std::vector<double> data;

  MC_OpenCL_PoreSizeDistribution(std::size_t numberOfBins);

  void run(const ForceField &forceField, const Framework &framework, double well_depth_factor,
           std::size_t number_of_iterations);
};
