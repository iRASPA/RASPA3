module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <print>
#include <string>
#include <vector>
#include <optional>
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

  std::size_t numberOfBins;
  std::vector<double> histogram;
  std::vector<double> histogram_cummulative;

  MC_OpenCL_PoreSizeDistribution(std::size_t numberOfBins);
  ~MC_OpenCL_PoreSizeDistribution();

  void run(const ForceField &forceField, const Framework &framework, double wellDepthFactor,
           std::optional<std::size_t> numberOfIterations, std::optional<std::size_t> numberOfInnerSteps,
           std::optional<double> maximumRange);
};
