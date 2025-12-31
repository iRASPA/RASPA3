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

export module mc_opencl_surface_area;

#ifdef USE_STD_IMPORT
import std;
#endif

import framework;
import forcefield;

export struct MC_OpenCL_SurfaceArea
{
  cl_program surfaceAreaProgram;
  cl_kernel surfaceAreaKernel;
  static const char* surfaceAreaKernelSource;
  size_t surfaceAreaWorkGroupSize;

  MC_OpenCL_SurfaceArea();

  void run(const ForceField &forceField, const Framework &framework, double wellDepthFactor,
           std::string probePseudoAtom, std::size_t numberOfSlices) const;
};
