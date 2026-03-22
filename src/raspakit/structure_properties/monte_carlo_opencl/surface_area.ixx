module;

#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

export module mc_opencl_surface_area;

import std;

import framework;
import forcefield;

export struct MC_OpenCL_SurfaceArea
{
  cl_program surfaceAreaProgram;
  cl_kernel surfaceAreaKernel;
  static const char* surfaceAreaKernelSource;
  size_t surfaceAreaWorkGroupSize;

  MC_OpenCL_SurfaceArea();

  void run(const ForceField &forceField, const Framework &framework, double wellDepthFactor, std::string probePseudoAtom,
           std::optional<std::size_t> numberOfIterations, std::optional<std::size_t> numberOfInnerSteps) const;
};
