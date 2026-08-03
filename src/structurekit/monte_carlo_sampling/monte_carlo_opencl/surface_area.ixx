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

import sampled_structure;

// The Monte Carlo accessible surface area on the GPU: one work item per atom, each sweeping the same set of
// random directions over its own contact sphere and counting the ones that reach no other atom. See
// `mc_surface_area` for what is being estimated; the directions are drawn on the processor so that all the
// atoms of a pass see the same ones, exactly as they do there.
export struct MC_OpenCL_SurfaceArea
{
  cl_program surfaceAreaProgram;
  cl_kernel surfaceAreaKernel;
  static const char* surfaceAreaKernelSource;
  size_t surfaceAreaWorkGroupSize;

  double surfaceArea{0.0};  // Å², averaged over the passes
  double seconds{0.0};

  MC_OpenCL_SurfaceArea();

  void run(const SampledStructure& structure, const SampledProbe& probe,
           std::optional<std::size_t> numberOfIterations, std::optional<std::size_t> numberOfInnerSteps);
};
