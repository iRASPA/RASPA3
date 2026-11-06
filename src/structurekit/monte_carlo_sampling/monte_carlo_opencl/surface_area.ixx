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
// atoms of a pass see the same ones, exactly as they do there. Each pass is one reading of the total; the
// sum and sum of squares of those readings give the mean and a 95% confidence interval.
export struct MC_OpenCL_SurfaceArea
{
  // Made once, when the route is, and let go of when the route goes. They are nothing where there is no
  // device, and a route is not copied, there being one program and one kernel to let go of once.
  cl_program surfaceAreaProgram{nullptr};
  cl_kernel surfaceAreaKernel{nullptr};
  static const char* surfaceAreaKernelSource;
  size_t surfaceAreaWorkGroupSize{0};

  double surfaceArea{0.0};       // Å², averaged over the passes
  double surfaceAreaError{0.0};  // Å², half-width of the 95% confidence interval of the mean
  double seconds{0.0};

  MC_OpenCL_SurfaceArea();
  ~MC_OpenCL_SurfaceArea();
  MC_OpenCL_SurfaceArea(const MC_OpenCL_SurfaceArea&) = delete;
  MC_OpenCL_SurfaceArea& operator=(const MC_OpenCL_SurfaceArea&) = delete;

  void run(const SampledStructure& structure, const SampledProbe& probe,
           std::optional<std::size_t> numberOfIterations, std::optional<std::size_t> numberOfInnerSteps);
};
