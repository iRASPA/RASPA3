module;

#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

export module integration_opencl_surface_area;

import std;

import sampled_structure;

// The gridded accessible surface area on the GPU: one work item per atom, each walking the same latitude
// and longitude grid over its own contact sphere. See `integration_surface_area` for what the grid is and
// what it costs in accuracy; the whole sweep is one kernel launch, there being nothing random to redraw.
export struct Integration_OpenCL_SurfaceArea
{
  // Made once, when the route is, and let go of when the route goes. They are nothing where there is no
  // device, and a route is not copied, there being one program and one kernel to let go of once.
  cl_program surfaceAreaProgram{nullptr};
  cl_kernel surfaceAreaKernel{nullptr};
  static const char* surfaceAreaKernelSource;
  size_t surfaceAreaWorkGroupSize{0};

  double surfaceArea{0.0};  // Å²
  double seconds{0.0};

  Integration_OpenCL_SurfaceArea();
  ~Integration_OpenCL_SurfaceArea();
  Integration_OpenCL_SurfaceArea(const Integration_OpenCL_SurfaceArea&) = delete;
  Integration_OpenCL_SurfaceArea& operator=(const Integration_OpenCL_SurfaceArea&) = delete;

  void run(const SampledStructure& structure, const SampledProbe& probe, std::optional<std::size_t> numberOfSlices);
};
