module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <array>
#include <vector>
#include <optional>
#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
  #include <OpenCL/cl.h>
#elif _WIN32
  #include <CL/cl.h>
#else
  #include <CL/opencl.h>
#endif
#endif

export module isosurface;

import int3;
import double2;
import double3;
import double4;
import float4;
import double3x3;


export struct Isosurface
{
  Isosurface();
  ~Isosurface();

  cl_program program;
  cl_kernel constructHPLevelKernel;
  cl_kernel classifyCubesKernel;
  cl_kernel traverseHPKernel[10];
  size_t constructHPLevelKernelWorkGroupSize;
  size_t classifyCubesKernelWorkGroupSize;
  size_t traverseHPKernelWorkGroupSize[10];
  static std::string marchingCubesKernel;

  std::vector<float4> computeIsosurface(int3 dimensions, std::vector<cl_float>* voxels, double isoValue) noexcept(false);
};
