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

export module energy_void_fraction;

import int3;
import double2;
import double3;
import double3x3;


export struct EnergyVoidFraction
{
  EnergyVoidFraction();
  ~EnergyVoidFraction();

  cl_program program;
  cl_kernel kernel;
  size_t workGroupSize;
  static const char* energyVoidFractionKernel;

  double computeVoidFraction(std::vector<cl_float> *voxels);
};
