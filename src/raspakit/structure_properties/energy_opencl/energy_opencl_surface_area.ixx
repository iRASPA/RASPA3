module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <array>
#include <vector>
#include <optional>
#include <string>
#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
  #include <OpenCL/cl.h>
#elif _WIN32
  #include <CL/cl.h>
#else
  #include <CL/opencl.h>
#endif
#endif

export module energy_opencl_surface_area;

import int3;
import double2;
import double3;
import double4;
import float4;
import double3x3;
import forcefield;
import framework;


export struct EnergyOpenCLSurfaceArea
{
  EnergyOpenCLSurfaceArea();
  ~EnergyOpenCLSurfaceArea();

  cl_program energyGridProgram;
  cl_kernel energyGridKernel;
  static const char* energyGridKernelSource;
  size_t energyGridWorkGroupSize;

  cl_program energyEnergyOpenCLSurfaceAreaProgram;
  cl_kernel constructHPLevelKernel;
  cl_kernel classifyCubesKernel;
  cl_kernel traverseHPKernel[10];
  size_t constructHPLevelKernelWorkGroupSize;
  size_t classifyCubesKernelWorkGroupSize;
  size_t traverseHPKernelWorkGroupSize[10];
  static std::string marchingCubesKernelSource;

  void run(const ForceField &forceField, const Framework &framework, int3 grid_size);
};
