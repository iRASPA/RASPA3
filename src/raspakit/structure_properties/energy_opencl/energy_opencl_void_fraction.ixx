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

export module energy_opencl_void_fraction;

import int3;
import double2;
import double3;
import double3x3;
import forcefield;
import framework;

export struct EnergyOpenCLVoidFraction
{
  EnergyOpenCLVoidFraction();
  ~EnergyOpenCLVoidFraction();

  cl_program energyGridProgram;
  cl_kernel energyGridKernel;
  static const char* energyGridKernelSource;
  size_t energyGridWorkGroupSize;

  cl_program energyVoidFractionProgram;
  cl_kernel energyVoidFractionKernel;
  static const char* energyVoidFractionKernelSource;
  size_t energyVoidFractionWorkGroupSize;

  void run(const ForceField &forceField, const Framework &framework);
};
