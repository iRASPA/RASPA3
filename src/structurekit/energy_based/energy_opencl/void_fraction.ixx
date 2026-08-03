module;

#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

export module energy_opencl_void_fraction;

import std;

import int3;
import uint3;
import double2;
import double3;
import double3x3;
import pair_interactions;
import crystal;

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

  // The void fraction the probe named sees, at the temperature given. Both matter: a Boltzmann average is
  // meaningless without saying how hot, and a probe is what the average is over.
  void run(const PairInteractions& interactions, const Crystal& framework, std::string probePseudoAtom, uint3 grid_size,
           double temperature = 298.0);
};
