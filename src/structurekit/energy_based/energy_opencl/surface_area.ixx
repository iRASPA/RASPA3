module;

#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

export module energy_opencl_surface_area;

import std;

import int3;
import uint3;
import double2;
import double3;
import double4;
import float4;
import double3x3;
import pair_interactions;
import crystal;

import energy_shared_isosurface;

export struct EnergyOpenCLSurfaceArea
{
  EnergyOpenCLSurfaceArea();
  ~EnergyOpenCLSurfaceArea();

  cl_program energyGridProgram;
  cl_kernel energyGridKernel;
  static const char *energyGridKernelSource;
  size_t energyGridWorkGroupSize;

  cl_program energyEnergyOpenCLSurfaceAreaProgram;
  cl_kernel constructHPLevelKernel;
  cl_kernel classifyCubesKernel;
  cl_kernel traverseHPKernel[10];
  size_t constructHPLevelKernelWorkGroupSize;
  size_t classifyCubesKernelWorkGroupSize;
  size_t traverseHPKernelWorkGroupSize[10];
  static std::string marchingCubesKernelSource;

  // Extracts the surface where an arbitrary field crosses 'isoValue'. Which field is handed in decides what
  // the area means: the energy field of a single probe atom gives the surface that atom sees, a molecular
  // free-energy field gives the surface a whole molecule sees once its orientations are averaged over.
  // The triangles themselves, three corners to a triangle, in fractional coordinates. Handing these back is
  // what lets a surface be divided among the atoms rather than only measured.
  std::vector<double3> trianglesOfIsosurface(std::span<const float> field, uint3 gridSize, double isoValue);

  IsosurfaceArea areaOfIsosurface(const Crystal &framework, std::span<const float> field, uint3 grid_size,
                                  double isoValue);

  void run(const PairInteractions &interactions, const Crystal &framework, double isoValue,
                         std::string probePseudoAtom, uint3 grid_size);
};
