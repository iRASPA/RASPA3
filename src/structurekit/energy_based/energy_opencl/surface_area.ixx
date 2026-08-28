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
  //
  // When `gradients` is given it is filled with the field's gradient at each of those same corners, held per
  // grid step and pointing towards larger values of the field, which on an energy field is into the wall. That
  // is the same sense the processor extractor uses, so the two are interchangeable and a consumer of either can
  // be written once; `FieldSense` on the far side is what turns a gradient into an outward normal.
  //
  // Nothing has to be recomputed or fetched to obtain it: the kernel already writes a gradient beside every
  // vertex it emits, and the host already reads the whole buffer back and was discarding those slots.
  //
  // The length is not comparable between the two extractors --- this one leaves the difference unscaled and the
  // processor one normalises --- and nothing may depend on it.
  std::vector<double3> trianglesOfIsosurface(std::span<const float> field, uint3 gridSize, double isoValue,
                                             std::vector<double3> *gradients = nullptr);

  IsosurfaceArea areaOfIsosurface(const Crystal &framework, std::span<const float> field, uint3 grid_size,
                                  double isoValue);

  void run(const PairInteractions &interactions, const Crystal &framework, double isoValue,
                         std::string probePseudoAtom, uint3 grid_size);
};
