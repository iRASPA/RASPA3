module;


#define CL_TARGET_OPENCL_VERSION 120
#define CL_SILENCE_DEPRECATION
#ifdef __APPLE__
  #include <OpenCL/cl.h>
#elif _WIN32
  #include <CL/cl.h>
#else
  #include <CL/opencl.h>
#endif

export module pore_size_distribution_ban_vlugt;

import std;


import int3;
import uint3;
import framework;
import forcefield;

// The pore size of Ban and Vlugt on a grid over the unit cell: at each grid point, the radius of the largest
// sphere that holds the point and no framework atom.
//
// It is found in two steps, as the method has it. A sphere is grown from every grid point out to the surface of
// the nearest atom, an atom being a hard sphere of half its Lennard-Jones size. That sphere is then laid over
// every grid point inside it, and a point keeps the widest sphere laid over it. A point inside an atom is in no
// sphere and keeps zero.
//
// The cell is periodic and the grid is of the cell, so point i of n along an axis stands at the fraction i / n
// of it and the point after the last is the first again: a sphere near a face is laid on round the cell, and one
// wider than half the cell meets copies of itself. The cell may be of any shape. Distances are taken in
// Ångström through the cell matrix, and how many grid points a sphere can reach along an axis is set by the
// distance between the faces of the cell across that axis, which in a skewed cell is less than the cell vector
// is long.
export struct BanVlugtPoreSizeDistribution
{
  uint3 grid_size;

  // The pore size at each grid point, in Ångström, with x varying fastest. Filled by `run`.
  std::vector<float> poreSizes;

  cl_program PSDProgram;
  cl_kernel PSDKernel;
  static const char* PSDKernelSource;
  size_t PSDWorkGroupSize;

  BanVlugtPoreSizeDistribution(uint3 grid_size);
  ~BanVlugtPoreSizeDistribution();

  void run(const ForceField &forceField, const Framework &framework);

};

