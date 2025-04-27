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

export module energy_grid;

import int3;
import double2;
import double3;
import double3x3;


export struct EnergyGrid
{
  enum class GridSizeType: size_t
  {
    custom = 0,
    size2x2x2 = 1,
    size4x4x4 = 2,
    size8x8x8 = 3,
    size16x16x16 = 4,
    size32x32x32 = 5,
    size64x64x64 = 6,
    size128x128x128 = 7,
    size256x256x256 = 8,
    size512x512x512 = 9,
    multiple_values = 10
  };

  EnergyGrid();
  ~EnergyGrid();

  cl_program program;
  cl_kernel kernel;
  size_t workGroupSize;
  static const char* energyGridKernel;

  std::vector<cl_float> computeEnergyGrid(int3 size, double2 probeParameter,
                                   std::vector<double3> positions, std::vector<double2> potentialParameters,
                                   double3x3 unitCell, int3 numberOfReplicas) noexcept(false);
  std::vector<cl_float> computeEnergyGridCPUImplementation(int3 size, double2 probeParameter,
                                                    std::vector<double3> positions, std::vector<double2> potentialParameters,
                                                    double3x3 unitCell, int3 numberOfReplicas) noexcept;

  std::vector<cl_float> computeEnergyGridGPUImplementation(int3 size, double2 probeParameter,
                                             std::vector<double3> positions, std::vector<double2> potentialParameters,
                                             double3x3 unitCell, int3 numberOfReplicas) noexcept(false);
};
