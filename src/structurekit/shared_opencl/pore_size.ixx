module;

export module grid_pore_size_opencl;

import std;

import uint3;
import unit_cell;

// The Gelb-Gubbins sphere sweep on the GPU: a sphere about every point of the void, each point keeping the
// largest sphere that reached it.
//
// This is the one part of a pore-size distribution that wants a GPU for its own sake rather than for the
// field beneath it. It is a sphere's worth of grid points per grid point of the void, which is far more
// arithmetic than filling the field in was, and it parallelises perfectly: the spheres do not interact, and
// a point being told a size is a single atomic maximum.
//
// It reads a distance field and knows nothing of where the distance came from, so the clearance route and
// the energy route run the same sweep and can only differ in the field they measured.
export std::vector<float> poreRadiusFieldOpenCL(uint3 gridSize, const UnitCell &unitCell,
                                                std::span<const float> distance, double slack);

export extern const char *poreSizeKernelSource;
