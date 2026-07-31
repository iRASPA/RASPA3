module;

export module opencl_pore_size_distribution;

import std;

import uint3;
import double3;
import framework;
import forcefield;
import opencl_clearance_grid;
import opencl_connected_components;

export struct GridPoreSizeDistributionPoint
{
  double diameter{0.0};

  // The share of the void lying in pores at least this wide, and the derivative of it, which is the
  // distribution itself.
  double cumulativeVoidFraction{0.0};
  double distribution{0.0};

  // The same two, over the part of the void a probe can reach.
  double cumulativeAccessibleFraction{0.0};
  double accessibleDistribution{0.0};
};


// The pore-size distribution in the sense of Gelb and Gubbins: the size at a point in the void is the diameter
// of the largest sphere that fits in the void and covers that point, and the distribution is how that size is
// spread over the void.
//
// Written out, the size at a point x is max{2 d(c) : |x - c| <= d(c)} over the centres c: a sphere is drawn
// about every point of the void, and every point a sphere covers is told how big it was and keeps the largest
// it hears. That is a great deal of work, one sphere's worth of grid points per grid point of the void, and it
// is the reason this is the one part of the route that wants a GPU for its own sake rather than for the
// clearance field: the spheres do not interact, and a point being told a size is one atomic maximum.
//
// It cannot be cut down by dropping spheres that sit inside larger ones, tempting as that looks. The field is a
// distance, so it grows by exactly one along its own gradient and by less in every other direction, and one
// ball therefore contains another only when the two centres and the nearest atom fall on a single line. Off
// that line no sphere is redundant, and on a grid the line is never hit exactly, so the test that would drop a
// sphere almost never fires.
//
// A sphere is counted as covering a grid point when it reaches the point's own voxel rather than the point
// itself. Without that, the shell of void lying within a voxel of an atom's surface reports the size of
// nothing at all: the largest sphere covering such a point touches it from inside, and a sphere centred on a
// grid point misses that tangency by the distance from the grid to where its centre should have been.
export struct GridPoreSizeDistribution
{
  double probeRadius{0.0};

  std::vector<GridPoreSizeDistributionPoint> points;

  // The largest sphere that fits in the void, which is the same as Di.
  double largestDiameter{0.0};

  double voidFraction{0.0};
  double accessibleVoidFraction{0.0};

  // How far a sphere may reach past its radius and still count as covering a grid point, which is what lets the
  // shell of void nearest the atoms be given a size at all.
  double slack{0.0};

  std::size_t numberOfCentres{0};
  std::size_t numberOfVoidVoxels{0};
  double gpuSeconds{0.0};
  double seconds{0.0};

  GridPoreSizeDistribution();
  ~GridPoreSizeDistribution();

  void run(const ForceField &forceField, const Framework &framework, std::string probePseudoAtom, uint3 gridSize,
           std::optional<double> maximumDiameter, std::optional<std::size_t> numberOfBins);
  void run(const ForceField &forceField, const Framework &framework, std::string probePseudoAtom,
           const ClearanceGrid &grid, std::optional<double> maximumDiameter, std::optional<std::size_t> numberOfBins);

  static const char *poreSizeKernelSource;
};
