module;

export module opencl_pore_analysis;

import std;

import int3;
import uint3;
import crystal;
import pair_interactions;
import opencl_clearance_grid;
import grid_connected_components;

// The three pore diameters, read off the clearance field.
//
// Di is the widest point of the void, so it is the largest value the field takes. Df and Dif come from a
// single sweep down through the field: the points are switched on from the widest downwards and joined to
// the neighbours already on, and the moment some pore closes on a periodic image of itself the sweep has
// found the widest probe that still gets through, which is Df. Dif is then the widest point of that pore,
// because a probe of radius Df/2 that travels along it can reach every point of it.
//
// The sweep costs one pass instead of a search over probe sizes with a fresh connectivity pass each time, and
// it yields the probe size at which the pore system stops running in three, two, and one directions as a
// by-product, since those are just later moments of the same sweep.
export struct GridPoreDiameters
{
  double includedSphereDiameter{0.0};
  double freeSphereDiameter{0.0};
  double includedAlongFreePathDiameter{0.0};
  bool percolates{false};

  // The dimensionality of the pore that first percolated, at the moment it did.
  int dimensionalityAtThreshold{0};

  // The largest free sphere that still runs in at least one, two, and three directions. Zero where the pore
  // system never runs that far.
  std::array<double, 3> freeSphereDiameterByDimension{0.0, 0.0, 0.0};

  std::size_t numberOfVoidVoxels{0};
  double seconds{0.0};

  static GridPoreDiameters compute(const ClearanceGrid &grid);
};


export struct GridPoreAnalysis
{
  GridPoreDiameters diameters;
  GridComponents channels;
  double probeRadius{0.0};

  GridPoreAnalysis();
  ~GridPoreAnalysis();

  void run(const PairInteractions &interactions, const Crystal &framework, std::string probePseudoAtom, uint3 gridSize);
  void run(const PairInteractions &interactions, const Crystal &framework, std::string probePseudoAtom,
           const ClearanceGrid &grid);
};
