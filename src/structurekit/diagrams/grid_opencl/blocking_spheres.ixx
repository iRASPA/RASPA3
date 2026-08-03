module;

export module opencl_blocking_spheres;

import std;

import int3;
import uint3;
import double3;
import crystal;
import pair_interactions;
import opencl_clearance_grid;
import grid_connected_components;
import grid_blocking_cover;

// Spheres covering the pockets a probe cannot reach, in the format a simulation reads back.
//
// The pockets are the pieces of the region where the probe fits that do not run anywhere: exactly what the
// connected components of the clearance field at the probe's radius report as having dimensionality zero. So
// nothing has to be decided here about what is sealed off and what is not, and the grid is doing what a grid
// is good at, which is knowing the shape of a region rather than guessing at it from samples.
//
// A sphere here is compared against the centres of the atoms a simulation tries to insert, so what it covers is
// the region the probe's centre can occupy in the pocket and not the pocket's whole room, which is larger by
// the probe's radius all round.
//
// The covering itself is the same for any route that can say which regions are sealed and where the probe is
// entitled to be, so it lives in `grid_blocking_cover` and is shared with the energy route, which decides
// what is sealed by a barrier against kT rather than by connectivity.
export struct GridBlockingSpheres
{
  double probeRadius{0.0};

  std::vector<GridBlockingSphere> spheres;

  std::size_t numberOfPockets{0};
  std::size_t numberOfSinglePointPockets{0};
  std::size_t numberOfPocketVoxels{0};
  std::size_t numberOfClippedSpheres{0};

  // Points of a pocket that lie nearer a channel than the grid can resolve, where no sphere can be placed
  // without taking channel with it. These are necks the grid has closed rather than pockets, and they hold no
  // volume worth blocking.
  std::size_t numberOfRefusedPoints{0};

  // The share of the cell the pockets hold, and the share the spheres cover, the second being the larger
  // because a sphere reaches to the pocket's wall and beyond it into the framework.
  double pocketFraction{0.0};

  double gridSeconds{0.0};
  double seconds{0.0};

  GridBlockingSpheres();
  ~GridBlockingSpheres();

  void run(const PairInteractions &interactions, const Crystal &framework, std::string probePseudoAtom, uint3 gridSize);
  void run(const PairInteractions &interactions, const Crystal &framework, std::string probePseudoAtom,
           const ClearanceGrid &grid);
};
