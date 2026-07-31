module;

export module opencl_blocking_spheres;

import std;

import int3;
import uint3;
import double3;
import framework;
import forcefield;
import opencl_clearance_grid;
import opencl_connected_components;

export struct GridBlockingSphere
{
  // What a simulation reads: where the sphere sits in the cell, and how far it reaches.
  double3 centreFractional{0.0, 0.0, 0.0};
  double radius{0.0};  // Å

  // Which pocket it was drawn for, which grid point it was centred on, and how much of the pocket it took.
  std::size_t pocket{0};
  int3 centreVoxel{0, 0, 0};
  std::size_t voxelsCovered{0};

  // Whether a channel, rather than the pocket's own extent, is what stopped it growing.
  bool clipped{false};

  // Whether there was no room to place it at all, the point being nearer a channel than the grid can resolve.
  // Such a sphere is counted but not written.
  bool refused{false};
};

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
// A sphere is placed at the middle of the pocket and grown until it reaches the pocket's furthest point, which
// is one sphere for a pocket of any ordinary shape and is where the Voronoi and Apollonius routes place theirs
// too. Growing that far takes in framework, and other pockets, and void too narrow for the probe, and none of
// that costs anything, there being nowhere in any of it that a centre was going to be allowed.
//
// What does cost something is reaching a channel. A sphere over that line blocks pore volume the probe is
// entitled to, and a simulation given such a file reports an uptake too low for a reason nothing in its output
// would explain. So a sphere is cut back at the nearest point of the region where the probe both fits and can
// get out, which the grid knows rather than estimates, and a pocket that a cut leaves partly uncovered is given
// further spheres until none of it is left.
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

  void run(const ForceField &forceField, const Framework &framework, std::string probePseudoAtom, uint3 gridSize);
  void run(const ForceField &forceField, const Framework &framework, std::string probePseudoAtom,
           const ClearanceGrid &grid);
};
