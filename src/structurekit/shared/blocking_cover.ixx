module;

export module grid_blocking_cover;

import std;

import int3;
import uint3;
import double3;
import unit_cell;
import grid_connected_components;

// Covering a set of pockets with spheres, in the format a simulation reads back.
//
// What counts as a pocket is not decided here, and the two routes that use this do not decide it the same
// way. The geometric one calls a pocket a piece of the void that does not run anywhere, so that a probe put
// into it can never leave; the energy one calls a pocket a region whose escape barrier is large against kT,
// so that a molecule put into it will not leave in any time worth waiting. Both hand the answer down as a
// list of regions to be covered and a mask of the points a sphere must keep out of, and from there the
// covering is the same arithmetic.
//
// A sphere here is compared against the centres of the atoms a simulation tries to insert, so what it covers
// is the region a centre could occupy and not the pocket's whole room, which is larger by the probe all
// round.
//
// A sphere is placed at the middle of the pocket and grown until it reaches the pocket's furthest point,
// which is one sphere for a pocket of any ordinary shape. Growing that far takes in framework, and other
// pockets, and void too narrow to be entered, and none of that costs anything, there being nowhere in any of
// it that a centre was going to be allowed.
//
// What does cost something is reaching somewhere the molecule can get to on its own. A sphere over that line
// blocks pore volume the molecule is entitled to, and a simulation given such a file reports an uptake too
// low for a reason nothing in its output would explain. So a sphere is cut back at the nearest reachable
// point, and a pocket that a cut leaves partly uncovered is given further spheres until none of it is left.
export struct GridBlockingSphere
{
  // What a simulation reads: where the sphere sits in the cell, and how far it reaches.
  double3 centreFractional{0.0, 0.0, 0.0};
  double radius{0.0};  // Å

  // Which region it was drawn for, which grid point it was centred on, and how much of the region it took.
  std::size_t pocket{0};
  int3 centreVoxel{0, 0, 0};
  std::size_t voxelsCovered{0};

  // Whether somewhere reachable, rather than the pocket's own extent, is what stopped it growing.
  bool clipped{false};
};

export struct GridBlockingCover
{
  std::vector<GridBlockingSphere> spheres;

  std::size_t numberOfClippedSpheres{0};

  // Points of a pocket that lie nearer somewhere reachable than the grid can resolve, where no sphere can be
  // placed without taking reachable room with it. These are necks the grid has closed rather than pockets,
  // and they hold no volume worth blocking.
  std::size_t numberOfRefusedPoints{0};

  std::size_t numberOfPocketVoxels{0};

  double seconds{0.0};
};

// Covers every region of `pores` that `needsCover` marks, using `voxelPore` to find its points.
//
// `width` says how much room there is at each point, and is only used to order them, so any field that is
// larger where the pocket is roomier will do: a clearance, or a distance to a level set.
//
// `reachable` marks the points a sphere must not grow into. It is the caller's statement of where the
// molecule is entitled to be, and it is not simply the complement of the pockets: the space between them is
// neither, being room the molecule cannot occupy at all.
export GridBlockingCover coverPockets(uint3 gridSize, const UnitCell &unitCell,
                                  std::span<const std::int32_t> voxelPore, const std::vector<GridPore> &pores,
                                  std::span<const std::uint8_t> needsCover, std::span<const float> width,
                                  std::span<const std::uint8_t> reachable);
