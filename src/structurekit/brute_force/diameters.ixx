module;

export module brute_force_diameters;

import std;

import double3;
import brute_force_structure;
import brute_force_voxels;

// The three pore diameters, worked out on a regular grid and refined by walking.
//
// The routes being checked read these off a pore network: the Voronoi or Apollonius diagram reduces the
// void to a graph of vertices and edges, and the diameters are then a maximum over vertices and a widest
// path over edges. That is fast and it is exact if the diagram is right, which is the thing in question.
//
// Here there is no diagram. The clearance is evaluated at the centre of every voxel of a fine grid, and:
//
//   Di is the largest of them, refined by walking uphill from the roomiest few until the step is smaller
//   than any tolerance that matters. The grid only has to get near the right basin; the walk does the rest,
//   so Di is limited by the walk and not by the spacing.
//
//   Df is the widest bottleneck of a path that leaves the cell and comes back to where it started as a
//   different image of itself, which is what percolation is. Every pair of neighbouring void voxels is an
//   edge, as wide as the room at its narrower end, and the edges are then added widest first until the piece
//   they are in repeats itself under some lattice vector. The width at which that happens is Df/2. This is
//   limited by the spacing, and both ways: the true bottleneck lies wherever it lies, and a grid line passes
//   near it but not through it, which makes Df come out small; and the room along an edge can be less than at
//   its ends, by half its length at most, which lets it come out large. `guaranteedFreeSphereDiameter` is the
//   same pass with that half length taken off every edge, a width every sphere is shown to manage, and the true
//   Df is at least that.
//
//   Dif is the roomiest point reachable without ever passing anywhere narrower than Df/2, found on the same
//   pass and then walked uphill from the roomiest few voxels of that part of the void.
//
// So Di should agree with a correct diagram to the last few digits, and Df and Dif to about the spacing. Di
// and Dif err one way only, a grid being able to miss a cavity but not to invent one; Df is bracketed.
export struct BruteForceDiameters
{
  double includedSphereDiameter{0.0};          // Di, Å
  double freeSphereDiameter{0.0};              // Df, Å
  double includedAlongFreePathDiameter{0.0};   // Dif, Å

  // Df again, with every hop credited only with the width a sphere is shown to manage along the whole of
  // it rather than the width at its ends. The clearance changes by no more than the distance moved, so
  // subtracting half a hop's length from its ends can only understate it. The true Df is at least this,
  // and the pair of them bracket where it lies.
  double guaranteedFreeSphereDiameter{0.0};    // Å

  double3 includedSphereCentre{};              // where Di was attained
  double3 includedAlongFreePathCentre{};

  bool percolates{false};
  std::size_t dimensionality{0};  // in how many directions the widest path runs away, the most of any where several tie

  double3 spacing{};              // Å, of the grid it was read off
  std::size_t numberOfVoidVoxels{0};
  std::size_t numberOfEdges{0};

  // How much the walk uphill improved on the best voxel, which is how much the grid alone was missing. A
  // large number here means the grid was too coarse for the walk to have started in the right basin.
  double walkGainedForDi{0.0};  // Å

  double seconds{0.0};

  // `voxels` must have been built from `structure` at threshold zero.
  static BruteForceDiameters compute(const BruteForceStructure &structure, const BruteForceVoxels &voxels);
};
