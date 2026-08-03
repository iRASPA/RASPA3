module;

export module brute_force_voxels;

import std;

import int3;
import double3;
import brute_force_structure;

// One connected piece of the void, as the voxels see it.
export struct BruteForceRegion
{
  std::size_t numberOfVoxels{0};
  double volume{0.0};  // Å³
  bool percolates{false};
  std::size_t dimensionality{0};  // in how many independent directions it runs away

  std::size_t deepestVoxel{0};  // the roomiest voxel of the piece, an index into the grid
  double deepestClearance{0.0};  // Å, the room there
};

// The void of a crystal on a regular grid, and how it falls apart into pieces.
//
// This is the brute force's answer to the pore classifier, and it is deliberately the stupidest thing that
// works: evaluate the clearance at the centre of every voxel of the cell, call a voxel void where that is
// not negative, and flood the void voxels. Nothing is solved for and no diagram is built, so it agrees with
// the classifier only if they are both right.
//
// Whether a piece of the void runs away through the crystal is decided the same blunt way. Flooding a
// periodic grid, a piece can arrive back at a voxel it has already reached by a different route, and if the
// two routes crossed the boundary a different number of times their difference is a lattice vector the
// piece is periodic under. Collecting those vectors and taking the rank of what they span gives the
// dimensionality, and a rank above zero is what makes the piece a channel rather than a pocket.
//
// What a grid on its own cannot do is see a neck narrower than a voxel, and in a zeolite that is not a
// corner case but the rule. The room a probe centre has pinches to nothing at every window the probe barely
// fits through, so the piece of space that gets it from one cage to the next can be a filament thinner than
// any grid one can afford, and those are precisely the windows that decide whether a cage is a pocket or
// part of a channel. A flood that only joins neighbouring voxels calls such a cage sealed.
//
// So the flood is followed by a second look that does not depend on the spacing at all. Where two pieces
// come close, the straight line from a voxel of one to a voxel of the other is offered as a passage, and
// whether a sphere fits along the whole of it is not sampled but worked out: the room along a segment is the
// distance from the nearest atom to that segment, less its radius. A join made that way is a join that was
// demonstrated. The grid is left deciding only which pairs are worth asking about, which is a far weaker
// thing to ask of it than deciding what is connected to what.
//
// What remains one-sided is separation. A pair not joined here is a pair no straight line got through,
// which is not the same as a pair nothing gets through, so the pieces can still be finer than the truth and
// never coarser. Refining the spacing and re-running is what says whether that matters on a given crystal.
export struct BruteForceVoxels
{
  int3 counts{0, 0, 0};
  double3 spacing{0.0, 0.0, 0.0};  // Å, along each cell axis
  double voxelVolume{0.0};         // Å³

  std::vector<float> clearance;       // Å, at the centre of each voxel
  std::vector<std::int32_t> regionOf;  // -1 inside an atom, else which piece of the void

  std::vector<BruteForceRegion> regions;

  std::size_t numberOfVoidVoxels{0};
  std::size_t numberOfChannels{0};
  std::size_t numberOfPockets{0};

  // Pieces the flood had apart and a straight line proved joined, and how many such lines were tried. A
  // count well above zero says the spacing is coarse against the necks of this crystal, and that the split
  // between channel and pocket would have been wrong without the second look.
  std::size_t necksProved{0};
  std::size_t necksTried{0};

  double voidFraction{0.0};         // of the cell
  double accessibleFraction{0.0};   // of the cell, the part in a channel
  double blockedFraction{0.0};      // of the cell, the part in a pocket

  double seconds{0.0};

  std::size_t indexOf(std::int32_t i, std::int32_t j, std::int32_t k) const
  {
    return (static_cast<std::size_t>(k) * static_cast<std::size_t>(counts.y) + static_cast<std::size_t>(j)) *
               static_cast<std::size_t>(counts.x) +
           static_cast<std::size_t>(i);
  }

  // The centre of a voxel, in fractional and in Cartesian coordinates.
  double3 fractionalCentre(std::size_t index) const;
  double3 centre(const BruteForceStructure &structure, std::size_t index) const;

  // Which voxel holds a point, wrapping into the home cell.
  std::size_t voxelOf(const BruteForceStructure &structure, const double3 &position) const;

  // Which piece of the void a point belongs to, for a point that is in the void but need not be at the
  // centre of anything.
  //
  // A voxel is labelled by the clearance at its centre, and a point can perfectly well have room around it
  // while the centre of the voxel it falls in does not: near a wall the two disagree over the last half
  // voxel. Taking the label of the voxel it falls in would then throw away a shell of the void one voxel
  // thick, which on a tight structure is most of it. So where that voxel is labelled solid the nearest
  // labelled one within a ring or two is taken instead, and -1 comes back only when there is genuinely no
  // void nearby.
  std::int32_t regionNear(const BruteForceStructure &structure, const double3 &position,
                          std::int32_t rings = 2) const;

  // Whether a point is in a piece of the void that runs away through the crystal. Points inside an atom, and
  // points in a pocket, are both false; `isInVoid` tells them apart.
  bool isAccessible(const BruteForceStructure &structure, const double3 &position) const;
  bool isInVoid(const BruteForceStructure &structure, const double3 &position) const;

  // The grid, at a spacing of about `targetSpacing` along each axis. `threshold` is how much room a voxel
  // must have to count as void, which is zero for the void itself and more when a probe of that radius has
  // to fit; the pieces are then the pieces that probe can be in.
  static BruteForceVoxels build(const BruteForceStructure &structure, double targetSpacing,
                                double threshold = 0.0);
};
