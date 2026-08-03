module;

export module brute_force_pore_volume;

import std;

import double3;
import brute_force_structure;
import brute_force_voxels;

// How much room a crystal has, and how much of it a molecule let in from outside can reach, by throwing
// points at the cell and by counting voxels.
//
// The route being checked gets the total from the volume of the union of the inflated atoms, which is a
// finite sum of closed-form pieces, and splits it by applying the divergence theorem to the surface around
// each sealed pocket. Neither step samples anything, and both depend on the same decomposition of the
// boundary that the surface area does.
//
// Here the total is a binomial proportion: draw a point uniformly in the cell, ask whether it is clear of
// every atom of every image, count. That is unbiased and its error is known exactly, so it says how far the
// exact total can be from the truth without any argument about arcs. The split is made twice over, once by
// the same drawn points -- looking up which piece of the flooded void each landed in -- and once by
// counting the voxels of the flood directly. The two disagree by the difference between a sample and a
// grid, which is worth seeing.
export struct BruteForcePoreVolume
{
  double voidFraction{0.0};       // of the cell, from the drawn points
  double voidFractionError{0.0};  // its standard error, the sample being binomial

  double accessibleFraction{0.0};    // of the cell, the part in a channel
  double inaccessibleFraction{0.0};  // of the cell, the part sealed in a pocket

  // Void the grid never found: a point with room around it that no voxel near it was labelled with,
  // which is void thinner than the spacing. It belongs to one of the two above and this says how much
  // is unaccounted for, which is the honest way to say how coarse the grid was for this structure.
  double unresolvedFraction{0.0};
  double unresolvedVolume{0.0};  // Å³

  double voidVolume{0.0};          // Å³, per unit cell
  double accessibleVolume{0.0};    // Å³
  double inaccessibleVolume{0.0};  // Å³

  // The same three from the voxels of the flood rather than from the drawn points. A grid is biased by its
  // own spacing where a sample is not, so these are the coarser of the two; where they part company by more
  // than the sample's error the grid is too coarse to be splitting anything.
  double voidFractionFromVoxels{0.0};
  double accessibleFractionFromVoxels{0.0};
  double inaccessibleFractionFromVoxels{0.0};

  std::size_t numberOfChannels{0};
  std::size_t numberOfPockets{0};

  // Carried over from the flood: how many pieces it had apart that a straight line then proved joined.
  std::size_t necksProved{0};
  std::size_t necksTried{0};

  std::size_t numberOfPoints{0};
  double seconds{0.0};

  static BruteForcePoreVolume compute(const BruteForceStructure &structure, const BruteForceVoxels &voxels,
                                      std::size_t numberOfPoints);
};
