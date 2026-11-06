module;

export module brute_force_pore_spikes;

import std;

import double3;
import brute_force_structure;
import brute_force_voxels;

// The last few pore-size spikes, from local clearance maxima and a Monte Carlo union volume.
//
// The exact PSD finds spikes as cliffs in the cumulative volume under refinement. That shares the
// Apollonius diagram and the solvent-excluded sweep with the rest of the exact routes. Here the same
// objects are read off the bare void alone:
//
//   A local maximum of the clearance is the centre of a maximal inscribed sphere. Walking uphill finds
//   those centres without a diagram: from every grid local maximum and from the roomiest voxels, as for
//   Di, once with a stride of a voxel and once with a short one, and from points thrown at the void.
//   Centres whose diameters agree are one family of pores — one spike. The walks are relied on for the
//   diameters. They are not relied on to end in every sphere of a family: where two maxima stand a
//   tenth of an Ångström apart, or one rises a few thousandths beside a slope, some are not found.
//
//   The spike weight is the fraction of the void whose pore size is that diameter exactly: the points
//   a sphere as wide as the family's holds and no wider sphere does. Either is any sphere that fits in
//   the void, wherever its centre is, and not only one of the maximal ones the walks ended in. Beside a
//   room a channel opens into there are spheres of every size between the channel's and the room's,
//   centred on the way in, and they hold part of what the channel's own maximal ball holds; that part
//   has a pore size above the channel's and is no part of its spike. Charging a point to the widest
//   maximal ball that holds it leaves it in, which is right for the widest family and for a family of
//   cages apart from all others, and too much for any other: on silicalite it gave the second spike
//   0.115 of the void where it is 0.107, and the third 0.045 where it is 0.017. And a point held only
//   by a sphere of the family that no walk ended in is the family's all the same: counted by the
//   centres found, the 1.87 Å spike of ITQ-29 came to 0.0129 of the void where it is 0.0140, half of
//   its spheres standing a tenth of an Ångström from the other half. It is estimated by throwing points
//   at the cell, the same binomial sample the pore-volume check uses.
//
// The last spike's diameter is Di. The families are reported largest first, as many as are asked for.
export struct BruteForceSpikeFamily
{
  double diameter{0.0};              // Å
  std::vector<double3> centres;      // Cartesian centres of the maximal spheres the walks ended in
  double weight{0.0};                // fraction of the void whose pore size is this diameter exactly
  double weightError{0.0};           // binomial standard error of that fraction
};

export struct BruteForcePoreSpikes
{
  std::vector<BruteForceSpikeFamily> families;  // as many as were asked for, largest diameter first

  std::size_t numberOfMaxima{0};    // distinct local clearance maxima before clustering
  std::size_t numberOfFamilies{0};  // distinct diameter families before taking the widest
  double voidVolume{0.0};           // Å³, from the same Monte Carlo draw
  double voidFraction{0.0};
  double seconds{0.0};

  // `voxels` must have been built from `structure` at threshold zero (bare radii).
  //
  // The weights are fractions of the void of `structure`, all of it. A void with pockets a probe cannot reach
  // is measured by handing in the structure with those pockets filled, which is how the exact route measures
  // it. Keeping only the points at which that probe's centre could stand is not the same thing: it drops the
  // whole of the void within a probe radius of a wall, in the channels as much as in the pockets, and a
  // weight over what is left is a fraction of another volume.
  static BruteForcePoreSpikes compute(const BruteForceStructure &structure, const BruteForceVoxels &voxels,
                                      std::size_t volumePoints, std::size_t maxFamilies = 3);
};
