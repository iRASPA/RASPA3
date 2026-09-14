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
//   A local maximum of the clearance is the centre of a maximal inscribed sphere. Walking uphill from
//   every grid local maximum (and from the roomiest voxels, as for Di) finds those centres without a
//   diagram. Centres whose diameters agree are one family of pores — one spike.
//
//   The spike weight is the fraction of the void covered by the union of those maximal balls. That is
//   the Gelb–Gubbins mass at that diameter for a cage-like family. It is estimated by throwing points
//   at the cell, the same binomial sample the pore-volume check uses.
//
// The last spike's diameter is Di. The check reports up to three distinct families, largest first.
export struct BruteForceSpikeFamily
{
  double diameter{0.0};              // Å
  std::vector<double3> centres;      // Cartesian centres of the maximal spheres
  double weight{0.0};                // fraction of the void covered by their union
  double weightError{0.0};           // binomial standard error of that fraction
};

export struct BruteForcePoreSpikes
{
  std::vector<BruteForceSpikeFamily> families;  // up to three, largest diameter first

  std::size_t numberOfMaxima{0};    // distinct local clearance maxima before clustering
  std::size_t numberOfFamilies{0};  // distinct diameter families before taking the top three
  double voidVolume{0.0};           // Å³, from the same Monte Carlo draw
  double voidFraction{0.0};
  double seconds{0.0};

  // `voxels` must have been built from `structure` at threshold zero (bare radii).
  static BruteForcePoreSpikes compute(const BruteForceStructure &structure, const BruteForceVoxels &voxels,
                                      std::size_t volumePoints, std::size_t maxFamilies = 3);
};
