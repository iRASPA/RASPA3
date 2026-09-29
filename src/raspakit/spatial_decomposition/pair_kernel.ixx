module;

export module spatial_decomposition_pair_kernel;

import std;

// The tabulated Ewald real-space term (EwaldRealSpaceTable) lives in the potentials layer, shared with
// Potentials::potentialCoulomb; re-exported here for the fast kernel.
export import potential_ewald_real_space_table;

/// Lennard-Jones parameters of one pair of pseudo-atom types, laid out for the fast kernel.
export struct LennardJonesPair
{
  double epsilon4{0.0};       ///< 4 epsilon [energy units].
  double inverseSigma2{1.0};  ///< 1 / sigma^2.
  double shift{0.0};          ///< Energy shift at the cutoff (0 for truncated potentials).
};
