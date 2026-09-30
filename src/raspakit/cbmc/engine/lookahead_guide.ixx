module;

export module cbmc_lookahead_guide;

import std;

import double3;
import bond_potential;
import bend_potential;
import torsion_potential;
import intra_molecular_potentials;

// The lookahead guidance of the torsion-spin selection of a flexible attach step.
//
// A growth step judges the spin of its grown beads by the terms that exist at that moment: its own
// torsions and the pairs to placed beads. Every interaction that reaches FORWARD in the growth order
// -- the torsion of the next bond, the 1-5 and 1-6 pairs and the charges of a group placed a step or
// two later -- is invisible when the spin is chosen and is paid for afterwards as a small weight. A
// bare torsion that was fitted together with such terms (the TraPPE-UA ester and ether torsions, whose
// bare form favours cis while the 1-5 pairs enforce trans) then steers the spin into a geometry the
// next step cannot repair, and it does so in one growth direction only, since the mirror step sees
// the same pairs as placed ones. The guide is the missing information: g(phi), the Boltzmann average
// over the beads grown in the next few steps ('Constants::lookaheadGuideDepth' bonds beyond the
// grown beads) of every term that couples them to the step, as a function of the dihedral of the
// first grown bead against a placed reference neighbour of the previous bead. It multiplies the
// trial selection and is divided out of the weight of the selected trial (cbmc_torsion_selection),
// so ANY positive g leaves the sampled distribution exact; a good g only makes the single growth
// land where the complete molecule wants to be, which is what insertion, reinsertion, creation and
// Widom sampling need. This is the tabulated form of the 'arbitrary trial distribution' of Martin and
// Frischknecht (Mol. Phys. 104, 2439 (2006)) with the trial distribution chosen as the ideal-chain
// potential of mean force of the spin, and in spirit the rotational-isomeric-state look-ahead of the
// Theodorou-Suter amorphous-cell construction.
//
// The model is an ideal sub-molecule around the step in the step's own frame: the previous and
// current bead on the axis; the placed neighbours of the previous bead (the 'backward' beads, the
// lowest-indexed one being the reference of the dihedral) drawn from their bonds and bends to the
// axis; the grown beads drawn as the step's base sampler draws them; and the future beads drawn in
// breadth-first order from their bond, anchor bend and torsion, with sibling bends imposed by
// rejection. The remaining terms among the positioned beads that involve at least one future bead
// ('cross' terms: torsions not used as a sampler, the non-bonded pairs) are averaged: the ones that
// involve no backward bead are invariant under the spin and are evaluated once per sample, the rest
// on every grid point after spinning the grown and future beads about the axis. Placed beads further
// back than the neighbours of the previous bead, and future beads in rigid bodies or rings, are not
// part of the model (their terms are simply absent from the guide). The table is frozen at setup,
// once per distinct model signature and temperature, so grow and retrace see one time-independent
// bias.
export namespace CBMC
{
/**
 * \brief One bead of the lookahead model, placed from its positioned 'parent': a bond length from
 * 'bond' (or the default length), a direction on the cone about parent -> bendReference with the angle
 * of 'bend' (or uniform), and a torsion about the bond with respect to 'torsionReference' imposed by
 * rejection together with the 'siblingBends' to the beads of the same parent positioned before it.
 */
struct LookaheadBead
{
  std::size_t bead{};
  std::size_t parent{};
  std::optional<std::size_t> bendReference{};
  std::optional<std::size_t> torsionReference{};
  std::optional<BondPotential> bond{};
  std::optional<BendPotential> bend{};
  std::optional<TorsionPotential> torsion{};
  std::vector<BendPotential> siblingBends{};
};

/**
 * \brief The ideal sub-molecule of a step (see the module description). Atom identifiers are those of
 * the component; the tabulation works on a scratch molecule of 'numberOfAtoms' positions.
 */
struct LookaheadGuideModel
{
  std::size_t numberOfAtoms{};
  std::size_t previousBead{};
  std::size_t currentBead{};
  std::optional<BondPotential> axisBond{};  ///< previous - current.
  std::size_t referenceBead{};              ///< The backward bead defining the dihedral.
  std::size_t firstNextBead{};              ///< The grown bead defining the dihedral.
  std::vector<LookaheadBead> backwardBeads{};  ///< Reference bead first.
  std::vector<LookaheadBead> nextBeads{};
  std::vector<LookaheadBead> futureBeads{};  ///< Breadth-first from the grown beads.
  Potentials::IntraMolecularPotentials invariantTerms{};  ///< Cross terms without a backward bead.
  Potentials::IntraMolecularPotentials spinTerms{};       ///< Cross terms with a backward bead.

  /// Temperature-independent memo key with atom identifiers mapped to step-local roles, so congruent
  /// steps (the two ends of a symmetric molecule, every repeat unit of a polymer) share one table.
  [[nodiscard]] std::string signature() const;
};

/**
 * \brief The tabulated guide log g(phi) on a uniform periodic grid of the dihedral, linearly
 * interpolated, floored at a small fraction of its maximum so the old configuration of a retrace
 * always keeps a finite weight.
 */
struct LookaheadGuideTable
{
  std::vector<double> logGuide{};  ///< Point k at dihedral -pi + k * 2 pi / size.

  [[nodiscard]] double logGuideAt(double dihedral) const;
};

/// The signed dihedral a-b-c-d in [-pi, pi] (the convention shared by the tabulation and the lookup).
[[nodiscard]] double dihedralAngle(const double3 &a, const double3 &b, const double3 &c, const double3 &d);

/// Builds the guide table of a model for inverse temperature 'beta' (see the module description).
[[nodiscard]] LookaheadGuideTable buildLookaheadGuideTable(double beta, const LookaheadGuideModel &model);
}  // namespace CBMC
