module;

export module cbmc_closure_guide;

import std;

import bond_potential;
import bend_potential;
import torsion_potential;

// The closure guidance of fixed-endpoint (bridging) CBMC regrowth.
//
// When a growth plan regrows an interior segment between two placed anchors, every bead of the
// segment but the last is placed by an ordinary attach step, and the last bead by a bridge-closure
// step (cbmc_bridge_closure). An attach step knows nothing about the far anchor: left alone it would
// grow the segment away from it, and the closure step would find no geometry that reaches the anchor.
// The guide is a bias on the torsion-spin selection of such a step, a positive function g(D) of the
// distance D between a trial position of the grown bead and the placed bead it must eventually reach
// (its 'closure target'). It is divided out of the step's Rosenbluth weight again (see
// cbmc_torsion_selection), so ANY positive g leaves the sampled distribution exact; the choice of g
// only decides how many grows reach a closable geometry. The g used here is the ideal-chain
// approximation of the partition function of the remaining segment: the Boltzmann-weighted density
// of the target at distance D from the grown bead, over the bonded terms (bonds, bends, torsions)
// along the shortest path from the bead to the target. It is tabulated once per distinct path (per
// temperature) at setup: exactly for a two-bond path (the geometry of the closure step itself), by
// sampling the ideal sub-chain for longer paths. This is the tabulated form of the 'self-adapting'
// closing bias of Wick and Siepmann (Macromolecules 33, 7207 (2000)), frozen at setup so grow and
// retrace see one time-independent bias.
export namespace CBMC
{
/**
 * \brief The bonded terms along the shortest path from a grown bead to its closure target: m bonds,
 * m-1 bends and m-2 torsions for a path of m bonds (m >= 2). An absent term means the engine's default
 * geometry (a bond at 'Constants::defaultBondLength', a free bend or torsion).
 */
struct ClosureGuidePath
{
  std::vector<std::optional<BondPotential>> bonds{};
  std::vector<std::optional<BendPotential>> bends{};
  std::vector<std::optional<TorsionPotential>> torsions{};

  [[nodiscard]] std::size_t numberOfBonds() const noexcept { return bonds.size(); }

  /// Temperature-independent memo key: the potential types and parameters along the path. Congruent
  /// paths (the same sequence of terms) share one table.
  [[nodiscard]] std::string signature() const;
};

/**
 * \brief The tabulated guide log g(D) on a uniform grid of the bead-target distance D, with linear
 * interpolation and a floor: g is clamped from below at a small fraction of its maximum (also beyond
 * the reach of the path), so the old configuration of a retrace always has a finite weight and every
 * trial keeps a non-zero, if negligible, probability.
 */
struct ClosureGuideTable
{
  double spacing{0.0};
  std::vector<double> logGuide{};

  [[nodiscard]] double logGuideAt(double distance) const;
};

/// Builds the guide table of a path for inverse temperature 'beta' (see the module description).
[[nodiscard]] ClosureGuideTable buildClosureGuideTable(double beta, const ClosureGuidePath &path);
}  // namespace CBMC
