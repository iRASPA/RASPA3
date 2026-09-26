module;

export module cbmc_bridge_closure;

import std;

import atom;
import double3;
import randomnumbers;
import cbmc_grow_step;
import cbmc_closure_guide;

// The base sampler of a bridge-closure step (fixed-endpoint CBMC regrowth), and the preparation of
// the closure guides of a plan.
//
// The closing bead n of a regrown interior segment is bonded to two placed beads: the anchor c the
// step grows from and the closure bead B. Given c and B at distance D, its position is described by
// bipolar coordinates (r1 = |n - c|, r2 = |n - B|, phi = the azimuth about the c-B axis), whose volume
// element is
//
//     d^3 n = (r1 r2 / D) dr1 dr2 dphi.
//
// The sampler draws r1 and r2 from the two bond samplers -- each the density r^2 exp(-beta u_bond(r))
// -- which fixes the bead to the circle at distance r1 from c and r2 from B (no circle exists when
// |r1 - r2| > D or r1 + r2 < D: an infeasible draw, a dead trial direction), and the torsion selection
// spins it about the axis. The generated density per unit volume is therefore
//
//     p(n) ~ r1^2 e^{-beta u1} r2^2 e^{-beta u2} / (2 pi) * D / (r1 r2),
//
// while the target is e^{-beta (u1 + u2 + u_mid + u_variant + u_ext)}; their ratio, up to constants
// that cancel between grow and retrace, is the base weight
//
//     w_base = exp(-beta u_mid(r1, r2, D)) / (r1 r2 D),
//
// with u_mid the bend c-n-B (a function of r1, r2 and D only, hence invariant under the spin). Every
// other bonded term of the step -- the bend at c, the bend at B, all torsions, all classically
// unsampled terms -- changes with phi and is Rosenbluth-weighted in the torsion selection. The factor
// 1/D is NOT a constant: D is set by the placement of c, which differs between the grown and the
// retraced configuration, so it must stay in the weight. Exact for any bond type including Fixed
// (delta-distributed lengths, where the bipolar element is the only freedom left).
export namespace CBMC
{
/// A sampled closure base: the closing bead on the c-B circle (random azimuth), and its base weight.
struct BridgeClosureBase
{
  Atom atom;
  double weight;
};

/// Draws the two bond lengths and places the closing bead; std::nullopt when the draw cannot span the
/// anchor-closure distance (a dead trial direction). The chain is not modified.
[[nodiscard]] std::optional<BridgeClosureBase> sampleBridgeClosureBase(RandomNumber &random, double beta,
                                                                       std::span<const Atom> chainAtoms,
                                                                       const GrowStep &step);

/// The base weight of an existing closing-bead position (the retrace's pinned old configuration).
[[nodiscard]] double bridgeClosureBaseWeight(double beta, std::span<const Atom> chainAtoms, const GrowStep &step);

/// The spin axis of a closure step: the unit vector from the anchor to the closure bead (a degenerate
/// axis falls back to z, like the junction axis of an attach step).
[[nodiscard]] double3 bridgeClosureAxis(std::span<const Atom> chainAtoms, const GrowStep &step);

/**
 * \brief Fills the closure-guide tables of the attach steps of a plan for inverse temperature 'beta'
 * ('GrowStep::SpinSelectionData::ClosureGuide::table'), memoised per path signature so congruent
 * paths share one table. Steps without guides are untouched. Called by 'Component::prepareGrowthPlans'
 * and by 'Component::growthPlan' for a plan built after the component was prepared.
 */
void prepareClosureGuides(double beta, std::span<GrowStep> plan,
                          std::map<std::string, std::shared_ptr<const ClosureGuideTable>> &memo);
}  // namespace CBMC
