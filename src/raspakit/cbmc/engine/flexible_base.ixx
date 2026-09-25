module;

export module cbmc_flexible_base;

import std;

import atom;
import randomnumbers;
import component;
import cbmc_growth_plan;

// The exact sampler of a flexible attach step's base conformation, and everything that belongs to its
// normalization: the base-coupling energy, the frozen per-signature coupling constants (estimated at
// setup for the system's temperature), the clamp-excess weight of a draw, and the log normalization
// of a plan that reptation needs. Ring-closure and rigid-body steps have their own samplers
// (cbmc_ring_closure, cbmc_rigid_tilt); the dispatch between them is in cbmc_operators.
export namespace CBMC
{
/// A sampled flexible base conformation: the positions of the step's next beads, and the clamp-excess
/// weight of the accepted draw (one unless the coupling energy fell below the reference constant).
struct FlexibleBase
{
  std::vector<Atom> nextBeadAtoms;
  double clampWeight;
};

/**
 * \brief Samples the base conformation of a flexible attach step exactly.
 *
 * Bond lengths are drawn from their one-dimensional Boltzmann densities, each direction from its
 * anchor-bend density on the cone about the previous-current axis (uniform azimuth), and the step's
 * internal couplings -- sibling-sibling bends and the base-routed unsampled terms (see the growth
 * plan's term classification) -- are imposed by clamped rejection sampling. The draw carries no
 * Rosenbluth weight apart from the rare clamp excess (see 'flexibleBaseClampWeight'), so this
 * distribution must be the exact bonded Boltzmann of its terms: the former internal Metropolis MC
 * (finite, reservoir-seeded, adaptive step sizes) only approximated it, and its deviations depended on
 * the structure of the growth step -- harmless for moves that pair grow and retrace on the same growth
 * plan, but a systematic bias for reptation, which pairs the grow weight of one chain end's plan
 * against the retrace weight of the other's.
 *
 * The base density is Boltzmann but its NORMALIZATION is plan-dependent (per-bead bond and
 * anchor-bend integrals, the coupling average, and a factor one half per determined chiral center);
 * 'logBaseSamplerNormalization' computes it, and reptation corrects its acceptance by the
 * normalization ratio of its two plans. Spin-variant bends (to placed atoms other than the previous
 * bead), all torsions, and the spin-routed unsampled terms must NOT shape the base: they are
 * Rosenbluth-weighted in the torsion-spin stage, identically on growth and retrace. Declared chiral
 * centers that are fully determined by this step are enforced by parity rejection: the exact
 * conditional distribution within the declared-parity sector, which carries probability exactly one
 * half by the reflection symmetry (through planes containing the previous-current axis) of the base
 * density -- the coupling terms preserve it, since distances, angles, and the cosine-even dihedral
 * forms are reflection-invariant. (The torsion spin afterwards is a proper rotation and preserves
 * that parity.)
 *
 * Every per-step ingredient (bonds, anchor bends, coupling terms, determined chiral centres, and the
 * frozen coupling constants) is precomputed in the growth plan; this function only samples. Throws
 * std::runtime_error when the rejection budget is exhausted (an inconsistent state, see the error
 * contract in the 'cbmc' module).
 *
 * 'chainAtoms' is the chain the step grows in; its next-beads are used as scratch for the draws and
 * are restored on return (see CBMC::ScratchBeads), so the chain is unchanged for the caller.
 */
FlexibleBase sampleExactFlexibleBase(RandomNumber &random, double beta, const Component &component,
                                     std::vector<Atom> &chainAtoms, const GrowStep &step);

/**
 * \brief The clamp-excess weight of a base conformation, max(1, e^{-beta (u - u_ref)}).
 *
 * One wherever the rejection acceptance was unclamped, the acceptance excess where the coupling
 * energy fell below the reference. It multiplies the step's trial weights on growth (freshly sampled
 * base) and retrace (the old positions ARE the base), restoring exactness of the density x weight
 * product for any reference constant. Returns one for steps without base coupling.
 */
double flexibleBaseClampWeight(double beta, const Component &component, const GrowStep &step,
                               std::span<const Atom> atoms);

/**
 * \brief The log of the base-sampler normalization of a growth plan's flexible steps.
 *
 * The exact flexible base sampler draws each bead from the normalized density
 * r^2 exp(-beta u_bond) x sin(theta) exp(-beta u_anchor) / Z_step, restricted by clamped rejection
 * to the coupling terms (sibling bends plus the base-routed unsampled terms, contributing
 * <min(1, e^{-beta(u-u_ref)})> e^{-beta u_ref} to Z_step) and optionally
 * to the declared chirality sector (probability one half). The Rosenbluth acceptance
 * W_grow / W_retrace of a CBMC move is therefore only correct up to the ratio of the two plans'
 * base normalizations: for moves that grow and retrace with the same plan the ratio is one and
 * cancels, but reptation pairs the grow weight of one chain end's plan against the retrace weight
 * of the other end's, and must correct its acceptance by exp(logZ_growPlan - logZ_retracePlan).
 *
 * Rigid-body and ring-closure steps are skipped: their internal-MC samplers have no closed-form
 * normalization, which is why reptation requires (at parse time) congruent end plans whenever a
 * repeat unit contains such steps -- congruent steps contribute equal factors that cancel.
 */
double logBaseSamplerNormalization(double beta, const Component &component, const std::vector<GrowStep> &plan);

/**
 * \brief The total coupling energy imposed on a flexible attach step's base conformation by rejection:
 * the sibling bends plus the base-routed unsampled terms (see the classification in the growth plan),
 * evaluated on real positions. Zero for steps without base coupling.
 */
double baseCouplingEnergy(const GrowStep &step, std::span<const Atom> atoms);

/**
 * \brief Estimates the frozen base-coupling constants of one flexible attach step at inverse
 * temperature 'beta' (see 'BaseCouplingConstants').
 *
 * A fixed-seed Monte Carlo that mirrors the base sampler's independent per-bead draws (Boltzmann bond
 * lengths and anchor-bend cone directions about a unit axis): the reference energy is the rigorous sum
 * of per-bend minima when only sibling bends couple, otherwise the sampled minimum minus a slack; the
 * clamped acceptance is averaged in batches until its standard error drops below a relative tolerance.
 * Deterministic for a given (signature, beta), so any two calls on congruent steps agree exactly.
 * 'numberOfAtoms' sizes the evaluation frame (the molecule's atom count). Returns zeros for a step
 * without base coupling.
 */
BaseCouplingConstants estimateBaseCouplingConstants(double beta, std::size_t numberOfAtoms, const GrowStep &step);

/**
 * \brief Fills 'baseCouplingConstants' of every step of 'plan' for inverse temperature 'beta'.
 *
 * Steps without base coupling are left without constants. Estimates are memoised in 'memo' per
 * 'baseCouplingSignature' so congruent steps -- within one plan and across the plans of one component
 * -- share one frozen number (required for their factors to cancel exactly in the reptation
 * acceptance). The memo must belong to one temperature; the caller clears it when 'beta' changes.
 */
void prepareBaseCouplingConstants(double beta, std::size_t numberOfAtoms, std::span<GrowStep> plan,
                                  std::map<std::string, BaseCouplingConstants> &memo);
}  // namespace CBMC
