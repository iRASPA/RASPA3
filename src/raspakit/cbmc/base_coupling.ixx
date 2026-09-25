module;

export module cbmc_base_coupling;

import std;

import atom;
import cbmc_growth_plan;

export namespace CBMC
{
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
