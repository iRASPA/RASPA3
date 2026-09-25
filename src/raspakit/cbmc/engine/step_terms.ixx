module;

export module cbmc_step_terms;

import std;

import chiral_center;
import cbmc_growth_plan;

export namespace CBMC
{
/**
 * \brief Fills the derived (temperature-independent) sampler data of a growth step from its topology
 * and filtered potentials; see the derived fields of 'GrowStep'.
 *
 * The topological classification of a step's internal terms: which bends are spin-variant (weighted
 * in the torsion selection) or sibling bends (imposed on the base by rejection), each next bead's bond
 * and anchor bend, the ring-internal terms of a ring-closure step, the split of a flexible attach
 * step's unsampled terms into base coupling and spin terms, the chiral centres the step determines,
 * and the memo signature of the base coupling. Everything depends only on the step (which beads are
 * previous / current / grown), so it is evaluated once per step when the plan is built
 * ('buildGrowthPlan'); the operator engine then only does lookups.
 */
void prepareStep(GrowStep &step, const std::vector<ChiralCenter> &chiralCenters);
}  // namespace CBMC
