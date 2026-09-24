module;

export module cbmc_operators;

import std;

import atom;
import randomnumbers;
import forcefield;
import component;
import cbmc_growth_plan;

export namespace CBMC
{
/// One generated trial direction of a growth step: candidate positions of the step's next-beads plus
/// the torsion Rosenbluth weight accumulated while selecting the torsion rotation about the junction.
struct StepTrial
{
  std::vector<Atom> positions{};
  double torsionWeight{1.0};
};

/**
 * \brief Generates the trial-direction set of a growth step (grow direction).
 *
 * The single place the three growth operators live (eliminating the former parallel copies in
 * insertion / deletion / recoil):
 *  - PlaceSeedFragment: a flexible bead is placed with a Boltzmann bond length in a uniformly random
 *    direction; a rigid body or ring seed gets uniformly random orientations about the anchor (no
 *    orientational reference exists yet, all torsion weights are 1).
 *  - AttachFragment: flexible beads sample a base conformation exactly (Boltzmann bond lengths and
 *    anchor-bend cone directions per bead, sibling bends between beads of the same step imposed by
 *    rejection, declared chiral centers enforced by parity rejection) carrying no weight; every
 *    other bonded term of the step (torsions and spin-variant bends) is Rosenbluth-weighted in the
 *    torsion (spin) selection about the junction bond. A rigid body samples its junction-bend tilt
 *    with a rigid-rotation Metropolis Monte-Carlo.
 *  - CloseRing: the internal conformation of the cyclic cluster is sampled from its Boltzmann
 *    distribution by an internal Monte-Carlo (the closure bonds keep every ring closed -- simple,
 *    fused, and bridged), and the spin about the junction bond is biased by the junction-crossing
 *    torsions and bends.
 *
 * All trial directions of one step share one freshly sampled base conformation (they differ only by
 * the torsion spin), matching the coupled-decoupled bookkeeping on grow and retrace.
 */
std::vector<StepTrial> generateGrowTrials(RandomNumber &random, const ForceField &forceField, double beta,
                                          const Component &component, const std::vector<Atom> &chainAtoms,
                                          const GrowStep &step, std::size_t numberOfTrialDirections);

/**
 * \brief Generates the trial-direction set of a growth step for the retrace direction.
 *
 * The existing (old) positions of the step's next-beads -- read from 'chainAtoms' -- are pinned as
 * trial direction 0, with their torsion Rosenbluth weight computed by a pinned torsion selection;
 * the remaining directions mirror the grow scheme.
 */
std::vector<StepTrial> generateRetraceTrials(RandomNumber &random, const ForceField &forceField, double beta,
                                             const Component &component, const std::vector<Atom> &chainAtoms,
                                             const GrowStep &step, std::size_t numberOfTrialDirections);

/**
 * \brief Generates one independent trial direction of a growth step (recoil growth).
 *
 * Bond lengths and bend angles are drawn once from their Boltzmann distributions. When 'biasTorsion'
 * is true the spin about the preceding bond is chosen with the usual torsion Rosenbluth weight; feeler
 * look-ahead passes false and takes one random spin, because a feeler only tests whether an open path
 * exists and does not enter the chain weight.
 */
StepTrial generateRecoilTrial(RandomNumber &random, const ForceField &forceField, double beta,
                              const Component &component, const std::vector<Atom> &contextAtoms, const GrowStep &step,
                              bool biasTorsion = true);

/**
 * \brief Torsion Rosenbluth weight of the existing (old) orientation of a step (recoil retrace).
 *
 * Mirrors the CBMC deletion scheme: the real orientation is trial 0 and the remaining torsion trials
 * are random rotations around the last bond vector.
 */
double oldConfigurationTorsionWeight(RandomNumber &random, const ForceField &forceField, double beta,
                                     const std::vector<Atom> &oldAtoms, const GrowStep &step);

/**
 * \brief The log of the base-sampler normalization of a growth plan's flexible steps.
 *
 * The exact flexible base sampler draws each bead from the normalized density
 * r^2 exp(-beta u_bond) x sin(theta) exp(-beta u_anchor) / Z_step, restricted by rejection to the
 * sibling-bend coupling (contributing the mean coupling Boltzmann factor to Z_step) and optionally
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
}  // namespace CBMC
