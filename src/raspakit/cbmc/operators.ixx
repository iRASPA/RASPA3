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
 *  - AttachFragment: flexible beads sample a base conformation (bond lengths, bend angles, and at
 *    branch points the coupled branch arrangement, with chiral centers protected against parity
 *    flips) with an internal Metropolis Monte-Carlo carrying no weight; a rigid body samples its
 *    junction-bend tilt with a rigid-rotation Metropolis Monte-Carlo. Each trial direction then gets
 *    its own coupled-decoupled torsion (spin) selection about the junction bond, weighted per
 *    direction.
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
 * Each call samples a fresh base conformation, so consecutive calls are i.i.d. -- required by the
 * lazy one-at-a-time trial generation of the recoil-growth search.
 */
StepTrial generateRecoilTrial(RandomNumber &random, const ForceField &forceField, double beta,
                              const Component &component, const std::vector<Atom> &contextAtoms, const GrowStep &step);

/**
 * \brief Torsion Rosenbluth weight of the existing (old) orientation of a step (recoil retrace).
 *
 * Mirrors the CBMC deletion scheme: the real orientation is trial 0 and the remaining torsion trials
 * are random rotations around the last bond vector.
 */
double oldConfigurationTorsionWeight(RandomNumber &random, const ForceField &forceField, double beta,
                                     const std::vector<Atom> &oldAtoms, const GrowStep &step);
}  // namespace CBMC
