module;

export module cbmc_operators;

import std;

import atom;
import randomnumbers;
import forcefield;
import component;
import cbmc_growth_plan;

// The operator engine's dispatch: per step kind it picks the sampler of the base conformation
// (cbmc_flexible_base, cbmc_rigid_tilt, cbmc_ring_closure) and the torsion-spin selection
// (cbmc_torsion_selection), and packs the result as trial directions. The growth schemes (CBMC and
// recoil growth) see only this interface and the 'StepTrial' it returns.
//
// Scratch contract: every function below takes the chain it grows in by non-const reference. The
// step's next-beads are used as scratch for the trial placements (so the energy routines see a
// complete, correctly indexed molecule without a per-trial copy of the whole chain) and are restored
// before the function returns, on every exit path; no other atom is touched. The chain is therefore
// unchanged for the caller.
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
 *    anchor-bend cone directions per bead; sibling bends and the step-local spin-invariant coupling
 *    terms -- Urey-Bradley, central inversion bends, bond/bond, bond/bend, bend/bend, improper
 *    torsions -- imposed by clamped rejection, whose excess weight max(1, e^{-beta(u-u_ref)}) rides
 *    on every trial; declared chiral centers enforced by parity rejection); every other bonded term
 *    of the step (torsions, spin-variant bends, and the remaining unsampled terms that touch placed
 *    geometry) is Rosenbluth-weighted in the torsion (spin) selection about the junction bond. A
 *    rigid body samples its junction-bend tilt with a rigid-rotation Metropolis Monte-Carlo.
 *  - CloseRing: the internal conformation of the cyclic cluster is sampled from its Boltzmann
 *    distribution by an internal Monte-Carlo (the closure bonds keep every ring closed -- simple,
 *    fused, and bridged), and the spin about the junction bond is biased by the junction-crossing
 *    torsions and bends.
 *
 * All trial directions of one step share one freshly sampled base conformation (they differ only by
 * the torsion spin), matching the coupled-decoupled bookkeeping on grow and retrace.
 */
std::vector<StepTrial> generateGrowTrials(RandomNumber &random, const ForceField &forceField, double beta,
                                          const Component &component, std::vector<Atom> &chainAtoms,
                                          const GrowStep &step, std::size_t numberOfTrialDirections);

/**
 * \brief Generates the trial-direction set of a growth step for the retrace direction.
 *
 * The existing (old) positions of the step's next-beads -- read from 'chainAtoms' -- are pinned as
 * trial direction 0, with their torsion Rosenbluth weight computed by a pinned torsion selection;
 * the remaining directions mirror the grow scheme.
 */
std::vector<StepTrial> generateRetraceTrials(RandomNumber &random, const ForceField &forceField, double beta,
                                             const Component &component, std::vector<Atom> &chainAtoms,
                                             const GrowStep &step, std::size_t numberOfTrialDirections);

/**
 * \brief Generates one independent trial direction of a growth step (recoil growth).
 *
 * Bond lengths and bend angles are drawn once from their Boltzmann distributions and the spin about
 * the preceding bond is chosen with the usual torsion Rosenbluth selection; the returned
 * 'torsionWeight' is that selection's weight.
 *
 * Recoil growth uses this one generator for every trial it ever draws: the directions of the growth
 * itself, the alternatives counted after the chain is complete, and every bead of a feeler. That is a
 * requirement, not a convenience: a trial direction counts as 'available' when it is open and a feeler
 * can be grown from it, and for the directions tried during growth the growth attempt itself IS the
 * feeler. The count of available directions m_i enters the acceptance ratio, so the probability that a
 * direction is found available must be the same random experiment whether it is decided by the growth
 * attempt or by an explicit feeler (on grow and on retrace). A feeler that took a plain random spin,
 * while the growth spins were torsion-selected, would probe the openness of a differently distributed
 * bead and bias m_i between the two directions of the move.
 */
StepTrial generateRecoilTrial(RandomNumber &random, const ForceField &forceField, double beta,
                              const Component &component, std::vector<Atom> &contextAtoms, const GrowStep &step);

/**
 * \brief Torsion Rosenbluth weight of the existing (old) orientation of a step (recoil retrace).
 *
 * Mirrors the CBMC deletion scheme: the real orientation is trial 0 and the remaining torsion trials
 * are random rotations around the last bond vector.
 */
double oldConfigurationTorsionWeight(RandomNumber &random, const ForceField &forceField, double beta,
                                     const Component &component, std::vector<Atom> &oldAtoms, const GrowStep &step);

/**
 * \brief Whether a step's growth machinery samples/weights the classically unsampled internal terms
 * itself (Urey-Bradley, inversion and out-of-plane bends, improper torsions, and the cross terms).
 *
 * True for flexible attach steps: their base sampler imposes the step-local ("base coupling") share
 * of these terms by rejection, and the torsion-spin selection Rosenbluth-weights the rest. Callers
 * must then NOT apply the post-selection e^{-beta u} factor for such steps (it would double count).
 * Rigid-body, ring-closure, and seed steps still rely on the caller's post-selection factor.
 */
bool stepHandlesUnsampledInternalTerms(const GrowStep &step);
}  // namespace CBMC
