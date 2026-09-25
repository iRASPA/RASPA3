module;

export module cbmc_torsion_selection;

import std;

import atom;
import double3;
import randomnumbers;
import cbmc_growth_plan;

export namespace CBMC
{
/// The outcome of the torsion (spin) selection of a growth step: the positions of the step's
/// next-beads after the selected spin, and the selection's Rosenbluth weight (the mean Boltzmann
/// factor over the torsion trials).
struct TorsionOrientation
{
  std::vector<Atom> positions;
  double rosenbluthWeight;
};

/**
 * \brief The coupled-decoupled torsion (spin) step, shared by all operators and both directions.
 *
 * The spin about the previous-current axis is selected among 'numberOfTorsionTrials' rotations of
 * 'baseOrientation' (the positions of the step's next-beads) by a Rosenbluth selection over the step's
 * precomputed torsion-selection data: its torsions, its spin-variant bends (bends to placed atoms other
 * than the previous bead, which are not invariant under the spin), and -- for a flexible attach step --
 * the spin-routed unsampled terms. All are Rosenbluth-weighted identically on growth and retrace.
 *
 * With 'pinFirstToBase' the base orientation itself is torsion trial 0 and is the selected one (the
 * retrace of an existing orientation); the remaining trials are random spins that only enter the
 * weight.
 *
 * 'chainAtoms' is the chain the step grows in; its next-beads are used as scratch for the trial spins
 * and are restored on return (see CBMC::ScratchBeads), so the chain is unchanged for the caller.
 */
TorsionOrientation selectTorsionOrientation(RandomNumber &random, std::size_t numberOfTorsionTrials, double beta,
                                            std::vector<Atom> &chainAtoms, const std::vector<Atom> &baseOrientation,
                                            const GrowStep &step, double3 lastBondVector, bool pinFirstToBase);
}  // namespace CBMC
