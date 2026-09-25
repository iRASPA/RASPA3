module;

export module cbmc_rigid_tilt;

import std;

import atom;
import randomnumbers;
import component;
import cbmc_grow_step;

export namespace CBMC
{
/**
 * \brief Samples the orientation of a rigid-body fragment hinged on the step's anchor.
 *
 * For a seed step (no previous bead) or a junction without bends every orientation is equally likely
 * and a uniform random rotation is returned. Otherwise the junction-bend tilt is sampled: the bend
 * angle previous-current-inner is drawn from its Boltzmann distribution, the body is aligned to that
 * cone direction, the roll about it is Rosenbluth-selected on a grid, and the whole is relaxed by a
 * rigid-rotation Metropolis Monte-Carlo on the step's bends ('numberOfTrialMovesPerOpenBead' moves
 * per bead, adaptive step size from the anchor bead's CBMC statistics). Carries no Rosenbluth weight.
 *
 * Returns the positions of the step's next-beads (in 'step.nextBeads' order); the internal geometry
 * of the body is that of the component's reference atoms, exactly. 'chainAtoms' is the chain the step
 * grows in; its next-beads are used as scratch and restored on return (see CBMC::ScratchBeads).
 */
std::vector<Atom> generateRigidTilt(RandomNumber &random, std::size_t numberOfTrialMovesPerOpenBead, double beta,
                                    const Component &component, std::vector<Atom> &chainAtoms, const GrowStep &step);
}  // namespace CBMC
