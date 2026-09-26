module;

export module mc_moves_bead_flip;

import std;

import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * @brief Performs a bead-flip move (kink jump / end rotation) on a flexible molecule.
 *
 * One bead is chosen uniformly from the precomputed flip beads of the component (see
 * Component::flipBeads): beads with one or two bonded neighbours, outside rigid fragments, whose
 * rotation does not change a FIXED or RIGID bend or torsion. The bead is rotated such that all its
 * bond lengths are preserved:
 *  - a bead with two neighbours rotates about the axis through the neighbours (kink jump; the
 *    bend centred on the bead is preserved as well, the bends and torsions at the neighbours
 *    change);
 *  - a terminal bead rotates about a random axis through its neighbour, i.e. its bond direction
 *    moves on the sphere (end rotation; the bend and torsions at the neighbour change).
 * Every other atom stays in place. Only the bonded terms at the junctions, the intramolecular
 * non-bonded energy, and the external (framework, intermolecular, Ewald) energy change. The
 * proposals are symmetric, so the move is accepted with the Metropolis criterion on the total
 * energy difference.
 *
 * The rotation is sampled with mixed step sizes: a fraction component.beadFlipRandomizationFraction
 * of the attempts randomizes the rotation completely (angle uniform in [-pi, pi]; for a terminal
 * bead a uniformly random direction on the sphere; statistics channel 1), the remainder perturbs
 * the angle within an adaptive window (channel 0) tuned by the standard move statistics machinery
 * (Vitalis & Pappu, Methods 46 (2009)).
 *
 * The move is the smallest bond-length-preserving local move; it complements the pivot (which
 * cannot move terminal beads at all) and the crankshaft. Molecules without a valid flip bead
 * (rigid molecules, molecules with FIXED bends on every bead) and simulations with polarization
 * are not supported and lead to rejection.
 *
 * @param random Random number generator instance.
 * @param system The current state of the simulation system.
 * @param selectedComponent Index of the component to be moved.
 * @param selectedMolecule Index of the selected molecule within the component.
 * @return An optional `RunningEnergy` containing the energy difference if the move is accepted;
 *         `std::nullopt` if the move is rejected.
 */
std::optional<RunningEnergy> beadFlipMove(RandomNumber &random, System &system, std::size_t selectedComponent,
                                          std::size_t selectedMolecule);
}  // namespace MC_Moves
