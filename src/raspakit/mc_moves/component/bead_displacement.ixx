module;

export module mc_moves_bead_displacement;

import std;

import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * @brief Performs a single-bead displacement move on a flexible molecule.
 *
 * One bead of the molecule is chosen uniformly from the precomputed displaceable beads of the
 * component (beads outside rigid fragments that take part in no FIXED or RIGID bond, bend or
 * torsion) and displaced along one randomly chosen Cartesian direction by a random amount drawn
 * uniformly from [-maxChange, +maxChange]; every other atom stays in place. The move changes the
 * bond lengths, bend angles and torsions the bead takes part in, the intramolecular non-bonded
 * energy, and the external (framework, intermolecular, Ewald) energy. The proposal is symmetric,
 * so the move is accepted with the Metropolis criterion on the total energy difference. The
 * maximum displacement is tuned per direction towards a 50% acceptance ratio (statistics channels
 * 0, 1, 2 for x, y, z, as for the translation move).
 *
 * The move is the elementary local relaxation of a flexible molecule: it is the only move that
 * samples the bond-length distribution of flexible bonds (all other conformational moves are
 * rotations that preserve bond lengths), and it relaxes bends locally. Molecules without any
 * displaceable bead (rigid molecules, molecules with FIXED bonds) and simulations with
 * polarization are not supported and lead to rejection.
 *
 * @param random Random number generator instance.
 * @param system The current state of the simulation system.
 * @param selectedComponent Index of the component to be moved.
 * @param selectedMolecule Index of the selected molecule within the component.
 * @return An optional `RunningEnergy` containing the energy difference if the move is accepted;
 *         `std::nullopt` if the move is rejected.
 */
std::optional<RunningEnergy> beadDisplacementMove(RandomNumber &random, System &system,
                                                  std::size_t selectedComponent, std::size_t selectedMolecule);
}  // namespace MC_Moves
