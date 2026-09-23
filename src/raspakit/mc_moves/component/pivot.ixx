module;

export module mc_moves_pivot;

import std;

import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * @brief Performs a pivot move for a given flexible molecule.
 *
 * Selects a pivot axis uniformly from the precomputed valid pivot bonds of the
 * component (bonds that are not part of a ring, not interior to a rigid fragment,
 * and have a non-empty part to rotate) and rigidly rotates the smaller of the two
 * chain parts hanging off that bond about the bond axis, by a random angle drawn
 * uniformly from [-maxChange, +maxChange]. The rotation preserves all bond
 * lengths and bend angles; only the torsions through the pivot bond, the
 * intramolecular non-bonded energy, and the external (framework, intermolecular,
 * Ewald) energy change. The proposal is symmetric, so the move is accepted with
 * the Metropolis criterion on the total energy difference.
 *
 * The angle is sampled with mixed step sizes: a fraction
 * component.pivotRandomizationFraction of the attempts randomizes the angle
 * completely (uniform in [-pi, pi], statistics channel 1), the remainder perturbs
 * it within an adaptive window (channel 0) tuned by the standard move statistics
 * machinery. Sampling both length scales is robustly more efficient than either
 * alone (Vitalis and Pappu, Methods 46 (2009)).
 *
 * The move directly samples the internal torsional degrees of freedom of chain
 * molecules, for which full-chain reinsertion becomes ineffective with increasing
 * chain length. Rigid fragments and rings lying entirely on one side of the pivot
 * bond rotate rigidly, which preserves their internal geometry exactly. Molecules
 * without any valid pivot bond (fully rigid molecules, pure rings, diatomics) and
 * simulations with polarization are not supported and lead to rejection.
 *
 * @param random Random number generator instance.
 * @param system The current state of the simulation system.
 * @param selectedComponent Index of the component to be moved.
 * @param selectedMolecule Index of the selected molecule within the component.
 * @return An optional `RunningEnergy` containing the energy difference if the move is accepted;
 *         `std::nullopt` if the move is rejected.
 */
std::optional<RunningEnergy> pivotMove(RandomNumber &random, System &system, std::size_t selectedComponent,
                                       std::size_t selectedMolecule);
}  // namespace MC_Moves
