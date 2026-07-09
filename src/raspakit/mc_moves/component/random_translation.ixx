module;

export module mc_moves_random_translation;

import std;

import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * @brief Performs a random translation move for a given molecule.
 *
 * Attempts to translate a molecule in a random direction and evaluates
 * the energy difference resulting from the move. The move is accepted or
 * rejected based on the Metropolis criterion.
 *
 * @param random Random number generator instance.
 * @param system The current state of the simulation system.
 * @param selectedComponent Index of the component to be moved.
 * @param selectedMolecule Index of the selected molecule within the component.
 * @return An optional `RunningEnergy` containing the energy difference if the move is accepted;
 *         `std::nullopt` if the move is rejected.
 */
std::optional<RunningEnergy> randomTranslationMove(RandomNumber &random, System &system, std::size_t selectedComponent,
                                                     std::size_t selectedMolecule);
}  // namespace MC_Moves
