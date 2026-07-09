module;

export module mc_moves_translation;

import std;

import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * \brief Attempts to translate a molecule within the system.
 *
 * This function performs a translation move on a selected molecule by applying
 * a random displacement along a randomly chosen axis. It computes the energy
 * difference resulting from the move and decides whether to accept or reject
 * the move based on the Metropolis criterion.
 *
 * \param random Reference to the random number generator.
 * \param system Reference to the simulation system.
 * \param selectedComponent Index of the selected component.
 * \param selectedMolecule Index of the selected molecule within the component.
 *
 * \return An optional RunningEnergy object containing the energy difference if
 *         the move is accepted; returns std::nullopt if the move is rejected.
 */
std::optional<RunningEnergy> translationMove(RandomNumber &random, System &system, std::size_t selectedComponent,
                                             std::size_t selectedMolecule);
}  // namespace MC_Moves
