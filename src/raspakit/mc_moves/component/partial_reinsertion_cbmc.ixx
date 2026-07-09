module;

export module mc_moves_partial_reinsertion;

import std;

import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * \brief Performs a CBMC reinsertion move for a selected molecule.
 *
 * Attempts to reinsert a molecule using Configurational Bias Monte Carlo (CBMC) method.
 * Updates the system if the move is accepted and returns the energy difference.
 *
 * \param random Random number generator instance.
 * \param system The simulation system.
 * \param selectedComponent Index of the selected component.
 * \param selectedMolecule Index of the selected molecule within the component.
 *
 * \return An optional RunningEnergy containing the energy difference if the move is accepted; std::nullopt otherwise.
 */
std::optional<RunningEnergy> partialReinsertionMove(RandomNumber &random, System &system, std::size_t selectedComponent,
                                                    std::size_t selectedMolecule);
}  // namespace MC_Moves
