module;

export module mc_moves_identity_change;

import std;

import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * \brief Performs a CBMC identity-change move for a selected molecule.
 *
 * Attempts to change the identity of a molecule from one component to another using
 * Configurational Bias Monte Carlo (CBMC). Updates the system if the move is accepted
 * and returns the energy difference.
 *
 * \param random Random number generator instance.
 * \param system The simulation system.
 * \param selectedComponent Index of the component selected for the move attempt.
 *
 * \return An optional RunningEnergy containing the energy difference if the move is accepted; std::nullopt otherwise.
 */
std::optional<RunningEnergy> identityChangeMove(RandomNumber &random, System &system,
                                                std::size_t selectedComponent);
}  // namespace MC_Moves
