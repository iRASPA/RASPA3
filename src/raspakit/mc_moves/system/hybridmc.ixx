module;

export module mc_moves_hybridmc;

import std;

import running_energy;
import randomnumbers;
import system;

export namespace MC_Moves
{
/**
 * \brief Performs a hybrid Monte Carlo move on the system.
 *
 * Short NVE Velocity-Verlet trajectory used as a collective Monte Carlo proposal.
 * Supports rigid and flexible molecules in a rigid framework, a flexible framework,
 * or no framework. Flexible-framework atom positions and velocities are included in
 * the trial state and restored only on acceptance.
 *
 * \param random A reference to the random number generator.
 * \param system A reference to the system on which the move is performed.
 * \return An optional RunningEnergy object containing the system's energy if the move is accepted, or std::nullopt if
 * rejected.
 */
std::optional<RunningEnergy> hybridMCMove(RandomNumber& random, System& system);
}  // namespace MC_Moves
