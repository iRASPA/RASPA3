module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#endif

export module mc_moves_hybridmc;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import running_energy;
import randomnumbers;
import system;

export namespace MC_Moves
{
/**
 * \brief Performs a hybrid Monte Carlo move on the system.
 *
 * This function executes a hybrid Monte Carlo (MC) move by performing the following steps:
 * - Copies necessary data such as molecule positions and thermostat settings.
 * - Initializes velocities according to the Boltzmann distribution and removes the center of mass velocity.
 * - Computes translational and rotational kinetic energies and scales velocities to match the system temperature.
 * - Updates gradients and computes the reference energy.
 * - Integrates the system's equations of motion for a specified number of hybrid MC steps.
 * - Measures the CPU time taken for integration.
 * - Calculates the energy drift and decides whether to accept or reject the move based on the energy difference.
 * - If accepted, updates the system's state with the new positions and velocities.
 *
 * \param random A reference to the random number generator.
 * \param system A reference to the system on which the move is performed.
 * \return An optional RunningEnergy object containing the system's energy if the move is accepted, or std::nullopt if
 * rejected.
 */
std::optional<RunningEnergy> hybridMCMove(RandomNumber& random, System& system);
}  // namespace MC_Moves
