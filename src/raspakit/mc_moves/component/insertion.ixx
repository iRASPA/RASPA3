module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#endif

export module mc_moves_insertion;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double3;
import randomnumbers;
import running_energy;
import atom;
import system;

export namespace MC_Moves
{
/**
 * \brief Performs an insertion move in the Monte Carlo simulation.
 *
 * Attempts to insert a molecule of the specified component into the system.
 * Calculates energy differences, acceptance probability, and updates the system state if the move is accepted.
 *
 * \param random Reference to the random number generator.
 * \param system Reference to the simulation system.
 * \param selectedComponent Index of the component to be inserted.
 * \return A pair containing an optional RunningEnergy and a double3 vector.
 *         - If the move is accepted, the RunningEnergy contains the energy difference.
 *         - The double3 vector contains acceptance statistics.
 */
std::pair<std::optional<RunningEnergy>, double3> insertionMove(RandomNumber& random, System& system,
                                                               std::size_t selectedComponent);
}  // namespace MC_Moves
