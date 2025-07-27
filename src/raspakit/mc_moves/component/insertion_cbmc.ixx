module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#endif

export module mc_moves_insertion_cbmc;

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
 * \brief Performs a Configurational-Bias Monte Carlo (CBMC) insertion move.
 *
 * Attempts to insert a molecule of the specified component into the system using CBMC.
 * Calculates the acceptance probability and updates the system accordingly.
 *
 * \param random Random number generator instance.
 * \param system Reference to the simulation system.
 * \param selectedComponent Index of the component to insert.
 * \return A pair containing an optional RunningEnergy and a double3 with acceptance statistics.
 *         The RunningEnergy is present if the move is accepted.
 */
std::pair<std::optional<RunningEnergy>, double3> insertionMoveCBMC(RandomNumber& random, System& system,
                                                                   std::size_t selectedComponent);
}  // namespace MC_Moves
