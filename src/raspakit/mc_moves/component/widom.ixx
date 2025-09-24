module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#endif

export module mc_moves_widom;

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
 * \brief Performs a Widom insertion move for the specified component.
 *
 * Attempts to insert a molecule of the selected component into the system using
 * Configurational Bias Monte Carlo (CBMC) method. Calculates the energy differences
 * and computes the insertion weight used for chemical potential estimation.
 *
 * \param random Reference to the random number generator.
 * \param system Reference to the simulation system.
 * \param selectedComponent Index of the component to perform the Widom move on.
 * \return The Widom insertion weight if successful, or std::nullopt if the move was rejected.
 */
double WidomMove(RandomNumber& random, System& system, std::size_t selectedComponent);
}  // namespace MC_Moves
