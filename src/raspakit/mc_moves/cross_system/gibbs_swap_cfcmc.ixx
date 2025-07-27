module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#endif

export module mc_moves_gibbs_swap_cfcmc;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double3;
import randomnumbers;
import running_energy;
import atom;
import system;

/**
 * \brief Performs a Gibbs swap move in a CFCMC simulation.
 *
 * Attempts to perform a Gibbs swap move between two systems in a continuous fractional component Monte Carlo (CFCMC)
 * simulation. Depending on the randomly selected move type, it may swap molecules between systems, move the fractional
 * molecule between systems, or change the lambda scaling parameter of the fractional molecule.
 *
 * All systems have a fractional molecule; only one of these is 'active', the others are switched off with 'lambda=0'.
 * This implementation ensures that the number of fractional molecules per system remains constant.
 *
 * \param random Reference to the random number generator.
 * \param systemA Reference to the first system (usually containing the fractional molecule).
 * \param systemB Reference to the second system.
 * \param selectedComponent The index of the selected component for the move.
 * \param fractionalMoleculeSystem The system index that currently contains the fractional molecule (unused in this
 * implementation).
 * \return An optional pair of RunningEnergy objects representing the energy differences for the move in systemA and
 * systemB. Returns std::nullopt if the move is rejected.
 */
export namespace MC_Moves
{
std::optional<std::pair<RunningEnergy, RunningEnergy>> GibbsSwapMove_CFCMC(
    RandomNumber& random, System& systemA, System& systemB, std::size_t selectedComponent,
    [[maybe_unused]] std::size_t& fractionalMoleculeSystem);
}
