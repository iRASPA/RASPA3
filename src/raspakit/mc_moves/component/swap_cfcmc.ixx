module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#endif

export module mc_moves_swap_cfcmc;

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
 * \brief Performs a swap move in Continuous Fractional Component Monte Carlo (CFCMC).
 *
 * Attempts to perform an insertion, deletion, or lambda-change move on a selected component and molecule in the system.
 * The type of move is determined based on a randomly selected new lambda bin and the allowed maximum change.
 *
 * The steps for each move are as follows:
 * - **Insertion move** (selectedNewBin >= lambda.numberOfSamplePoints):
 *   1. The existing fractional molecule is made integer (lambda = 1.0), and the energy difference is computed.
 *   2. A new fractional molecule is inserted with a small lambda value (epsilon), and the energy difference is
 * computed.
 * - **Deletion move** (selectedNewBin < 0):
 *   1. The existing fractional molecule is removed (lambda = 0.0), and the energy difference is computed.
 *   2. A new fractional molecule is selected with a high lambda value (1 - epsilon), and the energy difference is
 * computed.
 * - **Lambda-change move** (0 <= selectedNewBin < lambda.numberOfSamplePoints):
 *   - The lambda value of the fractional molecule is changed, and the energy difference is computed.
 *
 * \param random Reference to the random number generator.
 * \param system Reference to the simulation system.
 * \param selectedComponent Index of the selected component.
 * \param selectedMolecule Index of the selected molecule.
 * \param insertionDisabled If true, insertion moves are disabled.
 * \param deletionDisabled If true, deletion moves are disabled.
 * \return A pair containing:
 *   - An optional RunningEnergy representing the energy difference if the move is accepted.
 *   - A double3 representing the acceptance probabilities for deletion, lambda-change, and insertion moves
 * respectively.
 */
std::pair<std::optional<RunningEnergy>, double3> swapMove_CFCMC(RandomNumber& random, System& system,
                                                                std::size_t selectedComponent,
                                                                std::size_t selectedMolecule,
                                                                bool insertionDisabled = false,
                                                                bool deletionDisabled = false);
}  // namespace MC_Moves
