module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#endif

export module mc_moves_swap_cfcmc_cbmc;

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
 * \brief Performs a swap move in Configurational Bias Monte Carlo (CBMC) with Continuous Fractional Component Monte
 * Carlo (CFCMC).
 *
 * This function attempts to perform a swap move by either inserting, deleting, or changing the lambda value of a
 * molecule. The move type is determined based on the selected new bin. It calculates energy differences, updates the
 * system state, and applies acceptance criteria based on the Metropolis algorithm.
 *
 * \param random             Random number generator.
 * \param system             The simulation system containing all molecules and interactions.
 * \param selectedComponent  Index of the selected component for the move.
 * \param selectedMolecule   Index of the selected molecule for the move.
 * \param insertionDisabled  Flag to disable insertion moves.
 * \param deletionDisabled   Flag to disable deletion moves.
 *
 * \return A pair containing an optional RunningEnergy and a double3 vector. The RunningEnergy is present if the move is
 * accepted. The double3 vector contains acceptance probabilities for different move types.
 */
std::pair<std::optional<RunningEnergy>, double3> swapMove_CFCMC_CBMC(RandomNumber& random, System& system,
                                                                     std::size_t selectedComponent,
                                                                     std::size_t selectedMolecule,
                                                                     bool insertionDisabled = false,
                                                                     bool deletionDisabled = false);
}  // namespace MC_Moves
