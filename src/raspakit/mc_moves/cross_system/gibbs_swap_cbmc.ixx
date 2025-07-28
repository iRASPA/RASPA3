module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#endif

export module mc_moves_gibbs_swap_cbmc;

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
 * \brief Performs a Gibbs ensemble CBMC swap move between two systems.
 *
 * Attempts to swap a molecule of the selected component from system B to system A
 * using Configurational-Bias Monte Carlo (CBMC). If successful, the molecule is
 * deleted from system B and inserted into system A.
 *
 * \param random Random number generator.
 * \param systemA The first system (receiving the molecule).
 * \param systemB The second system (donating the molecule).
 * \param selectedComponent The index of the component to swap.
 * \return A pair of RunningEnergy differences for systems A and B if the move is accepted;
 *         std::nullopt otherwise.
 */
std::optional<std::pair<RunningEnergy, RunningEnergy>> GibbsSwapMove_CBMC(RandomNumber& random, System& systemA,
                                                                          System& systemB,
                                                                          std::size_t selectedComponent);
}  // namespace MC_Moves
