module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#endif

export module mc_moves_gibbs_volume;

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
 * \brief Performs a Gibbs ensemble volume move between two systems.
 *
 * Attempts to change the volumes of systemA and systemB while keeping the total volume constant,
 * as part of a Gibbs ensemble Monte Carlo simulation. The function scales the simulation boxes
 * and positions, computes the new energies, and decides whether to accept the move based on
 * the Metropolis criterion.
 *
 * \param random A random number generator.
 * \param systemA The first system involved in the volume move.
 * \param systemB The second system involved in the volume move.
 * \return An optional pair of RunningEnergy objects containing the new energies of systemA and systemB if the move is
 * accepted; std::nullopt otherwise.
 */
std::optional<std::pair<RunningEnergy, RunningEnergy>> GibbsVolumeMove(RandomNumber &random, System &systemA,
                                                                       System &systemB);
}  // namespace MC_Moves
