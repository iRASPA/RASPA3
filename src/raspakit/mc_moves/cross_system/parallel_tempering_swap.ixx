module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#endif

export module mc_moves_parallel_tempering_swap;

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
 * \brief Performs a parallel tempering swap between two systems.
 *
 * Attempts to swap the configurations of two systems in a parallel tempering Monte Carlo simulation.
 * The swap is accepted or rejected based on the Metropolis criterion, taking into account differences
 * in temperature, pressure, and potential energy between the two systems.
 *
 * Reference: "Hyper-parallel tempering Monte Carlo: Application to the Lennard-Jones fluid and the
 * restricted primitive model", G. Yan and J.J. de Pablo, JCP, 111(21): 9509-9516, 1999.
 *
 * \param random   Random number generator used for acceptance probability.
 * \param systemA  First system involved in the swap.
 * \param systemB  Second system involved in the swap.
 * \return         An optional pair of RunningEnergy objects if the swap is accepted; std::nullopt otherwise.
 */
std::optional<std::pair<RunningEnergy, RunningEnergy>> ParallelTemperingSwap(RandomNumber &random, System &systemA,
                                                                             System &systemB);
}  // namespace MC_Moves
