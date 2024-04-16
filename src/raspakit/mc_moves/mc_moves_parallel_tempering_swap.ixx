module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#endif

export module mc_moves_parallel_tempering_swap;

#ifndef USE_LEGACY_HEADERS
import <cstddef>;
import <optional>;
import <span>;
import <tuple>;
#endif

import double3;
import randomnumbers;
import running_energy;
import atom;
import system;

export namespace MC_Moves
{
std::optional<std::pair<RunningEnergy, RunningEnergy>> ParallelTemperingSwap(RandomNumber &random, System &systemA,
                                                                             System &systemB);
}
