module;

#ifdef USE_LEGACY_HEADERS
#include <optional>
#include <span>
#endif

export module mc_moves_hybridmc;

#ifndef USE_LEGACY_HEADERS
import <span>;
import <optional>;
#endif

import running_energy;
import randomnumbers;
import system;

export namespace MC_Moves
{
std::optional<RunningEnergy> hybridMCMove(RandomNumber& random, System& system);
}