module;

#ifdef USE_LEGACY_HEADERS
#include <optional>
#endif

export module mc_moves_ncmc;

#ifndef USE_LEGACY_HEADERS
import <optional>;
#endif

import running_energy;
import randomnumbers;
import system;


export namespace MC_Moves
{
    std::optional<RunningEnergy> nonEquilibriumCandidate(RandomNumber& random, System& system);
}