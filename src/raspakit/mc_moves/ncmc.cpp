module;

#ifdef USE_LEGACY_HEADERS
#include <optional>
#endif

module mc_moves_ncmc;

#ifndef USE_LEGACY_HEADERS
import <optional>;
#endif

import running_energy;
import randomnumbers;
import system;

std::optional<RunningEnergy> MC_Moves::nonEquilibriumCandidate(RandomNumber& random, System& system)
{
    return std::nullopt;
}
