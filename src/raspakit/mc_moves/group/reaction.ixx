module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#endif

export module mc_moves_reaction;

#ifdef USE_STD_IMPORT
import std;
#endif

import double3;
import randomnumbers;
import running_energy;
import atom;
import system;

export namespace MC_Moves
{
std::optional<RunningEnergy> reactionMove([[maybe_unused]] RandomNumber& random, System& system,
                                          [[maybe_unused]] const std::vector<std::size_t> reactantStoichiometry,
                                          [[maybe_unused]] const std::vector<std::size_t> productStoichiometry);
}
