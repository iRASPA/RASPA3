module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#endif

export module mc_moves_reaction;

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
  std::optional<RunningEnergy>
  reactionMove([[maybe_unused]] RandomNumber &random, System& system,
               [[maybe_unused]] const std::vector<size_t> reactantStoichiometry,
               [[maybe_unused]] const std::vector<size_t> productStoichiometry);
}