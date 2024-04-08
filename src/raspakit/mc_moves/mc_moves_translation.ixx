module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#endif

export module mc_moves_translation;

#ifndef USE_LEGACY_HEADERS
import <cstddef>;
import <optional>;
import <span>;
#endif

import randomnumbers;
import running_energy;
import atom;
import system;

export namespace MC_Moves
{
  std::optional<RunningEnergy>                                                                                            
  translationMove(RandomNumber &random, System &system, size_t selectedComponent, std::span<Atom> molecule);
}
