module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#endif

export module mc_moves_noneq_cbmc;

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


std::pair<std::optional<RunningEnergy>, double3> NonEqCBMC(RandomNumber& random, System& system,
                                                                   size_t selectedComponent);