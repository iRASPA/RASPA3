module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <print>
#include <string>
#include <vector>
#endif

export module mc_void_fraction;

#ifdef USE_STD_IMPORT
import std;
#endif

import framework;
import forcefield;

export struct MC_VoidFraction
{
  std::vector<double> data;

  MC_VoidFraction() {};

  void run(const ForceField &forceField, const Framework &framework, std::size_t number_of_iterations);
};
