module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <print>
#include <string>
#include <vector>
#endif

export module mc_void_fraction;

#ifndef USE_LEGACY_HEADERS
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
