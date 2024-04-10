module;

#ifdef USE_LEGACY_HEADERS
#include <fstream>
#if defined(__has_include) && __has_include(<format>)
#include <format>
#endif
#include <exception>
#include <source_location>
#include <complex>
#if defined(__has_include) && __has_include(<print>)
#include <print>
#endif
#endif

module energy_factor;

#ifndef USE_LEGACY_HEADERS
import <fstream>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;
#if defined(__has_include) && __has_include(<print>)
  import <print>;
#endif
#endif

#if !(defined(__has_include) && __has_include(<print>))
  import print;
#endif


import archive;


Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const EnergyFactor &e)
{
  archive << e.energy;
  archive << e.dUdlambda;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, EnergyFactor &e)
{
  archive >> e.energy;
  archive >> e.dUdlambda;

  return archive;
}
