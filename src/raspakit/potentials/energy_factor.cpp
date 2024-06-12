module;

#ifdef USE_LEGACY_HEADERS
#include <fstream>
#include <format>
#include <exception>
#include <source_location>
#include <complex>
#include <vector>
#include <array>
#include <map>
#include <algorithm>
#include <print>
#endif

module energy_factor;

#ifndef USE_LEGACY_HEADERS
import <fstream>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;
import <vector>;
import <array>;
import <map>;
import <algorithm>;
import <print>;
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
