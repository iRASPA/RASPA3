module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <exception>
#include <format>
#include <fstream>
#include <map>
#include <print>
#include <source_location>
#include <vector>
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

Archive<std::ofstream> &Potentials::operator<<(Archive<std::ofstream> &archive, const Potentials::EnergyFactor &e)
{
  archive << e.energy;
  archive << e.dUdlambda;

  return archive;
}

Archive<std::ifstream> &Potentials::operator>>(Archive<std::ifstream> &archive, Potentials::EnergyFactor &e)
{
  archive >> e.energy;
  archive >> e.dUdlambda;

  return archive;
}
