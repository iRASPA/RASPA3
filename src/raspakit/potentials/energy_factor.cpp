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

#if DEBUG_ARCHIVE
  archive << static_cast<uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &Potentials::operator>>(Archive<std::ifstream> &archive, Potentials::EnergyFactor &e)
{
  archive >> e.energy;
  archive >> e.dUdlambda;

#if DEBUG_ARCHIVE
  uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("Potentials::EnergyFactor: Error in binary restart\n"));
  }
#endif

  return archive;
}
