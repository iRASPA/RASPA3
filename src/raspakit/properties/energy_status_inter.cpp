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

module energy_status_inter;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import energy_factor;

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const EnergyInter &e)
{
  archive << e.versionNumber;

  archive << e.VanDerWaals;
  archive << e.VanDerWaalsTailCorrection;
  archive << e.CoulombicReal;
  archive << e.CoulombicFourier;
  archive << e.totalInter;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, EnergyInter &e)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > e.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'EnergyInter' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> e.VanDerWaals;
  archive >> e.VanDerWaalsTailCorrection;
  archive >> e.CoulombicReal;
  archive >> e.CoulombicFourier;
  archive >> e.totalInter;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("EnergyInter: Error in binary restart\n"));
  }
#endif

  return archive;
}
