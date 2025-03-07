module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <algorithm>
#include <array>
#include <complex>
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
import <fstream>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;
import <map>;
import <vector>;
import <array>;
import <map>;
import <algorithm>;
import <print>;
#endif

import archive;

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const EnergyInter &e)
{
  archive << e.versionNumber;

  archive << e.VanDerWaals;
  archive << e.VanDerWaalsTailCorrection;
  archive << e.CoulombicReal;
  archive << e.CoulombicFourier;
  archive << e.totalInter;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, EnergyInter &e)
{
  uint64_t versionNumber;
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

  return archive;
}
