module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <complex>
#include <exception>
#include <fstream>
#include <map>
#include <print>
#include <ranges>
#include <source_location>
#include <vector>
#endif

module property_simulationbox;

#ifndef USE_LEGACY_HEADERS
import <fstream>;
import <exception>;
import <source_location>;
import <complex>;
import <map>;
import <array>;
import <vector>;
import <algorithm>;
import <print>;
#endif

import archive;

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertySimulationBox &box)
{
  archive << box.versionNumber;
  archive << box.numberOfBlocks;
  archive << box.bookKeepingSimulationBox;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertySimulationBox &box)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > box.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'PropertySimulationBox' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> box.numberOfBlocks;
  archive >> box.bookKeepingSimulationBox;

  return archive;
}
