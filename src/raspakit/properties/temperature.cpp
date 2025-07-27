module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <exception>
#include <fstream>
#include <map>
#include <print>
#include <ranges>
#include <source_location>
#include <vector>
#endif

module property_temperature;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyTemperature &temp)
{
  archive << temp.versionNumber;
  archive << temp.numberOfBlocks;
  archive << temp.bookKeepingTemperature;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyTemperature &temp)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > temp.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'PropertyTemperature' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> temp.numberOfBlocks;
  archive >> temp.bookKeepingTemperature;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("PropertySimulationBox: Error in binary restart\n"));
  }
#endif

  return archive;
}
