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
#include <ostream>
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
import units;
import stringutils;

std::string PropertyTemperature::writeAveragesStatistics(const std::string tag)
{
  std::ostringstream stream;
  double conv = Units::TemperatureConversionFactor;

  std::pair<double, double> temperatureAverages = averageTemperature();
  for (std::size_t i = 0; i < bookKeepingTemperature.size(); ++i)
  {
    double blockAverage = averagedTemperature(i);
    std::print(stream, "    Block[ {:2d}] {: .6e}\n", i, conv * blockAverage);
  }
  std::print(stream, "    ---------------------------------------------------------------------------\n");
  std::print(stream, "    {} temperature  {: .6e} +/- {: .6e} [K]\n", tag, conv * temperatureAverages.first,
             temperatureAverages.second);
  std::print(stream, "\n");

  return stream.str();
}

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
