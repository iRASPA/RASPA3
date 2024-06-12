module;

#ifdef USE_LEGACY_HEADERS
#include <fstream>
#include <exception>
#include <source_location>
#include <complex>
#include <map>
#include <array>
#include <vector>
#include <ranges>
#include <algorithm>
#include <print>
#endif
  
module property_temperature;

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

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyTemperature &temp)
{
  archive << temp.versionNumber;
  archive << temp.numberOfBlocks;
  archive << temp.bookKeepingTemperature;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyTemperature &temp)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if(versionNumber > temp.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'PropertyTemperature' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> temp.numberOfBlocks;
  archive >> temp.bookKeepingTemperature;

  return archive;
}
