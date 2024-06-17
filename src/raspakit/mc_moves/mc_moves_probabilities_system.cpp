module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <complex>
#include <exception>
#include <format>
#include <fstream>
#include <map>
#include <print>
#include <source_location>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#endif

module mc_moves_probabilities_system;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <sstream>;
import <fstream>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;
import <vector>;
import <array>;
import <map>;
import <utility>;
import <algorithm>;
import <print>;
#endif

import archive;
import double3;
import stringutils;

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveProbabilitiesSystem &p)
{
  archive << p.versionNumber;

  archive << p.probabilityVolumeMove;
  archive << p.probabilityGibbsVolumeMove;
  archive << p.probabilityParallelTemperingSwap;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveProbabilitiesSystem &p)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > p.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'MCMoveProbabilitiesSystem' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> p.probabilityVolumeMove;
  archive >> p.probabilityGibbsVolumeMove;
  archive >> p.probabilityParallelTemperingSwap;

  return archive;
}
