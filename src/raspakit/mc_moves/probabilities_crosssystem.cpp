module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <complex>
#include <exception>
#include <format>
#include <fstream>
#include <iostream>
#include <map>
#include <print>
#include <source_location>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#endif

module mc_moves_probabilities_crosssystem;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <sstream>;
import <fstream>;
import <iostream>;
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

void MCMoveProbabilitiesCrossSystem::print() { std::cout << probability; }

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveProbabilitiesCrossSystem &p)
{
  archive << p.versionNumber;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveProbabilitiesCrossSystem &p)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > p.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'MCMoveProbabilitiesSystem' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  return archive;
}
