module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <numbers>
#include <print>
#include <sstream>
#include <string>
#endif

export module cbmc_move_statistics;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import move_statistics;

export struct CBMCMoveStatistics
{
  std::uint64_t versionNumber{1};

  MoveStatistics<double> bondLengthChange{
      .maxChange = double(0.3), .lowerLimit = double(0.01), .upperLimit = double(0.5)};
  MoveStatistics<double> bendAngleChange{
      .maxChange = double(0.3), .lowerLimit = double(0.01), .upperLimit = double(0.5)};
  MoveStatistics<double> conePositionChange{
      .maxChange = double(0.3), .lowerLimit = double(0.01), .upperLimit = double(0.5)};

  const std::string writeMCMoveStatistics() const;

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const CBMCMoveStatistics& p);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, CBMCMoveStatistics& p);
};
