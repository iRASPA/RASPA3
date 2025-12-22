module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <numbers>
#include <print>
#include <sstream>
#include <string>
#endif

export module cbmc_move_statistics;

#ifdef USE_STD_IMPORT
import std;
#endif

import archive;
import move_statistics;

export struct CBMCMoveStatistics
{
  std::uint64_t versionNumber{1};

  MoveStatistics<double> bondLengthChange{
      .maxChange = 0.3, .lowerLimit = 0.01, .upperLimit = 0.5};
  MoveStatistics<double> bendAngleChange{
      .maxChange = 0.3, .lowerLimit = 0.01, .upperLimit = 0.5};
  MoveStatistics<double> conePositionChange{
      .maxChange = 0.3, .lowerLimit = 0.01, .upperLimit = 0.5};

  const std::string writeMCMoveStatistics() const;

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const CBMCMoveStatistics& p);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, CBMCMoveStatistics& p);
};
