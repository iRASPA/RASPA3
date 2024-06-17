module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <fstream>
#include <string>
#endif

export module mc_moves_statistics_system;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <chrono>;
import <fstream>;
#endif

import archive;
import double3;
import move_statistics;
import json;

export struct MCMoveStatisticsSystem
{
  uint64_t versionNumber{1};

  bool operator==(MCMoveStatisticsSystem const &) const = default;

  MoveStatistics<double> volumeMove{.maxChange = 0.1};
  MoveStatistics<double> GibbsVolumeMove{.maxChange = 0.1};
  MoveStatistics<double> ParallelTemperingSwap{.maxChange = 0.1};

  void optimizeAcceptance();
  void clear();
  const std::string writeMCMoveStatistics() const;
  const nlohmann::json jsonMCMoveStatistics() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveStatisticsSystem &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveStatisticsSystem &p);
};
