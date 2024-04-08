module;

#ifdef USE_LEGACY_HEADERS
#include <string>
#include <chrono>
#include <fstream>
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

export struct MCMoveStatisticsSystem
{
  uint64_t versionNumber{ 1 };

  bool operator==(MCMoveStatisticsSystem const&) const = default;

  MoveStatistics<double> volumeMove{ .maxChange = 0.1 };
  MoveStatistics<double> GibbsVolumeMove{ .maxChange = 0.1 };

  void optimizeAcceptance();
  void clear();

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveStatisticsSystem &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveStatisticsSystem &p);
};
