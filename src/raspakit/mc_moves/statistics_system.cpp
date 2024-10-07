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
#include <vector>
#endif

module mc_moves_statistics_system;

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
import <algorithm>;
import <print>;
#endif

import archive;
import double3;
import move_statistics;
import stringutils;
import json;

void MCMoveStatisticsSystem::clear()
{
  volumeMove.clear();
  GibbsVolumeMove.clear();
  ParallelTemperingSwap.clear();
}

void MCMoveStatisticsSystem::optimizeAcceptance()
{
  volumeMove.optimizeAcceptance(0.01, 1.5);
  GibbsVolumeMove.optimizeAcceptance(0.01, 1.5);
  ParallelTemperingSwap.optimizeAcceptance(0.01, 1.5);
}

static std::string formatStatistics(const std::string name, const MoveStatistics<double> &move)
{
  std::ostringstream stream;
  std::print(stream, "    {} total:        {:10}\n", name, move.totalCounts);
  std::print(stream, "    {} constructed:  {:10}\n", name, move.totalConstructed);
  std::print(stream, "    {} accepted:     {:10}\n", name, move.totalAccepted);
  std::print(stream, "    {} fraction:     {:10f}\n", name,
             move.totalAccepted / std::max(1.0, double(move.totalCounts)));
  std::print(stream, "    {} max-change:   {:10f}\n\n", name, move.maxChange);
  return stream.str();
}

static nlohmann::json jsonStatistics(const MoveStatistics<double> &move)
{
  nlohmann::json status;
  status["all"] = move.allCounts;
  status["total"] = move.totalCounts;
  status["constructed"] = move.totalConstructed;
  status["accepted"] = move.totalAccepted;
  status["fraction"] = move.totalAccepted / std::max(1.0, double(move.totalCounts));
  status["maxChange"] = move.maxChange;
  return status;
}

const std::string MCMoveStatisticsSystem::writeMCMoveStatistics() const
{
  std::ostringstream stream;
  if (volumeMove.totalCounts > 0)
  {
    std::print(stream, "{}", formatStatistics("Volume Move", volumeMove));
  }
  if (GibbsVolumeMove.totalCounts > 0)
  {
    std::print(stream, "{}", formatStatistics("Gibbs Volume Move", GibbsVolumeMove));
  }
  if (ParallelTemperingSwap.totalCounts > 0)
  {
    std::print(stream, "{}", formatStatistics("Parallel Tempering Swap", ParallelTemperingSwap));
  }

  return stream.str();
}

const nlohmann::json MCMoveStatisticsSystem::jsonMCMoveStatistics() const
{
  nlohmann::json status;
  if (volumeMove.totalCounts > 0)
  {
    status["volumeMove"] = jsonStatistics(volumeMove);
  }
  if (GibbsVolumeMove.totalCounts > 0)
  {
    status["gibbsVolumeMove"] = jsonStatistics(GibbsVolumeMove);
  }
  if (ParallelTemperingSwap.totalCounts > 0)
  {
    status["parallelTemperingSwap"] = jsonStatistics(ParallelTemperingSwap);
  }

  return status;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveStatisticsSystem &p)
{
  archive << p.versionNumber;

  archive << p.volumeMove;
  archive << p.GibbsVolumeMove;
  archive << p.ParallelTemperingSwap;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveStatisticsSystem &p)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > p.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'MCMoveProbabilitiesSystem' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }

  archive >> p.volumeMove;
  archive >> p.GibbsVolumeMove;
  archive >> p.ParallelTemperingSwap;

  return archive;
}
