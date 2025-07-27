module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <exception>
#include <format>
#include <fstream>
#include <map>
#include <ostream>
#include <print>
#include <source_location>
#include <sstream>
#include <string>
#include <vector>

#endif

module mc_moves_statistics;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import mc_moves_move_types;

void MCMoveStatistics::clearMoveStatistics()
{
  for (auto& [moveType, statistics] : statsMapDouble)
  {
    statistics.clear();
  }
  for (auto& [moveType, statistics] : statsMapDouble3)
  {
    statistics.clear();
  }
}

void MCMoveStatistics::optimizeMCMoves()
{
  for (auto& [moveType, statistics] : statsMapDouble)
  {
    statistics.optimizeAcceptance();
  }
  for (auto& [moveType, statistics] : statsMapDouble3)
  {
    statistics.optimizeAcceptance();
  }
}

void MCMoveStatistics::addAllCounts(const MoveTypes& move)
{
  if (statsMapDouble.find(move) != statsMapDouble.end())
  {
    statsMapDouble[move].allCounts += 1;
  }
  else
  {
    statsMapDouble3[move].allCounts += 1;
  }
}

void MCMoveStatistics::addTrial(const MoveTypes& move)
{
  statsMapDouble[move].counts += 1;
  statsMapDouble[move].totalCounts += 1;
}

void MCMoveStatistics::addTrial(const MoveTypes& move, std::size_t direction)
{
  statsMapDouble3[move].counts[direction] += 1;
  statsMapDouble3[move].totalCounts[direction] += 1;
}

void MCMoveStatistics::addConstructed(const MoveTypes& move)
{
  statsMapDouble[move].constructed += 1;
  statsMapDouble[move].totalConstructed += 1;
}

void MCMoveStatistics::addConstructed(const MoveTypes& move, std::size_t direction)
{
  statsMapDouble3[move].constructed[direction] += 1;
  statsMapDouble3[move].totalConstructed[direction] += 1;
}

void MCMoveStatistics::addAccepted(const MoveTypes& move)
{
  statsMapDouble[move].accepted += 1;
  statsMapDouble[move].totalAccepted += 1;
}

void MCMoveStatistics::addAccepted(const MoveTypes& move, std::size_t direction)
{
  statsMapDouble3[move].accepted[direction] += 1;
  statsMapDouble3[move].totalAccepted[direction] += 1;
}

double MCMoveStatistics::getMaxChange(const MoveTypes& move) { return statsMapDouble[move].maxChange; }

void MCMoveStatistics::setMaxChange(const MoveTypes& move, double value) { statsMapDouble[move].maxChange = value; }

double MCMoveStatistics::getMaxChange(const MoveTypes& move, std::size_t direction)
{
  return statsMapDouble3[move].maxChange[direction];
}

static std::string formatStatistics(const std::string name, const MoveStatistics<double>& move)
{
  std::ostringstream stream;
  std::print(stream, "    {:20} all:          {:10}\n", name, move.allCounts);
  std::print(stream, "    {:20} total:        {:10}\n", name, move.totalCounts);
  std::print(stream, "    {:20} constructed:  {:10}\n", name, move.totalConstructed);
  std::print(stream, "    {:20} accepted:     {:10}\n", name, move.totalAccepted);
  std::print(stream, "    {:20} fraction:     {:10f}\n", name,
             move.totalAccepted / std::max(1.0, double(move.totalCounts)));
  std::print(stream, "    {:20} max-change:   {:10f}\n\n", name, move.maxChange);
  return stream.str();
}

static std::string formatStatistics(const std::string name, const MoveStatistics<double3>& move)
{
  std::ostringstream stream;
  std::print(stream, "    {:20} all:          {:10}\n", name, move.allCounts);
  std::print(stream, "    {:20} total:        {:10} {:10} {:10}\n", name, move.totalCounts.x, move.totalCounts.y,
             move.totalCounts.z);
  std::print(stream, "    {:20} constructed:  {:10} {:10} {:10}\n", name, move.totalConstructed.x,
             move.totalConstructed.y, move.totalConstructed.z);
  std::print(stream, "    {:20} accepted:     {:10} {:10} {:10}\n", name, move.totalAccepted.x, move.totalAccepted.y,
             move.totalAccepted.z);
  std::print(stream, "    {:20} fraction:     {:10f} {:10f} {:10f}\n", name,
             move.totalAccepted.x / std::max(1.0, double(move.totalCounts.x)),
             move.totalAccepted.y / std::max(1.0, double(move.totalCounts.y)),
             move.totalAccepted.z / std::max(1.0, double(move.totalCounts.z)));
  std::print(stream, "    {:20} max-change:   {:10f} {:10f} {:10f}\n\n", name, move.maxChange.x, move.maxChange.y,
             move.maxChange.z);
  return stream.str();
}

static nlohmann::json jsonStatistics(const MoveStatistics<double>& move)
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

static nlohmann::json jsonStatistics(const MoveStatistics<double3>& move)
{
  nlohmann::json status;
  status["all"] = move.allCounts;
  status["total"] = {move.totalCounts.x, move.totalCounts.y, move.totalCounts.z};
  status["constructed"] = {move.totalConstructed.x, move.totalConstructed.y, move.totalConstructed.z};
  status["accepted"] = {move.totalAccepted.x, move.totalAccepted.y, move.totalAccepted.z};
  status["fraction"] = {move.totalAccepted.x / std::max(1.0, double(move.totalCounts.x)),
                        move.totalAccepted.y / std::max(1.0, double(move.totalCounts.y)),
                        move.totalAccepted.z / std::max(1.0, double(move.totalCounts.z))};
  status["maxChange"] = {move.maxChange.x, move.maxChange.y, move.maxChange.z};
  return status;
}

const std::string MCMoveStatistics::writeMCMoveStatistics() const
{
  std::ostringstream stream;
  for (auto& [moveType, statistics] : statsMapDouble)
  {
    if (statistics.totalCounts > 0.0)
    {
      std::print(stream, "{}", formatStatistics(moveNames[moveType], statistics));
    }
  }
  for (auto& [moveType, statistics] : statsMapDouble3)
  {
    if (statistics.totalCounts.x > 0.0)
    {
      std::print(stream, "{}", formatStatistics(moveNames[moveType], statistics));
    }
  }
  return stream.str();
}

const std::string MCMoveStatistics::writeMCMoveStatistics(std::size_t countTotal) const
{
  std::ostringstream stream;

  std::size_t summed = 0;
  for (auto& [moveType, statistics] : statsMapDouble)
  {
    double moveCount = static_cast<double>(statistics.allCounts);
    if (moveCount > 0.0)
    {
      std::print(stream, "{:<29}{:14} ({:<6.4f} [%])\n", moveNames[moveType], moveCount,
                 100.0 * moveCount / static_cast<double>(countTotal));
      summed += statistics.allCounts;
    }
  }
  for (auto& [moveType, statistics] : statsMapDouble3)
  {
    double moveCount = static_cast<double>(statistics.allCounts);
    if (moveCount > 0.0)
    {
      std::print(stream, "{:<29}{:14} ({:<6.4f} [%])\n", moveNames[moveType], moveCount,
                 100.0 * moveCount / static_cast<double>(countTotal));
      summed += statistics.allCounts;
    }
  }

  std::print(stream, "\n");
  std::print(stream, "Production count MC-steps:   {:14d} [-]\n", countTotal);
  std::print(stream, "               All summed:   {:14d} [-]\n", summed);
  std::print(stream, "               difference:   {:14d} [-]\n", countTotal - summed);

  return stream.str();
}

const nlohmann::json MCMoveStatistics::jsonMCMoveStatistics() const
{
  nlohmann::json status;
  for (auto& [moveType, statistics] : statsMapDouble)
  {
    if (statistics.totalCounts > 0.0)
    {
      status[moveNames[moveType]] = jsonStatistics(statistics);
    }
  }
  for (auto& [moveType, statistics] : statsMapDouble3)
  {
    if (statistics.totalCounts.x > 0.0)
    {
      status[moveNames[moveType]] = jsonStatistics(statistics);
    }
  }
  return status;
}

Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const MCMoveStatistics& p)
{
  archive << p.versionNumber;
  archive << p.statsMapDouble;
  archive << p.statsMapDouble3;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, MCMoveStatistics& p)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > p.versionNumber)
  {
    const std::source_location& location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'MCMoveProbabilitiesSystem' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }
  archive >> p.statsMapDouble;
  archive >> p.statsMapDouble3;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("MCMoveStatistics: Error in binary restart\n"));
  }
#endif

  return archive;
}
