module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <format>
#include <map>
#include <print>
#include <string>
#include <variant>
#endif

export module mc_moves_statistics;

#ifndef USE_LEGACY_HEADERS
import <variant>;
import <string>;
import <map>;
import <print>;
import <format>;
import <string>;
#endif

import archive;
import double3;
import json;
import move_statistics;
import mc_moves_move_types;

export struct MCMoveStatistics
{
  uint64_t versionNumber{2};

  bool operator==(MCMoveStatistics const&) const = default;

  std::map<MoveTypes, MoveStatistics<double>> statsMapDouble;
  std::map<MoveTypes, MoveStatistics<double3>> statsMapDouble3;

  MCMoveStatistics()
  {
    statsMapDouble3[MoveTypes::Translation] =
        MoveStatistics<double3>{.maxChange = double3(1.0), .lowerLimit = double3(0.01), .upperLimit = double3(1.5)};
    statsMapDouble3[MoveTypes::RandomTranslation] = MoveStatistics<double3>{};
    statsMapDouble3[MoveTypes::Rotation] =
        MoveStatistics<double3>{.maxChange = double3(1.0), .lowerLimit = double3(0.01), .upperLimit = double3(1.5)};
    statsMapDouble3[MoveTypes::RandomRotation] = MoveStatistics<double3>{};
    statsMapDouble3[MoveTypes::Swap] = MoveStatistics<double3>{};
    statsMapDouble3[MoveTypes::SwapCBMC] = MoveStatistics<double3>{};
    statsMapDouble3[MoveTypes::SwapCFCMC] = MoveStatistics<double3>{
        .maxChange = double3(0.0, 0.0, 0.5), .lowerLimit = double3(0.0), .upperLimit = double3(1.0)};
    statsMapDouble3[MoveTypes::SwapCBCFCMC] = MoveStatistics<double3>{
        .maxChange = double3(0.0, 0.0, 0.5), .lowerLimit = double3(0.0), .upperLimit = double3(1.0)};
    statsMapDouble3[MoveTypes::GibbsSwapCFCMC] = MoveStatistics<double3>{
        .maxChange = double3(0.0, 0.0, 0.5), .lowerLimit = double3(0.0), .upperLimit = double3(1.0)};

    statsMapDouble[MoveTypes::VolumeChange] =
        MoveStatistics<double>{.maxChange = 0.1, .lowerLimit = 0.01, .upperLimit = 1.5};
    statsMapDouble[MoveTypes::ReinsertionCBMC] = MoveStatistics<double>{};
    statsMapDouble[MoveTypes::IdentityChangeCBMC] = MoveStatistics<double>{};
    statsMapDouble[MoveTypes::GibbsVolume] =
        MoveStatistics<double>{.maxChange = 0.1, .lowerLimit = 0.01, .upperLimit = 1.5};
    statsMapDouble[MoveTypes::GibbsSwapCBMC] = MoveStatistics<double>{};
    statsMapDouble[MoveTypes::Widom] = MoveStatistics<double>{};
    statsMapDouble[MoveTypes::WidomCFCMC] = MoveStatistics<double>{};
    statsMapDouble[MoveTypes::WidomCBCFCMC] = MoveStatistics<double>{};
    statsMapDouble[MoveTypes::ParallelTempering] = MoveStatistics<double>{};
    statsMapDouble[MoveTypes::HybridMC] =
        MoveStatistics<double>{.maxChange = 0.0005, .lowerLimit = 0.000001, .upperLimit = 0.01};
  };

  void clearMoveStatistics();
  void optimizeMCMoves();

  void addAllCounts(const MoveTypes& move);
  void addTrial(const MoveTypes& move);
  void addTrial(const MoveTypes& move, size_t direction);
  void addConstructed(const MoveTypes& move);
  void addConstructed(const MoveTypes& move, size_t direction);
  void addAccepted(const MoveTypes& move);
  void addAccepted(const MoveTypes& move, size_t direction);
  double getMaxChange(const MoveTypes& move);
  double getMaxChange(const MoveTypes& move, size_t direction);
  void setMaxChange(const MoveTypes& move, double value);

  const std::string writeMCMoveStatistics() const;
  const std::string writeMCMoveStatistics(size_t countTotal) const;
  const nlohmann::json jsonMCMoveStatistics() const;

  inline MCMoveStatistics& operator+=(const MCMoveStatistics& b)
  {
    for (auto& [moveType, statistics] : statsMapDouble)
    {
      statistics += b.statsMapDouble.at(moveType);
    }
    for (auto& [moveType, statistics] : statsMapDouble3)
    {
      statistics += b.statsMapDouble3.at(moveType);
    }
    return *this;
  }

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const MCMoveStatistics& p);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, MCMoveStatistics& p);
};

export inline MCMoveStatistics operator+(const MCMoveStatistics& a, const MCMoveStatistics& b)
{
  MCMoveStatistics c;
  for (auto& [moveType, statistics] : a.statsMapDouble)
  {
    c.statsMapDouble[moveType] = a.statsMapDouble.at(moveType) + b.statsMapDouble.at(moveType);
  }
  for (auto& [moveType, statistics] : a.statsMapDouble3)
  {
    c.statsMapDouble3[moveType] = a.statsMapDouble3.at(moveType) + b.statsMapDouble3.at(moveType);
  }
  return c;
}
