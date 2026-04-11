module;

export module mc_moves_statistics;

import std;

import archive;
import double3;
import json;
import move_statistics;
import mc_moves_move_types;


export struct MCMoveStatistics
{
  std::uint64_t versionNumber{2};

  bool operator==(MCMoveStatistics const&) const = default;

  std::array<std::variant<MoveStatistics<double>, MoveStatistics<double3>>, 
             std::to_underlying(MoveTypes::Count)> stats{};

  MCMoveStatistics()
  {
    stats[std::to_underlying(MoveTypes::Translation)] =
        MoveStatistics<double3>{.maxChange = double3(1.0), .lowerLimit = double3(0.01), .upperLimit = double3(1.5)};
    stats[std::to_underlying(MoveTypes::RandomTranslation)] = MoveStatistics<double3>{};
    stats[std::to_underlying(MoveTypes::Rotation)] =
        MoveStatistics<double3>{.maxChange = double3(1.0), .lowerLimit = double3(0.01), .upperLimit = double3(1.5)};
    stats[std::to_underlying(MoveTypes::RandomRotation)] = MoveStatistics<double3>{};
    stats[std::to_underlying(MoveTypes::Swap)] = MoveStatistics<double3>{};
    stats[std::to_underlying(MoveTypes::SwapCBMC)] = MoveStatistics<double3>{};
    stats[std::to_underlying(MoveTypes::SwapCFCMC)] = MoveStatistics<double3>{
        .maxChange = double3(0.0, 0.0, 0.5), .lowerLimit = double3(0.1), .upperLimit = double3(1.0)};
    stats[std::to_underlying(MoveTypes::SwapCBCFCMC)] = MoveStatistics<double3>{
        .maxChange = double3(0.0, 0.0, 0.5), .lowerLimit = double3(0.1), .upperLimit = double3(1.0)};
    stats[std::to_underlying(MoveTypes::GibbsSwapCFCMC)] = MoveStatistics<double3>{
        .maxChange = double3(0.0, 0.0, 0.5), .lowerLimit = double3(0.1), .upperLimit = double3(1.0)};

    stats[std::to_underlying(MoveTypes::VolumeChange)] =
        MoveStatistics<double>{.maxChange = 0.1, .lowerLimit = 0.01, .upperLimit = 1.5};
    stats[std::to_underlying(MoveTypes::ReinsertionCBMC)] = MoveStatistics<double>{};
    stats[std::to_underlying(MoveTypes::PartialReinsertionCBMC)] = MoveStatistics<double>{};
    stats[std::to_underlying(MoveTypes::IdentityChangeCBMC)] = MoveStatistics<double>{};
    stats[std::to_underlying(MoveTypes::GibbsVolume)] =
        MoveStatistics<double>{.maxChange = 0.1, .lowerLimit = 0.01, .upperLimit = 1.5};
    stats[std::to_underlying(MoveTypes::GibbsSwapCBMC)] = MoveStatistics<double>{};
    stats[std::to_underlying(MoveTypes::Widom)] = MoveStatistics<double>{};
    stats[std::to_underlying(MoveTypes::WidomCFCMC)] = MoveStatistics<double>{};
    stats[std::to_underlying(MoveTypes::WidomCBCFCMC)] = MoveStatistics<double>{};
    stats[std::to_underlying(MoveTypes::ParallelTempering)] = MoveStatistics<double>{};
    stats[std::to_underlying(MoveTypes::HybridMC)] =
        MoveStatistics<double>{.maxChange = 0.0005, .lowerLimit = 0.000001, .upperLimit = 0.01};
  };

  void clearMoveStatistics();
  void optimizeMCMoves();

  void addAllCounts(const MoveTypes& move);
  void addTrial(const MoveTypes& move);
  void addTrial(const MoveTypes& move, std::size_t direction);
  void addConstructed(const MoveTypes& move);
  void addConstructed(const MoveTypes& move, std::size_t direction);
  void addAccepted(const MoveTypes& move);
  void addAccepted(const MoveTypes& move, std::size_t direction);

  double getMaxChange(const MoveTypes& move)
  {
    return std::get<MoveStatistics<double>>(stats[std::to_underlying(move)]).maxChange;
  }

  double getMaxChange(const MoveTypes& move, std::size_t direction)
  {
    return std::get<MoveStatistics<double3>>(stats[std::to_underlying(move)]).maxChange[direction];
  }

  void setMaxChange(const MoveTypes& move, double value) 
  {
    std::get<MoveStatistics<double>>(stats[std::to_underlying(move)]).maxChange = value;
  };

  const std::string writeMCMoveStatistics() const;
  const std::string writeMCMoveStatistics(std::size_t countTotal) const;
  const nlohmann::json jsonMCMoveStatistics() const;

  inline MCMoveStatistics& operator+=(const MCMoveStatistics& b)
  {
    for (std::size_t i=0; i != stats.size(); ++i)
    {
      if(std::holds_alternative<MoveStatistics<double>>(stats[i]))
      {
        std::get<MoveStatistics<double>>(this->stats[i]) += std::get<MoveStatistics<double>>(b.stats[i]);
      }
      if(std::holds_alternative<MoveStatistics<double3>>(stats[i]))
      {
        std::get<MoveStatistics<double3>>(this->stats[i]) += std::get<MoveStatistics<double3>>(b.stats[i]);
      }
    }
    return *this;
  }

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const MCMoveStatistics& p);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, MCMoveStatistics& p);
};

export inline MCMoveStatistics operator+(const MCMoveStatistics& a, const MCMoveStatistics& b)
{
  MCMoveStatistics c;
  for (std::size_t i=0; i != a.stats.size(); ++i)
  {
    if(std::holds_alternative<MoveStatistics<double>>(a.stats[i]))
    {
      std::get<MoveStatistics<double>>(c.stats[i]) = std::get<MoveStatistics<double>>(a.stats[i]) + std::get<MoveStatistics<double>>(b.stats[i]);
    }
    if(std::holds_alternative<MoveStatistics<double3>>(a.stats[i]))
    {
      std::get<MoveStatistics<double3>>(c.stats[i]) = std::get<MoveStatistics<double3>>(a.stats[i]) + std::get<MoveStatistics<double3>>(b.stats[i]);
    }
  }
  return c;
}
