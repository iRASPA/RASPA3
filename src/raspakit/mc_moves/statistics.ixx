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
  std::uint64_t versionNumber{1};

  bool operator==(MCMoveStatistics const&) const = default;

  std::array<std::variant<MoveStatistics<double>, MoveStatistics<double3>>, 
             std::to_underlying(Move::Types::Count)> stats{};

  const std::variant<MoveStatistics<double>, MoveStatistics<double3>> &operator[](Move::Types i) const
  {
    return stats[std::to_underlying(i)];
  }

  MCMoveStatistics()
  {
    stats[std::to_underlying(Move::Types::Translation)] =
        MoveStatistics<double3>{.maxChange = double3(1.0), .lowerLimit = double3(0.01), .upperLimit = double3(1.5)};
    stats[std::to_underlying(Move::Types::RandomTranslation)] = MoveStatistics<double3>{};
    stats[std::to_underlying(Move::Types::Rotation)] =
        MoveStatistics<double3>{.maxChange = double3(1.0), .lowerLimit = double3(0.01), .upperLimit = double3(1.5)};
    stats[std::to_underlying(Move::Types::RandomRotation)] = MoveStatistics<double3>{};
    stats[std::to_underlying(Move::Types::Swap)] = MoveStatistics<double3>{};
    stats[std::to_underlying(Move::Types::SwapCBMC)] = MoveStatistics<double3>{};
    stats[std::to_underlying(Move::Types::PairSwapCBMC)] = MoveStatistics<double3>{};
    stats[std::to_underlying(Move::Types::PairSwap)] = MoveStatistics<double3>{};
    stats[std::to_underlying(Move::Types::GroupSwap)] = MoveStatistics<double3>{};
    stats[std::to_underlying(Move::Types::GroupSwapCBMC)] = MoveStatistics<double3>{};
    stats[std::to_underlying(Move::Types::GroupSwapCFCMC)] = MoveStatistics<double3>{
        .maxChange = double3(0.0, 0.0, 0.5), .lowerLimit = double3(0.1), .upperLimit = double3(1.0)};
    stats[std::to_underlying(Move::Types::GroupSwapCBCFCMC)] = MoveStatistics<double3>{
        .maxChange = double3(0.0, 0.0, 0.5), .lowerLimit = double3(0.1), .upperLimit = double3(1.0)};
    stats[std::to_underlying(Move::Types::PairSwapCFCMC)] = MoveStatistics<double3>{
        .maxChange = double3(0.0, 0.0, 0.5), .lowerLimit = double3(0.1), .upperLimit = double3(1.0)};
    stats[std::to_underlying(Move::Types::PairSwapCBCFCMC)] = MoveStatistics<double3>{
        .maxChange = double3(0.0, 0.0, 0.5), .lowerLimit = double3(0.1), .upperLimit = double3(1.0)};
    stats[std::to_underlying(Move::Types::SwapCFCMC)] = MoveStatistics<double3>{
        .maxChange = double3(0.0, 0.0, 0.5), .lowerLimit = double3(0.1), .upperLimit = double3(1.0)};
    stats[std::to_underlying(Move::Types::SwapCBCFCMC)] = MoveStatistics<double3>{
        .maxChange = double3(0.0, 0.0, 0.5), .lowerLimit = double3(0.1), .upperLimit = double3(1.0)};
    stats[std::to_underlying(Move::Types::GibbsSwapCFCMC)] = MoveStatistics<double3>{
        .maxChange = double3(0.0, 0.0, 0.5), .lowerLimit = double3(0.1), .upperLimit = double3(1.0)};
    stats[std::to_underlying(Move::Types::GibbsSwapCBCFCMC)] = MoveStatistics<double3>{
        .maxChange = double3(0.0, 0.0, 0.5), .lowerLimit = double3(0.1), .upperLimit = double3(1.0)};
    stats[std::to_underlying(Move::Types::GibbsConventionalCFCMC)] = MoveStatistics<double3>{
        .maxChange = double3(0.0, 0.0, 0.5), .lowerLimit = double3(0.1), .upperLimit = double3(1.0)};
    stats[std::to_underlying(Move::Types::GibbsConventionalCBCFCMC)] = MoveStatistics<double3>{
        .maxChange = double3(0.0, 0.0, 0.5), .lowerLimit = double3(0.1), .upperLimit = double3(1.0)};

    stats[std::to_underlying(Move::Types::VolumeChange)] =
        MoveStatistics<double>{.maxChange = 0.1, .lowerLimit = 0.01, .upperLimit = 1.5};
    stats[std::to_underlying(Move::Types::AnisotropicVolumeChange)] =
        MoveStatistics<double3>{.maxChange = double3(0.1), .lowerLimit = double3(0.01), .upperLimit = double3(1.5)};
    stats[std::to_underlying(Move::Types::ReinsertionCBMC)] = MoveStatistics<double>{};
    stats[std::to_underlying(Move::Types::PartialReinsertionCBMC)] = MoveStatistics<double>{};
    stats[std::to_underlying(Move::Types::IdentityChangeCBMC)] = MoveStatistics<double>{};
    stats[std::to_underlying(Move::Types::GibbsVolume)] =
        MoveStatistics<double>{.maxChange = 0.1, .lowerLimit = 0.01, .upperLimit = 1.5};
    stats[std::to_underlying(Move::Types::GibbsSwapCBMC)] = MoveStatistics<double>{};
    stats[std::to_underlying(Move::Types::GibbsIdentityChangeCBMC)] = MoveStatistics<double>{};
    stats[std::to_underlying(Move::Types::Widom)] = MoveStatistics<double>{};
    stats[std::to_underlying(Move::Types::WidomCFCMC)] = MoveStatistics<double>{};
    stats[std::to_underlying(Move::Types::WidomCBCFCMC)] = MoveStatistics<double>{};
    stats[std::to_underlying(Move::Types::ParallelTempering)] = MoveStatistics<double>{};
    stats[std::to_underlying(Move::Types::HybridMC)] =
        MoveStatistics<double>{.maxChange = 0.0005, .lowerLimit = 0.000001, .upperLimit = 0.01};
    // Translation smart MC: 'maxChange' is the standard deviation (in Angstrom) of the
    // Gaussian part of the trial displacement; the drift is derived from this via beta * sigma^2 / 2.
    stats[std::to_underlying(Move::Types::TranslationSmartMC)] =
        MoveStatistics<double>{.maxChange = 0.1, .lowerLimit = 0.001, .upperLimit = 1.5};
    stats[std::to_underlying(Move::Types::TranslationSmartMCAll)] =
        MoveStatistics<double>{.maxChange = 0.05, .lowerLimit = 0.001, .upperLimit = 1.5};
    // Rotation smart MC: 'maxChange' is the angular standard deviation (radians); the torque drift
    // uses the same b = beta * sigma^2 / 2 relation, applied as a quaternion rotation.
    stats[std::to_underlying(Move::Types::RotationSmartMC)] =
        MoveStatistics<double>{.maxChange = 0.1, .lowerLimit = 0.001, .upperLimit = 1.5};
    stats[std::to_underlying(Move::Types::RotationSmartMCAll)] =
        MoveStatistics<double>{.maxChange = 0.05, .lowerLimit = 0.001, .upperLimit = 1.5};
    // Combined translation-rotation smart MC: component 0 is the translational standard deviation
    // (Angstrom), component 1 the angular standard deviation (radians); component 2 is unused.
    // Each drift follows b = beta * sigma^2 / 2 with its own sigma.
    stats[std::to_underlying(Move::Types::TranslationRotationSmartMC)] = MoveStatistics<double3>{
        .maxChange = double3(0.1, 0.1, 0.0), .lowerLimit = double3(0.001, 0.001, 0.0),
        .upperLimit = double3(1.5, 1.5, 0.0)};
    stats[std::to_underlying(Move::Types::ReactionCBMC)] = MoveStatistics<double>{};
    stats[std::to_underlying(Move::Types::ReactionConventionalCFCMC)] = MoveStatistics<double>{};
    stats[std::to_underlying(Move::Types::ReactionConventionalCBCFCMC)] = MoveStatistics<double>{};
    stats[std::to_underlying(Move::Types::ReactionCFCMC)] = MoveStatistics<double>{};
    stats[std::to_underlying(Move::Types::ReactionCBCFCMC)] = MoveStatistics<double>{};
  };

  void clearMoveStatistics();
  void optimizeMCMoves();

  void addAllCounts(const Move::Types& move);
  void addTrial(const Move::Types& move);
  void addTrial(const Move::Types& move, std::size_t direction);
  void addConstructed(const Move::Types& move);
  void addConstructed(const Move::Types& move, std::size_t direction);
  void addAccepted(const Move::Types& move);
  void addAccepted(const Move::Types& move, std::size_t direction);

  double getMaxChange(const Move::Types& move)
  {
    return std::get<MoveStatistics<double>>(stats[std::to_underlying(move)]).maxChange;
  }

  double getMaxChange(const Move::Types& move, std::size_t direction)
  {
    return std::get<MoveStatistics<double3>>(stats[std::to_underlying(move)]).maxChange[direction];
  }

  void setMaxChange(const Move::Types& move, double value) 
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


  std::string repr() const;

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
