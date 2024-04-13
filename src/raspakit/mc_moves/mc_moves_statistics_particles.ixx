module;

#ifdef USE_LEGACY_HEADERS
#include <string>
#include <chrono>
#include <fstream>
#endif

export module mc_moves_statistics_particles;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <chrono>;
import <fstream>;
#endif

import double3;
import archive;
import move_statistics;

export struct MCMoveStatisticsParticles
{
  uint64_t versionNumber{ 1 };

  bool operator==(MCMoveStatisticsParticles const&) const = default;

  MoveStatistics<double3> translationMove{ .maxChange = double3(1.0,1.0,1.0) };
  MoveStatistics<double3> randomTranslationMove{};
  MoveStatistics<double3> rotationMove{ .maxChange = double3(1.0,1.0,1.0) };
  MoveStatistics<double3> randomRotationMove{};
  MoveStatistics<double> reinsertionMove_CBMC{};
  MoveStatistics<double> identityChangeMove_CBMC{};
  MoveStatistics<double> swapInsertionMove{};
  MoveStatistics<double> swapDeletionMove{};
  MoveStatistics<double> swapInsertionMove_CBMC{};
  MoveStatistics<double> swapDeletionMove_CBMC{};
  MoveStatistics<double3> swapMove_CFCMC{ .maxChange = double3(0.0,0.0,0.5) };
  MoveStatistics<double3> swapMove_CFCMC_CBMC{ .maxChange = double3(0.0,0.0,0.5) };
  MoveStatistics<double3> WidomMove_CBMC{};
  MoveStatistics<double3> WidomMove_CFCMC{};
  MoveStatistics<double3> WidomMove_CFCMC_CBMC{};

  MoveStatistics<double> GibbsSwapMove_CBMC{};
  MoveStatistics<double3> GibbsSwapMove_CFCMC{ .maxChange = double3(0.0,0.0,0.5) };
  MoveStatistics<double3> GibbsSwapMove_CFCMC_CBMC{ .maxChange = double3(0.0,0.0,0.5)};

  void clearMoveStatistics();
  void optimizeMCMoves();
  const std::string writeMCMoveStatistics() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveStatisticsParticles &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveStatisticsParticles &p);
};