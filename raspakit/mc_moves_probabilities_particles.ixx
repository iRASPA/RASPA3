export module mc_moves_probabilities_particles;

import <string>;
import <chrono>;
import <fstream>;

import double3;
import archive;
import move_statistics;


export struct MCMoveProbabilitiesParticles
{
  uint64_t versionNumber{ 1 };

  bool operator==(MCMoveProbabilitiesParticles const&) const = default;

  double probabilityTranslationMove{ 0.0 };
  double probabilityRandomTranslationMove{ 0.0 };
  double probabilityRotationMove{ 0.0 };
  double probabilityRandomRotationMove{ 0.0 };
  double probabilityVolumeMove{ 0.0 };
  double probabilityReinsertionMove_CBMC{ 0.0 };
  double probabilityIdentityChangeMove_CBMC{ 0.0 };
  double probabilitySwapMove_CBMC{ 0.0 };
  double probabilitySwapMove_CFCMC{ 0.0 };
  double probabilitySwapMove_CFCMC_CBMC{ 0.0 };
  double probabilityGibbsVolumeMove{ 0.0 };
  double probabilityGibbsSwapMove_CBMC{ 0.0 };
  double probabilityGibbsSwapMove_CFCMC{ 0.0 };
  double probabilityGibbsSwapMove_CFCMC_CBMC{ 0.0 };
  double probabilityWidomMove{ 0.0 };
  double probabilityWidomMove_CFCMC{ 0.0 };
  double probabilityWidomMove_CFCMC_CBMC{ 0.0 };

  double accumulatedProbabilityTranslationMove{ 0.0 };
  double accumulatedProbabilityRandomTranslationMove{ 0.0 };
  double accumulatedProbabilityRotationMove{ 0.0 };
  double accumulatedProbabilityRandomRotationMove{ 0.0 };
  double accumulatedProbabilityVolumeMove{ 0.0 };
  double accumulatedProbabilityReinsertionMove_CBMC{ 0.0 };
  double accumulatedProbabilityIdentityChangeMove_CBMC{ 0.0 };
  double accumulatedProbabilitySwapMove_CBMC{ 0.0 };
  double accumulatedProbabilitySwapMove_CFCMC{ 0.0 };
  double accumulatedProbabilitySwapMove_CFCMC_CBMC{ 0.0 };
  double accumulatedProbabilityGibbsVolumeMove{ 0.0 };
  double accumulatedProbabilityGibbsSwapMove_CBMC{ 0.0 };
  double accumulatedProbabilityGibbsSwapMove_CFCMC{ 0.0 };
  double accumulatedProbabilityGibbsSwapMove_CFCMC_CBMC{ 0.0 };
  double accumulatedProbabilityWidomMove{ 0.0 };
  double accumulatedProbabilityWidomMove_CFCMC{ 0.0 };
  double accumulatedProbabilityWidomMove_CFCMC_CBMC{ 0.0 };

  MoveStatistics<double3> statistics_TranslationMove{ .maxChange = double3(1.0,1.0,1.0) };
  MoveStatistics<double3> statistics_RandomTranslationMove{};
  MoveStatistics<double3> statistics_RotationMove{ .maxChange = double3(1.0,1.0,1.0) };
  MoveStatistics<double3> statistics_RandomRotationMove{};
  MoveStatistics<double> statistics_ReinsertionMove_CBMC{};
  MoveStatistics<double> statistics_IdentityChangeMove_CBMC{};
  MoveStatistics<double> statistics_SwapInsertionMove_CBMC{};
  MoveStatistics<double> statistics_SwapDeletionMove_CBMC{};
  MoveStatistics<double3> statistics_SwapMove_CFCMC{ .maxChange = double3(0.0,0.0,0.5) };
  MoveStatistics<double3> statistics_SwapMove_CFCMC_CBMC{ .maxChange = double3(0.0,0.0,0.5) };
  MoveStatistics<double3> statistics_WidomMove_CBMC{};
  MoveStatistics<double3> statistics_WidomMove_CFCMC{};
  MoveStatistics<double3> statistics_WidomMove_CFCMC_CBMC{};

  MoveStatistics<double> statistics_GibbsSwapMove_CBMC{};
  MoveStatistics<double3> statistics_GibbsSwapMove_CFCMC{ .maxChange = double3(0.0,0.0,0.5) };
  MoveStatistics<double3> statistics_GibbsSwapMove_CFCMC_CBMC{ .maxChange = double3(0.0,0.0,0.5)};

  void clearMoveStatistics();
  void optimizeMCMoves();

  void normalizeMoveProbabilties();
  const std::string writeMCMoveStatistics() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveProbabilitiesParticles &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveProbabilitiesParticles &p);
};
