export module mc_moves_probabilities_particles;

import <string>;
import <chrono>;

import double3;
import move_statistics;

export struct MCMoveProbabilitiesParticles
{
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
  MoveStatistics<double3> statistics_SwapMove_CFCMC{};
  MoveStatistics<double3> statistics_SwapMove_CFCMC_CBMC{};
  MoveStatistics<double3> statistics_WidomMove_CBMC{};
  MoveStatistics<double3> statistics_WidomMove_CFCMC{};
  MoveStatistics<double3> statistics_WidomMove_CFCMC_CBMC{};

  void clearMoveStatistics();
  void optimizeMCMoves();

  void normalizeMoveProbabilties();
  const std::string writeMCMoveStatistics() const;

  std::chrono::duration<double> cpuTime_TranslationMove{ 0.0 };
  std::chrono::duration<double> cpuTime_RandomTranslationMove{ 0.0 };
  std::chrono::duration<double> cpuTime_RotationMove{ 0.0 };
  std::chrono::duration<double> cpuTime_RandomRotationMove{ 0.0 };
  std::chrono::duration<double> cpuTime_ReinsertionMove_CBMC{ 0.0 };
  std::chrono::duration<double> cpuTime_IdentityChangeMove_CBMC{ 0.0 };
  std::chrono::duration<double> cpuTime_SwapInsertionMove_CBMC{ 0.0 };
  std::chrono::duration<double> cpuTime_SwapDeletionMove_CBMC{ 0.0 };
  std::chrono::duration<double> cpuTime_SwapMove_CFCMC{ 0.0 };
  std::chrono::duration<double> cpuTime_SwapMove_CFCMC_CBMC{ 0.0 };
  std::chrono::duration<double> cpuTime_WidomMove_CBMC{ 0.0 };
  std::chrono::duration<double> cpuTime_WidomMove_CFCMC{ 0.0 };
  std::chrono::duration<double> cpuTime_WidomMove_CFCMC_CBMC{ 0.0 };

  std::chrono::duration<double> cpuTime_TranslationMove_NonEwald{ 0.0 };
  std::chrono::duration<double> cpuTime_RandomTranslationMove_NonEwald{ 0.0 };
  std::chrono::duration<double> cpuTime_RotationMove_NonEwald{ 0.0 };
  std::chrono::duration<double> cpuTime_RandomRotationMove_NonEwald{ 0.0 };
  std::chrono::duration<double> cpuTime_ReinsertionGrowMove_CBMC_NonEwald{ 0.0 };
  std::chrono::duration<double> cpuTime_ReinsertionRetraceMove_CBMC_NonEwald{ 0.0 };
  std::chrono::duration<double> cpuTime_IdentityChangeMove_CBMC_NonEwald{ 0.0 };
  std::chrono::duration<double> cpuTime_SwapInsertionMove_CBMC_NonEwald{ 0.0 };
  std::chrono::duration<double> cpuTime_SwapDeletionMove_CBMC_NonEwald{ 0.0 };
  std::chrono::duration<double> cpuTime_SwapMove_CFCMC_NonEwald{ 0.0 };
  std::chrono::duration<double> cpuTime_SwapInsertionMove_CFCMC_CBMC_NonEwald{ 0.0 };
  std::chrono::duration<double> cpuTime_SwapInsertionGrowMove_CFCMC_CBMC_NonEwald{ 0.0 };
  std::chrono::duration<double> cpuTime_SwapDeletionMove_CFCMC_CBMC_NonEwald{ 0.0 };
  std::chrono::duration<double> cpuTime_SwapDeletionRetraceMove_CFCMC_CBMC_NonEwald{ 0.0 };
  std::chrono::duration<double> cpuTime_SwapLambdaMove_CFCMC_CBMC_NonEwald{ 0.0 };
  std::chrono::duration<double> cpuTime_WidomMove_CBMC_NonEwald{ 0.0 };
  std::chrono::duration<double> cpuTime_WidomMove_CFCMC_NonEwald{ 0.0 };
  std::chrono::duration<double> cpuTime_WidomMove_CFCMC_CBMC_NonEwald{ 0.0 };

  std::chrono::duration<double> cpuTime_TranslationMove_Ewald{ 0.0 };
  std::chrono::duration<double> cpuTime_RandomTranslationMove_Ewald{ 0.0 };
  std::chrono::duration<double> cpuTime_RotationMove_Ewald{ 0.0 };
  std::chrono::duration<double> cpuTime_RandomRotationMove_Ewald{ 0.0 };
  std::chrono::duration<double> cpuTime_ReinsertionMove_CBMC_Ewald{ 0.0 };
  std::chrono::duration<double> cpuTime_IdentityChangeMove_CBMC_Ewald{ 0.0 };
  std::chrono::duration<double> cpuTime_SwapInsertionMove_CBMC_Ewald{ 0.0 };
  std::chrono::duration<double> cpuTime_SwapDeletionMove_CBMC_Ewald{ 0.0 };
  std::chrono::duration<double> cpuTime_SwapMove_CFCMC_Ewald{ 0.0 };
  std::chrono::duration<double> cpuTime_SwapInsertionMove_CFCMC_CBMC_Ewald{ 0.0 };
  std::chrono::duration<double> cpuTime_SwapInsertionGrowMove_CFCMC_CBMC_Ewald{ 0.0 };
  std::chrono::duration<double> cpuTime_SwapDeletionMove_CFCMC_CBMC_Ewald{ 0.0 };
  std::chrono::duration<double> cpuTime_SwapDeletionRetraceMove_CFCMC_CBMC_Ewald{ 0.0 };
  std::chrono::duration<double> cpuTime_SwapLambdaMove_CFCMC_CBMC_Ewald{ 0.0 };
  std::chrono::duration<double> cpuTime_WidomMove_CBMC_Ewald{ 0.0 };
  std::chrono::duration<double> cpuTime_WidomMove_CFCMC_Ewald{ 0.0 };
  std::chrono::duration<double> cpuTime_WidomMove_CFCMC_CBMC_Ewald{ 0.0 };

  void clearTimingStatistics();
  const std::string writeMCMoveCPUTimeStatistics() const;
};
