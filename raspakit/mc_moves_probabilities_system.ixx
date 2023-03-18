export module mc_moves_probabilities_system;

import <string>;
import <chrono>;

import double3;
import move_statistics;

export struct MCMoveProbabilitiesSystem
{
  double probabilityVolumeMove{ 0.0 };
  double probabilityGibbsVolumeMove{ 0.0 };
  double probabilityGibbsSwapMove_CBMC{ 0.0 };
  double probabilityGibbsSwapMove_CFCMC{ 0.0 };
  double probabilityGibbsSwapMove_CFCMC_CBMC{ 0.0 };

  MoveStatistics<double> statistics_VolumeMove{ .maxChange = 0.1 };
  MoveStatistics<double> statistics_GibbsVolumeMove{ .maxChange = 0.1 };
  MoveStatistics<double3> statistics_GibbsSwapMove_CBMC{};
  MoveStatistics<double3> statistics_GibbsSwapMove_CFCMC{};
  MoveStatistics<double3> statistics_GibbsSwapMove_CFCMC_CBMC{};

  void clear();
  void optimizeAcceptance();

  std::chrono::duration<double> cpuTime_VolumeMove{ 0.0 };
  std::chrono::duration<double> cpuTime_GibbsVolumeMove{ 0.0 };
  std::chrono::duration<double> cpuTime_GibbsSwapMove_CBMC{ 0.0 };
  std::chrono::duration<double> cpuTime_GibbsSwapLambdaMove_CFCMC{ 0.0 };
  std::chrono::duration<double> cpuTime_GibbsSwapLambdaMove_CFCMC_CBMC{ 0.0 };

  std::chrono::duration<double> cpuTime_VolumeMove_NonEwald{ 0.0 };
  std::chrono::duration<double> cpuTime_GibbsVolumeMove_NonEwald{ 0.0 };
  std::chrono::duration<double> cpuTime_GibbsSwapMove_CBMC_NonEwald{ 0.0 };
  std::chrono::duration<double> cpuTime_GibbsSwapLambdaMove_CFCMC_NonEwald{ 0.0 };
  std::chrono::duration<double> cpuTime_GibbsSwapLambdaMove_CFCMC_CBMC_NonEwald{ 0.0 };

  std::chrono::duration<double> cpuTime_VolumeMove_Ewald{ 0.0 };
  std::chrono::duration<double> cpuTime_GibbsVolumeMove_Ewald{ 0.0 };
  std::chrono::duration<double> cpuTime_GibbsSwapMove_CBMC_Ewald{ 0.0 };
  std::chrono::duration<double> cpuTime_GibbsSwapLambdaMove_CFCMC_Ewald{ 0.0 };
  std::chrono::duration<double> cpuTime_GibbsSwapLambdaMove_CFCMC_CBMC_Ewald{ 0.0 };

  void clearTimingStatistics();
  const std::string writeMCMoveCPUTimeStatistics() const;
};