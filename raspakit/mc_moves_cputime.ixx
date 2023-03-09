export module mc_moves_cputime;

import <chrono>;
import <string>;

export struct MC_Move_CPUtimings
{
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