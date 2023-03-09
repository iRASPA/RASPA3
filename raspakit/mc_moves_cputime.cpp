module mc_moves_cputime;

import <string>;
import <sstream>;

import print;

void MC_Move_CPUtimings::clearTimingStatistics()
{
  cpuTime_TranslationMove = std::chrono::duration<double>(0.0);
  cpuTime_RandomTranslationMove = std::chrono::duration<double>(0.0);
  cpuTime_RotationMove = std::chrono::duration<double>(0.0);
  cpuTime_RandomRotationMove = std::chrono::duration<double>(0.0);
  cpuTime_ReinsertionMove_CBMC = std::chrono::duration<double>(0.0);
  cpuTime_IdentityChangeMove_CBMC = std::chrono::duration<double>(0.0);
  cpuTime_SwapInsertionMove_CBMC = std::chrono::duration<double>(0.0);
  cpuTime_SwapDeletionMove_CBMC = std::chrono::duration<double>(0.0);
  cpuTime_SwapMove_CFCMC = std::chrono::duration<double>(0.0);
  cpuTime_SwapMove_CFCMC_CBMC = std::chrono::duration<double>(0.0);
  cpuTime_WidomMove_CBMC = std::chrono::duration<double>(0.0);
  cpuTime_WidomMove_CFCMC = std::chrono::duration<double>(0.0);
  cpuTime_WidomMove_CFCMC_CBMC = std::chrono::duration<double>(0.0);

  cpuTime_TranslationMove_NonEwald = std::chrono::duration<double>(0.0);
  cpuTime_RandomTranslationMove_NonEwald = std::chrono::duration<double>(0.0);
  cpuTime_RotationMove_NonEwald = std::chrono::duration<double>(0.0);
  cpuTime_RandomRotationMove_NonEwald = std::chrono::duration<double>(0.0);
  cpuTime_ReinsertionGrowMove_CBMC_NonEwald = std::chrono::duration<double>(0.0);
  cpuTime_ReinsertionRetraceMove_CBMC_NonEwald = std::chrono::duration<double>(0.0);
  cpuTime_IdentityChangeMove_CBMC_NonEwald = std::chrono::duration<double>(0.0);
  cpuTime_SwapInsertionMove_CBMC_NonEwald = std::chrono::duration<double>(0.0);
  cpuTime_SwapDeletionMove_CBMC_NonEwald = std::chrono::duration<double>(0.0);
  cpuTime_SwapMove_CFCMC_NonEwald = std::chrono::duration<double>(0.0);
  cpuTime_SwapInsertionMove_CFCMC_CBMC_NonEwald = std::chrono::duration<double>(0.0);
  cpuTime_SwapInsertionGrowMove_CFCMC_CBMC_NonEwald = std::chrono::duration<double>(0.0);
  cpuTime_SwapDeletionMove_CFCMC_CBMC_NonEwald = std::chrono::duration<double>(0.0);
  cpuTime_SwapDeletionRetraceMove_CFCMC_CBMC_NonEwald = std::chrono::duration<double>(0.0);
  cpuTime_SwapLambdaMove_CFCMC_CBMC_NonEwald = std::chrono::duration<double>(0.0);
  cpuTime_WidomMove_CBMC_NonEwald = std::chrono::duration<double>(0.0);
  cpuTime_WidomMove_CFCMC_NonEwald = std::chrono::duration<double>(0.0);
  cpuTime_WidomMove_CFCMC_CBMC_NonEwald = std::chrono::duration<double>(0.0);

  cpuTime_TranslationMove_Ewald = std::chrono::duration<double>(0.0);
  cpuTime_RandomTranslationMove_Ewald = std::chrono::duration<double>(0.0);
  cpuTime_RotationMove_Ewald = std::chrono::duration<double>(0.0);
  cpuTime_RandomRotationMove_Ewald = std::chrono::duration<double>(0.0);
  cpuTime_ReinsertionMove_CBMC_Ewald = std::chrono::duration<double>(0.0);
  cpuTime_IdentityChangeMove_CBMC_Ewald = std::chrono::duration<double>(0.0);
  cpuTime_SwapInsertionMove_CBMC_Ewald = std::chrono::duration<double>(0.0);
  cpuTime_SwapDeletionMove_CBMC_Ewald = std::chrono::duration<double>(0.0);
  cpuTime_SwapMove_CFCMC_Ewald = std::chrono::duration<double>(0.0);
  cpuTime_SwapInsertionMove_CFCMC_CBMC_Ewald = std::chrono::duration<double>(0.0);
  cpuTime_SwapInsertionGrowMove_CFCMC_CBMC_Ewald = std::chrono::duration<double>(0.0);
  cpuTime_SwapDeletionMove_CFCMC_CBMC_Ewald = std::chrono::duration<double>(0.0);
  cpuTime_SwapDeletionRetraceMove_CFCMC_CBMC_Ewald = std::chrono::duration<double>(0.0);
  cpuTime_SwapLambdaMove_CFCMC_CBMC_Ewald = std::chrono::duration<double>(0.0);
  cpuTime_WidomMove_CBMC_Ewald = std::chrono::duration<double>(0.0);
  cpuTime_WidomMove_CFCMC_Ewald = std::chrono::duration<double>(0.0);
  cpuTime_WidomMove_CFCMC_CBMC_Ewald = std::chrono::duration<double>(0.0);
}

const std::string MC_Move_CPUtimings::writeMCMoveCPUTimeStatistics() const
{
  std::ostringstream stream;
  if (cpuTime_TranslationMove.count() > 0.0)
  {
    std::print(stream, "    Translation move:       {:14f} [s]\n", cpuTime_TranslationMove.count());
    std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_TranslationMove_NonEwald.count());
    std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_TranslationMove_Ewald.count());
    std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_TranslationMove.count() - cpuTime_TranslationMove_NonEwald.count() - cpuTime_TranslationMove_Ewald.count());
  }
  if (cpuTime_RandomTranslationMove.count() > 0.0)
  {
    std::print(stream, "    Random translation move:    {:14f} [s]\n", cpuTime_RandomTranslationMove.count());
    std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_RandomTranslationMove_NonEwald.count());
    std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_RandomTranslationMove_Ewald.count());
    std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_RandomTranslationMove.count() -
      cpuTime_RandomTranslationMove_NonEwald.count() - cpuTime_RandomTranslationMove_Ewald.count());
  }
  if (cpuTime_RotationMove.count() > 0.0)
  {
    std::print(stream, "    Rotation move:          {:14f} [s]\n", cpuTime_RotationMove.count());
    std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_RotationMove_NonEwald.count());
    std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_RotationMove_Ewald.count());
    std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_RotationMove.count() - cpuTime_RotationMove_NonEwald.count() - cpuTime_RotationMove_Ewald.count());
  }
  if (cpuTime_RandomRotationMove.count() > 0.0)
  {
    std::print(stream, "    Random rotation move:       {:14f} [s]\n", cpuTime_RandomRotationMove.count());
    std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_RandomRotationMove_NonEwald.count());
    std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_RandomRotationMove_Ewald.count());
    std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_RandomRotationMove.count() - cpuTime_RandomRotationMove_NonEwald.count() - cpuTime_RandomRotationMove_Ewald.count());
  }
  if (cpuTime_ReinsertionMove_CBMC.count() > 0.0)
  {
    std::print(stream, "    Reinsertion (CBMC):     {:14f} [s]\n", cpuTime_ReinsertionMove_CBMC.count());
    std::print(stream, "        Grow Non-Ewald:         {:14f} [s]\n", cpuTime_ReinsertionGrowMove_CBMC_NonEwald.count());
    std::print(stream, "        Retrace Non-Ewald:      {:14f} [s]\n", cpuTime_ReinsertionRetraceMove_CBMC_NonEwald.count());
    std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_ReinsertionMove_CBMC_Ewald.count());
    std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_ReinsertionMove_CBMC.count() - cpuTime_ReinsertionGrowMove_CBMC_NonEwald.count() -
      cpuTime_ReinsertionRetraceMove_CBMC_NonEwald.count() - cpuTime_ReinsertionMove_CBMC_Ewald.count());
  }
  if (cpuTime_IdentityChangeMove_CBMC.count() > 0.0)
  {
    std::print(stream, "    Identity-change (CBMC): {:14f} [s]\n", cpuTime_IdentityChangeMove_CBMC.count());
    std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_IdentityChangeMove_CBMC_NonEwald.count());
    std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_IdentityChangeMove_CBMC_Ewald.count());
    std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_IdentityChangeMove_CBMC.count() - cpuTime_IdentityChangeMove_CBMC_NonEwald.count() - cpuTime_IdentityChangeMove_CBMC_Ewald.count());
  }
  if (cpuTime_SwapInsertionMove_CBMC.count() > 0.0)
  {
    std::print(stream, "    Swap insert (CBMC):     {:14f} [s]\n", cpuTime_SwapInsertionMove_CBMC.count());
    std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_SwapInsertionMove_CBMC_NonEwald.count());
    std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_SwapInsertionMove_CBMC_Ewald.count());
    std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_SwapInsertionMove_CBMC.count() - cpuTime_SwapInsertionMove_CBMC_NonEwald.count() - cpuTime_SwapInsertionMove_CBMC_Ewald.count());
  }
  if (cpuTime_SwapDeletionMove_CBMC.count() > 0.0)
  {
    std::print(stream, "    Swap delete (CBMC):     {:14f} [s]\n", cpuTime_SwapDeletionMove_CBMC.count());
    std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_SwapDeletionMove_CBMC_NonEwald.count());
    std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_SwapDeletionMove_CBMC_Ewald.count());
    std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_SwapDeletionMove_CBMC.count() - cpuTime_SwapDeletionMove_CBMC_NonEwald.count() - cpuTime_SwapDeletionMove_CBMC_Ewald.count());
  }
  if (cpuTime_SwapMove_CFCMC.count() > 0.0)
  {
    std::print(stream, "    Swap (CFCMC):           {:14f} [s]\n", cpuTime_SwapMove_CFCMC.count());
    std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_SwapMove_CFCMC_NonEwald.count());
    std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_SwapMove_CFCMC_Ewald.count());
    std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_SwapMove_CFCMC.count() - cpuTime_SwapMove_CFCMC_NonEwald.count() - cpuTime_SwapMove_CFCMC_Ewald.count());
  }
  if (cpuTime_SwapMove_CFCMC_CBMC.count() > 0.0)
  {
    std::print(stream, "    Swap (CB/CFCMC):        {:14f} [s]\n", cpuTime_SwapMove_CFCMC_CBMC.count());
    std::print(stream, "        Ins. Non-Ewald:         {:14f} [s]\n", cpuTime_SwapInsertionMove_CFCMC_CBMC_NonEwald.count());
    std::print(stream, "        Ins. Ewald:             {:14f} [s]\n", cpuTime_SwapInsertionMove_CFCMC_CBMC_Ewald.count());
    std::print(stream, "        Ins. Grow Non-Ewald:    {:14f} [s]\n", cpuTime_SwapInsertionGrowMove_CFCMC_CBMC_NonEwald.count());
    std::print(stream, "        Ins. Grow Ewald:        {:14f} [s]\n", cpuTime_SwapInsertionGrowMove_CFCMC_CBMC_Ewald.count());
    std::print(stream, "        Del. Retrace Non-Ewald  {:14f} [s]\n", cpuTime_SwapDeletionRetraceMove_CFCMC_CBMC_NonEwald.count());
    std::print(stream, "        Del. Retrace Ewald:     {:14f} [s]\n", cpuTime_SwapDeletionRetraceMove_CFCMC_CBMC_Ewald.count());
    std::print(stream, "        Del. Non-Ewald:         {:14f} [s]\n", cpuTime_SwapDeletionMove_CFCMC_CBMC_NonEwald.count());
    std::print(stream, "        Del. Ewald:             {:14f} [s]\n", cpuTime_SwapDeletionMove_CFCMC_CBMC_Ewald.count());
    std::print(stream, "        Lambda Non-Ewald:       {:14f} [s]\n", cpuTime_SwapLambdaMove_CFCMC_CBMC_NonEwald.count());
    std::print(stream, "        Lambda Ewald:           {:14f} [s]\n", cpuTime_SwapLambdaMove_CFCMC_CBMC_Ewald.count());
    std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_SwapMove_CFCMC_CBMC.count() -
      cpuTime_SwapInsertionMove_CFCMC_CBMC_NonEwald.count() -
      cpuTime_SwapInsertionMove_CFCMC_CBMC_Ewald.count() -
      cpuTime_SwapInsertionGrowMove_CFCMC_CBMC_NonEwald.count() -
      cpuTime_SwapInsertionGrowMove_CFCMC_CBMC_Ewald.count() -
      cpuTime_SwapDeletionRetraceMove_CFCMC_CBMC_NonEwald.count() -
      cpuTime_SwapDeletionRetraceMove_CFCMC_CBMC_Ewald.count() -
      cpuTime_SwapDeletionMove_CFCMC_CBMC_NonEwald.count() -
      cpuTime_SwapDeletionMove_CFCMC_CBMC_Ewald.count() -
      cpuTime_SwapLambdaMove_CFCMC_CBMC_NonEwald.count() -
      cpuTime_SwapLambdaMove_CFCMC_CBMC_Ewald.count());
  }

  if (cpuTime_WidomMove_CBMC.count() > 0.0)
  {
    std::print(stream, "    Widom:                  {:14f} [s]\n", cpuTime_WidomMove_CBMC.count());
    std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_WidomMove_CBMC_NonEwald.count());
    std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_WidomMove_CBMC_Ewald.count());
    std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_WidomMove_CBMC.count() - cpuTime_WidomMove_CBMC_NonEwald.count() - cpuTime_WidomMove_CBMC_Ewald.count());
  }
  if (cpuTime_WidomMove_CFCMC.count() > 0.0)
  {
    std::print(stream, "    Widom (CFCMC):          {:14f} [s]\n", cpuTime_WidomMove_CFCMC.count());
    std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_WidomMove_CFCMC_NonEwald.count());
    std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_WidomMove_CFCMC_Ewald.count());
    std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_WidomMove_CFCMC.count() - cpuTime_WidomMove_CFCMC_NonEwald.count() - cpuTime_WidomMove_CFCMC_Ewald.count());
  }
  if (cpuTime_WidomMove_CFCMC_CBMC.count() > 0.0)
  {
    std::print(stream, "    Widom (CB/CFCMC):       {:14f} [s]\n", cpuTime_WidomMove_CFCMC_CBMC.count());
    std::print(stream, "        Non-Ewald:              {:14f} [s]\n", cpuTime_WidomMove_CFCMC_CBMC_NonEwald.count());
    std::print(stream, "        Ewald:                  {:14f} [s]\n", cpuTime_WidomMove_CFCMC_CBMC_Ewald.count());
    std::print(stream, "        Overhead:               {:14f} [s]\n", cpuTime_WidomMove_CFCMC_CBMC.count() - cpuTime_WidomMove_CFCMC_CBMC_NonEwald.count() - cpuTime_WidomMove_CFCMC_CBMC_Ewald.count());
  }
  std::print(stream, "\n");

  return stream.str();
}