module mc_moves_probabilities_particles;

import <string>;
import <sstream>;

import double3;

import print;

import move_statistics;

void MCMoveProbabilitiesParticles::clearMoveStatistics()
{
  statistics_TranslationMove.clear();
  statistics_RandomTranslationMove.clear();
  statistics_RotationMove.clear();
  statistics_RandomRotationMove.clear();
  statistics_ReinsertionMove_CBMC.clear();
  statistics_IdentityChangeMove_CBMC.clear();
  statistics_SwapInsertionMove_CBMC.clear();
  statistics_SwapDeletionMove_CBMC.clear();
  statistics_SwapMove_CFCMC.clear();
  statistics_SwapMove_CFCMC_CBMC.clear();
  statistics_WidomMove_CBMC.clear();
  statistics_WidomMove_CFCMC.clear();
  statistics_WidomMove_CFCMC_CBMC.clear();
}

void MCMoveProbabilitiesParticles::optimizeMCMoves()
{
  statistics_TranslationMove.optimizeAcceptance();
  statistics_RotationMove.optimizeAcceptance();
}

void MCMoveProbabilitiesParticles::normalizeMoveProbabilties()
{
  double totalProbability =
    probabilityTranslationMove +
    probabilityRandomTranslationMove +
    probabilityRotationMove +
    probabilityRandomRotationMove +
    probabilityVolumeMove +
    probabilityReinsertionMove_CBMC +
    probabilityIdentityChangeMove_CBMC +
    probabilitySwapMove_CBMC +
    probabilitySwapMove_CFCMC +
    probabilitySwapMove_CFCMC_CBMC +
    probabilityGibbsVolumeMove +
    probabilityGibbsSwapMove_CBMC +
    probabilityGibbsSwapMove_CFCMC +
    probabilityGibbsSwapMove_CFCMC_CBMC +
    probabilityWidomMove +
    probabilityWidomMove_CFCMC +
    probabilityWidomMove_CFCMC_CBMC;


  if (totalProbability > 1e-5)
  {
    probabilityTranslationMove /= totalProbability;
    probabilityRandomTranslationMove /= totalProbability;
    probabilityRotationMove /= totalProbability;
    probabilityRandomRotationMove /= totalProbability;
    probabilityVolumeMove /= totalProbability;
    probabilityReinsertionMove_CBMC /= totalProbability;
    probabilityIdentityChangeMove_CBMC /= totalProbability;
    probabilitySwapMove_CBMC /= totalProbability;
    probabilitySwapMove_CFCMC /= totalProbability;
    probabilitySwapMove_CFCMC_CBMC /= totalProbability;
    probabilityGibbsVolumeMove /= totalProbability;
    probabilityGibbsSwapMove_CBMC /= totalProbability;
    probabilityGibbsSwapMove_CFCMC /= totalProbability;
    probabilityGibbsSwapMove_CFCMC_CBMC /= totalProbability;
    probabilityWidomMove /= totalProbability;
    probabilityWidomMove_CFCMC /= totalProbability;
    probabilityWidomMove_CFCMC_CBMC /= totalProbability;
  }

  accumulatedProbabilityTranslationMove = probabilityTranslationMove;
  accumulatedProbabilityRandomTranslationMove = probabilityRandomTranslationMove;
  accumulatedProbabilityRotationMove = probabilityRotationMove;
  accumulatedProbabilityRandomRotationMove = probabilityRandomRotationMove;
  accumulatedProbabilityVolumeMove = probabilityVolumeMove;
  accumulatedProbabilityReinsertionMove_CBMC = probabilityReinsertionMove_CBMC;
  accumulatedProbabilityIdentityChangeMove_CBMC = probabilityIdentityChangeMove_CBMC;
  accumulatedProbabilitySwapMove_CBMC = probabilitySwapMove_CBMC;
  accumulatedProbabilitySwapMove_CFCMC = probabilitySwapMove_CFCMC;
  accumulatedProbabilitySwapMove_CFCMC_CBMC = probabilitySwapMove_CFCMC_CBMC;
  accumulatedProbabilityGibbsVolumeMove = probabilityGibbsVolumeMove;
  accumulatedProbabilityGibbsSwapMove_CBMC = probabilityGibbsSwapMove_CBMC;
  accumulatedProbabilityGibbsSwapMove_CFCMC = probabilityGibbsSwapMove_CFCMC;
  accumulatedProbabilityGibbsSwapMove_CFCMC_CBMC = probabilityGibbsSwapMove_CFCMC_CBMC;
  accumulatedProbabilityWidomMove = probabilityWidomMove;
  accumulatedProbabilityWidomMove_CFCMC = probabilityWidomMove_CFCMC;
  accumulatedProbabilityWidomMove_CFCMC_CBMC = probabilityWidomMove_CFCMC_CBMC;

  accumulatedProbabilityRandomTranslationMove += accumulatedProbabilityTranslationMove;
  accumulatedProbabilityRotationMove += accumulatedProbabilityRandomTranslationMove;
  accumulatedProbabilityRandomRotationMove += accumulatedProbabilityRotationMove;
  accumulatedProbabilityVolumeMove += accumulatedProbabilityRandomRotationMove;
  accumulatedProbabilityReinsertionMove_CBMC += accumulatedProbabilityVolumeMove;
  accumulatedProbabilityIdentityChangeMove_CBMC += accumulatedProbabilityReinsertionMove_CBMC;
  accumulatedProbabilitySwapMove_CBMC += accumulatedProbabilityIdentityChangeMove_CBMC;
  accumulatedProbabilitySwapMove_CFCMC += accumulatedProbabilitySwapMove_CBMC;
  accumulatedProbabilitySwapMove_CFCMC_CBMC += accumulatedProbabilitySwapMove_CFCMC;
  accumulatedProbabilityGibbsVolumeMove += accumulatedProbabilitySwapMove_CFCMC_CBMC;
  accumulatedProbabilityGibbsSwapMove_CBMC += accumulatedProbabilityGibbsVolumeMove;
  accumulatedProbabilityGibbsSwapMove_CFCMC += accumulatedProbabilityGibbsSwapMove_CBMC;
  accumulatedProbabilityGibbsSwapMove_CFCMC_CBMC += accumulatedProbabilityGibbsSwapMove_CFCMC;
  accumulatedProbabilityWidomMove += accumulatedProbabilityGibbsSwapMove_CFCMC_CBMC;
  accumulatedProbabilityWidomMove_CFCMC += accumulatedProbabilityWidomMove;
  accumulatedProbabilityWidomMove_CFCMC_CBMC += accumulatedProbabilityWidomMove_CFCMC;

}

std::string formatStatistics(const std::string name, const MoveStatistics<double>& move)
{
  std::ostringstream stream;
  std::print(stream, "    {} total:        {:10}\n", name, move.counts);
  std::print(stream, "    {} constructed:  {:10}\n", name, move.constructed);
  std::print(stream, "    {} accepted:     {:10}\n", name, move.accepted);
  std::print(stream, "    {} fraction:     {:10f}\n", name, move.accepted / std::max(1.0, double(move.counts)));
  std::print(stream, "    {} max-change:   {:10f}\n\n", name, move.maxChange);
  return stream.str();
}

std::string formatStatistics(const std::string name, const MoveStatistics<double3>& move)
{
  std::ostringstream stream;
  std::print(stream, "    {} total:        {:10} {:10} {:10}\n", name, move.counts.x, move.counts.y, move.counts.z);
  std::print(stream, "    {} constructed:  {:10} {:10} {:10}\n", name, move.constructed.x, move.constructed.y, move.constructed.z);
  std::print(stream, "    {} accepted:     {:10} {:10} {:10}\n", name, move.accepted.x, move.accepted.y, move.accepted.z);
  std::print(stream, "    {} fraction:     {:10f} {:10f} {:10f}\n", name, move.accepted.x / std::max(1.0, double(move.counts.x)),
    move.accepted.y / std::max(1.0, double(move.counts.y)), move.accepted.z / std::max(1.0, double(move.counts.z)));
  std::print(stream, "    {} max-change:   {:10f} {:10f} {:10f}\n\n", name, move.maxChange.x, move.maxChange.y, move.maxChange.z);
  return stream.str();
}

const std::string MCMoveProbabilitiesParticles::writeMCMoveStatistics() const
{
  std::ostringstream stream;
  if (probabilityTranslationMove > 0.0) std::print(stream, formatStatistics("Translation", statistics_TranslationMove));
  if (probabilityRandomTranslationMove > 0.0) std::print(stream, formatStatistics("Random translation", statistics_RandomTranslationMove));
  if (probabilityRotationMove > 0.0) std::print(stream, formatStatistics("Rotation", statistics_RotationMove));
  if (probabilityRandomRotationMove > 0.0) std::print(stream, formatStatistics("Random rotation", statistics_RandomRotationMove));
  if (probabilityReinsertionMove_CBMC > 0.0) std::print(stream, formatStatistics("Reinsertion(CBMC)", statistics_ReinsertionMove_CBMC));
  if (probabilityIdentityChangeMove_CBMC > 0.0) std::print(stream, formatStatistics("Identity Swap (CBMC)", statistics_IdentityChangeMove_CBMC));
  if (probabilitySwapMove_CBMC > 0.0) std::print(stream, formatStatistics("Swap Insertion (CBMC)", statistics_SwapInsertionMove_CBMC));
  if (probabilitySwapMove_CBMC > 0.0) std::print(stream, formatStatistics("Swap Deletion (CBMC)", statistics_SwapDeletionMove_CBMC));
  if (probabilitySwapMove_CFCMC > 0.0) std::print(stream, formatStatistics("Swap (CFCMC)", statistics_SwapMove_CFCMC));
  if (probabilitySwapMove_CFCMC_CBMC > 0.0) std::print(stream, formatStatistics("Swap (CB/CFCMC)", statistics_SwapMove_CFCMC_CBMC));
  if (probabilityWidomMove > 0.0) std::print(stream, formatStatistics("Widom (CBMC)", statistics_WidomMove_CBMC));
  if (probabilityWidomMove_CFCMC > 0.0) std::print(stream, formatStatistics("Widom (CFCMC)", statistics_WidomMove_CFCMC));
  if (probabilityWidomMove_CFCMC_CBMC > 0.0) std::print(stream, formatStatistics("Widom (CB/CFCMC)", statistics_WidomMove_CFCMC_CBMC));

  return stream.str();
}


void MCMoveProbabilitiesParticles::clearTimingStatistics()
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

const std::string MCMoveProbabilitiesParticles::writeMCMoveCPUTimeStatistics() const
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