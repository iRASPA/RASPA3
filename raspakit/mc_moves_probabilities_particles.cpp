module mc_moves_probabilities_particles;

import <string>;
import <sstream>;
import <print>;

import double3;
import move_statistics;
import stringutils;

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

  statistics_GibbsSwapMove_CBMC.clear();
  statistics_GibbsSwapMove_CFCMC.clear();
}
void MCMoveProbabilitiesParticles::optimizeMCMoves()
{
  statistics_TranslationMove.optimizeAcceptance(0.01, 1.5);
  statistics_RotationMove.optimizeAcceptance(0.01, 1.5);

  statistics_SwapMove_CFCMC_CBMC.optimizeAcceptance(0.0, 1.0);
  statistics_GibbsSwapMove_CFCMC.optimizeAcceptance(0.0, 1.0);
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
  accumulatedProbabilityWidomMove += accumulatedProbabilityGibbsSwapMove_CFCMC;
  accumulatedProbabilityWidomMove_CFCMC += accumulatedProbabilityWidomMove;
  accumulatedProbabilityWidomMove_CFCMC_CBMC += accumulatedProbabilityWidomMove_CFCMC;

}

std::string formatStatistics(const std::string name, const MoveStatistics<double>& move)
{
  std::ostringstream stream;
  std::print(stream, "    {} total:        {:10}\n", name, move.totalCounts);
  std::print(stream, "    {} constructed:  {:10}\n", name, move.totalConstructed);
  std::print(stream, "    {} accepted:     {:10}\n", name, move.totalAccepted);
  std::print(stream, "    {} fraction:     {:10f}\n", name, move.totalAccepted / std::max(1.0, double(move.totalCounts)));
  std::print(stream, "    {} max-change:   {:10f}\n\n", name, move.maxChange);
  return stream.str();
}

std::string formatStatistics(const std::string name, const MoveStatistics<double3>& move)
{
  std::ostringstream stream;
  std::print(stream, "    {} total:        {:10} {:10} {:10}\n", name, move.totalCounts.x, move.totalCounts.y, move.totalCounts.z);
  std::print(stream, "    {} constructed:  {:10} {:10} {:10}\n", name, move.totalConstructed.x, move.totalConstructed.y, move.totalConstructed.z);
  std::print(stream, "    {} accepted:     {:10} {:10} {:10}\n", name, move.totalAccepted.x, move.totalAccepted.y, move.totalAccepted.z);
  std::print(stream, "    {} fraction:     {:10f} {:10f} {:10f}\n", name, move.totalAccepted.x / std::max(1.0, double(move.totalCounts.x)),
    move.totalAccepted.y / std::max(1.0, double(move.totalCounts.y)), move.totalAccepted.z / std::max(1.0, double(move.totalCounts.z)));
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

  if (probabilityGibbsSwapMove_CBMC > 0.0) std::print(stream, formatStatistics("Gibbs Swap (CBMC)", statistics_GibbsSwapMove_CBMC));
  if (probabilityGibbsSwapMove_CFCMC > 0.0) std::print(stream, formatStatistics("Gibbs Swap (CFCMC)", statistics_GibbsSwapMove_CFCMC));

  return stream.str();
}
