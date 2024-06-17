module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <complex>
#include <exception>
#include <format>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <print>
#include <source_location>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#endif

module mc_moves_probabilities_particles;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <iostream>;
import <sstream>;
import <fstream>;
import <format>;
import <exception>;
import <source_location>;
import <complex>;
import <vector>;
import <array>;
import <map>;
import <utility>;
import <algorithm>;
import <functional>;
import <print>;
#endif

import archive;
import double3;
import stringutils;

void MCMoveProbabilitiesParticles::normalizeMoveProbabilties()
{
  double totalProbability = probabilityTranslationMove + probabilityRandomTranslationMove + probabilityRotationMove +
                            probabilityRandomRotationMove + probabilityVolumeMove + probabilityReinsertionMove_CBMC +
                            probabilityIdentityChangeMove_CBMC + probabilitySwapMove + probabilitySwapMove_CBMC +
                            probabilitySwapMove_CFCMC + probabilitySwapMove_CFCMC_CBMC + probabilityGibbsVolumeMove +
                            probabilityGibbsSwapMove_CBMC + probabilityGibbsSwapMove_CFCMC + probabilityWidomMove +
                            probabilityWidomMove_CFCMC + probabilityWidomMove_CFCMC_CBMC +
                            probabilityParallelTemperingSwap;

  if (totalProbability > 1e-5)
  {
    probabilityTranslationMove /= totalProbability;
    probabilityRandomTranslationMove /= totalProbability;
    probabilityRotationMove /= totalProbability;
    probabilityRandomRotationMove /= totalProbability;
    probabilityVolumeMove /= totalProbability;
    probabilityReinsertionMove_CBMC /= totalProbability;
    probabilityIdentityChangeMove_CBMC /= totalProbability;
    probabilitySwapMove /= totalProbability;
    probabilitySwapMove_CBMC /= totalProbability;
    probabilitySwapMove_CFCMC /= totalProbability;
    probabilitySwapMove_CFCMC_CBMC /= totalProbability;
    probabilityGibbsVolumeMove /= totalProbability;
    probabilityGibbsSwapMove_CBMC /= totalProbability;
    probabilityGibbsSwapMove_CFCMC /= totalProbability;
    probabilityWidomMove /= totalProbability;
    probabilityWidomMove_CFCMC /= totalProbability;
    probabilityWidomMove_CFCMC_CBMC /= totalProbability;
    probabilityParallelTemperingSwap /= totalProbability;
  }

  accumulatedProbabilityTranslationMove = probabilityTranslationMove;
  accumulatedProbabilityRandomTranslationMove = probabilityRandomTranslationMove;
  accumulatedProbabilityRotationMove = probabilityRotationMove;
  accumulatedProbabilityRandomRotationMove = probabilityRandomRotationMove;
  accumulatedProbabilityVolumeMove = probabilityVolumeMove;
  accumulatedProbabilityReinsertionMove_CBMC = probabilityReinsertionMove_CBMC;
  accumulatedProbabilityIdentityChangeMove_CBMC = probabilityIdentityChangeMove_CBMC;
  accumulatedProbabilitySwapMove = probabilitySwapMove;
  accumulatedProbabilitySwapMove_CBMC = probabilitySwapMove_CBMC;
  accumulatedProbabilitySwapMove_CFCMC = probabilitySwapMove_CFCMC;
  accumulatedProbabilitySwapMove_CFCMC_CBMC = probabilitySwapMove_CFCMC_CBMC;
  accumulatedProbabilityGibbsVolumeMove = probabilityGibbsVolumeMove;
  accumulatedProbabilityGibbsSwapMove_CBMC = probabilityGibbsSwapMove_CBMC;
  accumulatedProbabilityGibbsSwapMove_CFCMC = probabilityGibbsSwapMove_CFCMC;
  accumulatedProbabilityWidomMove = probabilityWidomMove;
  accumulatedProbabilityWidomMove_CFCMC = probabilityWidomMove_CFCMC;
  accumulatedProbabilityWidomMove_CFCMC_CBMC = probabilityWidomMove_CFCMC_CBMC;
  accumulatedProbabilityParallelTemperingSwap = probabilityParallelTemperingSwap;

  accumulatedProbabilityRandomTranslationMove += accumulatedProbabilityTranslationMove;
  accumulatedProbabilityRotationMove += accumulatedProbabilityRandomTranslationMove;
  accumulatedProbabilityRandomRotationMove += accumulatedProbabilityRotationMove;
  accumulatedProbabilityVolumeMove += accumulatedProbabilityRandomRotationMove;
  accumulatedProbabilityReinsertionMove_CBMC += accumulatedProbabilityVolumeMove;
  accumulatedProbabilityIdentityChangeMove_CBMC += accumulatedProbabilityReinsertionMove_CBMC;
  accumulatedProbabilitySwapMove += accumulatedProbabilityIdentityChangeMove_CBMC;
  accumulatedProbabilitySwapMove_CBMC += accumulatedProbabilitySwapMove;
  accumulatedProbabilitySwapMove_CFCMC += accumulatedProbabilitySwapMove_CBMC;
  accumulatedProbabilitySwapMove_CFCMC_CBMC += accumulatedProbabilitySwapMove_CFCMC;
  accumulatedProbabilityGibbsVolumeMove += accumulatedProbabilitySwapMove_CFCMC_CBMC;
  accumulatedProbabilityGibbsSwapMove_CBMC += accumulatedProbabilityGibbsVolumeMove;
  accumulatedProbabilityGibbsSwapMove_CFCMC += accumulatedProbabilityGibbsSwapMove_CBMC;
  accumulatedProbabilityWidomMove += accumulatedProbabilityGibbsSwapMove_CFCMC;
  accumulatedProbabilityWidomMove_CFCMC += accumulatedProbabilityWidomMove;
  accumulatedProbabilityWidomMove_CFCMC_CBMC += accumulatedProbabilityWidomMove_CFCMC;
  accumulatedProbabilityParallelTemperingSwap += accumulatedProbabilityWidomMove_CFCMC_CBMC;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveProbabilitiesParticles &p)
{
  archive << p.versionNumber;

  archive << p.probabilityTranslationMove;
  archive << p.probabilityRandomTranslationMove;
  archive << p.probabilityRotationMove;
  archive << p.probabilityRandomRotationMove;
  archive << p.probabilityVolumeMove;
  archive << p.probabilityReinsertionMove_CBMC;
  archive << p.probabilityIdentityChangeMove_CBMC;
  archive << p.probabilitySwapMove;
  archive << p.probabilitySwapMove_CBMC;
  archive << p.probabilitySwapMove_CFCMC;
  archive << p.probabilitySwapMove_CFCMC_CBMC;
  archive << p.probabilityGibbsVolumeMove;
  archive << p.probabilityGibbsSwapMove_CBMC;
  archive << p.probabilityGibbsSwapMove_CFCMC;
  archive << p.probabilityGibbsSwapMove_CFCMC_CBMC;
  archive << p.probabilityWidomMove;
  archive << p.probabilityWidomMove_CFCMC;
  archive << p.probabilityWidomMove_CFCMC_CBMC;
  archive << p.probabilityParallelTemperingSwap;

  archive << p.accumulatedProbabilityTranslationMove;
  archive << p.accumulatedProbabilityRandomTranslationMove;
  archive << p.accumulatedProbabilityRotationMove;
  archive << p.accumulatedProbabilityRandomRotationMove;
  archive << p.accumulatedProbabilityVolumeMove;
  archive << p.accumulatedProbabilityReinsertionMove_CBMC;
  archive << p.accumulatedProbabilityIdentityChangeMove_CBMC;
  archive << p.accumulatedProbabilitySwapMove;
  archive << p.accumulatedProbabilitySwapMove_CBMC;
  archive << p.accumulatedProbabilitySwapMove_CFCMC;
  archive << p.accumulatedProbabilitySwapMove_CFCMC_CBMC;
  archive << p.accumulatedProbabilityGibbsVolumeMove;
  archive << p.accumulatedProbabilityGibbsSwapMove_CBMC;
  archive << p.accumulatedProbabilityGibbsSwapMove_CFCMC;
  archive << p.accumulatedProbabilityGibbsSwapMove_CFCMC_CBMC;
  archive << p.accumulatedProbabilityWidomMove;
  archive << p.accumulatedProbabilityWidomMove_CFCMC;
  archive << p.accumulatedProbabilityWidomMove_CFCMC_CBMC;
  archive << p.accumulatedProbabilityParallelTemperingSwap;

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveProbabilitiesParticles &p)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > p.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(
        std::format("Invalid version reading 'MCMoveProbabilitiesParticles' "
                    "at line {} in file {}\n",
                    location.line(), location.file_name()));
  }

  archive >> p.probabilityTranslationMove;
  archive >> p.probabilityRandomTranslationMove;
  archive >> p.probabilityRotationMove;
  archive >> p.probabilityRandomRotationMove;
  archive >> p.probabilityVolumeMove;
  archive >> p.probabilityReinsertionMove_CBMC;
  archive >> p.probabilityIdentityChangeMove_CBMC;
  archive >> p.probabilitySwapMove;
  archive >> p.probabilitySwapMove_CBMC;
  archive >> p.probabilitySwapMove_CFCMC;
  archive >> p.probabilitySwapMove_CFCMC_CBMC;
  archive >> p.probabilityGibbsVolumeMove;
  archive >> p.probabilityGibbsSwapMove_CBMC;
  archive >> p.probabilityGibbsSwapMove_CFCMC;
  archive >> p.probabilityGibbsSwapMove_CFCMC_CBMC;
  archive >> p.probabilityWidomMove;
  archive >> p.probabilityWidomMove_CFCMC;
  archive >> p.probabilityWidomMove_CFCMC_CBMC;
  archive >> p.probabilityParallelTemperingSwap;

  archive >> p.accumulatedProbabilityTranslationMove;
  archive >> p.accumulatedProbabilityRandomTranslationMove;
  archive >> p.accumulatedProbabilityRotationMove;
  archive >> p.accumulatedProbabilityRandomRotationMove;
  archive >> p.accumulatedProbabilityVolumeMove;
  archive >> p.accumulatedProbabilityReinsertionMove_CBMC;
  archive >> p.accumulatedProbabilityIdentityChangeMove_CBMC;
  archive >> p.accumulatedProbabilitySwapMove;
  archive >> p.accumulatedProbabilitySwapMove_CBMC;
  archive >> p.accumulatedProbabilitySwapMove_CFCMC;
  archive >> p.accumulatedProbabilitySwapMove_CFCMC_CBMC;
  archive >> p.accumulatedProbabilityGibbsVolumeMove;
  archive >> p.accumulatedProbabilityGibbsSwapMove_CBMC;
  archive >> p.accumulatedProbabilityGibbsSwapMove_CFCMC;
  archive >> p.accumulatedProbabilityGibbsSwapMove_CFCMC_CBMC;
  archive >> p.accumulatedProbabilityWidomMove;
  archive >> p.accumulatedProbabilityWidomMove_CFCMC;
  archive >> p.accumulatedProbabilityWidomMove_CFCMC_CBMC;
  archive >> p.accumulatedProbabilityParallelTemperingSwap;

  return archive;
}
