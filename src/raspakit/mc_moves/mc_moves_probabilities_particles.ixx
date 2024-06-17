module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <fstream>
#include <string>
#endif

export module mc_moves_probabilities_particles;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <chrono>;
import <fstream>;
#endif

import double3;
import archive;

export struct MCMoveProbabilitiesParticles
{
  uint64_t versionNumber{1};

  bool operator==(MCMoveProbabilitiesParticles const &) const = default;

  double probabilityTranslationMove{0.0};
  double probabilityRandomTranslationMove{0.0};
  double probabilityRotationMove{0.0};
  double probabilityRandomRotationMove{0.0};
  double probabilityVolumeMove{0.0};
  double probabilityReinsertionMove_CBMC{0.0};
  double probabilityIdentityChangeMove_CBMC{0.0};
  double probabilitySwapMove{0.0};
  double probabilitySwapMove_CBMC{0.0};
  double probabilitySwapMove_CFCMC{0.0};
  double probabilitySwapMove_CFCMC_CBMC{0.0};
  double probabilityGibbsVolumeMove{0.0};
  double probabilityGibbsSwapMove_CBMC{0.0};
  double probabilityGibbsSwapMove_CFCMC{0.0};
  double probabilityGibbsSwapMove_CFCMC_CBMC{0.0};
  double probabilityWidomMove{0.0};
  double probabilityWidomMove_CFCMC{0.0};
  double probabilityWidomMove_CFCMC_CBMC{0.0};
  double probabilityParallelTemperingSwap{0.0};

  double accumulatedProbabilityTranslationMove{0.0};
  double accumulatedProbabilityRandomTranslationMove{0.0};
  double accumulatedProbabilityRotationMove{0.0};
  double accumulatedProbabilityRandomRotationMove{0.0};
  double accumulatedProbabilityVolumeMove{0.0};
  double accumulatedProbabilityReinsertionMove_CBMC{0.0};
  double accumulatedProbabilityIdentityChangeMove_CBMC{0.0};
  double accumulatedProbabilitySwapMove{0.0};
  double accumulatedProbabilitySwapMove_CBMC{0.0};
  double accumulatedProbabilitySwapMove_CFCMC{0.0};
  double accumulatedProbabilitySwapMove_CFCMC_CBMC{0.0};
  double accumulatedProbabilityGibbsVolumeMove{0.0};
  double accumulatedProbabilityGibbsSwapMove_CBMC{0.0};
  double accumulatedProbabilityGibbsSwapMove_CFCMC{0.0};
  double accumulatedProbabilityGibbsSwapMove_CFCMC_CBMC{0.0};
  double accumulatedProbabilityWidomMove{0.0};
  double accumulatedProbabilityWidomMove_CFCMC{0.0};
  double accumulatedProbabilityWidomMove_CFCMC_CBMC{0.0};
  double accumulatedProbabilityParallelTemperingSwap{0.0};

  void normalizeMoveProbabilties();

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveProbabilitiesParticles &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveProbabilitiesParticles &p);
};
