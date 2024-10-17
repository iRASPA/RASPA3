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

  MCMoveProbabilitiesParticles(double translationProbability = 0.0, double randomTranslationProbability = 0.0,
                               double rotationProbability = 0.0, double randomRotationProbability = 0.0,
                               double volumeChangeProbability = 0.0, double reinsertionCBMCProbability = 0.0,
                               double identityChangeCBMCProbability = 0.0, double swapProbability = 0.0,
                               double swapCBMCProbability = 0.0, double swapCFCMCProbability = 0.0,
                               double swapCBCFCMCProbability = 0.0, double gibbsVolumeChangeProbability = 0.0,
                               double gibbsSwapCBMCProbability = 0.0, double gibbsSwapCFCMCProbability = 0.0,
                               double gibbsSwapCBCFCMCProbability = 0.0, double widomProbability = 0.0,
                               double widomCFCMCProbability = 0.0, double widomCBCFCMCProbability = 0.0,
                               double parallelTemperingProbability = 0.0, double hybridMCProbability = 0.0);

  double translationProbability;
  double randomTranslationProbability;
  double rotationProbability;
  double randomRotationProbability;
  double volumeChangeProbability;
  double reinsertionCBMCProbability;
  double identityChangeCBMCProbability;
  double swapProbability;
  double swapCBMCProbability;
  double swapCFCMCProbability;
  double swapCBCFCMCProbability;
  double gibbsVolumeChangeProbability;
  double gibbsSwapCBMCProbability;
  double gibbsSwapCFCMCProbability;
  double gibbsSwapCBCFCMCProbability;
  double widomProbability;
  double widomCFCMCProbability;
  double widomCBCFCMCProbability;
  double parallelTemperingProbability;
  double hybridMCProbability;

  double accumulatedTranslationProbability{0.0};
  double accumulatedRandomTranslationProbability{0.0};
  double accumulatedRotationProbability{0.0};
  double accumulatedRandomRotationProbability{0.0};
  double accumulatedVolumeChangeProbability{0.0};
  double accumulatedReinsertionCBMCProbability{0.0};
  double accumulatedIdentityChangeCBMCProbability{0.0};
  double accumulatedSwapProbability{0.0};
  double accumulatedSwapCBMCProbability{0.0};
  double accumulatedSwapCFCMCProbability{0.0};
  double accumulatedSwapCBCFCMCProbability{0.0};
  double accumulatedGibbsVolumeChangeProbability{0.0};
  double accumulatedGibbsSwapCBMCProbability{0.0};
  double accumulatedGibbsSwapCFCMCProbability{0.0};
  double accumulatedGibbsSwapCBCFCMCProbability{0.0};
  double accumulatedWidomProbability{0.0};
  double accumulatedWidomCFCMCProbability{0.0};
  double accumulatedWidomCBCFCMCProbability{0.0};
  double accumulatedParallelTemperingProbability{0.0};
  double accumulatedHybridMCProbability{0.0};

  void normalizeMoveProbabilities();

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveProbabilitiesParticles &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveProbabilitiesParticles &p);
};
