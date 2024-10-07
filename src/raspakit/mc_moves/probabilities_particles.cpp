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

MCMoveProbabilitiesParticles::MCMoveProbabilitiesParticles(
    double translationProbability, double randomTranslationProbability, double rotationProbability,
    double randomRotationProbability, double volumeChangeProbability, double reinsertionCBMCProbability,
    double identityChangeCBMCProbability, double swapProbability, double swapCBMCProbability,
    double swapCFCMCProbability, double swapCBCFCMCProbability, double gibbsVolumeChangeProbability,
    double gibbsSwapCBMCProbability, double gibbsSwapCFCMCProbability,
    double gibbsSwapCBCFCMCProbability, double widomProbability, double widomCFCMCProbability,
    double widomCBCFCMCProbability, double parallelTemperingProbability)
    : translationProbability(translationProbability),
      randomTranslationProbability(randomTranslationProbability),
      rotationProbability(rotationProbability),
      randomRotationProbability(randomRotationProbability),
      volumeChangeProbability(volumeChangeProbability),
      reinsertionCBMCProbability(reinsertionCBMCProbability),
      identityChangeCBMCProbability(identityChangeCBMCProbability),
      swapProbability(swapProbability),
      swapCBMCProbability(swapCBMCProbability),
      swapCFCMCProbability(swapCFCMCProbability),
      swapCBCFCMCProbability(swapCBCFCMCProbability),
      gibbsVolumeChangeProbability(gibbsVolumeChangeProbability),
      gibbsSwapCBMCProbability(gibbsSwapCBMCProbability),
      gibbsSwapCFCMCProbability(gibbsSwapCFCMCProbability),
      gibbsSwapCBCFCMCProbability(gibbsSwapCBCFCMCProbability),
      widomProbability(widomProbability),
      widomCFCMCProbability(widomCFCMCProbability),
      widomCBCFCMCProbability(widomCBCFCMCProbability),
      parallelTemperingProbability(parallelTemperingProbability)
{
  normalizeMoveProbabilities();
}

void MCMoveProbabilitiesParticles::normalizeMoveProbabilities()
{
  double totalProbability = translationProbability + randomTranslationProbability + rotationProbability +
                            randomRotationProbability + volumeChangeProbability + reinsertionCBMCProbability +
                            identityChangeCBMCProbability + swapProbability + swapCBMCProbability +
                            swapCFCMCProbability + swapCBCFCMCProbability + gibbsVolumeChangeProbability +
                            gibbsSwapCBMCProbability + gibbsSwapCFCMCProbability + widomProbability +
                            widomCFCMCProbability + widomCBCFCMCProbability +
                            parallelTemperingProbability;

  if (totalProbability > 1e-5)
  {
    translationProbability /= totalProbability;
    randomTranslationProbability /= totalProbability;
    rotationProbability /= totalProbability;
    randomRotationProbability /= totalProbability;
    volumeChangeProbability /= totalProbability;
    reinsertionCBMCProbability /= totalProbability;
    identityChangeCBMCProbability /= totalProbability;
    swapProbability /= totalProbability;
    swapCBMCProbability /= totalProbability;
    swapCFCMCProbability /= totalProbability;
    swapCBCFCMCProbability /= totalProbability;
    gibbsVolumeChangeProbability /= totalProbability;
    gibbsSwapCBMCProbability /= totalProbability;
    gibbsSwapCFCMCProbability /= totalProbability;
    widomProbability /= totalProbability;
    widomCFCMCProbability /= totalProbability;
    widomCBCFCMCProbability /= totalProbability;
    parallelTemperingProbability /= totalProbability;
  }

  accumulatedTranslationProbability = translationProbability;
  accumulatedRandomTranslationProbability = randomTranslationProbability;
  accumulatedRotationProbability = rotationProbability;
  accumulatedRandomRotationProbability = randomRotationProbability;
  accumulatedVolumeChangeProbability = volumeChangeProbability;
  accumulatedReinsertionCBMCProbability = reinsertionCBMCProbability;
  accumulatedIdentityChangeCBMCProbability = identityChangeCBMCProbability;
  accumulatedSwapProbability = swapProbability;
  accumulatedSwapCBMCProbability = swapCBMCProbability;
  accumulatedSwapCFCMCProbability = swapCFCMCProbability;
  accumulatedSwapCBCFCMCProbability = swapCBCFCMCProbability;
  accumulatedGibbsVolumeChangeProbability = gibbsVolumeChangeProbability;
  accumulatedGibbsSwapCBMCProbability = gibbsSwapCBMCProbability;
  accumulatedGibbsSwapCFCMCProbability = gibbsSwapCFCMCProbability;
  accumulatedWidomProbability = widomProbability;
  accumulatedWidomCFCMCProbability = widomCFCMCProbability;
  accumulatedWidomCBCFCMCProbability = widomCBCFCMCProbability;
  accumulatedParallelTemperingProbability = parallelTemperingProbability;

  accumulatedRandomTranslationProbability += accumulatedTranslationProbability;
  accumulatedRotationProbability += accumulatedRandomTranslationProbability;
  accumulatedRandomRotationProbability += accumulatedRotationProbability;
  accumulatedVolumeChangeProbability += accumulatedRandomRotationProbability;
  accumulatedReinsertionCBMCProbability += accumulatedVolumeChangeProbability;
  accumulatedIdentityChangeCBMCProbability += accumulatedReinsertionCBMCProbability;
  accumulatedSwapProbability += accumulatedIdentityChangeCBMCProbability;
  accumulatedSwapCBMCProbability += accumulatedSwapProbability;
  accumulatedSwapCFCMCProbability += accumulatedSwapCBMCProbability;
  accumulatedSwapCBCFCMCProbability += accumulatedSwapCFCMCProbability;
  accumulatedGibbsVolumeChangeProbability += accumulatedSwapCBCFCMCProbability;
  accumulatedGibbsSwapCBMCProbability += accumulatedGibbsVolumeChangeProbability;
  accumulatedGibbsSwapCFCMCProbability += accumulatedGibbsSwapCBMCProbability;
  accumulatedWidomProbability += accumulatedGibbsSwapCFCMCProbability;
  accumulatedWidomCFCMCProbability += accumulatedWidomProbability;
  accumulatedWidomCBCFCMCProbability += accumulatedWidomCFCMCProbability;
  accumulatedParallelTemperingProbability += accumulatedWidomCBCFCMCProbability;
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveProbabilitiesParticles &p)
{
  archive << p.versionNumber;

  archive << p.translationProbability;
  archive << p.randomTranslationProbability;
  archive << p.rotationProbability;
  archive << p.randomRotationProbability;
  archive << p.volumeChangeProbability;
  archive << p.reinsertionCBMCProbability;
  archive << p.identityChangeCBMCProbability;
  archive << p.swapProbability;
  archive << p.swapCBMCProbability;
  archive << p.swapCFCMCProbability;
  archive << p.swapCBCFCMCProbability;
  archive << p.gibbsVolumeChangeProbability;
  archive << p.gibbsSwapCBMCProbability;
  archive << p.gibbsSwapCFCMCProbability;
  archive << p.gibbsSwapCBCFCMCProbability;
  archive << p.widomProbability;
  archive << p.widomCFCMCProbability;
  archive << p.widomCBCFCMCProbability;
  archive << p.parallelTemperingProbability;

  archive << p.accumulatedTranslationProbability;
  archive << p.accumulatedRandomTranslationProbability;
  archive << p.accumulatedRotationProbability;
  archive << p.accumulatedRandomRotationProbability;
  archive << p.accumulatedVolumeChangeProbability;
  archive << p.accumulatedReinsertionCBMCProbability;
  archive << p.accumulatedIdentityChangeCBMCProbability;
  archive << p.accumulatedSwapProbability;
  archive << p.accumulatedSwapCBMCProbability;
  archive << p.accumulatedSwapCFCMCProbability;
  archive << p.accumulatedSwapCBCFCMCProbability;
  archive << p.accumulatedGibbsVolumeChangeProbability;
  archive << p.accumulatedGibbsSwapCBMCProbability;
  archive << p.accumulatedGibbsSwapCFCMCProbability;
  archive << p.accumulatedGibbsSwapCBCFCMCProbability;
  archive << p.accumulatedWidomProbability;
  archive << p.accumulatedWidomCFCMCProbability;
  archive << p.accumulatedWidomCBCFCMCProbability;
  archive << p.accumulatedParallelTemperingProbability;

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

  archive >> p.translationProbability;
  archive >> p.randomTranslationProbability;
  archive >> p.rotationProbability;
  archive >> p.randomRotationProbability;
  archive >> p.volumeChangeProbability;
  archive >> p.reinsertionCBMCProbability;
  archive >> p.identityChangeCBMCProbability;
  archive >> p.swapProbability;
  archive >> p.swapCBMCProbability;
  archive >> p.swapCFCMCProbability;
  archive >> p.swapCBCFCMCProbability;
  archive >> p.gibbsVolumeChangeProbability;
  archive >> p.gibbsSwapCBMCProbability;
  archive >> p.gibbsSwapCFCMCProbability;
  archive >> p.gibbsSwapCBCFCMCProbability;
  archive >> p.widomProbability;
  archive >> p.widomCFCMCProbability;
  archive >> p.widomCBCFCMCProbability;
  archive >> p.parallelTemperingProbability;

  archive >> p.accumulatedTranslationProbability;
  archive >> p.accumulatedRandomTranslationProbability;
  archive >> p.accumulatedRotationProbability;
  archive >> p.accumulatedRandomRotationProbability;
  archive >> p.accumulatedVolumeChangeProbability;
  archive >> p.accumulatedReinsertionCBMCProbability;
  archive >> p.accumulatedIdentityChangeCBMCProbability;
  archive >> p.accumulatedSwapProbability;
  archive >> p.accumulatedSwapCBMCProbability;
  archive >> p.accumulatedSwapCFCMCProbability;
  archive >> p.accumulatedSwapCBCFCMCProbability;
  archive >> p.accumulatedGibbsVolumeChangeProbability;
  archive >> p.accumulatedGibbsSwapCBMCProbability;
  archive >> p.accumulatedGibbsSwapCFCMCProbability;
  archive >> p.accumulatedGibbsSwapCBCFCMCProbability;
  archive >> p.accumulatedWidomProbability;
  archive >> p.accumulatedWidomCFCMCProbability;
  archive >> p.accumulatedWidomCBCFCMCProbability;
  archive >> p.accumulatedParallelTemperingProbability;

  return archive;
}
