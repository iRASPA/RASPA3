module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <fstream>
#include <map>
#include <random>
#include <source_location>
#include <sstream>
#include <vector>
#endif

module mc_moves_probabilities;

#ifndef USE_LEGACY_HEADERS
import <random>;
import <map>;
import <vector>;
import <fstream>;
import <sstream>;
import <source_location>;
#endif

import randomnumbers;
import mc_moves_move_types;

MCMoveProbabilities::MCMoveProbabilities(double translationProbability, double randomTranslationProbability,
                                         double rotationProbability, double randomRotationProbability,
                                         double volumeChangeProbability, double reinsertionCBMCProbability,
                                         double identityChangeCBMCProbability, double swapProbability,
                                         double swapCBMCProbability, double swapCFCMCProbability,
                                         double swapCBCFCMCProbability, double gibbsVolumeChangeProbability,
                                         double gibbsSwapCBMCProbability, double gibbsSwapCFCMCProbability,
                                         double widomProbability, double widomCFCMCProbability,
                                         double widomCBCFCMCProbability, double parallelTemperingProbability,
                                         double hybridMCProbability)
{
  probabilities[MoveTypes::Translation] = translationProbability;
  probabilities[MoveTypes::RandomTranslation] = randomTranslationProbability;
  probabilities[MoveTypes::Rotation] = rotationProbability;
  probabilities[MoveTypes::RandomRotation] = randomRotationProbability;
  probabilities[MoveTypes::VolumeChange] = volumeChangeProbability;
  probabilities[MoveTypes::ReinsertionCBMC] = reinsertionCBMCProbability;
  probabilities[MoveTypes::IdentityChangeCBMC] = identityChangeCBMCProbability;
  probabilities[MoveTypes::Swap] = swapProbability;
  probabilities[MoveTypes::SwapCBMC] = swapCBMCProbability;
  probabilities[MoveTypes::SwapCFCMC] = swapCFCMCProbability;
  probabilities[MoveTypes::SwapCBCFCMC] = swapCBCFCMCProbability;
  probabilities[MoveTypes::GibbsVolume] = gibbsVolumeChangeProbability;
  probabilities[MoveTypes::GibbsSwapCBMC] = gibbsSwapCBMCProbability;
  probabilities[MoveTypes::GibbsSwapCFCMC] = gibbsSwapCFCMCProbability;
  probabilities[MoveTypes::Widom] = widomProbability;
  probabilities[MoveTypes::WidomCFCMC] = widomCFCMCProbability;
  probabilities[MoveTypes::WidomCBCFCMC] = widomCBCFCMCProbability;
  probabilities[MoveTypes::ParallelTempering] = parallelTemperingProbability;
  probabilities[MoveTypes::HybridMC] = hybridMCProbability;
}

const std::map<MoveTypes, double> MCMoveProbabilities::normalizedMap() const
{
  double totalProbability = 0.0;
  std::map<MoveTypes, double> normalizedMap(probabilities);
  for (auto &[moveType, probability] : normalizedMap)
  {
    totalProbability += probability;
  }
  for (auto &[moveType, probability] : normalizedMap)
  {
    probability /= totalProbability;
  }
  return normalizedMap;
}

void MCMoveProbabilities::removeRedundantMoves()
{
  if (probabilities[MoveTypes::WidomCFCMC] > 0.0 && probabilities[MoveTypes::SwapCFCMC] > 0.0)
  {
    setProbability((MoveTypes::WidomCFCMC), 0.0);
  }
  if (probabilities[MoveTypes::WidomCBCFCMC] > 0.0 && probabilities[MoveTypes::SwapCBCFCMC] > 0.0)
  {
    setProbability((MoveTypes::WidomCBCFCMC), 0.0);
  }
}

void MCMoveProbabilities::join(const MCMoveProbabilities &other)
{
  for (size_t i = 0; i < static_cast<size_t>(MoveTypes::Count); i++)
  {
    if (probabilities[static_cast<MoveTypes>(i)] > 0.0 && other.getProbability(static_cast<MoveTypes>(i)) > 0.0)
    {
      throw std::runtime_error("Adding MCMoveProbabilities that define the same value is not permitted.");
    }
    else if (other.getProbability(static_cast<MoveTypes>(i)) > 0.0)
    {
      probabilities[static_cast<MoveTypes>(i)] = other.getProbability(static_cast<MoveTypes>(i));
    }
  }
}

MoveTypes MCMoveProbabilities::sample(RandomNumber &random)
{
  std::vector<double> vectorProbabilities(static_cast<size_t>(MoveTypes::Count));
  for (size_t i = 0; i < static_cast<size_t>(MoveTypes::Count); i++)
  {
    vectorProbabilities[i] = probabilities[static_cast<MoveTypes>(i)];
  }
  return static_cast<MoveTypes>(random.categoricalDistribution(vectorProbabilities));
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveProbabilities &p)
{
  archive << p.versionNumber;
  archive << p.probabilities;
  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveProbabilities &p)
{
  uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > p.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'MCMoveProbabilitiesSystem' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }
  archive >> p.probabilities;

  return archive;
}
