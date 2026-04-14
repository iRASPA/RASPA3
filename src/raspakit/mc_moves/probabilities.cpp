module;

module mc_moves_probabilities;

import std;

import randomnumbers;
import mc_moves_move_types;

MCMoveProbabilities::MCMoveProbabilities(double translationProbability, double randomTranslationProbability,
                                         double rotationProbability, double randomRotationProbability,
                                         double volumeChangeProbability, double reinsertionCBMCProbability,
                                         double partialReinsertionCBMCProbability, double identityChangeCBMCProbability,
                                         double swapProbability, double swapCBMCProbability,
                                         double swapCFCMCProbability, double swapCBCFCMCProbability,
                                         double gibbsVolumeChangeProbability, double gibbsSwapCBMCProbability,
                                         double gibbsSwapCFCMCProbability, double widomProbability,
                                         double widomCFCMCProbability, double widomCBCFCMCProbability,
                                         double parallelTemperingProbability, double hybridMCProbability):
                                                       probabilities(std::to_underlying(Move::Types::Count))
{
  probabilities[std::to_underlying(Move::Types::Translation)] = translationProbability;
  probabilities[std::to_underlying(Move::Types::RandomTranslation)] = randomTranslationProbability;
  probabilities[std::to_underlying(Move::Types::Rotation)] = rotationProbability;
  probabilities[std::to_underlying(Move::Types::RandomRotation)] = randomRotationProbability;
  probabilities[std::to_underlying(Move::Types::VolumeChange)] = volumeChangeProbability;
  probabilities[std::to_underlying(Move::Types::ReinsertionCBMC)] = reinsertionCBMCProbability;
  probabilities[std::to_underlying(Move::Types::PartialReinsertionCBMC)] = partialReinsertionCBMCProbability;
  probabilities[std::to_underlying(Move::Types::IdentityChangeCBMC)] = identityChangeCBMCProbability;
  probabilities[std::to_underlying(Move::Types::Swap)] = swapProbability;
  probabilities[std::to_underlying(Move::Types::SwapCBMC)] = swapCBMCProbability;
  probabilities[std::to_underlying(Move::Types::SwapCFCMC)] = swapCFCMCProbability;
  probabilities[std::to_underlying(Move::Types::SwapCBCFCMC)] = swapCBCFCMCProbability;
  probabilities[std::to_underlying(Move::Types::GibbsVolume)] = gibbsVolumeChangeProbability;
  probabilities[std::to_underlying(Move::Types::GibbsSwapCBMC)] = gibbsSwapCBMCProbability;
  probabilities[std::to_underlying(Move::Types::GibbsSwapCFCMC)] = gibbsSwapCFCMCProbability;
  probabilities[std::to_underlying(Move::Types::Widom)] = widomProbability;
  probabilities[std::to_underlying(Move::Types::WidomCFCMC)] = widomCFCMCProbability;
  probabilities[std::to_underlying(Move::Types::WidomCBCFCMC)] = widomCBCFCMCProbability;
  probabilities[std::to_underlying(Move::Types::ParallelTempering)] = parallelTemperingProbability;
  probabilities[std::to_underlying(Move::Types::HybridMC)] = hybridMCProbability;
}

const std::vector<double> MCMoveProbabilities::normalizedMap() const
{
  double totalProbability = 0.0;
  std::vector<double> normalizedMap(probabilities);
  for (auto &probability : normalizedMap)
  {
    totalProbability += probability;
  }
  for (auto &probability : normalizedMap)
  {
    probability /= totalProbability;
  }
  return normalizedMap;
}

void MCMoveProbabilities::removeRedundantMoves()
{
  if (probabilities[std::to_underlying(Move::Types::WidomCFCMC)] > 0.0 && 
      probabilities[std::to_underlying(Move::Types::SwapCFCMC)] > 0.0)
  {
    setProbability((Move::Types::WidomCFCMC), 0.0);
  }
  if (probabilities[std::to_underlying(Move::Types::WidomCBCFCMC)] > 0.0 && 
      probabilities[std::to_underlying(Move::Types::SwapCBCFCMC)] > 0.0)
  {
    setProbability((Move::Types::WidomCBCFCMC), 0.0);
  }
}

void MCMoveProbabilities::join(const MCMoveProbabilities &other)
{
  if (*this == other)
  {
    return;
  }
  for (std::size_t i = 0; i < static_cast<std::size_t>(Move::Types::Count); i++)
  {
    if (probabilities[i] > 0.0 && other.probabilities[i] > 0.0)
    {
      throw std::runtime_error("Adding MCMoveProbabilities that define the same value is not permitted.");
    }
    else if (other.probabilities[i] > 0.0)
    {
      probabilities[i] = other.probabilities[i];
    }
  }
}

Move::Types MCMoveProbabilities::sample(RandomNumber &random)
{
  return static_cast<Move::Types>(random.categoricalDistribution(probabilities));
}

std::string MCMoveProbabilities::repr()
{
  std::ostringstream stream;

  std::print(stream, "translationProbability: {}\n", probabilities[std::to_underlying(Move::Types::Translation)]);
  std::print(stream, "randomTranslationProbability: {}\n", probabilities[std::to_underlying(Move::Types::RandomTranslation)]);
  std::print(stream, "rotationProbability: {}\n", probabilities[std::to_underlying(Move::Types::Rotation)]);
  std::print(stream, "randomRotationProbability: {}\n", probabilities[std::to_underlying(Move::Types::RandomRotation)]);
  std::print(stream, "volumeChangeProbability: {}\n", probabilities[std::to_underlying(Move::Types::VolumeChange)]);
  std::print(stream, "reinsertionCBMCProbability: {}\n", probabilities[std::to_underlying(Move::Types::ReinsertionCBMC)]);
  std::print(stream, "partialReinsertionCBMCProbability: {}\n", probabilities[std::to_underlying(Move::Types::PartialReinsertionCBMC)]);
  std::print(stream, "identityChangeCBMCProbability: {}\n", probabilities[std::to_underlying(Move::Types::IdentityChangeCBMC)]);
  std::print(stream, "swapProbability: {}\n", probabilities[std::to_underlying(Move::Types::Swap)]);
  std::print(stream, "swapCBMCProbability: {}\n", probabilities[std::to_underlying(Move::Types::SwapCBMC)]);
  std::print(stream, "swapCFCMCProbability: {}\n", probabilities[std::to_underlying(Move::Types::SwapCFCMC)]);
  std::print(stream, "swapCBCFCMCProbability: {}\n", probabilities[std::to_underlying(Move::Types::SwapCBCFCMC)]);
  std::print(stream, "gibbsVolumeChangeProbability: {}\n", probabilities[std::to_underlying(Move::Types::GibbsVolume)]);
  std::print(stream, "gibbsSwapCBMCProbability: {}\n", probabilities[std::to_underlying(Move::Types::GibbsSwapCBMC)]);
  std::print(stream, "gibbsSwapCFCMCProbability: {}\n", probabilities[std::to_underlying(Move::Types::GibbsSwapCFCMC)]);
  std::print(stream, "widomProbability: {}\n", probabilities[std::to_underlying(Move::Types::Widom)]);
  std::print(stream, "widomCFCMCProbability: {}\n", probabilities[std::to_underlying(Move::Types::WidomCFCMC)]);
  std::print(stream, "widomCBCFCMCProbability: {}\n", probabilities[std::to_underlying(Move::Types::WidomCBCFCMC)]);
  std::print(stream, "parallelTemperingProbability: {}\n", probabilities[std::to_underlying(Move::Types::ParallelTempering)]);
  std::print(stream, "hybridMCProbability: {}\n", probabilities[std::to_underlying(Move::Types::HybridMC)]);

  return stream.str();
}

Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveProbabilities &p)
{
  archive << p.versionNumber;
  archive << p.probabilities;

#if DEBUG_ARCHIVE
  archive << static_cast<std::uint64_t>(0x6f6b6179);  // magic number 'okay' in hex
#endif

  return archive;
}

Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveProbabilities &p)
{
  std::uint64_t versionNumber;
  archive >> versionNumber;
  if (versionNumber > p.versionNumber)
  {
    const std::source_location &location = std::source_location::current();
    throw std::runtime_error(std::format("Invalid version reading 'MCMoveProbabilitiesSystem' at line {} in file {}\n",
                                         location.line(), location.file_name()));
  }
  archive >> p.probabilities;

#if DEBUG_ARCHIVE
  std::uint64_t magicNumber;
  archive >> magicNumber;
  if (magicNumber != static_cast<std::uint64_t>(0x6f6b6179))
  {
    throw std::runtime_error(std::format("MCMoveProbabilities: Error in binary restart\n"));
  }
#endif

  return archive;
}
