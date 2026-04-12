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
                                         double parallelTemperingProbability, double hybridMCProbability)
{
  probabilities[Move::Types::Translation] = translationProbability;
  probabilities[Move::Types::RandomTranslation] = randomTranslationProbability;
  probabilities[Move::Types::Rotation] = rotationProbability;
  probabilities[Move::Types::RandomRotation] = randomRotationProbability;
  probabilities[Move::Types::VolumeChange] = volumeChangeProbability;
  probabilities[Move::Types::ReinsertionCBMC] = reinsertionCBMCProbability;
  probabilities[Move::Types::PartialReinsertionCBMC] = partialReinsertionCBMCProbability;
  probabilities[Move::Types::IdentityChangeCBMC] = identityChangeCBMCProbability;
  probabilities[Move::Types::Swap] = swapProbability;
  probabilities[Move::Types::SwapCBMC] = swapCBMCProbability;
  probabilities[Move::Types::SwapCFCMC] = swapCFCMCProbability;
  probabilities[Move::Types::SwapCBCFCMC] = swapCBCFCMCProbability;
  probabilities[Move::Types::GibbsVolume] = gibbsVolumeChangeProbability;
  probabilities[Move::Types::GibbsSwapCBMC] = gibbsSwapCBMCProbability;
  probabilities[Move::Types::GibbsSwapCFCMC] = gibbsSwapCFCMCProbability;
  probabilities[Move::Types::Widom] = widomProbability;
  probabilities[Move::Types::WidomCFCMC] = widomCFCMCProbability;
  probabilities[Move::Types::WidomCBCFCMC] = widomCBCFCMCProbability;
  probabilities[Move::Types::ParallelTempering] = parallelTemperingProbability;
  probabilities[Move::Types::HybridMC] = hybridMCProbability;
}

const std::unordered_map<Move::Types, double> MCMoveProbabilities::normalizedMap() const
{
  double totalProbability = 0.0;
  std::unordered_map<Move::Types, double> normalizedMap(probabilities);
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
  if (probabilities[Move::Types::WidomCFCMC] > 0.0 && probabilities[Move::Types::SwapCFCMC] > 0.0)
  {
    setProbability((Move::Types::WidomCFCMC), 0.0);
  }
  if (probabilities[Move::Types::WidomCBCFCMC] > 0.0 && probabilities[Move::Types::SwapCBCFCMC] > 0.0)
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
    if (probabilities[static_cast<Move::Types>(i)] > 0.0 && other.getProbability(static_cast<Move::Types>(i)) > 0.0)
    {
      throw std::runtime_error("Adding MCMoveProbabilities that define the same value is not permitted.");
    }
    else if (other.getProbability(static_cast<Move::Types>(i)) > 0.0)
    {
      probabilities[static_cast<Move::Types>(i)] = other.getProbability(static_cast<Move::Types>(i));
    }
  }
}

Move::Types MCMoveProbabilities::sample(RandomNumber &random)
{
  std::vector<double> vectorProbabilities(static_cast<std::size_t>(Move::Types::Count));
  for (std::size_t i = 0; i < static_cast<std::size_t>(Move::Types::Count); i++)
  {
    vectorProbabilities[i] = probabilities[static_cast<Move::Types>(i)];
  }
  return static_cast<Move::Types>(random.categoricalDistribution(vectorProbabilities));
}

std::string MCMoveProbabilities::repr()
{
  std::ostringstream stream;

  std::print(stream, "translationProbability: {}\n", probabilities[Move::Types::Translation]);
  std::print(stream, "randomTranslationProbability: {}\n", probabilities[Move::Types::RandomTranslation]);
  std::print(stream, "rotationProbability: {}\n", probabilities[Move::Types::Rotation]);
  std::print(stream, "randomRotationProbability: {}\n", probabilities[Move::Types::RandomRotation]);
  std::print(stream, "volumeChangeProbability: {}\n", probabilities[Move::Types::VolumeChange]);
  std::print(stream, "reinsertionCBMCProbability: {}\n", probabilities[Move::Types::ReinsertionCBMC]);
  std::print(stream, "partialReinsertionCBMCProbability: {}\n", probabilities[Move::Types::PartialReinsertionCBMC]);
  std::print(stream, "identityChangeCBMCProbability: {}\n", probabilities[Move::Types::IdentityChangeCBMC]);
  std::print(stream, "swapProbability: {}\n", probabilities[Move::Types::Swap]);
  std::print(stream, "swapCBMCProbability: {}\n", probabilities[Move::Types::SwapCBMC]);
  std::print(stream, "swapCFCMCProbability: {}\n", probabilities[Move::Types::SwapCFCMC]);
  std::print(stream, "swapCBCFCMCProbability: {}\n", probabilities[Move::Types::SwapCBCFCMC]);
  std::print(stream, "gibbsVolumeChangeProbability: {}\n", probabilities[Move::Types::GibbsVolume]);
  std::print(stream, "gibbsSwapCBMCProbability: {}\n", probabilities[Move::Types::GibbsSwapCBMC]);
  std::print(stream, "gibbsSwapCFCMCProbability: {}\n", probabilities[Move::Types::GibbsSwapCFCMC]);
  std::print(stream, "widomProbability: {}\n", probabilities[Move::Types::Widom]);
  std::print(stream, "widomCFCMCProbability: {}\n", probabilities[Move::Types::WidomCFCMC]);
  std::print(stream, "widomCBCFCMCProbability: {}\n", probabilities[Move::Types::WidomCBCFCMC]);
  std::print(stream, "parallelTemperingProbability: {}\n", probabilities[Move::Types::ParallelTempering]);
  std::print(stream, "hybridMCProbability: {}\n", probabilities[Move::Types::HybridMC]);

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
