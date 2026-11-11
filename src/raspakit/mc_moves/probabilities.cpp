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
                                         double gibbsSwapCFCMCProbability, double gibbsIdentityChangeCBMCProbability,
                                         double widomProbability,
                                         double widomCFCMCProbability, double widomCBCFCMCProbability,
                                         double parallelTemperingProbability, double hybridMCProbability,
                                         double reactionCBMCProbability, double reactionConventionalCFCMCProbability,
                                         double reactionConventionalCBCFCMCProbability):
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
  probabilities[std::to_underlying(Move::Types::GibbsIdentityChangeCBMC)] = gibbsIdentityChangeCBMCProbability;
  probabilities[std::to_underlying(Move::Types::Widom)] = widomProbability;
  probabilities[std::to_underlying(Move::Types::WidomCFCMC)] = widomCFCMCProbability;
  probabilities[std::to_underlying(Move::Types::WidomCBCFCMC)] = widomCBCFCMCProbability;
  probabilities[std::to_underlying(Move::Types::ParallelTempering)] = parallelTemperingProbability;
  probabilities[std::to_underlying(Move::Types::HybridMC)] = hybridMCProbability;
  probabilities[std::to_underlying(Move::Types::ReactionCBMC)] = reactionCBMCProbability;
  probabilities[std::to_underlying(Move::Types::ReactionConventionalCFCMC)] = reactionConventionalCFCMCProbability;
  probabilities[std::to_underlying(Move::Types::ReactionConventionalCBCFCMC)] = reactionConventionalCBCFCMCProbability;
}

const std::vector<double> MCMoveProbabilities::normalizedMap() const
{
  double totalProbability = 0.0;
  std::vector<double> normalizedMap(probabilities);
  for (auto &probability : normalizedMap)
  {
    totalProbability += probability;
  }

  // a component without any move (all probabilities zero) has no distribution to normalize: all zeros, not 0/0
  if (!(totalProbability > 0.0)) return normalizedMap;

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
  for (std::size_t i = 0; i != static_cast<std::size_t>(Move::Types::Count); ++i)
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
  // A table without any move (an immobile species in a system without system moves) has no distribution to draw
  // from: 'std::discrete_distribution' over all-zero weights draws one fixed, meaningless index. No move.
  if (std::ranges::none_of(probabilities, [](double probability) { return probability > 0.0; }))
  {
    return Move::Types::Count;
  }
  return static_cast<Move::Types>(random.categoricalDistribution(probabilities));
}

std::string MCMoveProbabilities::repr()
{
  std::ostringstream stream;

  std::vector<double> normalized_probabilities = normalizedMap();

  std::print(stream, "translationProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::Translation)]);
  std::print(stream, "randomTranslationProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::RandomTranslation)]);
  std::print(stream, "translationSmartMCProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::TranslationSmartMC)]);
  std::print(stream, "translationSmartMCAllProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::TranslationSmartMCAll)]);
  std::print(stream, "rotationProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::Rotation)]);
  std::print(stream, "randomRotationProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::RandomRotation)]);
  std::print(stream, "rotationSmartMCProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::RotationSmartMC)]);
  std::print(stream, "rotationSmartMCAllProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::RotationSmartMCAll)]);
  std::print(stream, "translationRotationSmartMCProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::TranslationRotationSmartMC)]);
  std::print(stream, "volumeChangeProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::VolumeChange)]);
  std::print(stream, "reinsertionCBMCProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::ReinsertionCBMC)]);
  std::print(stream, "partialReinsertionCBMCProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::PartialReinsertionCBMC)]);
  std::print(stream, "identityChangeCBMCProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::IdentityChangeCBMC)]);
  std::print(stream, "identitySwitchCBMCProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::IdentitySwitchCBMC)]);
  std::print(stream, "swapProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::Swap)]);
  std::print(stream, "swapCBMCProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::SwapCBMC)]);
  std::print(stream, "swapMultiTrialCBMCProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::SwapMultiTrialCBMC)]);
  std::print(stream, "reinsertionMultiTrialCBMCProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::ReinsertionMultiTrialCBMC)]);
  std::print(stream, "swapCFCMCProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::SwapCFCMC)]);
  std::print(stream, "swapCBCFCMCProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::SwapCBCFCMC)]);
  std::print(stream, "gibbsVolumeChangeProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::GibbsVolume)]);
  std::print(stream, "gibbsSwapCBMCProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::GibbsSwapCBMC)]);
  std::print(stream, "gibbsSwapMultiTrialCBMCProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::GibbsSwapMultiTrialCBMC)]);
  std::print(stream, "gibbsSwapCFCMCProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::GibbsSwapCFCMC)]);
  std::print(stream, "gibbsIdentityChangeCBMCProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::GibbsIdentityChangeCBMC)]);
  std::print(stream, "widomProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::Widom)]);
  std::print(stream, "widomCFCMCProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::WidomCFCMC)]);
  std::print(stream, "widomCBCFCMCProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::WidomCBCFCMC)]);
  std::print(stream, "parallelTemperingProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::ParallelTempering)]);
  std::print(stream, "hybridMCProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::HybridMC)]);
  std::print(stream, "crossLinkSwapProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::CrossLinkSwap)]);
  std::print(stream, "crossLinkFormationProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::CrossLinkFormationScission)]);
  std::print(stream, "crossLinkExchangeProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::CrossLinkExchange)]);
  std::print(stream, "reactionCBMCProbability: {}\n",
             normalized_probabilities[std::to_underlying(Move::Types::ReactionCBMC)]);
  std::print(stream, "reactionConventionalCFCMCProbability: {}\n",
             normalized_probabilities[std::to_underlying(Move::Types::ReactionConventionalCFCMC)]);
  std::print(stream, "reactionConventionalCBCFCMCProbability: {}\n",
             normalized_probabilities[std::to_underlying(Move::Types::ReactionConventionalCBCFCMC)]);
  std::print(stream, "reactionCFCMCProbability: {}\n",
             normalized_probabilities[std::to_underlying(Move::Types::ReactionCFCMC)]);
  std::print(stream, "reactionCBCFCMCProbability: {}\n",
             normalized_probabilities[std::to_underlying(Move::Types::ReactionCBCFCMC)]);
  std::print(stream, "anisotropicVolumeChangeProbability: {}\n",
             normalized_probabilities[std::to_underlying(Move::Types::AnisotropicVolumeChange)]);
  std::print(stream, "gibbsConventionalCFCMCProbability: {}\n",
             normalized_probabilities[std::to_underlying(Move::Types::GibbsConventionalCFCMC)]);
  std::print(stream, "gibbsConventionalCBCFCMCProbability: {}\n",
             normalized_probabilities[std::to_underlying(Move::Types::GibbsConventionalCBCFCMC)]);
  std::print(stream, "gibbsSwapCBCFCMCProbability: {}\n",
             normalized_probabilities[std::to_underlying(Move::Types::GibbsSwapCBCFCMC)]);
  std::print(stream, "pairSwapProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::PairSwap)]);
  std::print(stream, "pairSwapCBMCProbability: {}\n",
             normalized_probabilities[std::to_underlying(Move::Types::PairSwapCBMC)]);
  std::print(stream, "pairSwapCFCMCProbability: {}\n",
             normalized_probabilities[std::to_underlying(Move::Types::PairSwapCFCMC)]);
  std::print(stream, "pairSwapCBCFCMCProbability: {}\n",
             normalized_probabilities[std::to_underlying(Move::Types::PairSwapCBCFCMC)]);
  std::print(stream, "groupSwapProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::GroupSwap)]);
  std::print(stream, "groupSwapCBMCProbability: {}\n",
             normalized_probabilities[std::to_underlying(Move::Types::GroupSwapCBMC)]);
  std::print(stream, "groupSwapCFCMCProbability: {}\n",
             normalized_probabilities[std::to_underlying(Move::Types::GroupSwapCFCMC)]);
  std::print(stream, "groupSwapCBCFCMCProbability: {}\n",
             normalized_probabilities[std::to_underlying(Move::Types::GroupSwapCBCFCMC)]);
  std::print(stream, "tetheredProtonHopProbability: {}\n",
             normalized_probabilities[std::to_underlying(Move::Types::TetheredProtonHop)]);
  std::print(stream, "pivotProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::Pivot)]);
  std::print(stream, "crankshaftProbability: {}\n",
             normalized_probabilities[std::to_underlying(Move::Types::Crankshaft)]);
  std::print(stream, "reptationProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::Reptation)]);
  std::print(stream, "concertedRotationProbability: {}\n",
             normalized_probabilities[std::to_underlying(Move::Types::ConcertedRotation)]);
  std::print(stream, "beadDisplacementProbability: {}\n",
             normalized_probabilities[std::to_underlying(Move::Types::BeadDisplacement)]);
  std::print(stream, "beadFlipProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::BeadFlip)]);
  std::print(stream, "doubleBridgingProbability: {}\n",
             normalized_probabilities[std::to_underlying(Move::Types::DoubleBridging)]);
  std::print(stream, "intramolecularDoubleRebridgingProbability: {}\n",
             normalized_probabilities[std::to_underlying(Move::Types::IntramolecularDoubleRebridging)]);
  std::print(stream, "pivotCBMCProbability: {}\n", normalized_probabilities[std::to_underlying(Move::Types::PivotCBMC)]);

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
