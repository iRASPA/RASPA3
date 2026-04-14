module;

export module mc_moves_probabilities;

import std;

import archive;
import randomnumbers;
import mc_moves_move_types;

export struct MCMoveProbabilities
{
  std::uint64_t versionNumber{2};

  bool operator==(MCMoveProbabilities const &) const = default;

  MCMoveProbabilities(double translationProbability = 0.0, double randomTranslationProbability = 0.0,
                      double rotationProbability = 0.0, double randomRotationProbability = 0.0,
                      double volumeChangeProbability = 0.0, double reinsertionCBMCProbability = 0.0,
                      double partialReinsertionCBMCProbability = 0.0, double identityChangeCBMCProbability = 0.0,
                      double swapProbability = 0.0, double swapCBMCProbability = 0.0, double swapCFCMCProbability = 0.0,
                      double swapCBCFCMCProbability = 0.0, double gibbsVolumeChangeProbability = 0.0,
                      double gibbsSwapCBMCProbability = 0.0, double gibbsSwapCFCMCProbability = 0.0,
                      double widomProbability = 0.0, double widomCFCMCProbability = 0.0,
                      double widomCBCFCMCProbability = 0.0, double parallelTemperingProbability = 0.0,
                      double hybridMCProbability = 0.0);

  // use .at such that access is const
  double getProbability(const Move::Types &move) const { return probabilities[std::to_underlying(move)]; };
  void setProbability(const Move::Types &move, double probability) { probabilities[std::to_underlying(move)] = probability; };

  const std::vector<double> normalizedMap() const;
  void removeRedundantMoves();
  Move::Types sample(RandomNumber &random);

  void join(const MCMoveProbabilities &other);

  // vector of unnormalized probabilities (not necessary due to std::discrete_distribution)
  std::vector<double> probabilities{};

  std::string repr();

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveProbabilities &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveProbabilities &p);
};
