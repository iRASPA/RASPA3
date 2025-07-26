module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <map>
#include <random>
#include <vector>
#endif

export module mc_moves_probabilities;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import randomnumbers;
import mc_moves_move_types;

export struct MCMoveProbabilities
{
  uint64_t versionNumber{2};

  bool operator==(MCMoveProbabilities const &) const = default;

  // vector of unnormalized probabilities (not necessary due to std::discrete_distribution)
  std::map<MoveTypes, double> probabilities;

  MCMoveProbabilities(double translationProbability = 0.0, double randomTranslationProbability = 0.0,
                      double rotationProbability = 0.0, double randomRotationProbability = 0.0,
                      double volumeChangeProbability = 0.0, double reinsertionCBMCProbability = 0.0,
                      double identityChangeCBMCProbability = 0.0, double swapProbability = 0.0,
                      double swapCBMCProbability = 0.0, double swapCFCMCProbability = 0.0,
                      double swapCBCFCMCProbability = 0.0, double gibbsVolumeChangeProbability = 0.0,
                      double gibbsSwapCBMCProbability = 0.0, double gibbsSwapCFCMCProbability = 0.0,
                      double widomProbability = 0.0, double widomCFCMCProbability = 0.0,
                      double widomCBCFCMCProbability = 0.0, double parallelTemperingProbability = 0.0,
                      double hybridMCProbability = 0.0);

  // use .at such that access is const
  double getProbability(const MoveTypes &move) const { return probabilities.at(move); };
  void setProbability(const MoveTypes &move, double probability) { probabilities[move] = probability; };

  const std::map<MoveTypes, double> normalizedMap() const;
  void removeRedundantMoves();
  MoveTypes sample(RandomNumber &random);

  void join(const MCMoveProbabilities &other);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveProbabilities &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveProbabilities &p);
};
