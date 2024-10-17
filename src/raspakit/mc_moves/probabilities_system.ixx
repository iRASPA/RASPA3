module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <fstream>
#include <string>
#endif

export module mc_moves_probabilities_system;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <chrono>;
import <fstream>;
#endif

import archive;
import double3;

export struct MCMoveProbabilitiesSystem
{
  uint64_t versionNumber{1};

  bool operator==(MCMoveProbabilitiesSystem const &) const = default;

  MCMoveProbabilitiesSystem(double volumeChangeProbability = 0.0, double gibbsVolumeChangeProbability = 0.0,
                            double parallelTemperingProbability = 0.0, double hybridMCProbability = 0.0);

  double volumeChangeProbability;
  double gibbsVolumeChangeProbability;
  double parallelTemperingProbability;
  double hybridMCProbability;

  void optimizeAcceptance();

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveProbabilitiesSystem &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveProbabilitiesSystem &p);
};
