module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <fstream>
#include <string>
#endif

export module mc_moves_probabilities_crosssystem;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <chrono>;
import <fstream>;
#endif

import archive;

export struct MCMoveProbabilitiesCrossSystem
{
  uint64_t versionNumber{1};

  MCMoveProbabilitiesCrossSystem() : probability(0.0) {};

  void print();

  double probability;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveProbabilitiesCrossSystem &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveProbabilitiesCrossSystem &p);
};
