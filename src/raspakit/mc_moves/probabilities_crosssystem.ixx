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

/**
 * \brief Represents the Monte Carlo move probabilities across systems.
 *
 * The MCMoveProbabilitiesCrossSystem struct holds the probability data used for cross-system
 * Monte Carlo moves. It includes methods for printing and serialization, and keeps track of the
 * version number for compatibility.
 */
export struct MCMoveProbabilitiesCrossSystem
{
  uint64_t versionNumber{1};  ///< Version number for serialization compatibility.

  /**
   * \brief Default constructor initializing the probability to zero.
   *
   * Initializes an MCMoveProbabilitiesCrossSystem object with probability set to 0.0.
   */
  MCMoveProbabilitiesCrossSystem() : probability(0.0) {};

  /**
   * \brief Prints the probability to the standard output.
   *
   * Outputs the probability value to the console.
   */
  void print();

  double probability;  ///< The probability of the Monte Carlo move.

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveProbabilitiesCrossSystem &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveProbabilitiesCrossSystem &p);
};
