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
  /**
   * \brief Manages Monte Carlo move probabilities for the simulation system.
   *
   * The MCMoveProbabilitiesSystem struct holds the probabilities for various Monte Carlo moves at the system level.
   * It includes probabilities for volume change, Gibbs ensemble volume change, and parallel tempering.
   * It provides methods to initialize and optimize these probabilities to improve acceptance rates during simulations.
   */
  uint64_t versionNumber{1};  ///< Version number for serialization purposes.

  bool operator==(MCMoveProbabilitiesSystem const &) const = default;

  /**
   * \brief Constructs an MCMoveProbabilitiesSystem with specified move probabilities.
   *
   * Initializes the system-level Monte Carlo move probabilities for volume change, Gibbs ensemble volume change,
   * and parallel tempering moves.
   *
   * \param volumeChangeProbability Probability for volume change moves.
   * \param gibbsVolumeChangeProbability Probability for Gibbs ensemble volume change moves.
   * \param parallelTemperingProbability Probability for parallel tempering moves.
   */
  MCMoveProbabilitiesSystem(double volumeChangeProbability = 0.0, double gibbsVolumeChangeProbability = 0.0,
                            double parallelTemperingProbability = 0.0);

  double volumeChangeProbability;       ///< Probability for volume change moves.
  double gibbsVolumeChangeProbability;  ///< Probability for Gibbs ensemble volume change moves.
  double parallelTemperingProbability;  ///< Probability for parallel tempering moves.

  /**
   * \brief Optimizes move probabilities based on acceptance rates.
   *
   * Adjusts the move probabilities to improve acceptance rates during simulations.
   */
  void optimizeAcceptance();

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveProbabilitiesSystem &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveProbabilitiesSystem &p);
};
