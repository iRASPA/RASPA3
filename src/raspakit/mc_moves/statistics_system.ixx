module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <fstream>
#include <string>
#endif

export module mc_moves_statistics_system;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <chrono>;
import <fstream>;
#endif

import archive;
import double3;
import move_statistics;
import json;

/**
 * \brief Holds statistics for Monte Carlo moves within a simulation system.
 *
 * The MCMoveStatisticsSystem struct keeps track of the counts and acceptance statistics
 * for different types of Monte Carlo moves, including volume changes, Gibbs volume changes,
 * and parallel tempering swaps. It provides methods to optimize acceptance rates, clear
 * statistics, and output statistics in string or JSON format.
 */
export struct MCMoveStatisticsSystem
{
  uint64_t versionNumber{1};  ///< Version number for serialization.

  bool operator==(MCMoveStatisticsSystem const &) const = default;

  MoveStatistics<double> volumeMove{.maxChange = 0.1};             ///< Statistics for volume move.
  MoveStatistics<double> GibbsVolumeMove{.maxChange = 0.1};        ///< Statistics for Gibbs volume move.
  MoveStatistics<double> ParallelTemperingSwap{.maxChange = 0.1};  ///< Statistics for parallel tempering swap.

  /**
   * \brief Optimizes the acceptance rates of the moves.
   *
   * Adjusts the maximum change parameters for each move type to optimize their acceptance rates.
   * The maxChange parameter is adjusted within a specified range to achieve desired acceptance ratios.
   */
  void optimizeAcceptance();

  /**
   * \brief Clears the move statistics.
   *
   * Resets the counts and acceptance statistics for all move types.
   */
  void clear();

  /**
   * \brief Writes the move statistics to a string.
   *
   * Generates a formatted string containing the statistics of all the moves.
   *
   * \return A string containing move statistics.
   */
  const std::string writeMCMoveStatistics() const;

  /**
   * \brief Outputs the move statistics in JSON format.
   *
   * Generates a JSON object containing the statistics of all the moves.
   *
   * \return A nlohmann::json object containing move statistics.
   */
  const nlohmann::json jsonMCMoveStatistics() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveStatisticsSystem &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveStatisticsSystem &p);
};
