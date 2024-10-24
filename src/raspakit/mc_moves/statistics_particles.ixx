module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <fstream>
#include <string>
#endif

export module mc_moves_statistics_particles;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <chrono>;
import <fstream>;
#endif

import double3;
import archive;
import move_statistics;
import json;

/**
 * \brief Stores statistics for various Monte Carlo moves involving particles.
 *
 * This struct keeps track of counts, acceptances, and other statistics for different types of Monte Carlo
 * moves used in particle simulations, such as translations, rotations, insertions, deletions, and swaps.
 */
export struct MCMoveStatisticsParticles
{
  uint64_t versionNumber{1};  ///< Version number for serialization compatibility.

  bool operator==(MCMoveStatisticsParticles const &) const = default;

  MoveStatistics<double3> translationMove{.maxChange = double3(1.0, 1.0, 1.0)};  ///< Statistics for translation moves.
  MoveStatistics<double3> randomTranslationMove{};  ///< Statistics for random translation moves.
  MoveStatistics<double3> rotationMove{.maxChange = double3(1.0, 1.0, 1.0)};  ///< Statistics for rotation moves.
  MoveStatistics<double3> randomRotationMove{};                               ///< Statistics for random rotation moves.
  MoveStatistics<double> reinsertionMove_CBMC{};     ///< Statistics for CBMC reinsertion moves.
  MoveStatistics<double> identityChangeMove_CBMC{};  ///< Statistics for CBMC identity change moves.
  MoveStatistics<double> swapInsertionMove{};        ///< Statistics for swap insertion moves.
  MoveStatistics<double> swapDeletionMove{};         ///< Statistics for swap deletion moves.
  MoveStatistics<double> swapInsertionMove_CBMC{};   ///< Statistics for CBMC swap insertion moves.
  MoveStatistics<double> swapDeletionMove_CBMC{};    ///< Statistics for CBMC swap deletion moves.
  MoveStatistics<double3> swapMove_CFCMC{.maxChange = double3(0.0, 0.0, 0.5)};  ///< Statistics for CFCMC swap moves.
  MoveStatistics<double3> swapMove_CFCMC_CBMC{.maxChange =
                                                  double3(0.0, 0.0, 0.5)};  ///< Statistics for CBMC/CFCMC swap moves.
  MoveStatistics<double> WidomMove_CBMC{};                                  ///< Statistics for CBMC Widom moves.
  MoveStatistics<double> WidomMove_CFCMC{};                                 ///< Statistics for CFCMC Widom moves.
  MoveStatistics<double> WidomMove_CFCMC_CBMC{};                            ///< Statistics for CBMC/CFCMC Widom moves.

  MoveStatistics<double> GibbsSwapMove_CBMC{};  ///< Statistics for CBMC Gibbs swap moves.
  MoveStatistics<double3> GibbsSwapMove_CFCMC{.maxChange =
                                                  double3(0.0, 0.0, 0.5)};  ///< Statistics for CFCMC Gibbs swap moves.
  MoveStatistics<double3> GibbsSwapMove_CFCMC_CBMC{
      .maxChange = double3(0.0, 0.0, 0.5)};  ///< Statistics for CBMC/CFCMC Gibbs swap moves.

  /**
   * \brief Clears the move statistics by resetting all counts to zero.
   */
  void clearMoveStatistics();

  /**
   * \brief Optimizes move parameters based on acceptance rates.
   *
   * Adjusts parameters of certain Monte Carlo moves to optimize their acceptance rates.
   */
  void optimizeMCMoves();

  /**
   * \brief Generates a formatted string of the move statistics.
   *
   * \return A string representing the statistics of the Monte Carlo moves.
   */
  const std::string writeMCMoveStatistics() const;

  /**
   * \brief Returns the move statistics in JSON format.
   *
   * \return A JSON object containing the statistics of the Monte Carlo moves.
   */
  const nlohmann::json jsonMCMoveStatistics() const;

  /**
   * \brief Generates a formatted string of the move statistics for a specific component.
   *
   * \param countTotal The total count of moves.
   * \param componentId The ID of the component.
   * \param componentName The name of the component.
   * \return A string representing the statistics of the Monte Carlo moves for the specified component.
   */
  const std::string writeMCMoveStatistics(size_t countTotal, size_t componentId,
                                          const std::string &componentName) const;

  /**
   * \brief Returns the move statistics in JSON format for a specific component.
   *
   * \param countTotal The total count of moves.
   * \return A JSON object containing the statistics of the Monte Carlo moves for the specified component.
   */
  const nlohmann::json jsonMCMoveStatistics(size_t countTotal) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MCMoveStatisticsParticles &p);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MCMoveStatisticsParticles &p);
};
