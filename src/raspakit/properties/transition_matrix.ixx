module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <fstream>
#include <vector>
#endif

export module transition_matrix;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import double3;

/**
 * \brief Represents the transition matrix for TMMC simulations.
 *
 * The TransitionMatrix struct encapsulates the data and methods required for performing
 * Transition Matrix Monte Carlo (TMMC) simulations. It maintains the collection matrix,
 * bias factors, natural logarithm of probability distributions, and histograms to
 * compute macrostates and adjust biases during the simulation. It provides methods
 * to initialize data structures, update matrices and histograms, calculate bias factors,
 * adjust biases, clear statistics, and write statistics to files.
 */
export struct TransitionMatrix
{
  std::uint64_t versionNumber{1};  ///< Version number for serialization compatibility.

  bool operator==(TransitionMatrix const &) const = default;

  std::vector<double3> cmatrix;        ///< Collection matrix: x=deletion, y=no change, z=insertion.
  std::vector<double> bias;            ///< Bias factors for each macrostate.
  std::vector<double> lnpi;            ///< Natural logarithm of the probability distribution.
  std::vector<double> forward_lnpi;    ///< Forward ln(pi) for debugging purposes.
  std::vector<double> reverse_lnpi;    ///< Reverse ln(pi) for debugging purposes.
  std::vector<std::size_t> histogram;  ///< Histogram of macrostate visits.

  std::size_t numberOfSteps = {0};        ///< Number of steps performed in the TMMC simulation.
  std::size_t minMacrostate = {0};        ///< Minimum macrostate value.
  std::size_t maxMacrostate = {100};      ///< Maximum macrostate value.
  std::size_t updateTMEvery = {1000000};  ///< Number of steps between bias updates.
  std::size_t numberOfUpdates = {0};      ///< Number of times the bias has been updated.

  bool doTMMC = {false};                     ///< Flag indicating whether to perform TMMC simulation.
  bool useBias = {false};                    ///< Flag indicating whether to use bias for changing macrostates.
  bool useTMBias = {true};                   ///< Flag indicating whether to use Transition Matrix bias.
  bool rejectOutofBound = {true};            ///< Flag indicating whether to reject moves outside macrostate bounds.
  bool rezeroAfterInitialization = {false};  ///< Flag indicating whether to reset statistics after initialization.

  /**
   * \brief Initializes the transition matrix and related data structures.
   *
   * Sets up the collection matrix, bias vectors, and histograms based on the specified
   * macrostate range. Initializes bias factors and probability distributions to default values.
   * Should be called before starting the TMMC simulation.
   */
  void initialize();

  /**
   * \brief Updates the collection matrix with acceptance probabilities.
   *
   * Records the probabilities of transitions between macrostates in the collection matrix.
   * The acceptance probabilities vector Pacc contains probabilities for deletion (x),
   * no change (y), and insertion (z) moves.
   *
   * \param Pacc A double3 vector containing acceptance probabilities for deletion (x),
   *             no change (y), and insertion (z) moves.
   * \param oldN The macrostate value before the move.
   */
  void updateMatrix(double3 Pacc, std::size_t oldN);

  /**
   * \brief Updates the histogram of macrostate visits.
   *
   * Increments the histogram count for the given macrostate N, tracking how often each
   * macrostate is visited during the simulation.
   *
   * \param N The current macrostate value.
   */
  void updateHistogram(std::size_t N);

  /**
   * \brief Calculates the bias factor between two macrostates.
   *
   * Computes the bias factor used to adjust the acceptance probability of moves
   * between macrostates. The bias factor is calculated as the exponential of the
   * difference in bias between the new and old macrostate.
   *
   * \param newN The macrostate value after the move.
   * \param oldN The macrostate value before the move.
   * \return The bias factor for the transition from oldN to newN.
   */
  double biasFactor(std::size_t newN, std::size_t oldN);

  /**
   * \brief Adjusts the bias factors based on collected statistics.
   *
   * Updates the bias vector and the natural logarithm of the probability distribution
   * (lnpi) using the current state of the collection matrix and histogram. This method
   * recalculates biases to improve sampling efficiency.
   */
  void adjustBias();

  /**
   * \brief Clears the collection matrix and resets counters.
   *
   * Resets the collection matrix, histogram, bias factors, and related counters to
   * their initial values. Typically used after initialization cycles to reset statistics.
   */
  void clearCMatrix();

  /**
   * \brief Writes the transition matrix statistics to a file.
   *
   * Outputs the current state of the collection matrix, bias factors, probability
   * distributions, and histogram data to a text file for analysis and debugging.
   */
  void writeStatistics();

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const TransitionMatrix &m);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, TransitionMatrix &m);
};
