module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>
#endif

export module monte_carlo_transition_matrix;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <iostream>;
import <fstream>;
import <chrono>;
#endif

import randomnumbers;
import threadpool;

import averages;
import system;
import mc_moves;
import input_reader;
import energy_status;

/**
 * \brief Class that implements the Monte Carlo Transition Matrix method.
 *
 * Manages the simulation workflow, including initialization, equilibration,
 * production, output, and cleanup phases of a Monte Carlo simulation using
 * the Transition Matrix method.
 */
export struct MonteCarloTransitionMatrix
{
  /**
   * \enum WeightingMethod
   * \brief Specifies the weighting method for the Transition Matrix.
   */
  enum class WeightingMethod : size_t
  {
    LambdaZero = 0,  ///< Use only the zero value of lambda for weighting.
    AllLambdas = 1   ///< Use all lambda values for weighting.
  };

  /**
   * \brief Constructs a MonteCarloTransitionMatrix object.
   *
   * Initializes simulation parameters using the provided InputReader.
   *
   * \param reader Reference to an InputReader object containing simulation parameters.
   */
  MonteCarloTransitionMatrix(InputReader& reader) noexcept;

  /**
   * \brief Runs the entire Monte Carlo Transition Matrix simulation.
   *
   * Executes the simulation by calling initialization, equilibration,
   * production, output, and cleanup methods in sequence.
   */
  void run();

  /**
   * \brief Initializes the simulation systems and parameters.
   *
   * Sets up the systems, initializes histograms, creates output directories,
   * and performs initial Monte Carlo moves for system preparation.
   */
  void initialize();

  /**
   * \brief Equilibrates the simulation systems.
   *
   * Performs equilibration cycles to adjust the systems before production.
   * Updates histograms, adjusts biasing factors, and optimizes Monte Carlo moves.
   */
  void equilibrate();

  /**
   * \brief Runs the production phase of the simulation.
   *
   * Performs the main Monte Carlo sampling, collects statistics,
   * and computes averages during the production cycles.
   */
  void production();

  /**
   * \brief Outputs the simulation results and statistics.
   *
   * Writes energy drifts, move statistics, timing information,
   * and other relevant data to the output files.
   */
  void output();

  /**
   * \brief Cleans up resources after the simulation.
   *
   * Performs any necessary cleanup operations after the simulation is complete.
   */
  void cleanup();

  /**
   * \brief Selects and returns a random system from the list of systems.
   *
   * \return Reference to a randomly selected System object.
   */
  System& randomSystem();

  size_t numberOfCycles;                ///< Number of cycles in the production phase.
  size_t numberOfInitializationCycles;  ///< Number of cycles in the initialization phase.
  size_t numberOfEquilibrationCycles;   ///< Number of cycles in the equilibration phase.
  size_t optimizeMCMovesEvery{5000};    ///< Frequency at which Monte Carlo moves are optimized.
  size_t printEvery;                    ///< Frequency at which status is printed to output.

  std::vector<System> systems;         ///< Vector of simulation systems.
  RandomNumber random;                 ///< Random number generator used in the simulation.
  size_t fractionalMoleculeSystem{0};  ///< The system where the fractional molecule is located

  std::vector<std::ofstream> streams;  ///< Output file streams for writing simulation data.

  BlockErrorEstimation estimation;  ///< Object for block error estimation during the simulation.

  std::chrono::system_clock::time_point t1;  ///< Time point marking the start of the production phase.
  std::chrono::system_clock::time_point t2;  ///< Time point marking the end of the production phase.
};
