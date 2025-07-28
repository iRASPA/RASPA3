module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <optional>
#include <vector>
#endif

export module monte_carlo_transition_matrix;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import randomnumbers;
import threadpool;

import averages;
import system;
import mc_moves;
import input_reader;
import energy_status;
import archive;
import json;

/**
 * \brief Performs Monte Carlo simulations for molecular systems.
 *
 * The MonteCarloTransitionMatrix struct orchestrates the execution of Monte Carlo simulations,
 * including initialization, equilibration, and production stages. It manages multiple systems,
 * random number generation, simulation parameters, and output generation.
 */
export struct MonteCarloTransitionMatrix
{
  /**
   * \brief Enumeration of weighting methods for sampling.
   */
  enum class WeightingMethod : std::size_t
  {
    LambdaZero = 0,  ///< Use only lambda = 0 in sampling.
    AllLambdas = 1   ///< Use all lambda values in sampling.
  };

  /**
   * \brief Enumeration of simulation stages.
   */
  enum class SimulationStage : std::size_t
  {
    Uninitialized = 0,   ///< Simulation not initialized.
    Initialization = 1,  ///< Initialization stage.
    Equilibration = 2,   ///< Equilibration stage.
    Production = 3       ///< Production stage.
  };

  /**
   * \brief Default constructor for the MonteCarloTransitionMatrix class.
   *
   * Initializes a MonteCarloTransitionMatrix object with default parameters.
   */
  MonteCarloTransitionMatrix();

  /**
   * \brief Constructs a MonteCarloTransitionMatrix object with simulation parameters from an InputReader.
   *
   * \param reader InputReader containing simulation parameters.
   */
  MonteCarloTransitionMatrix(InputReader &reader) noexcept;

  /**
   * \brief Constructs a MonteCarloTransitionMatrix object with specified simulation parameters.
   *
   * \param numberOfCycles Number of production cycles.
   * \param numberOfInitializationCycles Number of initialization cycles.
   * \param numberOfEquilibrationCycles Number of equilibration cycles.
   * \param printEvery Frequency of printing status reports.
   * \param writeBinaryRestartEvery Frequency of writing binary restart files.
   * \param rescaleWangLandauEvery Frequency of rescaling Wang-Landau factors.
   * \param optimizeMCMovesEvery Frequency of optimizing MC moves.
   * \param systems Vector of System objects to simulate.
   * \param randomSeed Random number generator seed.
   * \param numberOfBlocks Number of blocks for error estimation.
   */
  MonteCarloTransitionMatrix(std::size_t numberOfCycles, std::size_t numberOfInitializationCycles,
                             std::size_t numberOfEquilibrationCycles, std::size_t printEvery,
                             std::size_t writeBinaryRestartEvery, std::size_t rescaleWangLandauEvery,
                             std::size_t optimizeMCMovesEvery, std::vector<System> &systems, RandomNumber &randomSeed,
                             std::size_t numberOfBlocks);

  std::uint64_t versionNumber{1};  ///< Version number for serialization.

  std::size_t numberOfCycles;                ///< Number of production cycles.
  std::size_t numberOfSteps;                 ///< Total number of steps performed.
  std::size_t numberOfInitializationCycles;  ///< Number of initialization cycles.
  std::size_t numberOfEquilibrationCycles;   ///< Number of equilibration cycles.
  std::size_t printEvery;                    ///< Frequency of printing status reports.
  std::size_t writeBinaryRestartEvery;       ///< Frequency of writing binary restart files.
  std::size_t rescaleWangLandauEvery;        ///< Frequency of rescaling Wang-Landau factors.
  std::size_t optimizeMCMovesEvery;          ///< Frequency of optimizing MC moves.

  std::size_t currentCycle{0};                                      ///< Current cycle number.
  SimulationStage simulationStage{SimulationStage::Uninitialized};  ///< Current simulation stage.

  std::vector<System> systems;              ///< Vector of systems to simulate.
  RandomNumber random;                      ///< Random number generator.
  std::size_t fractionalMoleculeSystem{0};  // the system where the fractional molecule is located

  std::vector<std::ofstream> streams;            ///< Output streams for writing data.
  std::vector<std::string> outputJsonFileNames;  ///< Filenames for output JSON files.
  std::vector<nlohmann::json> outputJsons;       ///< Output data in JSON format.

  BlockErrorEstimation estimation;  ///< Block error estimation object.

  std::chrono::duration<double> totalInitializationSimulationTime{0};  ///< Total time for initialization stage.
  std::chrono::duration<double> totalEquilibrationSimulationTime{0};   ///< Total time for equilibration stage.
  std::chrono::duration<double> totalProductionSimulationTime{0};      ///< Total time for production stage.
  std::chrono::duration<double> totalSimulationTime{0};                ///< Total simulation time.

  /**
   * \brief Creates output files for writing simulation data.
   */
  void createOutputFiles();

  /**
   * \brief Runs the Monte Carlo simulation.
   *
   * Orchestrates the simulation by executing initialization, equilibration,
   * and production stages sequentially.
   */
  void run();

  /**
   * \brief Performs a single Monte Carlo cycle.
   *
   * Executes a number of Monte Carlo steps, updates system properties,
   * and handles output generation and restart file writing.
   */
  void performCycle();

  /**
   * \brief Performs the initialization stage of the simulation.
   *
   * Sets up the simulation systems, writes initial output, and runs
   * the specified number of initialization cycles.
   */
  void initialize();

  /**
   * \brief Performs the equilibration stage of the simulation.
   *
   * Equilibrates the systems by running the specified number of equilibration cycles,
   * and adjusts Wang-Landau biasing factors.
   */
  void equilibrate();

  /**
   * \brief Performs the production stage of the simulation.
   *
   * Runs the main simulation cycles, collects statistics, and samples properties.
   */
  void production();

  /**
   * \brief Generates the final output of the simulation.
   *
   * Writes energy statistics, move statistics, CPU timings, and
   * other properties to output files.
   */
  void output();

  /**
   * \brief Selects a random system from the list of systems.
   *
   * \return Reference to a randomly selected System object.
   */
  System &randomSystem();

  /**
   * \brief Returns a string representation of the MonteCarloTransitionMatrix object.
   *
   * \return A string describing the MonteCarloTransitionMatrix object.
   */
  std::string repr() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MonteCarloTransitionMatrix &mc);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MonteCarloTransitionMatrix &mc);
};
