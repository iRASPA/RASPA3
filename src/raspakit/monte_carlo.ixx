module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <optional>
#include <vector>
#endif

export module monte_carlo;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <iostream>;
import <fstream>;
import <chrono>;
import <optional>;
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
 * The MonteCarlo struct orchestrates the execution of Monte Carlo simulations,
 * including initialization, equilibration, and production stages. It manages multiple systems,
 * random number generation, simulation parameters, and output generation.
 */
export struct MonteCarlo
{
  /**
   * \brief Enumeration of weighting methods for sampling.
   */
  enum class WeightingMethod : size_t
  {
    LambdaZero = 0,  ///< Use only lambda = 0 in sampling.
    AllLambdas = 1   ///< Use all lambda values in sampling.
  };

  /**
   * \brief Enumeration of simulation stages.
   */
  enum class SimulationStage : size_t
  {
    Uninitialized = 0,   ///< Simulation not initialized.
    Initialization = 1,  ///< Initialization stage.
    Equilibration = 2,   ///< Equilibration stage.
    Production = 3       ///< Production stage.
  };

  /**
   * \brief Default constructor for the MonteCarlo class.
   *
   * Initializes a MonteCarlo object with default parameters.
   */
  MonteCarlo();

  /**
   * \brief Constructs a MonteCarlo object with simulation parameters from an InputReader.
   *
   * \param reader InputReader containing simulation parameters.
   */
  MonteCarlo(InputReader &reader) noexcept;

  /**
   * \brief Constructs a MonteCarlo object with specified simulation parameters.
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
  MonteCarlo(size_t numberOfCycles, size_t numberOfInitializationCycles, size_t numberOfEquilibrationCycles,
             size_t printEvery, size_t writeBinaryRestartEvery, size_t rescaleWangLandauEvery,
             size_t optimizeMCMovesEvery, std::vector<System> &systems, RandomNumber &randomSeed, size_t numberOfBlocks,
             bool outputToFiles = false);

  uint64_t versionNumber{1};  ///< Version number for serialization.

  bool outputToFiles{true};
  RandomNumber random;  ///< Random number generator.

  size_t numberOfCycles;                ///< Number of production cycles.
  size_t numberOfSteps;                 ///< Total number of steps performed.
  size_t numberOfInitializationCycles;  ///< Number of initialization cycles.
  size_t numberOfEquilibrationCycles;   ///< Number of equilibration cycles.
  size_t printEvery;                    ///< Frequency of printing status reports.
  size_t writeRestartEvery;             ///< Frequency of writing restart files.
  size_t writeBinaryRestartEvery;       ///< Frequency of writing binary restart files.
  size_t rescaleWangLandauEvery;        ///< Frequency of rescaling Wang-Landau factors.
  size_t optimizeMCMovesEvery;          ///< Frequency of optimizing MC moves.

  size_t currentCycle{0};                                           ///< Current cycle number.
  SimulationStage simulationStage{SimulationStage::Uninitialized};  ///< Current simulation stage.

  std::vector<System> systems;         ///< Vector of systems to simulate.
  size_t fractionalMoleculeSystem{0};  // the system where the fractional molecule is located

  std::vector<std::ofstream> streams;            ///< Output streams for writing data.
  std::vector<std::string> outputJsonFileNames;  ///< Filenames for output JSON files.
  std::vector<nlohmann::json> outputJsons;       ///< Output data in JSON format.

  BlockErrorEstimation estimation;  ///< Block error estimation object.

  std::chrono::duration<double> totalGridCreationTime{0};  ///< Total time for calculating the interpolation grid.
  std::chrono::duration<double> totalInitializationSimulationTime{0};  ///< Total time for initialization stage.
  std::chrono::duration<double> totalEquilibrationSimulationTime{0};   ///< Total time for equilibration stage.
  std::chrono::duration<double> totalProductionSimulationTime{0};      ///< Total time for production stage.
  std::chrono::duration<double> totalSimulationTime{0};                ///< Total simulation time.

  /**
   * \brief Creates output files for writing simulation data.
   */
  void createOutputFiles();

  /**
   * \brief Write the output header
   */
  void writeOutputHeader();

  /**
   * \brief Creates energy interpolation grids
   */
  void createInterpolationGrids();

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
   * \brief Returns a string representation of the MonteCarlo object.
   *
   * \return A string describing the MonteCarlo object.
   */
  std::string repr() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MonteCarlo &mc);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MonteCarlo &mc);
};
