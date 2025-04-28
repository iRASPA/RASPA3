module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <optional>
#include <vector>
#endif

export module molecular_dynamics;

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

/**
 * \brief Represents a molecular dynamics simulation.
 *
 * The MolecularDynamics struct manages the settings and execution of a molecular dynamics simulation,
 * including initialization, equilibration, and production stages.
 * It contains parameters controlling the simulation cycles, output, systems involved, and random number generation.
 */
export struct MolecularDynamics
{
  /**
   * \brief Enumeration for the weighting methods used in the simulation.
   */
  enum class WeightingMethod : size_t
  {
    LambdaZero = 0,  ///< Use weighting at lambda = 0.
    AllLambdas = 1   ///< Use weighting across all lambda values.
  };

  /**
   * \brief Enumeration representing the stages of the simulation.
   */
  enum class SimulationStage : size_t
  {
    Uninitialized = 0,   ///< The simulation has not been initialized.
    Initialization = 1,  ///< The simulation is in the initialization stage.
    Equilibration = 2,   ///< The simulation is in the equilibration stage.
    Production = 3       ///< The simulation is in the production stage.
  };

  /**
   * \brief Default constructor for MolecularDynamics.
   *
   * Initializes a MolecularDynamics object with default values.
   */
  MolecularDynamics();

  /**
   * \brief Constructs a MolecularDynamics object with specified input parameters.
   *
   * Initializes the MolecularDynamics simulation with parameters read from an InputReader.
   *
   * \param reader The InputReader object containing simulation parameters.
   */
  MolecularDynamics(InputReader &reader) noexcept;

  uint64_t versionNumber{1};  ///< Version number for serialization purposes.

  size_t numberOfCycles;                ///< Total number of production cycles.
  size_t numberOfSteps;                 ///< Total number of steps performed.
  size_t numberOfInitializationCycles;  ///< Number of initialization cycles.
  size_t numberOfEquilibrationCycles;   ///< Number of equilibration cycles.
  size_t printEvery;                    ///< Frequency of printing output.
  size_t writeBinaryRestartEvery;       ///< Frequency of writing binary restart files.
  size_t rescaleWangLandauEvery;        ///< Frequency of rescaling Wang-Landau factors.
  size_t optimizeMCMovesEvery;          ///< Frequency of optimizing Monte Carlo moves.

  size_t currentCycle{0};                                           ///< Current simulation cycle.
  SimulationStage simulationStage{SimulationStage::Uninitialized};  ///< Current stage of the simulation.

  std::vector<System> systems;         ///< Vector of systems in the simulation.
  RandomNumber random;                 ///< Random number generator.
  size_t fractionalMoleculeSystem{0};  ///< Index of the system where the fractional molecule is located.

  std::vector<std::ofstream> streams;  ///< Output file streams for each system.

  BlockErrorEstimation estimation;  ///< Object for block error estimation.

  std::chrono::duration<double> totalSimulationTime{0};  ///< Total simulation time.

  /**
   * \brief Creates output files for each system in the simulation.
   *
   * Initializes output streams for writing simulation data.
   */
  void createOutputFiles();


  void createInterpolationGrids();

  /**
   * \brief Runs the molecular dynamics simulation.
   *
   * Manages the simulation stages: initialization, equilibration, and production.
   */
  void run();

  /**
   * \brief Initializes the simulation.
   *
   * Performs the initialization stage, setting up initial conditions and variables.
   */
  void initialize();

  /**
   * \brief Equilibrates the simulation.
   *
   * Performs the equilibration stage, allowing the system to reach equilibrium.
   */
  void equilibrate();

  /**
   * \brief Runs the production stage of the simulation.
   *
   * Performs the main simulation cycles and collects data.
   */
  void production();

  /**
   * \brief Outputs the simulation results.
   *
   * Writes final statistics and results to the output files.
   */
  void output();

  /**
   * \brief Selects a random system from the simulation.
   *
   * \return A reference to a randomly selected System object.
   */
  System &randomSystem();

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MolecularDynamics &mc);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MolecularDynamics &mc);
};
