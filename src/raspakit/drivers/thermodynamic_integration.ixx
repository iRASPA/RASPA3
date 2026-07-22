module;

export module thermodynamic_integration;

import std;

import randomnumbers;
import averages;
import system;
import input_reader;
import json;

/**
 * \brief Fixed-lambda thermodynamic-integration driver.
 *
 * Runs a Monte Carlo simulation in which the fractional molecule of each component with a
 * 'LambdaBinIndex' is pinned at a constant lambda (lambda = binIndex * delta). No lambda-changing
 * CFCMC moves and no Wang-Landau biasing are involved; the configuration space is sampled with the
 * regular MC moves (translation, rotation, reinsertion, volume, ...) and during production the
 * ensemble average <dU/dlambda> is accumulated in the single fixed lambda-bin.
 *
 * The result is one point of the <dU/dlambda>(lambda) curve. Running one simulation per lambda-bin
 * and integrating <dU/dlambda> from lambda=0 to lambda=1 yields the excess chemical potential.
 */
export struct ThermodynamicIntegration
{
  enum class SimulationStage : std::size_t
  {
    Uninitialized = 0,   ///< Simulation not initialized.
    Initialization = 1,  ///< Initialization stage.
    Equilibration = 2,   ///< Equilibration stage.
    Production = 3       ///< Production stage.
  };

  ThermodynamicIntegration() = delete;
  ThermodynamicIntegration(const ThermodynamicIntegration&) = delete;
  ThermodynamicIntegration& operator=(const ThermodynamicIntegration&) = delete;

  /**
   * \brief Constructs a ThermodynamicIntegration object with simulation parameters from an InputReader.
   *
   * \param reader InputReader containing simulation parameters.
   */
  ThermodynamicIntegration(InputReader& reader) noexcept;

  RandomNumber random;  ///< Random number generator.

  std::size_t numberOfProductionCycles;      ///< Number of production cycles.
  std::size_t numberOfInitializationCycles;  ///< Number of initialization cycles.
  std::size_t numberOfEquilibrationCycles;   ///< Number of equilibration cycles.
  std::size_t numberOfSteps{0};              ///< Total number of production steps performed.

  std::size_t printEvery;               ///< Frequency of printing status reports.
  std::size_t writeRestartEvery{5000};  ///< Frequency of writing restart files.
  std::size_t optimizeMCMovesEvery;     ///< Frequency of optimizing MC moves.

  std::size_t currentCycle{0};                                      ///< Current cycle number.
  SimulationStage simulationStage{SimulationStage::Uninitialized};  ///< Current simulation stage.

  std::vector<System> systems;              ///< Vector of systems to simulate.
  std::size_t fractionalMoleculeSystem{0};  ///< The system where the fractional molecule is located.

  std::vector<std::ofstream> streams;            ///< Output streams for writing data.
  std::vector<std::string> outputJsonFileNames;  ///< Filenames for output JSON files.
  std::vector<nlohmann::json> outputJsons;       ///< Output data in JSON format.

  BlockErrorEstimation estimation{};  ///< Block error estimation object.

  std::chrono::duration<double> totalInitializationSimulationTime{0};  ///< Total time for initialization stage.
  std::chrono::duration<double> totalEquilibrationSimulationTime{0};   ///< Total time for equilibration stage.
  std::chrono::duration<double> totalProductionSimulationTime{0};      ///< Total time for production stage.
  std::chrono::duration<double> totalSimulationTime{0};                ///< Total simulation time.

  /**
   * \brief Runs the fixed-lambda thermodynamic-integration simulation.
   *
   * Orchestrates the simulation by executing initialization, equilibration,
   * and production stages sequentially, followed by the final output.
   */
  void run();

  /**
   * \brief Sets up the systems: fixed-lambda fractional molecules, output files, interpolation grids.
   */
  void setup();

  /**
   * \brief Performs a single Monte Carlo cycle (no lambda-changing moves, no Wang-Landau).
   */
  void performCycle();

  void initialize();
  void equilibrate();
  void production();

  /**
   * \brief Generates the final output: statistics, timings and the fixed-lambda <dU/dlambda> point.
   */
  void output();

  /**
   * \brief Writes the fixed-lambda <dU/dlambda> results of a system as human-readable text.
   */
  std::string writeThermodynamicIntegrationPoint(const System& system) const;

  /**
   * \brief The fixed-lambda <dU/dlambda> results of a system in JSON format.
   */
  nlohmann::json jsonThermodynamicIntegrationPoint(const System& system) const;
};
