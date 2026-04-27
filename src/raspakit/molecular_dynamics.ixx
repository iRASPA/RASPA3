module;

export module molecular_dynamics;

import std;

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
  enum class WeightingMethod : std::size_t
  {
    LambdaZero = 0,  ///< Use weighting at lambda = 0.
    AllLambdas = 1   ///< Use weighting across all lambda values.
  };

  /**
   * \brief Enumeration representing the stages of the simulation.
   */
  enum class SimulationStage : std::size_t
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

  /**
   * \brief Constructs a MolecularDynamicso object with specified simulation parameters.
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
  MolecularDynamics(std::size_t numberOfCycles, std::size_t numberOfInitializationCycles,
                    std::size_t numberOfEquilibrationCycles, std::size_t printEvery, std::size_t writeBinaryRestartEvery,
                    std::size_t rescaleWangLandauEvery, std::size_t optimizeMCMovesEvery, std::vector<System> &systems,
                    std::optional<std::size_t> randomSeed, std::size_t numberOfBlocks, bool outputToFiles = false);

  std::uint64_t versionNumber{1};  ///< Version number for serialization purposes.
  
  bool outputToFiles{true};
  RandomNumber random;             ///< Random number generator.

  std::size_t numberOfCycles;                ///< Total number of production cycles.
  std::size_t numberOfSteps;                 ///< Total number of steps performed.
  std::size_t numberOfInitializationCycles;  ///< Number of initialization cycles.
  std::size_t numberOfEquilibrationCycles;   ///< Number of equilibration cycles.
  
  std::size_t printEvery;                    ///< Frequency of printing output.
  std::size_t writeRestartEvery;             ///< Frequency of writing restart files. 
  std::size_t writeBinaryRestartEvery;       ///< Frequency of writing binary restart files.
  std::size_t rescaleWangLandauEvery;        ///< Frequency of rescaling Wang-Landau factors.
  std::size_t optimizeMCMovesEvery;          ///< Frequency of optimizing Monte Carlo moves.

  std::size_t currentCycle{};                                       ///< Current simulation cycle.
  std::size_t absoluteCurrentCycle{};
  SimulationStage simulationStage{SimulationStage::Uninitialized};  ///< Current stage of the simulation.

  std::vector<System> systems;              ///< Vector of systems in the simulation.
  std::size_t fractionalMoleculeSystem{0};  ///< Index of the system where the fractional molecule is located.

  std::vector<std::ofstream> streams;  ///< Output file streams for each system.
  std::vector<std::string> outputJsonFileNames;  ///< Filenames for output JSON files.
  std::vector<nlohmann::json> outputJsons;       ///< Output data in JSON format.
                                                 ///
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

  void setup();

  void tearDown();

  /**
   * \brief Initializes the simulation.
   *
   * Performs the initialization stage, setting up initial conditions and variables.
   */
  void initialize(std::function<void()> call_back_function = []{}, std::size_t callBackEvery = 100);

  /**
   * \brief Equilibrates the simulation.
   *
   * Performs the equilibration stage, allowing the system to reach equilibrium.
   */
  void equilibrate(std::function<void()> call_back_function = []{}, std::size_t callBackEvery = 100);

  /**
   * \brief Runs the production stage of the simulation.
   *
   * Performs the main simulation cycles and collects data.
   */
  void production(std::function<void()> call_back_function = []{}, std::size_t callBackEvery = 100);

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
