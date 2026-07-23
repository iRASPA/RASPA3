module;

export module hyper_parallel_tempering;

import std;

import randomnumbers;
import averages;
import system;
import input_reader;
import running_energy;
import json;

/**
 * \brief Multithreaded hyper-parallel-tempering driver for Monte Carlo.
 *
 * Ref: "Hyper-parallel tempering Monte Carlo: Application to the Lennard-Jones fluid and the
 * restricted primitive model", G. Yan and J.J. de Pablo, JCP, 111(21): 9509-9516, 1999.
 *
 * The single declared system is replicated into one replica per point of the (temperature,
 * pressure) grid spanned by the ladders 'ExternalTemperatures' x 'ExternalPressures'. The
 * pressures are converted to per-component fugacities internally through the Peng-Robinson
 * equation of state evaluated at each grid point. Every replica runs in its own thread with its
 * own random-number stream (no OpenMP; plain std::jthread worker threads). The threads only
 * synchronize on a std::barrier every 'ParallelTemperingSwapEvery' cycles, where configuration
 * swaps between replicas at neighboring grid points are attempted with the Yan & de Pablo
 * acceptance rule
 *
 *     acc = min(1, exp[(beta_B - beta_A)(U_B - U_A)]
 *                  x prod_i [(beta_A f_A,i)/(beta_B f_B,i)]^(N_B,i - N_A,i))
 *
 * The sweeps alternate between the temperature direction (neighboring temperatures at the same
 * pressure) and the pressure direction (neighboring pressures at the same temperature), each with
 * an alternating pairing offset, so a configuration can traverse the whole grid. The replicas keep
 * their (temperature, pressure) state points; only the configurations migrate through the grid.
 *
 * Each replica writes its own output file with the standard status reports and final averages;
 * a combined output file holds the swap statistics for both directions.
 */
export struct HyperParallelTempering
{
  enum class SimulationStage : std::size_t
  {
    Uninitialized = 0,      ///< Simulation not initialized.
    PreInitialization = 1,  ///< Pre-initialization stage (translation/rotation only).
    Initialization = 2,     ///< Initialization stage.
    Equilibration = 3,      ///< Equilibration stage (Wang-Landau biasing for CFCMC moves).
    Production = 4          ///< Production stage.
  };

  HyperParallelTempering() = delete;
  HyperParallelTempering(const HyperParallelTempering&) = delete;
  HyperParallelTempering& operator=(const HyperParallelTempering&) = delete;

  /**
   * \brief Constructs the driver and replicates the single declared system into one replica per
   *        (temperature, pressure) grid point, recomputing the per-replica fugacity coefficients
   *        with the Peng-Robinson equation of state.
   */
  HyperParallelTempering(InputReader& reader);

  RandomNumber random;  ///< Random number generator (seeding + swap acceptance).

  std::size_t numberOfProductionCycles;         ///< Number of production cycles.
  std::size_t numberOfPreInitializationCycles;  ///< Number of pre-initialization cycles.
  std::size_t numberOfInitializationCycles;     ///< Number of initialization cycles.
  std::size_t numberOfEquilibrationCycles;      ///< Number of equilibration cycles.

  std::size_t printEvery;              ///< Frequency of printing status reports.
  std::size_t optimizeMCMovesEvery;    ///< Frequency of optimizing MC moves.
  std::size_t rescaleWangLandauEvery;  ///< Frequency of adjusting the Wang-Landau biasing factors.

  std::size_t numberOfBlocks;              ///< Number of blocks for the block-error estimation.
  std::size_t parallelTemperingSwapEvery;  ///< Attempt a swap sweep every this many cycles (0 disables).

  SimulationStage simulationStage{SimulationStage::Uninitialized};  ///< Current simulation stage.

  std::vector<double> temperatures;  ///< The temperature ladder [K].
  std::vector<double> pressures;     ///< The pressure ladder [Pa].
  std::size_t numberOfTemperatures;  ///< Number of temperatures in the ladder.
  std::size_t numberOfPressures;     ///< Number of pressures in the ladder.
  std::size_t numberOfReplicas;      ///< numberOfTemperatures x numberOfPressures == number of threads.
  std::vector<System> systems;       ///< One replica per (temperature, pressure) grid point.
  std::vector<RandomNumber> randoms; ///< Independent random-number stream per replica.

  std::vector<std::size_t> stepsPerReplica;  ///< Production MC steps performed per replica.

  /// Cycles completed in the previous stages; the time-evolution properties (number of molecules,
  /// volume) are indexed by the absolute cycle number counted over all stages.
  std::size_t absoluteCycleOffset{0};

  std::ofstream stream;            ///< The combined output stream (swap statistics, timings).
  std::string outputJsonFileName;  ///< Filename for the combined output JSON file.
  nlohmann::json outputJson;       ///< Combined output data in JSON format.

  std::mutex outputMutex;  ///< Guards the combined output stream for in-run progress lines.

  // Per-replica output: each worker thread writes exclusively to its own stream, so no locking is
  // needed for the periodic status reports.
  std::vector<std::ofstream> replicaStreams;      ///< One output stream per replica.
  std::vector<std::string> replicaJsonFileNames;  ///< Filename of the JSON output file per replica.
  std::vector<nlohmann::json> replicaJsons;       ///< JSON output data per replica.

  std::size_t swapSweeps{0};       ///< Number of swap sweeps performed (all stages).
  std::size_t sweepsThisStage{0};  ///< Number of swap sweeps in the current stage.
  std::size_t swapAttempts{0};     ///< Number of pairwise swap attempts (both directions).
  std::size_t swapAccepted{0};     ///< Number of accepted pairwise swaps (both directions).

  /// Attempts/acceptances per neighboring temperature-pair (t, t+1), aggregated over the pressures.
  std::vector<std::size_t> swapAttemptsPerTemperaturePair;
  std::vector<std::size_t> swapAcceptedPerTemperaturePair;
  /// Attempts/acceptances per neighboring pressure-pair (p, p+1), aggregated over the temperatures.
  std::vector<std::size_t> swapAttemptsPerPressurePair;
  std::vector<std::size_t> swapAcceptedPerPressurePair;

  std::chrono::duration<double> totalPreInitializationSimulationTime{0};  ///< Total time for pre-initialization.
  std::chrono::duration<double> totalInitializationSimulationTime{0};    ///< Total time for initialization stage.
  std::chrono::duration<double> totalEquilibrationSimulationTime{0};     ///< Total time for equilibration stage.
  std::chrono::duration<double> totalProductionSimulationTime{0};        ///< Total time for production stage.
  std::chrono::duration<double> totalSimulationTime{0};                  ///< Total simulation time.

  /**
   * \brief The replica at grid point (temperature index, pressure index).
   */
  std::size_t replicaIndex(std::size_t temperatureIndex, std::size_t pressureIndex) const
  {
    return temperatureIndex * numberOfPressures + pressureIndex;
  }

  /**
   * \brief Runs the hyper-parallel-tempering simulation.
   */
  void run();

  /**
   * \brief Sets up the replicas, output files and interpolation grids.
   */
  void setup();

  /**
   * \brief Runs one simulation stage: every replica in its own thread, synchronized on a barrier
   *        every 'parallelTemperingSwapEvery' cycles for the swap sweeps.
   */
  void runStage(SimulationStage stage, std::size_t numberOfCycles);

  /**
   * \brief One Monte Carlo cycle of a single replica.
   */
  void performReplicaCycle(std::size_t replicaId, SimulationStage stage, std::size_t currentBlock);

  /**
   * \brief One sweep of configuration-swap attempts, alternating between the temperature and the
   *        pressure direction of the grid. Runs single-threaded inside the barrier completion
   *        (all worker threads are blocked).
   */
  void performSwapSweep(SimulationStage stage, std::size_t numberOfCycles) noexcept;

  /**
   * \brief Generates the final combined output: swap statistics and timings.
   */
  void output();

  /**
   * \brief Assembles the adsorption-isotherm data per temperature and writes one gnuplot-friendly
   *        file per temperature (rows ordered by increasing pressure; columns: fugacity, loading
   *        with confidence-interval error in molecules/cell, molecules/unit-cell, mol/kg-framework
   *        and mg/g-framework; one block per component). Overwritten on every write: called during
   *        production from the barrier completion (all worker threads are parked, so reading the
   *        replica averages is race-free) and once more at the end of the run.
   */
  void writeIsothermSnapshot() const;

  /**
   * \brief Writes the final per-replica reports (energy drift, move statistics and averages) to
   *        the per-replica output files.
   */
  void writeReplicaFinalReports(std::vector<RunningEnergy>& recomputed);
};
