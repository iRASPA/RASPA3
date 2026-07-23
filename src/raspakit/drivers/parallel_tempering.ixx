module;

export module parallel_tempering;

import std;

import randomnumbers;
import averages;
import system;
import input_reader;
import running_energy;
import json;

/**
 * \brief Multithreaded parallel-tempering driver for Monte Carlo.
 *
 * The single declared system is replicated into one replica per temperature of the ladder
 * ('ExternalTemperatures'); replica k runs at temperature T_k. Every replica runs in its own
 * thread with its own random-number stream (no OpenMP; plain std::jthread worker threads). The
 * threads only synchronize on a std::barrier every 'ParallelTemperingSwapEvery' cycles, where
 * configuration swaps between replicas at neighboring temperatures are attempted with the standard
 * parallel-tempering acceptance rule
 *
 *     acc = min(1, exp[(beta_B - beta_A) (U_B - U_A)])
 *
 * (extended with the Yan & de Pablo factor when the pressures differ). The replicas keep their
 * temperatures; only the configurations migrate through the ladder. Swaps alternate between the
 * (0,1),(2,3),... and (1,2),(3,4),... pairings so a configuration can traverse the whole ladder.
 *
 * Each replica writes its own output file with the standard status reports and final averages;
 * a combined output file holds the swap statistics.
 */
export struct ParallelTempering
{
  enum class SimulationStage : std::size_t
  {
    Uninitialized = 0,      ///< Simulation not initialized.
    PreInitialization = 1,  ///< Pre-initialization stage (translation/rotation only).
    Initialization = 2,     ///< Initialization stage.
    Equilibration = 3,      ///< Equilibration stage (Wang-Landau biasing for CFCMC moves).
    Production = 4          ///< Production stage.
  };

  ParallelTempering() = delete;
  ParallelTempering(const ParallelTempering&) = delete;
  ParallelTempering& operator=(const ParallelTempering&) = delete;

  /**
   * \brief Constructs the driver and replicates the single declared system into one replica per
   *        temperature of the ladder (replica k pinned at temperature k).
   */
  ParallelTempering(InputReader& reader);

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

  std::vector<double> temperatures;   ///< The temperature ladder (one replica per entry).
  std::size_t numberOfReplicas;       ///< Number of replicas == number of temperatures == number of threads.
  std::vector<System> systems;        ///< One replica per temperature.
  std::vector<RandomNumber> randoms;  ///< Independent random-number stream per replica.

  std::vector<std::size_t> stepsPerReplica;  ///< Production MC steps performed per replica.

  std::ofstream stream;            ///< The combined output stream (swap statistics, timings).
  std::string outputJsonFileName;  ///< Filename for the combined output JSON file.
  nlohmann::json outputJson;       ///< Combined output data in JSON format.

  std::mutex outputMutex;  ///< Guards the combined output stream for in-run progress lines.

  // Per-replica output: each worker thread writes exclusively to its own stream, so no locking is
  // needed for the periodic status reports.
  std::vector<std::ofstream> replicaStreams;      ///< One output stream per replica.
  std::vector<std::string> replicaJsonFileNames;  ///< Filename of the JSON output file per replica.
  std::vector<nlohmann::json> replicaJsons;       ///< JSON output data per replica.

  std::size_t swapSweeps{0};      ///< Number of swap sweeps performed (all stages).
  std::size_t sweepsThisStage{0};  ///< Number of swap sweeps in the current stage.
  std::size_t swapAttempts{0};    ///< Number of pairwise swap attempts.
  std::size_t swapAccepted{0};    ///< Number of accepted pairwise swaps.

  std::vector<std::size_t> swapAttemptsPerPair;  ///< Attempts per neighboring temperature-pair (k, k+1).
  std::vector<std::size_t> swapAcceptedPerPair;  ///< Acceptances per neighboring temperature-pair (k, k+1).

  std::chrono::duration<double> totalPreInitializationSimulationTime{0};  ///< Total time for pre-initialization.
  std::chrono::duration<double> totalInitializationSimulationTime{0};    ///< Total time for initialization stage.
  std::chrono::duration<double> totalEquilibrationSimulationTime{0};     ///< Total time for equilibration stage.
  std::chrono::duration<double> totalProductionSimulationTime{0};        ///< Total time for production stage.
  std::chrono::duration<double> totalSimulationTime{0};                  ///< Total simulation time.

  /**
   * \brief Runs the parallel-tempering simulation.
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
   * \brief One sweep of configuration-swap attempts between replicas at neighboring temperatures.
   *        Runs single-threaded inside the barrier completion (all worker threads are blocked).
   */
  void performSwapSweep(SimulationStage stage, std::size_t numberOfCycles) noexcept;

  /**
   * \brief Generates the final combined output: swap statistics and timings.
   */
  void output();

  /**
   * \brief Writes the final per-replica reports (energy drift, move statistics and averages) to
   *        the per-replica output files.
   */
  void writeReplicaFinalReports(std::vector<RunningEnergy>& recomputed);
};
