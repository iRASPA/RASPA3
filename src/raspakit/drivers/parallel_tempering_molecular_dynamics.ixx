module;

export module parallel_tempering_molecular_dynamics;

import std;

import randomnumbers;
import averages;
import system;
import input_reader;
import running_energy;
import integrators_cputime;
import archive;
import json;
import replica_round_trips;

/**
 * \brief Multithreaded replica-exchange (parallel-tempering) driver for molecular dynamics.
 *
 * The single declared system is replicated into one replica per temperature of the ladder
 * ('ExternalTemperatures'); replica k is integrated at temperature T_k. Every replica runs in its
 * own thread with its own random-number stream (plain std::jthread worker threads). The
 * pre-initialization and initialization stages are Monte Carlo (as in the plain MD driver); the
 * equilibration and production stages integrate the equations of motion, one time step per cycle.
 *
 * The threads only synchronize on a std::barrier every 'ParallelTemperingSwapEvery' cycles, where
 * configuration swaps between replicas at neighboring temperatures are attempted with the standard
 * parallel-tempering acceptance rule on the potential energies
 *
 *     acc = min(1, exp[(beta_B - beta_A) (U_B - U_A)])
 *
 * After an accepted swap the momenta that travelled with the configuration are rescaled by
 * sqrt(T_new / T_old) (Sugita & Okamoto, Chem. Phys. Lett. 314, 141-151, 1999), which makes the
 * kinetic contributions to the acceptance rule cancel; the thermostat chain is a heat-bath property
 * and stays with the replica. The replicas keep their temperatures; only the configurations migrate
 * through the ladder. Swaps alternate between the (0,1),(2,3),... and (1,2),(3,4),... pairings so a
 * configuration can traverse the whole ladder.
 *
 * Each replica writes its own output file with the standard status reports and final averages;
 * a combined output file holds the swap statistics.
 */
export struct ParallelTemperingMolecularDynamics
{
  enum class SimulationStage : std::size_t
  {
    Uninitialized = 0,      ///< Simulation not initialized.
    PreInitialization = 1,  ///< Pre-initialization stage (Monte Carlo: translation/rotation/reinsertion).
    Initialization = 2,     ///< Initialization stage (Monte Carlo).
    Equilibration = 3,      ///< Equilibration stage (molecular dynamics).
    Production = 4          ///< Production stage (molecular dynamics).
  };

  ParallelTemperingMolecularDynamics() = delete;
  ParallelTemperingMolecularDynamics(const ParallelTemperingMolecularDynamics&) = delete;
  ParallelTemperingMolecularDynamics& operator=(const ParallelTemperingMolecularDynamics&) = delete;

  /**
   * \brief Constructs the driver and replicates the single declared system into one replica per
   *        temperature of the ladder (replica k pinned at temperature k, thermostat included).
   */
  ParallelTemperingMolecularDynamics(InputReader& reader);

  std::uint64_t versionNumber{1};  ///< Version number for serialization.

  RandomNumber random;  ///< Random number generator (seeding + swap acceptance).

  std::size_t numberOfProductionCycles;         ///< Number of production cycles (MD steps).
  std::size_t numberOfPreInitializationCycles;  ///< Number of pre-initialization cycles (MC).
  std::size_t numberOfInitializationCycles;     ///< Number of initialization cycles (MC).
  std::size_t numberOfEquilibrationCycles;      ///< Number of equilibration cycles (MD steps).

  std::size_t printEvery;               ///< Frequency of printing status reports.
  std::size_t optimizeMCMovesEvery;     ///< Frequency of optimizing the MC moves (MC stages).
  std::size_t writeBinaryRestartEvery;  ///< Frequency of writing the binary restart file (0 disables).

  std::size_t numberOfBlocks;              ///< Number of blocks for the block-error estimation.
  std::size_t parallelTemperingSwapEvery;  ///< Attempt a swap sweep every this many cycles (0 disables).

  SimulationStage simulationStage{SimulationStage::Uninitialized};  ///< Current simulation stage.

  /// Completed cycles of the stage the binary restart file was written in; a restarted run
  /// continues that stage from this cycle.
  std::size_t cyclesCompletedThisStage{0};

  /// The end-of-cycle count communicated by the worker threads to the barrier completion
  /// (written by replica 0 before arriving; read while all threads are parked). Not serialized.
  std::atomic<std::size_t> checkpointCycle{0};

  std::vector<double> temperatures;   ///< The temperature ladder (one replica per entry).
  std::size_t numberOfReplicas;       ///< Number of replicas == number of temperatures == number of threads.
  std::vector<System> systems;        ///< One replica per temperature.
  std::vector<RandomNumber> randoms;  ///< Independent random-number stream per replica.

  std::vector<std::size_t> stepsPerReplica;  ///< Production MD steps performed per replica.

  /// Production-stage integrator timings per replica (gathered from the thread-local accumulators
  /// of the worker threads at the end of the stage).
  std::vector<IntegratorsCPUTime> integratorsCPUTimePerReplica;

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
  std::size_t swapAttempts{0};     ///< Number of pairwise swap attempts.
  std::size_t swapAccepted{0};     ///< Number of accepted pairwise swaps.

  std::vector<std::size_t> swapAttemptsPerPair;  ///< Attempts per neighboring temperature-pair (k, k+1).
  std::vector<std::size_t> swapAcceptedPerPair;  ///< Acceptances per neighboring temperature-pair (k, k+1).

  ReplicaRoundTrips roundTrips;  ///< Round-trip / up-fraction diagnostic of the configuration diffusion.

  std::chrono::duration<double> totalPreInitializationSimulationTime{0};  ///< Total time for pre-initialization.
  std::chrono::duration<double> totalInitializationSimulationTime{0};    ///< Total time for initialization stage.
  std::chrono::duration<double> totalEquilibrationSimulationTime{0};     ///< Total time for equilibration stage.
  std::chrono::duration<double> totalProductionSimulationTime{0};        ///< Total time for production stage.
  std::chrono::duration<double> totalSimulationTime{0};                  ///< Total simulation time.

  /**
   * \brief Runs the replica-exchange molecular-dynamics simulation.
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
   * \brief One Monte Carlo cycle of a single replica (pre-initialization and initialization stages).
   */
  void performReplicaMonteCarloCycle(std::size_t replicaId, SimulationStage stage);

  /**
   * \brief Per-replica preparation of a molecular-dynamics stage (velocities, thermostat, gradients,
   *        conserved-energy reference). Runs inside the replica's worker thread.
   */
  void prepareReplicaMolecularDynamicsStage(std::size_t replicaId, SimulationStage stage);

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
   * \brief Writes the final per-replica reports (energy drift, timings and averages) to the
   *        per-replica output files.
   */
  void writeReplicaFinalReports(std::vector<RunningEnergy>& recomputed);

  /**
   * \brief Writes the binary restart file. Called from the barrier completion, where all worker
   *        threads are parked, so the full driver state is consistent.
   */
  void writeBinaryRestartFile(std::size_t cyclesCompleted) noexcept;

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive,
                                            const ParallelTemperingMolecularDynamics& pt);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, ParallelTemperingMolecularDynamics& pt);
};
