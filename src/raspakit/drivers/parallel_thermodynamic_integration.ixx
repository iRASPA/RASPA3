module;

export module parallel_thermodynamic_integration;

import std;

import randomnumbers;
import averages;
import system;
import input_reader;
import property_lambda_probability_histogram;
import running_energy;
import json;

/**
 * \brief Parallel thermodynamic-integration driver with lambda-exchange.
 *
 * The single declared system is replicated into one replica per lambda-bin; replica k starts with
 * its fractional molecule(s) pinned at lambda-bin k. Every replica runs in its own thread with its
 * own random-number stream (no OpenMP; plain std::jthread worker threads). The threads only
 * synchronize on a std::barrier every 'LambdaExchangeEvery' cycles, where Hamiltonian
 * parallel-tempering exchanges of the lambda values between replicas at neighboring lambda-bins
 * are attempted (and at the start/end of each simulation stage).
 *
 * Because exchanges permute the lambda values over the replicas, the set of occupied bins is
 * always exactly {0, ..., numberOfLambdaBins-1}: every lambda-bin is sampled during every cycle.
 * During production each replica accumulates <dU/dlambda> in whatever bin it currently occupies.
 *
 * At the end the per-bin dU/dlambda book-keeping of all replicas is stitched together into a
 * single histogram; the excess chemical potential is obtained by integrating a natural cubic
 * spline through the <dU/dlambda>(lambda) curve from lambda=0 to lambda=1 (a Simpson estimate is
 * reported alongside), with the error estimated from the block-wise spline integrals.
 */
export struct ParallelThermodynamicIntegration
{
  enum class SimulationStage : std::size_t
  {
    Uninitialized = 0,   ///< Simulation not initialized.
    Initialization = 1,  ///< Initialization stage.
    Equilibration = 2,   ///< Equilibration stage.
    Production = 3       ///< Production stage.
  };

  ParallelThermodynamicIntegration() = delete;
  ParallelThermodynamicIntegration(const ParallelThermodynamicIntegration&) = delete;
  ParallelThermodynamicIntegration& operator=(const ParallelThermodynamicIntegration&) = delete;

  /**
   * \brief Constructs the driver and replicates the single declared system into one replica per
   *        lambda-bin (replica k pinned at bin k).
   */
  ParallelThermodynamicIntegration(InputReader& reader);

  RandomNumber random;  ///< Random number generator (seeding + lambda-exchange acceptance).

  std::size_t numberOfProductionCycles;      ///< Number of production cycles.
  std::size_t numberOfInitializationCycles;  ///< Number of initialization cycles.
  std::size_t numberOfEquilibrationCycles;   ///< Number of equilibration cycles.

  std::size_t printEvery;            ///< Frequency of printing status reports.
  std::size_t optimizeMCMovesEvery;  ///< Frequency of optimizing MC moves.

  std::size_t numberOfBlocks;       ///< Number of blocks for the block-error estimation.
  std::size_t numberOfLambdaBins;   ///< Number of lambda-bins == number of replicas == number of threads.
  std::size_t lambdaExchangeEvery;  ///< Attempt a lambda-exchange sweep every this many cycles (0 disables).

  SimulationStage simulationStage{SimulationStage::Uninitialized};  ///< Current simulation stage.

  std::vector<System> systems;        ///< One replica per lambda-bin.
  std::vector<RandomNumber> randoms;  ///< Independent random-number stream per replica.
  std::size_t tiComponentId{0};       ///< The component whose lambda is pinned.

  std::vector<std::size_t> stepsPerReplica;  ///< Production MC steps performed per replica.

  std::ofstream stream;             ///< The combined output stream (stitched results, exchange statistics).
  std::string outputJsonFileName;   ///< Filename for the combined output JSON file.
  std::string dudlambdaFileName;    ///< Filename for the periodic <dU/dlambda> snapshot file.
  nlohmann::json outputJson;        ///< Combined output data in JSON format.

  // Per-replica output: each worker thread writes exclusively to its own stream, so no locking is
  // needed for the periodic status reports.
  std::vector<std::ofstream> replicaStreams;        ///< One output stream per replica.
  std::vector<std::string> replicaJsonFileNames;    ///< Filename of the JSON output file per replica.
  std::vector<nlohmann::json> replicaJsons;         ///< JSON output data per replica.

  std::mutex outputMutex;  ///< Guards the output stream for in-run progress lines.

  std::size_t exchangeSweeps{0};      ///< Number of lambda-exchange sweeps performed (all stages).
  std::size_t sweepsThisStage{0};     ///< Number of lambda-exchange sweeps in the current stage.
  std::size_t exchangeAttempts{0};  ///< Number of pairwise lambda-exchange attempts.
  std::size_t exchangeAccepted{0};  ///< Number of accepted pairwise lambda-exchanges.

  std::vector<std::size_t> exchangeAttemptsPerPair;  ///< Attempts per neighboring bin-pair (bin, bin+1).
  std::vector<std::size_t> exchangeAcceptedPerPair;  ///< Acceptances per neighboring bin-pair (bin, bin+1).

  std::chrono::duration<double> totalInitializationSimulationTime{0};  ///< Total time for initialization stage.
  std::chrono::duration<double> totalEquilibrationSimulationTime{0};   ///< Total time for equilibration stage.
  std::chrono::duration<double> totalProductionSimulationTime{0};      ///< Total time for production stage.
  std::chrono::duration<double> totalSimulationTime{0};                ///< Total simulation time.

  /**
   * \brief Runs the parallel thermodynamic-integration simulation.
   */
  void run();

  /**
   * \brief Sets up the replicas, output file and interpolation grids.
   */
  void setup();

  /**
   * \brief Runs one simulation stage: every replica in its own thread, synchronized on a barrier
   *        every 'lambdaExchangeEvery' cycles for the lambda-exchange sweeps.
   */
  void runStage(SimulationStage stage, std::size_t numberOfCycles);

  /**
   * \brief One Monte Carlo cycle of a single replica (no lambda-changing moves, no Wang-Landau).
   */
  void performReplicaCycle(std::size_t replicaId, SimulationStage stage, std::size_t currentBlock);

  /**
   * \brief One sweep of lambda-exchange attempts between replicas at neighboring lambda-bins.
   *        Runs single-threaded inside the barrier completion (all worker threads are blocked).
   */
  void performLambdaExchangeSweep(SimulationStage stage, std::size_t numberOfCycles) noexcept;

  /**
   * \brief Writes the current stitched <dU/dlambda>(lambda) curve and its running spline integral
   *        to a gnuplot-friendly snapshot file (overwritten on every write). Called during
   *        production from the barrier completion, where all worker threads are parked.
   */
  void writeDUdlambdaSnapshot(std::size_t cycle, std::size_t numberOfCycles) const;

  /**
   * \brief Generates the final output: statistics, timings and the stitched <dU/dlambda> curve.
   */
  void output();

  /**
   * \brief Writes the final per-replica reports (energy drift, move statistics, averages and the
   *        replica's own per-bin dU/dlambda book-keeping) to the per-replica output files.
   */
  void writeReplicaFinalReports(std::vector<RunningEnergy>& recomputed);

  /**
   * \brief The per-bin dU/dlambda book-keeping of a single replica as human-readable text: which
   *        lambda-bins the replica visited and the samples it collected there.
   */
  std::string writeReplicaThermodynamicIntegration(std::size_t replicaId) const;

  /**
   * \brief Sums the per-bin dU/dlambda book-keeping of all replicas into one histogram.
   */
  PropertyLambdaProbabilityHistogram stitchedHistogram() const;

  /**
   * \brief The stitched <dU/dlambda>(lambda) curve and its spline-integrated excess chemical
   *        potential as human-readable text.
   */
  std::string writeStitchedThermodynamicIntegration() const;

  /**
   * \brief The stitched thermodynamic-integration results in JSON format.
   */
  nlohmann::json jsonStitchedThermodynamicIntegration() const;

  /**
   * \brief Integral over [0,1] of the natural cubic spline through equidistant data points.
   */
  static double cubicSplineIntegral(const std::vector<double>& data, double delta);
};
