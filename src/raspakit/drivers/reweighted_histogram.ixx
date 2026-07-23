module;

export module reweighted_histogram;

import std;

import randomnumbers;
import averages;
import system;
import input_reader;
import running_energy;
import json;

/**
 * \brief Multithreaded replica-grid driver with multiple-histogram reweighting.
 *
 * Refs: A.M. Ferrenberg and R.H. Swendsen, "Optimized Monte Carlo data analysis",
 * PRL 63(12), 1195-1198, 1989; S. Kumar et al., "The weighted histogram analysis method
 * for free-energy calculations on biomolecules", J. Comput. Chem. 13(8), 1011-1021, 1992.
 *
 * The simulation part is a hyper-parallel-tempering run: the single declared system is
 * replicated into one replica per point of the (temperature, pressure) grid spanned by the
 * ladders 'ExternalTemperatures' x 'ExternalPressures', every replica runs in its own thread,
 * and configuration swaps between neighboring grid points are attempted with the Yan & de Pablo
 * acceptance rule every 'ParallelTemperingSwapEvery' cycles.
 *
 * On top of that, every 'SampleReweightingEvery' production cycles each replica records a raw
 * (N, U) sample: the molecule count of the single adsorbate component and the potential energy.
 * At the end of the run the samples of all replicas are combined by solving the WHAM
 * self-consistent equations in log space for the grand-canonical density of states Omega(N, U),
 *
 *     ln Omega(N,U) = ln M(N,U) - ln sum_i n_i exp(g_i - beta_i U + N ln(beta_i f_i))
 *     g_i           = -ln sum_{N,U} Omega(N,U) exp(-beta_i U + N ln(beta_i f_i))
 *
 * (M is the total count in a bin, n_i the number of samples of state i, f_i the fugacity).
 * The density of states is then reweighted to arbitrary (temperature, fugacity) points, giving
 * continuous adsorption isotherms in between and beyond the simulated grid points. The error
 * bars are estimated by re-solving the WHAM equations per block. The single-simulation grid
 * points must overlap in (N, U) space for the reweighting to be reliable; the effective sample
 * size written with every reweighted point diagnoses the overlap.
 *
 * For bulk boxes (no framework) the analysis additionally locates the vapor-liquid coexistence
 * at every requested subcritical temperature with the equal-weight criterion (Wilding): the
 * fugacity is bisected until the two peaks of the bimodal reweighted P(N) carry equal
 * probability weight, yielding the coexistence fugacity, the saturated vapor and liquid
 * densities, and the saturation pressure (from beta p V = ln Xi, normalized by the empty-box
 * state).
 *
 * Because a single sampled U(x) is reused at all temperatures, the force field must not depend
 * on the temperature (temperature-dependent potentials such as Feynman-Hibbs are rejected), and
 * a single adsorbate component is required (the histograms are collected over its molecule count).
 */
export struct ReweightedHistogram
{
  enum class SimulationStage : std::size_t
  {
    Uninitialized = 0,      ///< Simulation not initialized.
    PreInitialization = 1,  ///< Pre-initialization stage (translation/rotation only).
    Initialization = 2,     ///< Initialization stage.
    Equilibration = 3,      ///< Equilibration stage (Wang-Landau biasing for CFCMC moves).
    Production = 4          ///< Production stage.
  };

  /// One raw reweighting sample of a replica: the potential energy (internal units) and the
  /// molecule count of the adsorbate component, tagged with the block it was collected in.
  struct Sample
  {
    double energy;                    ///< Potential energy U [internal energy units].
    std::uint32_t numberOfMolecules;  ///< Molecule count N of the single adsorbate component.
    std::uint32_t block;              ///< Block index (for the per-block error estimation).
  };

  ReweightedHistogram() = delete;
  ReweightedHistogram(const ReweightedHistogram&) = delete;
  ReweightedHistogram& operator=(const ReweightedHistogram&) = delete;

  /**
   * \brief Constructs the driver and replicates the single declared system into one replica per
   *        (temperature, pressure) grid point, recomputing the per-replica fugacity coefficients
   *        with the Peng-Robinson equation of state. Throws when the force field is
   *        temperature-dependent (the reweighting reuses one U(x) at all temperatures).
   */
  ReweightedHistogram(InputReader& reader);

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

  std::size_t sampleReweightingEvery;  ///< Record an (N, U) sample every this many production cycles.
  std::vector<double> reweightingTemperatures;  ///< Temperatures the reweighted isotherms are written at [K].
  std::pair<double, double> reweightingPressureRange;  ///< Pressure range of the reweighted isotherms [Pa].
  std::size_t reweightingNumberOfPressures;  ///< Number of log-spaced pressures of the reweighted isotherms.

  SimulationStage simulationStage{SimulationStage::Uninitialized};  ///< Current simulation stage.

  std::vector<double> temperatures;  ///< The temperature ladder [K].
  std::vector<double> pressures;     ///< The pressure ladder [Pa].
  std::size_t numberOfTemperatures;  ///< Number of temperatures in the ladder.
  std::size_t numberOfPressures;     ///< Number of pressures in the ladder.
  std::size_t numberOfReplicas;      ///< numberOfTemperatures x numberOfPressures == number of threads.
  std::vector<System> systems;       ///< One replica per (temperature, pressure) grid point.
  std::vector<RandomNumber> randoms; ///< Independent random-number stream per replica.

  /// Raw (N, U) samples per replica; each worker thread appends exclusively to its own vector.
  std::vector<std::vector<Sample>> reweightingSamples;

  std::vector<std::size_t> stepsPerReplica;  ///< Production MC steps performed per replica.

  /// Cycles completed in the previous stages; the time-evolution properties (number of molecules,
  /// volume) are indexed by the absolute cycle number counted over all stages.
  std::size_t absoluteCycleOffset{0};

  std::ofstream stream;            ///< The combined output stream (swap statistics, WHAM analysis, timings).
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
  std::chrono::duration<double> totalReweightingAnalysisTime{0};         ///< Total time for the WHAM analysis.
  std::chrono::duration<double> totalSimulationTime{0};                  ///< Total simulation time.

  /**
   * \brief The replica at grid point (temperature index, pressure index).
   */
  std::size_t replicaIndex(std::size_t temperatureIndex, std::size_t pressureIndex) const
  {
    return temperatureIndex * numberOfPressures + pressureIndex;
  }

  /**
   * \brief Runs the simulation and the reweighting analysis.
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
   * \brief Generates the final combined output: swap statistics, WHAM analysis and timings.
   */
  void output();

  /**
   * \brief Assembles the directly-measured adsorption-isotherm data per temperature and writes
   *        one gnuplot-friendly file per temperature (identical to the hyper-parallel-tempering
   *        isotherm files; the reweighted isotherms are written separately by the analysis).
   */
  void writeIsothermSnapshot() const;

  /**
   * \brief Solves the WHAM self-consistent equations over the pooled (N, U) samples of all
   *        replicas and writes the reweighted isotherms (one file per requested temperature,
   *        evaluated on a fine log-spaced pressure grid), the per-state free energies with a
   *        self-consistency check against the directly-measured loadings, and the convergence
   *        report to the combined output.
   */
  void performReweightingAnalysis();

  /**
   * \brief Writes the final per-replica reports (energy drift, move statistics and averages) to
   *        the per-replica output files.
   */
  void writeReplicaFinalReports(std::vector<RunningEnergy>& recomputed);
};
