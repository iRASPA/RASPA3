module;

export module parallel_tmmc;

import std;

import randomnumbers;
import averages;
import system;
import input_reader;
import running_energy;
import double3;
import json;

/**
 * \brief Multithreaded transition-matrix Monte Carlo (TMMC) with windowed macrostate walkers.
 *
 * Refs: J.R. Errington, "Direct calculation of liquid-vapor phase equilibria from transition
 * matrix Monte Carlo", J. Chem. Phys. 118(22), 9915-9925, 2003; V.K. Shen and J.R. Errington,
 * "Determination of fluid-phase behavior using transition-matrix Monte Carlo", J. Chem. Phys.
 * 122, 064508, 2005.
 *
 * The macrostate range [MacroStateMinimumNumberOfMolecules, MacroStateMaximumNumberOfMolecules]
 * of the single adsorbate component is split into 'NumberOfWindows' contiguous windows that share
 * their endpoint macrostates. The single declared system is replicated into one walker per
 * (temperature, window) pair of the ladder 'ExternalTemperatures' x windows (a single
 * 'ExternalTemperature' gives a one-temperature run); every walker runs grand-canonical Monte
 * Carlo in its own thread, confined to its window and flattened by the transition-matrix bias
 * (updated every 'TMMCUpdateEvery' steps). The collection matrix records the unbiased acceptance
 * probabilities of all attempted insertions/deletions - also those rejected at the window bounds -
 * so the collection matrices of the windows of one temperature simply add up, and the macrostate
 * probability distribution ln Pi(N) over the full range follows from detailed balance:
 *
 *     ln Pi(N+1) = ln Pi(N) + ln P(N -> N+1) - ln P(N+1 -> N)
 *
 * ln Pi(N) is exact at the reference fugacity f_ref (from 'ExternalPressure' through the
 * Peng-Robinson equation of state) and reweights exactly to any other fugacity,
 * ln Pi(N; f) = ln Pi(N; f_ref) + N ln(f / f_ref), giving continuous adsorption isotherms
 * <N>(f) per temperature. The error bars come from re-deriving ln Pi from the per-block
 * increments of the collection matrices.
 *
 * For bulk boxes (no framework) the analysis additionally locates the vapor-liquid coexistence
 * at every simulated subcritical temperature with the equal-weight criterion (Wilding): the
 * fugacity is bisected until the two peaks of the bimodal Pi(N) carry equal probability weight,
 * yielding the coexistence fugacity, the saturated vapor and liquid densities, and the
 * saturation pressure from beta p V = ln Xi (exactly normalized by the empty-box state N = 0
 * when the macrostate range starts at zero).
 *
 * Every temperature is solved from its own walkers only, so temperature-dependent potentials are
 * allowed; the ladder gives one isotherm and one coexistence point per simulated temperature.
 */
export struct ParallelTMMC
{
  enum class SimulationStage : std::size_t
  {
    Uninitialized = 0,      ///< Simulation not initialized.
    PreInitialization = 1,  ///< Pre-initialization stage (translation/rotation only).
    Initialization = 2,     ///< Initialization stage.
    Equilibration = 3,      ///< Equilibration stage (transition-matrix bias develops).
    Production = 4          ///< Production stage (collection matrices accumulate per block).
  };

  ParallelTMMC() = delete;
  ParallelTMMC(const ParallelTMMC&) = delete;
  ParallelTMMC& operator=(const ParallelTMMC&) = delete;

  /**
   * \brief Constructs the driver and replicates the single declared system into one walker per
   *        (temperature, window) pair, pinning the temperatures and window bounds and recomputing
   *        the per-temperature fugacity coefficients with the Peng-Robinson equation of state.
   */
  ParallelTMMC(InputReader& reader);

  RandomNumber random;  ///< Random number generator (seeding).

  std::size_t numberOfProductionCycles;         ///< Number of production cycles.
  std::size_t numberOfPreInitializationCycles;  ///< Number of pre-initialization cycles.
  std::size_t numberOfInitializationCycles;     ///< Number of initialization cycles.
  std::size_t numberOfEquilibrationCycles;      ///< Number of equilibration cycles.

  std::size_t printEvery;              ///< Frequency of printing status reports.
  std::size_t optimizeMCMovesEvery;    ///< Frequency of optimizing MC moves.
  std::size_t rescaleWangLandauEvery;  ///< Frequency of adjusting the Wang-Landau biasing factors.

  std::size_t numberOfBlocks;  ///< Number of blocks for the block-error estimation.

  std::pair<double, double> reweightingPressureRange;  ///< Pressure range of the reweighted isotherms [Pa].
  std::size_t reweightingNumberOfPressures;            ///< Number of log-spaced pressures of the reweighted isotherms.

  SimulationStage simulationStage{SimulationStage::Uninitialized};  ///< Current simulation stage.

  std::vector<double> temperatures;  ///< The temperature ladder [K].
  double referencePressure;          ///< The reference pressure ('ExternalPressure') [Pa].
  std::size_t numberOfTemperatures;  ///< Number of temperatures in the ladder.
  std::size_t numberOfWindows;       ///< Number of macrostate windows per temperature.
  std::size_t numberOfWalkers;       ///< numberOfTemperatures x numberOfWindows == number of threads.

  std::size_t minMacrostate;  ///< Global minimum macrostate (molecule count).
  std::size_t maxMacrostate;  ///< Global maximum macrostate (molecule count).
  /// Window boundaries: window w spans the macrostates [windowBoundaries[w], windowBoundaries[w+1]]
  /// (neighboring windows share their endpoint macrostate, which stitches the collection matrices).
  std::vector<std::size_t> windowBoundaries;

  std::vector<System> systems;        ///< One walker per (temperature, window) pair.
  std::vector<RandomNumber> randoms;  ///< Independent random-number stream per walker.

  /// Cumulative production-only collection-matrix snapshots per walker at the production block
  /// boundaries (the per-block increments give the block errors); each worker thread appends
  /// exclusively to its own vector.
  std::vector<std::vector<std::vector<double3>>> blockCollectionMatrices;

  /// The collection matrix and visit histogram of every walker at the start of the production
  /// stage. The collection matrix itself is never reset (its entries are unbiased acceptance
  /// probabilities, valid across bias updates, so all statistics accumulate); these snapshots
  /// are subtracted to obtain the production-only block increments and coverage diagnostics.
  std::vector<std::vector<double3>> productionStartCollectionMatrices;
  std::vector<std::vector<std::size_t>> productionStartHistograms;

  std::vector<std::size_t> stepsPerWalker;  ///< Production MC steps performed per walker.

  /// Cycles completed in the previous stages; the time-evolution properties (number of molecules,
  /// volume) are indexed by the absolute cycle number counted over all stages.
  std::size_t absoluteCycleOffset{0};

  std::ofstream stream;            ///< The combined output stream (analysis, timings).
  std::string outputJsonFileName;  ///< Filename for the combined output JSON file.
  nlohmann::json outputJson;       ///< Combined output data in JSON format.

  std::mutex outputMutex;  ///< Guards the combined output stream for in-run progress lines.

  // Per-walker output: each worker thread writes exclusively to its own stream, so no locking is
  // needed for the periodic status reports.
  std::vector<std::ofstream> walkerStreams;      ///< One output stream per walker.
  std::vector<std::string> walkerJsonFileNames;  ///< Filename of the JSON output file per walker.
  std::vector<nlohmann::json> walkerJsons;       ///< JSON output data per walker.

  std::chrono::duration<double> totalPreInitializationSimulationTime{0};  ///< Total time for pre-initialization.
  std::chrono::duration<double> totalInitializationSimulationTime{0};     ///< Total time for initialization stage.
  std::chrono::duration<double> totalEquilibrationSimulationTime{0};      ///< Total time for equilibration stage.
  std::chrono::duration<double> totalProductionSimulationTime{0};         ///< Total time for production stage.
  std::chrono::duration<double> totalAnalysisTime{0};                     ///< Total time for the TMMC analysis.
  std::chrono::duration<double> totalSimulationTime{0};                   ///< Total simulation time.

  /**
   * \brief The walker at grid point (temperature index, window index).
   */
  std::size_t walkerIndex(std::size_t temperatureIndex, std::size_t windowIndex) const
  {
    return temperatureIndex * numberOfWindows + windowIndex;
  }

  /**
   * \brief Runs the simulation and the transition-matrix analysis.
   */
  void run();

  /**
   * \brief Sets up the walkers (grows the initial configurations into their windows in parallel),
   *        the output files and the interpolation grids.
   */
  void setup();

  /**
   * \brief Runs one simulation stage: every walker in its own thread (the walkers are fully
   *        independent; the threads only join at the end of the stage).
   */
  void runStage(SimulationStage stage, std::size_t numberOfCycles);

  /**
   * \brief One Monte Carlo cycle of a single walker: random moves plus the TMMC state sampling
   *        (histogram update and periodic bias adjustment).
   */
  void performWalkerCycle(std::size_t walkerId, SimulationStage stage, std::size_t currentBlock);

  /**
   * \brief The cumulative collection matrix of a walker with the production-start snapshot
   *        subtracted: the statistics collected during the production stage only.
   */
  std::vector<double3> productionCollectionMatrix(std::size_t walkerId) const;

  /**
   * \brief Generates the final combined output: coverage statistics, the transition-matrix
   *        analysis and the timings.
   */
  void output();

  /**
   * \brief Combines the collection matrices of the windows per temperature into the macrostate
   *        probability distribution ln Pi(N) over the full range, writes it together with the
   *        reweighted isotherms (exact in the fugacity): the equilibrium isotherm plus the
   *        adsorption and desorption branches (conditional averages over the low-/high-density
   *        basin of Pi(N) where it is bimodal - the hysteresis loop). For bulk boxes the
   *        vapor-liquid coexistence is additionally located with the equal-weight criterion.
   *        The error bars come from re-deriving ln Pi(N) from the per-block collection-matrix
   *        increments.
   */
  void performTransitionMatrixAnalysis();

  /**
   * \brief Writes the final per-walker reports (energy drift, move statistics and averages) to
   *        the per-walker output files.
   */
  void writeWalkerFinalReports(std::vector<RunningEnergy>& recomputed);
};
