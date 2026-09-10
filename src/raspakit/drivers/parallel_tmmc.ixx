module;

export module parallel_tmmc;

import std;

import randomnumbers;
import averages;
import system;
import input_reader;
import running_energy;
import double3;
import archive;
import json;
import isotherm_bet;

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
 * At the low-N end that estimator is badly conditioned in a strongly binding adsorbent: the two
 * transition probabilities differ by many orders of magnitude (some 19 natural logs per molecule
 * for nitrogen in ferrierite at 77 K), the deletion side is an average with a heavy upper tail,
 * and every increment error propagates into ln Pi of all higher macrostates, which moves the
 * filling pressure. The same increment also follows exactly from a Widom test insertion,
 *
 *     ln Pi(N+1) - ln Pi(N) = ln(beta f V <W>_N / (N+1)),
 *
 * with <W>_N the mean Rosenbluth weight of a test insertion into the N-molecule system (the
 * N = 0 case is the Henry coefficient). Test insertions are sampled per macrostate during
 * production and replace the collection-matrix increments over the contiguous low-N block.
 * The block ends where the test-insertion mean becomes hopeless (relative error of order 1,
 * one lucky insertion in a full pore), not where it is merely noisy: a 0.35 relative-error
 * cut closed the ferrierite (N, λ) block at N = 1 after three increments at 0.35–0.39. Every
 * counted increment inside the block is used. Insertions become hopeless as the pore fills,
 * exactly where the collection matrix of a 1D-N walk is well conditioned, so the two
 * estimators are complementary and the handover is automatic.
 *
 * ln Pi(N) is exact at the reference fugacity f_ref (from 'ExternalPressure' through the
 * Peng-Robinson equation of state) and reweights exactly to any other fugacity,
 * ln Pi(N; f) = ln Pi(N; f_ref) + N ln(f / f_ref), giving continuous adsorption isotherms
 * <N>(f) per temperature. The error bars come from re-deriving ln Pi from the per-block
 * increments of the collection matrices. With ComputeBET the Rouquerol window is taken from
 * the full-data fit and slope/intercept are refit per block for a CI on the BET area.
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
/// Cycle counts and TMMC analysis controls for a programmatic (no input-file) construction.
export struct ParallelTMMCParameters
{
  std::size_t numberOfProductionCycles{10000};
  std::size_t numberOfPreInitializationCycles{0};
  std::size_t numberOfInitializationCycles{5000};
  std::size_t numberOfEquilibrationCycles{5000};
  std::size_t printEvery{5000};
  std::size_t optimizeMCMovesEvery{5000};
  std::size_t rescaleWangLandauEvery{5000};
  std::size_t writeBinaryRestartEvery{0};
  std::size_t numberOfBlocks{5};
  std::size_t numberOfWindows{8};
  std::size_t tmmcUpdateEvery{10000};
  std::pair<double, double> reweightingPressureRange{1.0, 101325.0};
  std::size_t reweightingNumberOfPressures{100};
  bool computeBET{false};
  std::optional<std::size_t> randomSeed{};
};

/// One point of a TMMC-reweighted isotherm, after the analysis has run.
export struct TMMCIsothermPoint
{
  double pressure{0.0};               ///< Pressure [Pa].
  double fugacity{0.0};               ///< Fugacity [Pa].
  double moleculesPerCell{0.0};       ///< Absolute loading [molecules / simulation cell].
  double moleculesPerCellError{0.0};  ///< Block error on the loading.
  bool bimodal{false};                ///< True when Pi(N) has two basins at this pressure.
};

/// The three TMMC isotherm branches at one temperature.
export struct TMMCIsotherm
{
  double temperature{0.0};
  std::vector<TMMCIsothermPoint> equilibrium;
  std::vector<TMMCIsothermPoint> adsorption;
  std::vector<TMMCIsothermPoint> desorption;
};

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

  /**
   * \brief Constructs the driver from a single template system and a temperature ladder, without
   *        an input file. The system is replicated into one walker per (temperature, window) pair
   *        as in the InputReader constructor. The template system's tmmc min/max macrostates must
   *        already be set (minimum < maximum); TMMC flags and the bias-update interval are applied
   *        from \p parameters.
   */
  ParallelTMMC(System templateSystem, std::vector<double> temperatures,
               ParallelTMMCParameters parameters = {});

  std::uint64_t versionNumber{1};  ///< Version number for serialization.

  RandomNumber random;  ///< Random number generator (seeding).

  std::size_t numberOfProductionCycles;         ///< Number of production cycles.
  std::size_t numberOfPreInitializationCycles;  ///< Number of pre-initialization cycles.
  std::size_t numberOfInitializationCycles;     ///< Number of initialization cycles.
  std::size_t numberOfEquilibrationCycles;      ///< Number of equilibration cycles.

  std::size_t printEvery;               ///< Frequency of printing status reports.
  std::size_t optimizeMCMovesEvery;     ///< Frequency of optimizing MC moves.
  std::size_t rescaleWangLandauEvery;   ///< Frequency of adjusting the Wang-Landau biasing factors.
  std::size_t writeBinaryRestartEvery;  ///< Frequency of writing the binary restart file (0 disables).

  std::size_t numberOfBlocks;  ///< Number of blocks for the block-error estimation.

  std::pair<double, double> reweightingPressureRange;  ///< Pressure range of the reweighted isotherms [Pa].
  std::size_t reweightingNumberOfPressures;            ///< Number of log-spaced pressures of the reweighted isotherms.

  /// Extract a nitrogen BET area from each reweighted isotherm (Rouquerol, P0 = 101325 Pa).
  bool computeBET{false};
  /// 'ReweightingPressureRange': 'auto' (or omitted under ComputeBET): Henry-to-P0 isotherm grid.
  bool autoReweightingPressureRange{false};
  /// 'MacroStateMaximumNumberOfMolecules': 'auto' (or omitted under ComputeBET): P0 occupancy scout.
  bool autoMacroStateMaximum{false};
  /// Filled when the pressure span was placed from a Widom Henry coefficient.
  std::optional<NitrogenBETPressurePlan> nitrogenBETPressurePlan;
  /// Filled when N_max was placed from a P0 occupancy scout.
  std::optional<NitrogenBETFillingCeiling> nitrogenBETFillingCeiling;
  /// MC moves per cycle: the global maxMacrostate, not the current occupancy or the window width.
  std::size_t numberOfStepsPerCycle{20};

  SimulationStage simulationStage{SimulationStage::Uninitialized};  ///< Current simulation stage.

  /// Completed cycles of the stage the binary restart file was written in; a restarted run
  /// continues that stage from this cycle.
  std::size_t cyclesCompletedThisStage{0};

  /// The end-of-cycle count communicated by the worker threads to the barrier completion
  /// (written by walker 0 before arriving; read while all threads are parked). Not serialized.
  std::atomic<std::size_t> checkpointCycle{0};

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

  /// Widom test-insertion statistics per walker on the global macrostate grid, accumulated during
  /// production only: the sum and the sum of squares of the Rosenbluth weight and the number of
  /// insertions, per macrostate. They give the exact increment
  /// ln Pi(N+1) - ln Pi(N) = ln(beta f V <W>_N / (N+1)) at the reference fugacity, which anchors
  /// the low-N end of ln Pi where the collection-matrix estimate is worst conditioned.
  std::vector<std::vector<double>> widomWeightSums;
  std::vector<std::vector<double>> widomWeightSquaredSums;
  std::vector<std::vector<std::size_t>> widomInsertions;

  std::vector<std::size_t> stepsPerWalker;  ///< Production MC steps performed per walker.

  /// Reweighted isotherms (equilibrium, adsorption and desorption branches), one per temperature,
  /// filled by performTransitionMatrixAnalysis.
  std::vector<TMMCIsotherm> reweightedIsotherms;

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
   * \brief Combines the windows per temperature into the macrostate probability distribution
   *        ln Pi(N) over the full range by reconstructing ln Pi in each window and matching
   *        the gauge at the shared endpoint. Writes it together with the reweighted isotherms
   *        (exact in the fugacity): the equilibrium isotherm plus the
   *        adsorption and desorption branches (conditional averages over the low-/high-density
   *        basin of Pi(N) where it is bimodal - the hysteresis loop). For bulk boxes the
   *        vapor-liquid coexistence is additionally located with the equal-weight criterion.
   *        The error bars come from re-deriving ln Pi(N) from the per-block collection-matrix
   *        increments.
   */
  void performTransitionMatrixAnalysis();

  /**
   * \brief Replicates \p templateSystem into one walker per (temperature, window) pair and pins
   *        each walker's temperature, Peng-Robinson fugacity and macrostate window.
   */
  void initializeWalkers(System templateSystem);

  /**
   * \brief Writes the final per-walker reports (energy drift, move statistics and averages) to
   *        the per-walker output files.
   */
  void writeWalkerFinalReports(std::vector<RunningEnergy>& recomputed);

  /**
   * \brief Writes the binary restart file. Called from the barrier completion, where all worker
   *        threads are parked, so the full driver state is consistent.
   */
  void writeBinaryRestartFile(std::size_t cyclesCompleted) noexcept;

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const ParallelTMMC& ptmmc);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, ParallelTMMC& ptmmc);
};
