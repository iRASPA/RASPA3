module;

export module dctst_free_energy;

import std;

import randomnumbers;
import double3;
import system;
import input_reader;
import running_energy;
import json;
import dctst_free_energy_settings;
import dctst_barrier;

/**
 * \brief A uniform grid of fine bins over the reaction coordinate.
 *
 * Unconditioned mode: x is the reduced reaction coordinate over one period [-lambda/2, lambda/2), the barrier
 * at x = 0 and the cage center at x = +-lambda/2. Occupancy-resolved mode: x is the coordinate relative to
 * the center of cage A along the normal, [-lambda/2, 3 lambda/4), the window to cage B at x = lambda/2.
 */
export struct DcTSTGrid
{
  double begin{0.0};
  double end{1.0};
  double width{1.0};
  std::size_t numberOfBins{1};
  double barrier{0.0};            ///< Coordinate of the dividing plane.
  double normalizationBegin{0.0};  ///< Range over which exp(-beta F) is integrated (the cage).
  double normalizationEnd{1.0};

  [[nodiscard]] bool contains(double x) const { return x >= begin && x < end; }
  [[nodiscard]] std::size_t bin(double x) const
  {
    const double f = (x - begin) / width;
    return std::min(numberOfBins - 1uz, static_cast<std::size_t>(std::max(0.0, std::floor(f))));
  }
  [[nodiscard]] double center(std::size_t b) const { return begin + (static_cast<double>(b) + 0.5) * width; }
};

/**
 * \brief One umbrella window: a harmonic bias on the reaction coordinate of the tagged molecule plus the
 * optional flattening bias, tabulated per fine bin (piecewise constant, which makes the reweighting exact).
 */
export struct DcTSTWindow
{
  double center{0.0};        ///< q_k [A]
  double stiffness{0.0};     ///< beta kappa [1/A^2]
  double supportHalfWidth{0.0};  ///< Half-width of the reinsertion support, 4 sigma [A].
  std::size_t replica{0};    ///< The replica (thread) that samples this window.
  /// With more replicas than window centers the centers are repeated so that every replica samples; this is
  /// the index of the first window with the same center (the window itself when it is not a repeat).
  std::size_t source{0};
  std::vector<double> bias;  ///< Reduced bias per fine bin [kT].
};

/**
 * \brief Per-window accumulators.
 */
export struct DcTSTWindowAccumulator
{
  std::size_t numberOfBlocks{};
  std::size_t numberOfBins{};
  std::vector<std::vector<double>> histogram;  ///< [block][fine bin], the tagged molecule.
  std::vector<double> count;                   ///< [block]
  std::vector<double> sumX;                    ///< [block]
  std::vector<double> sumXSquared;             ///< [block]
  std::vector<double> sumBias;                 ///< [block], the reduced bias w_k(x) of the samples [kT].

  std::size_t translationTrials{0};
  std::size_t translationAccepted{0};
  std::size_t translationTotalTrials{0};
  std::size_t translationTotalAccepted{0};
  double maximumTranslation{0.5};
  std::size_t reinsertionTrials{0};
  std::size_t reinsertionAccepted{0};
  std::size_t reinsertionTotalTrials{0};
  std::size_t reinsertionTotalAccepted{0};
  std::size_t rejectedByRange{0};      ///< Tagged moves rejected by the coordinate range (conditioned mode).
  std::size_t rejectedByOccupancy{0};  ///< Moves of other molecules reverted by the occupancy constraint.
  std::size_t cageChanges{0};
  double startCoordinate{0.0};  ///< x of the tagged molecule when the window started.
  double endCoordinate{0.0};
  /// Joint occupancy histogram (n_A, n_B) of the two cages sharing the window while the tagged molecule is
  /// within 'barrierOccupancyHalfWidth' of the dividing plane, the tagged molecule counted in A and the two
  /// assignments (A, B) = (behind, ahead) and (ahead, behind) each with weight 1/2. The umbrella bias acts on
  /// the coordinate of the tagged molecule only, so conditioned on it the sample is unbiased: this is
  /// P(n_A, n_B | q*), which with the pilot's P(n_A, n_B) gives the conditioned transition-state rates
  /// k^TST(n_A -> n_B) = k^TST P(n_A, n_B | q*) / P(n_A, n_B) without constrained runs. Index as in the pilot.
  std::vector<double> barrierOccupancyPairs;
  double barrierOccupancySamples{0.0};
  static constexpr double barrierOccupancyHalfWidth = 0.25;  ///< [A]

  DcTSTWindowAccumulator() = default;
  DcTSTWindowAccumulator(std::size_t numberOfBlocks, std::size_t numberOfBins, double maximumTranslation);

  void optimize(double targetAcceptance = 0.5);
};

/**
 * \brief Per-replica accumulators of the pilot stage.
 */
export struct DcTSTPilotAccumulator
{
  static constexpr std::size_t maximumOccupancy = 31;

  std::size_t numberOfBlocks{};
  std::size_t numberOfBins{};
  /// histogram[block][fine bin]: reaction beads of all molecules of the tagged component (unconditioned grid).
  std::vector<std::vector<double>> histogram;
  /// occupancyPairs[block][n_A * (maximumOccupancy + 1) + n_B], see the transmission driver.
  std::vector<std::vector<double>> occupancyPairs;
  /// conditionedHistogram[block][fine bin]: beads of molecules whose cage pair has the requested occupancies,
  /// on the conditioned grid (occupancy-resolved mode only).
  std::vector<std::vector<double>> conditionedHistogram;
  /// Thermal speed of the reaction bead along the normal, averaged over the sampled orientations: the sum
  /// over the sampled molecules of sqrt(2 kT sigma^2 / pi) with sigma^2 the (orientation-dependent) inverse
  /// effective mass of the bead along the normal, and the number of terms. For a single-bead or
  /// center-of-mass bead this is sqrt(2 kT / pi M) exactly.
  double beadSpeedSum{0.0};
  double beadSpeedCount{0.0};

  DcTSTPilotAccumulator() = default;
  DcTSTPilotAccumulator(std::size_t numberOfBlocks, std::size_t numberOfBins, std::size_t numberOfConditionedBins);

  [[nodiscard]] static std::size_t occupancyIndex(std::size_t nA, std::size_t nB)
  {
    return std::min(nA, maximumOccupancy) * (maximumOccupancy + 1uz) + std::min(nB, maximumOccupancy);
  }
};

/**
 * \brief Per-replica state of the window stage.
 */
export struct DcTSTWindowState
{
  std::size_t taggedLocalIndex{0};  ///< Index of the tagged molecule within the tagged component.
  // occupancy-resolved mode
  std::vector<std::size_t> cageOfMolecule;
  std::vector<std::size_t> occupancy;  ///< Per cage, the tagged molecule excluded.
  std::size_t cageA{0};
  std::size_t cageB{0};
};

/**
 * \brief dcTST free-energy driver: the free-energy profile F(q) of a tagged component along the barrier
 * normal, the probability density P(q*) at the dividing plane and the transition-state rate k^TST, from
 * umbrella sampling of one tagged molecule combined with the unbiased histogram of all molecules by binned
 * multistate reweighting.
 */
export struct DynamicallyCorrectedTSTFreeEnergy
{
  enum class SimulationStage : std::size_t
  {
    Initialization = 0,
    Pilot = 1,
    Windows = 2
  };

  DynamicallyCorrectedTSTFreeEnergy() = delete;
  DynamicallyCorrectedTSTFreeEnergy(const DynamicallyCorrectedTSTFreeEnergy&) = delete;
  DynamicallyCorrectedTSTFreeEnergy& operator=(const DynamicallyCorrectedTSTFreeEnergy&) = delete;

  DynamicallyCorrectedTSTFreeEnergy(InputReader& reader);

  RandomNumber random;
  DcTSTFreeEnergySettings settings;
  DcTSTBarrier barrier;
  DcTSTCageLattice cageLattice;
  DcTSTGrid grid;  ///< The grid of the reweighting (unconditioned or conditioned, see 'DcTSTGrid').
  DcTSTGrid pilotGrid;  ///< The unconditioned period grid of the pilot histogram.
  bool occupancyMode{false};

  std::size_t numberOfInitializationCycles;
  std::size_t printEvery;
  std::size_t optimizeMCMovesEvery;
  std::size_t numberOfBlocks;

  std::size_t numberOfReplicas;
  std::vector<System> systems;
  std::vector<RandomNumber> randoms;

  std::vector<DcTSTPilotAccumulator> pilot;        ///< Per replica.
  std::vector<DcTSTWindow> windows;                ///< All windows.
  std::vector<DcTSTWindowAccumulator> windowData;  ///< Per window.
  std::vector<DcTSTWindowState> windowStates;      ///< Per replica.
  std::vector<double> flatteningBias;              ///< w_0 per fine bin [kT] (zero when disabled).

  std::ofstream stream;
  std::string outputJsonFileName;
  nlohmann::json outputJson;
  std::vector<std::ofstream> replicaStreams;

  std::chrono::duration<double> totalInitializationTime{0};
  std::chrono::duration<double> totalPilotTime{0};
  std::chrono::duration<double> totalWindowTime{0};
  std::chrono::duration<double> totalSimulationTime{0};

  void run();
  void setup();
  void runInitialization();
  void runPilot();
  void buildWindows();
  void runWindows();
  void output();

  // --- per replica -------------------------------------------------------------------------------

  /// One unconstrained MC cycle (initialization and pilot).
  void performUnconstrainedCycle(std::size_t replicaId, SimulationStage stage);
  /// One MC cycle of the window stage: biased moves of the tagged molecule, free moves of the others.
  void performWindowCycle(std::size_t replicaId, std::size_t windowIndex);
  std::optional<RunningEnergy> biasedTranslationMove(std::size_t replicaId, std::size_t windowIndex);
  std::optional<RunningEnergy> biasedReinsertionMove(std::size_t replicaId, std::size_t windowIndex);
  /// Samples the pilot histograms.
  void samplePilot(std::size_t replicaId, std::size_t block);
  /// Samples the window histogram.
  void sampleWindow(std::size_t replicaId, std::size_t windowIndex, std::size_t block);
  /// Chooses the tagged molecule for a window (the molecule of the tagged component nearest to the window
  /// center) and, in the occupancy-resolved mode, prepares the cage occupancies.
  void startWindow(std::size_t replicaId, std::size_t windowIndex);

  /// Reaction coordinate x of a bead position on the grid of the reweighting (reduced; conditioned frame in
  /// the occupancy-resolved mode).
  [[nodiscard]] double coordinate(const System& system, const DcTSTWindowState& state, double3 position) const;
  /// Folded coordinate (|x| in mirror mode).
  [[nodiscard]] double foldedCoordinate(double x) const;
  /// The bias of window 'windowIndex' at coordinate x [kT]; a large value outside the grid.
  [[nodiscard]] double windowBias(std::size_t windowIndex, double x) const;
  /// Global index (into moleculeData) of the first molecule of the tagged component.
  [[nodiscard]] std::size_t taggedComponentOffset(const System& system) const;
  /// Occupancy-resolved mode: all molecules' cages.
  [[nodiscard]] std::vector<std::size_t> countOccupancies(const System& system) const;
  /// Occupancy-resolved mode: a point with coordinate x in the frame of cage A, uniform laterally over cage A.
  [[nodiscard]] std::optional<double3> pointInCageFrame(RandomNumber& rng, const System& system,
                                                        const DcTSTWindowState& state, double x) const;

  // --- analysis ----------------------------------------------------------------------------------

  struct ProfileResult
  {
    bool available{false};
    bool converged{false};
    std::size_t iterations{0};
    double residual{0.0};
    std::vector<double> stateFreeEnergies;  ///< g_k of the states (0 = pilot).
    std::vector<double> fineFreeEnergy;     ///< beta F per fine bin (minimum zero), 'defined' flags.
    std::vector<bool> fineDefined;
    std::vector<double> fineFreeEnergyError;
    std::vector<double> q;  ///< Reporting grid.
    std::vector<double> freeEnergy;
    std::vector<double> freeEnergyError;
    std::vector<bool> defined;
    std::vector<double> pilotFreeEnergy;  ///< beta F from the pilot histogram alone, reporting grid.
    std::vector<bool> pilotDefined;
    double barrierProbabilityBin{0.0};  ///< P(q*) from the fine bin containing the barrier [1/A].
    double barrierProbabilityBinError{0.0};
    double barrierProbability{0.0};  ///< P(q*) from the parabolic fit of the top [1/A].
    double barrierProbabilityError{0.0};
    double barrierProbabilityErrorHalfBlocks{0.0};  ///< The same with half the number of blocks (block-size check).
    double barrierFreeEnergy{0.0};                  ///< beta F(q*) - beta F(min) from the fit.
    double barrierFreeEnergyError{0.0};
    double barrierFreeEnergyBin{0.0};
    double barrierOffset{0.0};  ///< Position of the fitted maximum relative to the dividing plane [A].
    double barrierOffsetError{0.0};
    double curvature{0.0};  ///< beta F'' at the maximum [1/A^2] (negative for a maximum).
    bool fitValid{false};
    std::size_t fitPoints{0};
    double pilotBarrierProbability{0.0};  ///< P(q*) from the pilot histogram alone (bin value).
    double pilotBarrierProbabilityError{0.0};
    bool pilotBarrierDefined{false};
    std::vector<double> overlap;  ///< Overlap of consecutive windows (and pilot) in the order of the centers.
    std::vector<double> windowMean;
    std::vector<double> windowSigma;
    std::vector<double> windowCounts;
    std::vector<double> windowMeanBias;  ///< <w_k> over the samples of window k [kT].
    std::size_t numberOfHalfBlocks{0};   ///< Number of merged blocks behind barrierProbabilityErrorHalfBlocks.
    /// overlapMatrix[i * (K + 1) + j]: histogram overlap sum_b min(p_i(b), p_j(b)) of the states, state 0 the pilot.
    std::vector<double> overlapMatrix;
    /// Consistency of the pilot histogram with the reweighted profile where the pilot is well sampled (fine
    /// bins with at least 'pilotConsistencyMinimumCounts' pilot counts): rms of beta F_pilot - beta F after
    /// removing the mean offset, and the number of bins.
    double pilotConsistencyRms{0.0};
    std::size_t pilotConsistencyBins{0};
    static constexpr double pilotConsistencyMinimumCounts = 100.0;
  };

  struct OccupancyDistribution
  {
    std::size_t size{0};
    std::vector<double> joint;
    std::vector<double> marginalA;
    std::vector<double> marginalB;
    double meanA{0.0};
    double meanB{0.0};
    double total{0.0};
    double requestedProbability{0.0};
    double requestedProbabilityError{0.0};
    bool available{false};
    /// P(n_A, n_B | q*) from the window stage (see DcTSTWindowAccumulator::barrierOccupancyPairs).
    std::vector<double> jointAtBarrier;
    double totalAtBarrier{0.0};
    bool availableAtBarrier{false};
  };

  [[nodiscard]] ProfileResult analyzeProfile() const;
  [[nodiscard]] OccupancyDistribution analyzeOccupancies() const;
  /// Thermal speed <|qdot|> of the reaction bead along the normal from the pilot orientations (falls back to
  /// sqrt(2 kT / pi M) of the molecular mass when the pilot did not sample).
  [[nodiscard]] double beadThermalSpeed() const;
  /// The pilot histogram pooled over replicas and blocks on the grid of the reweighting (the conditioned
  /// histogram in the occupancy-resolved mode).
  [[nodiscard]] std::vector<double> pooledPilotHistogram(std::optional<std::size_t> block = std::nullopt) const;
};
