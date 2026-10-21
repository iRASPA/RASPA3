module;

export module dctst;

import std;

import randomnumbers;
import double3;
import system;
import input_reader;
import running_energy;
import json;
import dctst_settings;
import dctst_barrier;

/**
 * \brief Per-replica bookkeeping of the cage occupancies in the occupancy-resolved transmission stage.
 */
export struct DcTSTOccupancyState
{
  bool active{false};  ///< The constraint is enforced (from 'prepareOccupancies' at the start of the transmission stage).
  std::vector<std::size_t> cageOfMolecule;  ///< Cage (lattice index) of the reaction bead of every molecule.
  std::vector<std::size_t> occupancy;       ///< Number of molecules per cage, the tagged molecule excluded.
  std::size_t cageA{0};                     ///< The cage the tagged molecule hops from.
  std::size_t cageB{0};                     ///< The cage the tagged molecule hops to.
  std::size_t cageChanges{0};               ///< Accepted moves of other molecules that changed their cage.
  std::size_t rejectedMoves{0};             ///< Of these, the moves reverted by the constraint (A or B involved).
  std::size_t rejectedTaggedMoves{0};       ///< In-plane moves of the tagged molecule rejected by the constraint.
};

/**
 * \brief Per-replica accumulators of the self and collective transmission functions.
 *
 * For every launched trajectory pair (forward and time-reversed, initial reaction-bead speed
 * |qdot|) and every stored time t the accumulators receive
 *
 *   self[block][t]                 += |qdot| (Delta_S^+(t) - Delta_S^-(t))
 *   collective[block][s][c][t]     += |qdot| sum_{j != tagged, shell_j <= s, component_j = c}
 *                                             (Delta_j^+(t) - Delta_j^-(t))
 *   speed[block]                   += |qdot|
 *
 * with Delta_S = c_tagged(t) + 1/2 (+1/2 in cage B, -1/2 in cage A) and Delta_j = c_j(t) - c_j(0)
 * the net number of cages molecule j advanced along the normal. The transmission functions are
 * kappa_S(t) = self / speed and kappa_C^(s)(t) = kappa_S(t) + sum_c collective[s][c] / speed; both
 * are exactly one at t = 0+. Shell index 'maximumShell + 1' holds the sum over all molecules.
 */
export struct DcTSTTransmissionAccumulator
{
  std::size_t numberOfBlocks{};
  std::size_t numberOfShells{};  ///< maximumShell + 2 (shells 0..maximumShell, then "all").
  std::size_t numberOfComponents{};
  std::size_t numberOfTimes{};

  std::vector<double> speed;                      ///< [block]
  std::vector<std::vector<double>> self;          ///< [block][time]
  std::vector<std::vector<double>> selfTwoState;  ///< [block][time], side-of-the-plane variant (Bennett-Chandler)
  std::vector<std::vector<double>>
      selfCommitted;  ///< [block][time], side of the plane frozen once |q| >= commit distance (RASPA2 convention)
  /// [block][time], the committor-smoothed cage coordinate C(q) of the tagged molecule (Ruiz-Montero, Frenkel
  /// and Brey): the step of the net cage count at every window replaced by the committor p_B of the free-energy
  /// profile, which removes the variance of the A/B decision from kappa_S; zero when no committor is available.
  std::vector<std::vector<double>> selfSmooth;
  std::vector<std::vector<double>>
      collective;  ///< [block][(shell * numberOfComponents + component) * numberOfTimes + time]
  std::vector<std::size_t> numberOfTrajectoryPairs;  ///< [block]

  double sumSpeedSquared{0.0};         ///< Sum of |qdot|^2 (velocity-distribution check).
  double sumRelativeEnergyDrift{0.0};  ///< Sum over trajectories of the mean relative drift of the conserved energy.
  double maximumRelativeEnergyDrift{0.0};
  std::size_t numberOfTrajectories{0};
  /// Trajectories in which the tagged molecule left the pair of cages A/B (reached a second window):
  /// the signature of a molecule that was not thermalized in the cage.
  std::size_t numberOfLongJumpTrajectories{0};
  /// Trajectories in which the tagged molecule reached the commit distance (in either cage).
  std::size_t numberOfCommittedTrajectories{0};

  // --- event-terminated estimator -------------------------------------------------------------------
  /// Channels of the local axial hop events, in the frame of the tagged hop (for the time-reversed
  /// trajectory the cages A <-> B, A' <-> B' and the direction are mirrored so that "forward" is always
  /// the direction of the tagged molecule's commitment).
  enum class EventChannel : std::size_t
  {
    TaggedBackHop = 0,     ///< Tagged molecule returns from B to A.
    TaggedForwardHop = 1,  ///< Tagged molecule runs on from B to B'.
    OtherAtoB = 2,         ///< Another molecule follows through the same window A -> B.
    OtherBtoA = 3,         ///< Another molecule leaves B into A (the excess in B relaxes backwards).
    OtherAtoAPrime = 4,    ///< Another molecule leaves A backwards into A'.
    OtherAPrimeToA = 5,    ///< Another molecule refills A from A'.
    OtherBtoBp = 6,        ///< Another molecule leaves B forwards into B' (the excess relaxes forwards).
    OtherBptoB = 7,        ///< Another molecule enters B from B'.
    Count = 8
  };
  /// selfEvent[block] += |qdot| (Delta_S^+ - Delta_S^-) with Delta_S the net cage count of the tagged molecule at
  /// termination; collectiveEvent[block] likewise for Delta_S + Q_local, Q_local the signed sum of the axial
  /// window crossings of the other molecules in the column of cages A', A, B, B'.
  std::vector<double> selfEvent;        ///< [block]
  std::vector<double> collectiveEvent;  ///< [block]
  /// The same estimator restricted to the event-terminated trajectories, per trajectory (not per pair): the
  /// forward (s = +1) and time-reversed (s = -1) members contribute |qdot| s Delta separately, and
  /// speedEventTerminated[block] sums |qdot| over the contributing members, so that
  /// kappa^cond = 2 selfEventConditioned / speedEventTerminated (identical to kappa^event when no trajectory
  /// is capped). Only the identity of the first channel enters the first-order estimate, not its timing, so
  /// the restriction removes the bias of the capped trajectories (which carry no event) without introducing
  /// a new one.
  std::vector<double> speedEventTerminated;        ///< [block]
  std::vector<double> selfEventConditioned;        ///< [block]
  std::vector<double> collectiveEventConditioned;  ///< [block]
  /// Decomposition of Q_local (in the frame of the tagged hop, |qdot|-weighted, forward and time-reversed
  /// members added) into the crossings before and after the commitment of the tagged molecule, so that
  /// collectiveEvent - selfEvent = collectiveEventBeforeCommitment + collectiveEventAfterCommitment per block.
  std::vector<double> collectiveEventBeforeCommitment;  ///< [block]
  std::vector<double> collectiveEventAfterCommitment;   ///< [block]
  /// The same per channel (phase 0 = before commitment, 1 = after): localQChannel is the |qdot|-weighted signed
  /// contribution to Q_local (zero for the tagged channels, whose effect is in the self part), channelEvents the
  /// number of crossings and channelTimeSum the sum of their times (since launch before commitment, since
  /// commitment after it) [ps].
  std::array<std::array<double, static_cast<std::size_t>(EventChannel::Count)>, 2> localQChannel{};
  std::array<std::array<double, static_cast<std::size_t>(EventChannel::Count)>, 2> channelTimeSum{};
  std::array<std::array<std::size_t, static_cast<std::size_t>(EventChannel::Count)>, 2> channelEvents{};
  /// terminatingEvents[direction][channel]: number of trajectories terminated by that channel; direction 0 is the
  /// forward and 1 the time-reversed trajectory of the pair.
  std::array<std::array<std::size_t, static_cast<std::size_t>(EventChannel::Count)>, 2> terminatingEvents{};
  /// allEvents[direction][channel]: all local axial events counted after commitment (several per trajectory
  /// when 'TerminationEvents' > 1).
  std::array<std::array<std::size_t, static_cast<std::size_t>(EventChannel::Count)>, 2> allEvents{};
  std::size_t numberOfEventTerminatedTrajectories{0};  ///< Ended by the n-th event.
  std::size_t numberOfCappedTrajectories{0};           ///< Reached the cap after commitment.
  std::size_t numberOfUncommittedTrajectories{0};      ///< Reached the cap without ever committing.
  double terminationTimeSum{0.0};                      ///< Over the event-terminated trajectories [ps].
  double terminationTimeMaximum{0.0};                  ///< [ps]
  double integratedTimeSum{0.0};                       ///< Over all trajectories (for the cost per trajectory) [ps].

  // --- hybrid estimator (handover times) --------------------------------------------------------------
  /// Occupancy state of the column (all molecules, the tagged one included) at a handover time, in the frame of
  /// the tagged hop: the tagged molecule is in B.
  struct ColumnState
  {
    std::uint16_t nAPrime{0};
    std::uint16_t nA{0};
    std::uint16_t nB{0};
    std::uint16_t nBPrime{0};
    auto operator<=>(const ColumnState&) const = default;
  };
  /// Per handover time k and block, the hybrid analogues of selfEvent and collectiveEvent - selfEvent: for a
  /// member whose first post-commitment event occurred before tau the values at that event (as in the first-event
  /// estimator), otherwise the committor-smoothed cage coordinate C of the tagged molecule at tau, resp. the
  /// smoothed axial displacement sum_j delta_j of the other molecules of the column accumulated up to tau.
  /// Accumulated as |qdot| (value^+ - value^-) over the pair; members that did not reach tau enter with their
  /// final values.
  std::vector<std::vector<double>> hybridTagged;  ///< [handover][block]
  std::vector<std::vector<double>> hybridOthers;  ///< [handover][block]
  /// hybridOccupancy[k][state] = (sum of |qdot|, number of members) over the members that had no event before tau
  /// and were in that state at tau; the Markov closure acts on this histogram.
  std::vector<std::map<ColumnState, std::pair<double, std::size_t>>> hybridOccupancy;  ///< [handover]
  /// The full neighbourhood of the column: the column itself and the occupancies of the lateral neighbours of A
  /// and of B (sorted, so that states that differ only by a permutation coincide); the Markov closure is evaluated
  /// per state with the actual lateral occupancies instead of the mean one.
  struct NeighbourhoodState
  {
    ColumnState column{};
    std::array<std::uint16_t, 4> lateralA{};
    std::array<std::uint16_t, 4> lateralB{};
    auto operator<=>(const NeighbourhoodState&) const = default;
  };
  std::vector<std::map<NeighbourhoodState, std::pair<double, std::size_t>>> hybridNeighbourhood;  ///< [handover]
  std::vector<std::size_t> hybridEventBefore;  ///< [handover] members whose first event occurred before tau
  std::vector<std::size_t> hybridTruncated;    ///< [handover] members that did not reach tau

  // --- occupancy classes (post-stratification) ----------------------------------------------------------
  /// Per trajectory member, classified by the occupancies (n_from, n_to) of its departure cage (the tagged
  /// molecule included) and its arrival cage at launch: the forward member of a pair hops A -> B with
  /// (m_A + 1, m_B), the time-reversed member B -> A with (m_B + 1, m_A), m the numbers of other molecules at
  /// launch. kappa(n, n') = 2 sum s |qdot| Delta / sum |qdot| over the members of the class (s = +-1 the launch
  /// direction, Delta the plateau value of the member) is the occupancy-resolved recrossing factor without a
  /// constrained run, and with the conditioned k^TST(n -> n') of the free-energy input gives k(n -> n').
  struct ClassSums
  {
    std::vector<double> speed;          ///< [block] sum |qdot|
    std::vector<double> self;           ///< [block] sum s |qdot| Delta (net cage count, plateau)
    std::vector<double> selfSmooth;     ///< [block] committor-smoothed
    std::vector<double> selfCommitted;  ///< [block]
    std::size_t members{0};
  };
  std::map<std::pair<std::size_t, std::size_t>, ClassSums> classes;

  DcTSTTransmissionAccumulator() = default;
  DcTSTTransmissionAccumulator(std::size_t numberOfBlocks, std::size_t numberOfShells, std::size_t numberOfComponents,
                               std::size_t numberOfTimes, std::size_t numberOfHandovers);

  [[nodiscard]] std::size_t collectiveIndex(std::size_t shell, std::size_t component, std::size_t time) const
  {
    return (shell * numberOfComponents + component) * numberOfTimes + time;
  }
};

/**
 * \brief Statistics of the constrained Monte Carlo moves of the tagged molecule.
 */
export struct DcTSTConstrainedMoveStatistics
{
  std::size_t translationTrials{0};
  std::size_t translationAccepted{0};
  std::size_t translationTotalTrials{0};
  std::size_t translationTotalAccepted{0};
  double maximumTranslation{1.0};
  std::size_t reinsertionTrials{0};
  std::size_t reinsertionAccepted{0};
  std::size_t reinsertionTotalTrials{0};
  std::size_t reinsertionTotalAccepted{0};

  void optimize(double targetAcceptance = 0.5);
};

/**
 * \brief Dynamically corrected transition state theory (dcTST) driver.
 *
 * Computes the hopping rate of a tagged component across a dividing surface (a plane in the
 * framework, e.g. through the window between two cages) and from it the self-diffusion and the
 * collective diffusion coefficient:
 *
 *   k_{A->B} = kappa k^TST,   k^TST = <|qdot|>/2 exp(-beta F(q*)) / int_A exp(-beta F(q)) dq
 *   D_S = kappa_S k^TST lambda^2,   D_C = kappa_C k^TST lambda^2
 *
 * The free-energy profile F(q) and the density P(q*) at the dividing plane come from the separate
 * free-energy driver ('DynamicallyCorrectedTSTFreeEnergy'); its JSON output is read through
 * 'FreeEnergyInput' (or P(q*) is given directly with 'BarrierProbability').
 *
 * Stages (every replica runs all of them in its own thread, on its own copy of the system):
 *  1. Initialization: unconstrained Monte Carlo at the requested loading.
 *  2. Transmission: the reaction bead of one tagged molecule is constrained to the barrier plane
 *     (in-plane translations and in-plane CBMC reinsertions), all other molecules move freely.
 *     Every 'LaunchTrajectoryEvery' cycles the current configuration is copied, Maxwell-Boltzmann
 *     velocities are drawn, the tagged bead is oriented towards B, and the forward and the
 *     time-reversed NVE trajectories are integrated. Along both the net cage changes of every
 *     molecule are tracked and the self and collective transmission functions accumulated. The MC
 *     chain continues from the untouched configuration.
 *
 * The replicas are independent; their accumulators are summed at the end and the block-to-block
 * (and replica-to-replica) scatter gives the error estimate. Rigid frameworks only.
 */
export struct DynamicallyCorrectedTST
{
  enum class SimulationStage : std::size_t
  {
    Uninitialized = 0,
    Initialization = 1,
    Transmission = 2
  };

  DynamicallyCorrectedTST() = delete;
  DynamicallyCorrectedTST(const DynamicallyCorrectedTST&) = delete;
  DynamicallyCorrectedTST& operator=(const DynamicallyCorrectedTST&) = delete;

  DynamicallyCorrectedTST(InputReader& reader);

  RandomNumber random;
  DcTSTSettings settings;
  DcTSTBarrier barrier;
  DcTSTCageLattice cageLattice;

  std::size_t numberOfInitializationCycles;
  std::size_t printEvery;
  std::size_t optimizeMCMovesEvery;
  std::size_t numberOfBlocks;

  std::size_t numberOfReplicas;
  std::vector<System> systems;
  std::vector<RandomNumber> randoms;

  std::size_t numberOfTrajectorySteps{};  ///< Integration steps per trajectory.
  std::size_t numberOfTimes{};            ///< Stored time points per trajectory (including t = 0).

  std::vector<DcTSTTransmissionAccumulator> transmission;        ///< Per replica.
  std::vector<DcTSTConstrainedMoveStatistics> constrainedMoves;  ///< Per replica.
  std::vector<DcTSTOccupancyState> occupancyStates;              ///< Per replica (occupancy-resolved mode).

  std::ofstream stream;
  std::string outputJsonFileName;
  nlohmann::json outputJson;
  std::vector<std::ofstream> replicaStreams;
  std::vector<std::ofstream> pairStreams;  ///< Per replica: the per-pair records ('WritePairRecords').

  /// Stored time indices of the plateau window ('PlateauWindow', default the second half of the trajectory).
  std::vector<std::size_t> plateauIndices;

  // --- early termination and recovery ------------------------------------------------------------
  /// A termination signal (SIGTERM/SIGINT/SIGUSR1, e.g. the scheduler's warning before the wall-time kill) was
  /// received: the replicas stop at the next cycle boundary (a trajectory pair in progress is abandoned), and the
  /// analysis and output run over the completed pairs.
  std::atomic<bool> stoppedEarly{false};
  std::vector<std::size_t> completedTransmissionCycles;  ///< Per replica: production cycles completed.
  std::vector<std::chrono::steady_clock::time_point> lastSnapshotTime;  ///< Per replica ('SnapshotInterval').

  std::chrono::duration<double> totalInitializationTime{0};
  std::chrono::duration<double> totalTransmissionTime{0};
  std::chrono::duration<double> totalSimulationTime{0};

  /// P(q*) [1/A] at the dividing plane from 'FreeEnergyInput' or 'BarrierProbability', with its error and a
  /// description of its provenance; absent when neither was given.
  std::optional<double> barrierProbability{};
  double barrierProbabilityError{0.0};
  std::string barrierProbabilitySource{};

  // --- hybrid estimator ----------------------------------------------------------------------------
  /// Committor p_B(q) of the diffusive limit, p_B(q) = int_{-lambda/2}^{q} e^{beta F} / int_{-lambda/2}^{lambda/2}
  /// e^{beta F}, tabulated on 'committorTable.size()' points of q in [-lambda/2, lambda/2) (barrier at 0) from the
  /// profile of 'FreeEnergyInput'; empty when no profile is available (the hybrid estimator is then disabled).
  std::vector<double> committorTable{};
  /// The handover times in integration steps (sorted), empty when the hybrid estimator is off.
  std::vector<std::size_t> handoverSteps{};
  [[nodiscard]] bool hybridEnabled() const { return !handoverSteps.empty(); }
  /// Smooth cage coordinate of a reaction coordinate q: the net cage count floor(q / lambda) + 1/2 with the
  /// step at every window replaced by the committor, so that it equals the expected final cage count of a
  /// molecule found at q.
  [[nodiscard]] double smoothCageCoordinate(double q) const;
  /// Builds 'committorTable' from the free-energy profile (q grid and beta F) of the free-energy driver's JSON,
  /// 'barrierOffset' being the q of the barrier on that grid (lambda/2 for a conditioned profile over cage A).
  void buildCommittorTable(const std::vector<double>& q, const std::vector<double>& betaF, double barrierOffset);

  /// Occupancy-resolved hopping rates k(n -> n') [1/ps] (n the occupancy of the departure cage, the hopping
  /// molecule included; n' that of the arrival cage): from the occupancy classes of this run combined with the
  /// conditioned k^TST(n -> n') of the free-energy input (post-stratification), overridden by the conditioned
  /// runs of 'OccupancyRateInputs' where given, with a weighted fit log k = c + a n - b n' through them for the
  /// pairs that were not measured.
  struct OccupancyRates
  {
    struct Entry
    {
      double rate{0.0};
      double weight{1.0};  ///< Fit weight (number of trajectory members behind the value).
      double kappa{0.0};   ///< The recrossing factor of the class (post-stratified entries).
      double tstRatio{0.0};  ///< k^TST(n -> n') / k^TST (post-stratified entries).
      std::string source;  ///< "classes" or the file of a conditioned run.
    };
    std::map<std::pair<std::size_t, std::size_t>, Entry> measured;
    bool fitValid{false};
    double fitConstant{0.0};
    double fitDeparture{0.0};
    double fitArrival{0.0};
    double fitRms{0.0};  ///< Weighted root-mean-square residual of log k over the measured pairs.
    /// Rate from the table when measured (and 'tableFirst'), else from the fit (zero when neither is available).
    [[nodiscard]] double rate(std::size_t n, std::size_t nPrime, bool tableFirst) const;
    [[nodiscard]] bool isMeasured(std::size_t n, std::size_t nPrime) const { return measured.contains({n, nPrime}); }
  };
  /// Rates of conditioned runs ('OccupancyRateInputs'), read before the run.
  std::map<std::pair<std::size_t, std::size_t>, OccupancyRates::Entry> externalRates{};
  void loadOccupancyRates();
  /// P(n_A, n_B) and P(n_A, n_B | q*) of the free-energy input (keys "nA,nB"); their ratio is k^TST(n_A -> n_B)/k^TST.
  std::map<std::pair<std::size_t, std::size_t>, double> occupancyWeights{};
  std::map<std::pair<std::size_t, std::size_t>, double> barrierOccupancyWeights{};
  /// Minimum number of trajectory members for an occupancy class to supply a rate.
  static constexpr std::size_t minimumClassMembers = 20;

  void run();
  void setup();
  /// Reads P(q*) from the JSON output of the free-energy driver and checks that it describes the same
  /// barrier (tagged component, lattice distance, normal, occupancies).
  void loadFreeEnergyInput();
  void runStage(SimulationStage stage, std::size_t numberOfCycles);
  void output();

  // --- per replica, executed inside the replica's thread ---------------------------------------

  /// One Monte Carlo cycle in the constrained ensemble: the tagged molecule (molecule 0 of the tagged
  /// component) is restricted to its in-plane moves in both stages, all other molecules move with the
  /// single-molecule moves of their component. The occupancy constraint is enforced once the occupancy
  /// state is active (from the start of the transmission stage).
  void performMonteCarloCycle(std::size_t replicaId);

  /// In-plane translation of the tagged molecule (reaction bead stays on the barrier plane).
  std::optional<RunningEnergy> translationInPlaneMove(std::size_t replicaId);

  /// In-plane CBMC reinsertion of the tagged molecule (reaction bead at a random point of the plane).
  std::optional<RunningEnergy> reinsertionInPlaneMove(std::size_t replicaId);

  /// Builds the starting configuration of the constrained ensemble at the start of the initialization
  /// stage: the adsorbates created with the system are removed, the tagged molecule is grown first with
  /// its reaction bead on the barrier plane (into the empty framework, where this always succeeds; a
  /// window centre is tried first), and the remaining molecules of every component are then grown around
  /// it with CBMC. The tagged molecule never leaves the plane afterwards, so no re-placement into a
  /// crowded system is needed at the start of the transmission stage.
  void createConstrainedConfiguration(std::size_t replicaId);

  /// Occupancy-resolved mode: the cages A and B of the tagged molecule's window from the position of
  /// its reaction bead.
  [[nodiscard]] std::pair<std::size_t, std::size_t> cagesOfWindow(double3 beadPosition) const;

  /// Occupancy-resolved mode: whether the occupancies of 'cageA'/'cageB' are the requested ones.
  [[nodiscard]] bool occupanciesSatisfied(const DcTSTOccupancyState& state, std::size_t cageA, std::size_t cageB) const;

  /// Occupancy-resolved mode: rebuilds the cage bookkeeping from the current configuration and moves
  /// molecules in or out of cages A and B until the requested occupancies hold.
  void prepareOccupancies(std::size_t replicaId);

  /// Occupancy-resolved mode: a (per cage) occupancy count of all molecules of the system.
  [[nodiscard]] std::vector<std::size_t> countOccupancies(const System& system) const;

  /// Launches 'TrajectoriesPerConfiguration' forward/time-reversed trajectory pairs from the
  /// current configuration and accumulates the transmission functions. 'productionCycle' is the cycle of the
  /// transmission stage (recorded with the pair). A pair interrupted by a termination signal is discarded.
  void launchTrajectories(std::size_t replicaId, std::size_t block, std::size_t productionCycle);

  /// Writes the header of the per-pair record file of a replica ('WritePairRecords').
  void writePairRecordHeader(std::size_t replicaId);

  /// Writes the accumulators of a replica to 'dctst/accumulator_{T}_{P}.r{k}.json' (written to a temporary file
  /// and renamed, so the file is always complete). 'final' marks the snapshot at the end of the transmission stage.
  void writeAccumulatorSnapshot(std::size_t replicaId, bool final) const;

  /// The stored time indices inside the plateau window.
  [[nodiscard]] std::vector<std::size_t> computePlateauIndices() const;

  /// Global index (into moleculeData) of the tagged molecule.
  [[nodiscard]] std::size_t taggedMoleculeIndex(const System& system) const;

  // --- analysis ----------------------------------------------------------------------------------

  struct TransmissionResult
  {
    std::vector<double> time;
    std::vector<double> self;
    std::vector<double> selfError;
    std::vector<double> selfTwoState;
    std::vector<double> selfCommitted;
    std::vector<double> selfSmooth;
    std::vector<double> selfSmoothError;
    double selfTwoStatePlateau{0.0};
    double selfCommittedPlateau{0.0};
    double selfCommittedPlateauError{0.0};
    bool smoothAvailable{false};
    double selfSmoothPlateau{0.0};
    double selfSmoothPlateauError{0.0};
    double longJumpFraction{0.0};   ///< Fraction of trajectories in which the tagged molecule left cages A/B.
    double committedFraction{0.0};  ///< Fraction of trajectories in which the tagged molecule got committed.
    double plateauBegin{0.0};
    double plateauEnd{0.0};
    std::vector<std::vector<double>> collective;  ///< [shell][time]
    std::vector<std::vector<double>> collectiveError;
    std::vector<std::vector<std::vector<double>>> collectivePerComponent;  ///< [shell][component][time] (kappa_X part)
    double meanSpeed{0.0};
    double meanSpeedError{0.0};
    double selfPlateau{0.0};
    double selfPlateauError{0.0};
    std::vector<double> collectivePlateau;
    std::vector<double> collectivePlateauError;
    std::size_t numberOfTrajectoryPairs{0};
    double meanRelativeEnergyDrift{0.0};
    double maximumRelativeEnergyDrift{0.0};

    // event-terminated estimator
    double selfEvent{0.0};  ///< kappa_S^event
    double selfEventError{0.0};
    double collectiveEvent{0.0};  ///< kappa_C^event
    double collectiveEventError{0.0};
    double backCorrelationFactor{0.0};  ///< f = kappa_S^event / kappa_S^committed
    double backCorrelationFactorError{0.0};
    double collectiveExcess{0.0};  ///< epsilon = kappa_C^event / kappa_S^committed - 1
    double collectiveExcessError{0.0};
    /// Restricted to the event-terminated trajectories (unbiased by the cap at first-event order).
    double selfEventConditioned{0.0};
    double selfEventConditionedError{0.0};
    double collectiveEventConditioned{0.0};
    double collectiveEventConditionedError{0.0};
    double backCorrelationFactorConditioned{0.0};
    double collectiveExcessConditioned{0.0};
    /// Cap-corrected: the deviations from kappa_S^committed divided by the fraction of trajectories that
    /// ended on an event (a trajectory without an event carries no induced jump and no back-hop).
    double cappedFraction{0.0};  ///< (capped + never committed) / all trajectories
    double backCorrelationFactorCapCorrected{0.0};
    double collectiveExcessCapCorrected{0.0};
    /// Decomposition epsilon = (f - 1) + epsilon_before + epsilon_after (Q_local crossings before and after the
    /// commitment of the tagged molecule), and per channel and phase the contribution to epsilon, the number of
    /// crossings and their mean time (since launch before commitment, since commitment after it) [ps].
    double collectiveExcessBeforeCommitment{0.0};
    double collectiveExcessBeforeCommitmentError{0.0};
    double collectiveExcessAfterCommitment{0.0};
    double collectiveExcessAfterCommitmentError{0.0};
    std::array<std::array<double, static_cast<std::size_t>(DcTSTTransmissionAccumulator::EventChannel::Count)>, 2>
        collectiveExcessChannel{};
    std::array<std::array<double, static_cast<std::size_t>(DcTSTTransmissionAccumulator::EventChannel::Count)>, 2>
        channelMeanTime{};
    std::array<std::array<std::size_t, static_cast<std::size_t>(DcTSTTransmissionAccumulator::EventChannel::Count)>, 2>
        channelEvents{};
    std::array<std::array<std::size_t, static_cast<std::size_t>(DcTSTTransmissionAccumulator::EventChannel::Count)>, 2>
        terminatingEvents{};
    std::array<std::array<std::size_t, static_cast<std::size_t>(DcTSTTransmissionAccumulator::EventChannel::Count)>, 2>
        allEvents{};
    std::size_t numberOfEventTerminatedTrajectories{0};
    std::size_t numberOfCappedTrajectories{0};
    std::size_t numberOfUncommittedTrajectories{0};
    double meanTerminationTime{0.0};  ///< [ps], over the event-terminated trajectories
    double maximumTerminationTime{0.0};
    double meanIntegratedTime{0.0};  ///< [ps], over all trajectories

    /// Hybrid estimator at one handover time tau: epsilon_hyb(tau) = epsilon_conc(tau) + epsilon_Markov(tau), all
    /// relative to kappa_S^committed like epsilon. 'tagged' and 'others' split epsilon_conc into the committor-smoothed
    /// tagged part (the smooth analogue of f - 1) and the smoothed displacement of the other molecules.
    struct Hybrid
    {
      double time{0.0};  ///< tau [ps]
      double concerted{0.0};
      double concertedError{0.0};
      double concertedTagged{0.0};
      double concertedOthers{0.0};
      double fractionEventBefore{0.0};  ///< Members whose first event occurred before tau (they enter with their event value).
      double fractionTruncated{0.0};    ///< Members that did not reach tau (trajectory capped or never committed).
      std::size_t numberOfStates{0};    ///< Distinct occupancy states of the column at tau.
      /// Markov closure over the occupancy histogram at tau with the measured/fitted rates: lateral channels of
      /// A and B ending the race with s = 0 (the first-order closure of the article), and axial channels only.
      bool markovAvailable{false};
      double markovTable{0.0};      ///< Measured rates where available, fit elsewhere.
      double markovTableAxial{0.0};
      double markovFit{0.0};        ///< Fit for every pair.
      double markovFitAxial{0.0};
      double unmeasuredWeight{0.0};  ///< Weight fraction of the states that needed the fit for at least one channel.
      double hybrid{0.0};           ///< concerted + markovTable
      double hybridFit{0.0};        ///< concerted + markovFit
      /// Occupancy histogram at tau (weight fraction of |qdot|, frame of the tagged hop, tagged molecule counted in B).
      std::vector<std::pair<DcTSTTransmissionAccumulator::ColumnState, double>> occupancy;
      std::vector<std::size_t> occupancyCounts;
      std::size_t numberOfNeighbourhoodStates{0};
    };
    std::vector<Hybrid> hybrid;

    /// Occupancy classes of the trajectory members (post-stratification), sorted by weight.
    struct OccupancyClass
    {
      std::size_t nFrom{0};
      std::size_t nTo{0};
      std::size_t members{0};
      double weightFraction{0.0};  ///< Share of sum |qdot| over all members.
      double kappa{0.0};           ///< Net cage count, plateau.
      double kappaError{0.0};
      double kappaSmooth{0.0};
      double kappaSmoothError{0.0};
      double kappaCommitted{0.0};
      double pilotWeight{0.0};  ///< P(n_from, n_to) of the free-energy input (0 when absent).
      double tstRatio{0.0};     ///< k^TST(n_from -> n_to) / k^TST (0 when absent).
    };
    std::vector<OccupancyClass> classes;
    double classesReweightedKappa{0.0};  ///< sum_classes P_FE(n,n') kappa(n,n') tstRatio(n,n') over the measured classes.
    double classesReweightedCoverage{0.0};  ///< sum of P_FE over the measured classes.

    /// The rate table of the Markov closure (empty when no k^TST is available and no conditioned runs are given).
    std::optional<OccupancyRates> rates;
  };

  [[nodiscard]] TransmissionResult analyzeTransmission() const;
};
