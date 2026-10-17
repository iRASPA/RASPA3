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

  DcTSTTransmissionAccumulator() = default;
  DcTSTTransmissionAccumulator(std::size_t numberOfBlocks, std::size_t numberOfShells, std::size_t numberOfComponents,
                               std::size_t numberOfTimes);

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

  std::chrono::duration<double> totalInitializationTime{0};
  std::chrono::duration<double> totalTransmissionTime{0};
  std::chrono::duration<double> totalSimulationTime{0};

  /// P(q*) [1/A] at the dividing plane from 'FreeEnergyInput' or 'BarrierProbability', with its error and a
  /// description of its provenance; absent when neither was given.
  std::optional<double> barrierProbability{};
  double barrierProbabilityError{0.0};
  std::string barrierProbabilitySource{};

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
  /// current configuration and accumulates the transmission functions.
  void launchTrajectories(std::size_t replicaId, std::size_t block);

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
    double selfTwoStatePlateau{0.0};
    double selfCommittedPlateau{0.0};
    double selfCommittedPlateauError{0.0};
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
  };

  [[nodiscard]] TransmissionResult analyzeTransmission() const;
};
