module;

export module polymerization;

import std;

import archive;
import json;
import input_reader;
import molecular_dynamics_spatial_decomposition;
import polymerization_settings;

/**
 * \brief Stochastic free-radical chain-growth polymerization with the spatial-decomposition MD engine.
 *
 * Selected with 'SimulationType' : 'Polymerization'. The driver is the MolecularDynamicsSpatialDecomposition
 * driver with a reaction step: the Monte Carlo (pre-)initialization and the MD equilibration run on the
 * unreacted monomers; the production stage is the curing stage, in which every 'ScanEvery' MD steps the
 * reactive sites ('PolymerizationSites' of the molecule files) are scanned and reactions are accepted with
 * their probabilities (see PolymerizationSettings):
 *
 *   initiation   M -> R*                  (monomer site becomes a radical)
 *   propagation  R* + M -> R-R*           (bond radical tail - monomer head; the monomer is the new radical)
 *   termination  R* + R* -> R-R           (bond tail - tail, combination)
 *
 * A reaction retypes the head and the tail of the reacted site to the product pseudo-atoms and adds the
 * bond with its bends, torsions and non-bonded exclusions as an overlay on the component topology (see
 * PolymerLinkTable); the new terms are ramped in over 'RelaxationSteps' steps. The molecules of the system
 * never change: the polymer is the graph of bonds between them, written at the end to
 * 'output/polymer_topology.s<i>.json' (and every reaction to 'output/polymerization_bonds.s<i>.txt').
 *
 * Event frames ('WriteEventFrames'): the configuration at the moment of every reaction scan that accepted a
 * reaction is written as a LAMMPS data file to 'polymerization_events/s<i>/step_<step>.data' (with the companion
 * input script), with the polymer bonds and their bends and torsions merged into its topology and the reacted
 * atoms carrying their product types; 'EventFrameEvery' thins the data files. The reactions of every frame (the
 * atoms with their types before and after, the distance, the configuration, the new bonded terms) are appended
 * as one JSON object per line to 'output/polymerization_events.s<i>.jsonl', and the topology file is refreshed.
 *
 * Radical trapping is measured with the PropertyRadicalWaitingTimes of the system ('ComputeRadicalWaitingTimes'):
 * the driver adds the waiting time of a radical to it at the reaction (the only moment the closed interval is
 * known) and writes the property every 'WriteRadicalWaitingTimesEvery' cycles and at the end.
 *
 * Tacticity ('Tacticity' of the settings): the tail of a reacted site is a stereocentre whose configuration
 * is the sign of its chiral volume (see ChiralVolume). When a radical adds to a monomer the configuration of
 * its tail is fixed by the face from which the monomer approaches; the configuration a radical has to take is
 * set when it is created, relative to the centre before it in the chain (isotactic: the same sign, a meso
 * dyad; syndiotactic: the opposite sign, a racemo dyad; atactic: meso with 'MesoProbability'). Candidates
 * with the wrong approach face are skipped ('StereoSelective') or accepted and driven into the required
 * configuration by the restraint; every reacted centre is held by the restraint afterwards.
 *
 * Engine: with the device bonded kernels the overlay is evaluated on the device as well (DevicePolymerOverlay,
 * rebuilt after every reaction, its ramp factors written per step), so the resident (fully on-device) integrator
 * can be used; the scans then fetch the device state first and a reaction re-evaluates the forces on the host,
 * after which the state is uploaded again. Without the device bonded kernels the overlay is evaluated on the host
 * after the step. 'ReactiveSites' / 'CrossLinkBonds' (the equilibrium cross-link moves) are a separate mechanism
 * and are not combined with it.
 */
export struct Polymerization
{
  Polymerization(InputReader& reader);

  std::uint64_t versionNumber{5};

  MolecularDynamicsSpatialDecomposition md;
  PolymerizationSettings settings{};

  enum class ReactionKind : std::uint8_t
  {
    Initiation,
    Propagation,
    Termination
  };
  /// One accepted reaction of a scan (the bond log and the event frames are written from these records).
  struct ReactionEvent
  {
    ReactionKind kind{ReactionKind::Initiation};
    std::size_t siteA{};                              ///< The initiated site, or the attacking radical.
    std::size_t siteB{std::numeric_limits<std::size_t>::max()};  ///< The attacked monomer / the other radical.
    std::uint32_t atomA{};                            ///< Global atom of the bond (the tail of 'siteA').
    std::uint32_t atomB{std::numeric_limits<std::uint32_t>::max()};  ///< Its partner (none for an initiation).
    double distance{0.0};                             ///< [Angstrom] at the reaction.
    std::array<std::int8_t, 2> chirality{0, 0};       ///< Configuration taken by the tail(s) of a / b.
    std::string stereo{"-"};                          ///< The log's stereo field.
    std::optional<double> strain{};                   ///< Angular strain [K] with an 'OrientationTemperature'.
    std::array<std::optional<double>, 2> waitingTime{};  ///< [ps] of the reacting radical(s).
    std::size_t bond{std::numeric_limits<std::size_t>::max()};  ///< Index of the new bond in the table.
    /// The atoms whose pseudo-atom type changed: global atom, type before and after.
    std::vector<std::tuple<std::uint32_t, std::uint16_t, std::uint16_t>> retyped{};

    static constexpr std::string_view name(ReactionKind kind)
    {
      switch (kind)
      {
        case ReactionKind::Initiation:
          return "initiation";
        case ReactionKind::Propagation:
          return "propagation";
        case ReactionKind::Termination:
          return "termination";
      }
      return "unknown";
    }
  };

  /// Per system: whether the reactions have stopped (target conversion reached, or nothing left to react).
  std::vector<std::uint8_t> reactionsFinished{};
  /// Per system: the number of reaction scans performed.
  std::vector<std::size_t> numberOfScans{};
  /// Per system: the production cycle at which the reactions were activated.
  std::vector<std::size_t> activationCycle{};
  /// Per system: the candidates skipped because the monomer approached from the wrong face.
  std::vector<std::size_t> stereoRejections{};
  /// Per system: the candidates rejected by the orientation weight (with an 'OrientationTemperature').
  std::vector<std::size_t> orientationRejections{};
  /// Per system: the angular strain [K] of the accepted reactions, summed, and their number (the mean strain).
  std::vector<double> acceptedStrain{};
  std::vector<std::size_t> acceptedStrainCount{};
  /// Per system: the step at which the reactions stopped (0: still running); the ages of the living radicals of
  /// the waiting-time property are measured up to this step.
  std::vector<std::uint64_t> reactionsStoppedStep{};
  /// Per system: the number of reaction frames so far (the stride of the event data files counts these).
  std::vector<std::size_t> eventFrames{};

  std::vector<std::ofstream> bondLogs{};   ///< Not part of the restart (reopened in append mode).
  std::vector<std::ofstream> eventLogs{};  ///< The JSON-lines records of the event frames; reopened in append mode.

  void run();
  /// Reopens the output files of a restarted run (the MD driver's streams and the logs, in append mode).
  void createOutputFiles();
  /// Opens the bond logs and the event logs (append mode on a restart).
  void openBondLogs();

  /// The hook of the MD driver: called after every MD step of a system.
  void afterStep(std::size_t systemId);
  /// The status block of the MD driver's reports.
  std::string status(std::size_t systemId);
  /// The orientation-weight line of the status (empty without an 'OrientationTemperature').
  std::string orientationText(std::size_t systemId) const;
  /// The step up to which the ages of the living radicals are measured (the stop of the reactions, or now).
  std::uint64_t waitingTimeReferenceStep(std::size_t systemId) const;
  /// The waiting-time lines of the status (empty without 'ComputeRadicalWaitingTimes').
  std::string waitingTimeText(std::size_t systemId) const;
  /// Writes the waiting-time property of a system ('WriteRadicalWaitingTimesEvery', or forced at the end).
  void writeWaitingTimes(std::size_t systemId, bool force) const;

  /// Activates the table of a system (site layout, initial radicals) at the start of the reactive stage.
  void activate(std::size_t systemId);
  /// The configuration the next radical of a chain has to take, given the configuration of the centre before it.
  std::int8_t nextChirality(std::int8_t previous);
  /// One reaction scan of a system; returns the accepted reactions (logged to the bond log).
  std::vector<ReactionEvent> scan(std::size_t systemId);
  /// Appends an event to the bond log of a system.
  void logEvent(std::size_t systemId, const ReactionEvent& event);
  /// After accepted reactions: tells the engine, refreshes the forces and the energies of the system.
  void topologyChanged(std::size_t systemId);
  /// The event frame of the reactions just accepted (after topologyChanged(): the LAMMPS data file of the
  /// configuration with the polymer topology, the JSON record, the refreshed topology file).
  void writeEventFrame(std::size_t systemId, std::span<const ReactionEvent> events);
  /// The JSON record of an event frame ('dataFile': the data file written for it, empty when thinned out).
  nlohmann::json eventRecord(std::size_t systemId, std::span<const ReactionEvent> events,
                             std::string_view dataFile) const;
  void writeSettings();
  /// Writes 'output/polymer_topology.s<i>.json' of a system (the current site states and bonds).
  void writeTopology(std::size_t systemId);

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const Polymerization& p);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, Polymerization& p);
};
