module;

export module polymerization;

import std;

import archive;
import json;
import input_reader;
import molecular_dynamics_spatial_decomposition;
import polymerization_settings;
import polymer_reactions;
import polymer_links;

/**
 * \brief Stochastic reactive molecular dynamics with user-defined reaction rules (polymerization) on the
 * spatial-decomposition MD engine.
 *
 * Selected with 'SimulationType' : 'Polymerization'. The driver is the MolecularDynamicsSpatialDecomposition
 * driver with a reaction step: the Monte Carlo (pre-)initialization and the MD equilibration run on the
 * unreacted monomers; the production stage is the curing stage, in which every 'ScanEvery' MD steps the
 * reactive groups ('ReactiveGroups' of the molecule files) are scanned for the reactions of the 'Reactions' list
 * of the 'Polymerization' block (see ReactionRule, ReactionScheme and PolymerizationSettings). For example, the
 * free-radical chain growth of an acrylate:
 *
 *   initiation   M[monomer] -> M[radical]                      (unimolecular; 'InitialEvents' at the start)
 *   propagation  R[radical] + M[monomer] -> R[dead] + M[radical]  (bond R.tail - M.head)
 *   termination  A[radical] + B[radical] -> A[dead] + B[dead]     (bond A.tail - B.tail)
 *
 * A reaction retypes the atoms of the groups as the rule prescribes, changes their states and adds the bonds
 * with their bends, torsions and non-bonded exclusions as an overlay on the component topology (see
 * PolymerLinkTable); the new terms are ramped in over 'RelaxationSteps' steps. The molecules of the system never
 * change: the polymer is the graph of bonds between them, written at the end to
 * 'output/polymer_topology.s<i>.json' (and every reaction to 'output/polymerization_bonds.s<i>.txt').
 *
 * Every living chain end carries a radical id (see RadicalRecord): created at an initiation, handed on by every
 * propagation, closed by a termination. The bond log, the event records and the topology file name the radical
 * of every reaction, and the topology file holds the table of the radicals (birth, death, propagations), from
 * which the reactions of one radical can be followed.
 *
 * Event frames ('WriteEventFrames'): the configuration at the moment of every reaction scan that accepted a
 * reaction is written as a LAMMPS data file to 'polymerization_events/s<i>/step_<step>.data' (with the companion
 * input script), with the polymer bonds and their bends and torsions merged into its topology and the reacted
 * atoms carrying their product types; 'EventFrameEvery' thins the data files. The reactions of every frame (the
 * groups and atoms with their types before and after, the distance, the configuration, the new bonded terms) are
 * appended as one JSON object per line to 'output/polymerization_events.s<i>.jsonl', and the topology file is
 * refreshed.
 *
 * Radical trapping is measured with the PropertyRadicalWaitingTimes of the system ('ComputeRadicalWaitingTimes'):
 * the driver adds the waiting time of a reacting active group to it at the reaction (the only moment the closed
 * interval is known) and writes the property every 'WriteRadicalWaitingTimesEvery' cycles and at the end.
 *
 * Tacticity ('Tacticity' of the settings): a stereo entry of a rule creates a stereocentre whose configuration
 * is the sign of its chiral volume (see ChiralVolume) in the geometry of the reaction. The configuration a
 * group's centre has to take is set when the group is named 'Next' by a centre, relative to that centre
 * (isotactic: the same sign, a meso dyad; syndiotactic: the opposite sign, a racemo dyad; atactic: meso with
 * 'MesoProbability'). Candidates with the wrong approach face are skipped ('StereoSelective') or accepted and
 * driven into the required configuration by the restraint; every created centre is held by the restraint
 * afterwards.
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

  std::uint64_t versionNumber{1};

  MolecularDynamicsSpatialDecomposition md;
  PolymerizationSettings settings{};

  static constexpr std::size_t noGroup = std::numeric_limits<std::size_t>::max();
  static constexpr std::uint32_t noAtom = std::numeric_limits<std::uint32_t>::max();

  /// One accepted reaction of a scan (the bond log and the event frames are written from these records).
  struct ReactionEvent
  {
    std::uint16_t rule{0};                             ///< The reaction rule (index into the scheme).
    std::array<std::size_t, 2> groups{noGroup, noGroup};  ///< The reactant groups, in the order of the rule.
    std::array<std::uint32_t, 2> atoms{noAtom, noAtom};   ///< The capture atoms (bimolecular).
    double distance{0.0};                              ///< [Angstrom] of the capture atoms at the reaction.
    std::vector<std::int8_t> chirality{};              ///< Configuration taken by every stereocentre of the rule.
    std::string stereo{"-"};                           ///< The log's stereo field.
    std::optional<double> strain{};                    ///< Angular strain [K] with an 'OrientationTemperature'.
    std::array<std::optional<double>, 2> waitingTime{};  ///< [ps] of the reacting active group(s).
    /// Per reactant: the radical (see RadicalRecord) on it before and after the reaction (noRadical: none), and
    /// whether that radical was born or died in this reaction.
    std::array<std::size_t, 2> radicalBefore{PolymerLinkTable::noRadical, PolymerLinkTable::noRadical};
    std::array<std::size_t, 2> radicalAfter{PolymerLinkTable::noRadical, PolymerLinkTable::noRadical};
    std::array<bool, 2> radicalBorn{false, false};
    std::array<bool, 2> radicalDied{false, false};
    std::vector<std::size_t> bonds{};                  ///< Indices of the new bonds in the table.
    std::vector<std::size_t> brokenBonds{};            ///< Indices of the bonds broken by the event in the table.
    /// The atoms whose pseudo-atom type changed: global atom, type before and after.
    std::vector<std::tuple<std::uint32_t, std::uint16_t, std::uint16_t>> retyped{};

    bool bimolecular() const { return groups[1] != noGroup; }
  };

  /// Per system: whether the reactions have stopped (target conversion reached, or nothing left to react).
  std::vector<std::uint8_t> reactionsFinished{};
  /// Per system: the number of reaction scans performed.
  std::vector<std::size_t> numberOfScans{};
  /// Per system: the production cycle at which the reactions were activated.
  std::vector<std::size_t> activationCycle{};
  /// Per system: the candidates skipped because the approach gave the wrong configuration.
  std::vector<std::size_t> stereoRejections{};
  /// Per system: the candidates rejected by the orientation weight (with an 'OrientationTemperature').
  std::vector<std::size_t> orientationRejections{};
  /// Per system: the angular strain [K] of the accepted reactions, summed, and their number (the mean strain).
  std::vector<double> acceptedStrain{};
  std::vector<std::size_t> acceptedStrainCount{};
  /// Per system: the step at which the reactions stopped (0: still running); the ages of the active groups of
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
  /// The step up to which the ages of the active groups are measured (the stop of the reactions, or now).
  std::uint64_t waitingTimeReferenceStep(std::size_t systemId) const;
  /// The waiting-time lines of the status (empty without 'ComputeRadicalWaitingTimes').
  std::string waitingTimeText(std::size_t systemId) const;
  /// Writes the waiting-time property of a system ('WriteRadicalWaitingTimesEvery', or forced at the end).
  void writeWaitingTimes(std::size_t systemId, bool force) const;

  /// Activates the table of a system (group layout, initial events of the unimolecular rules) at the start of
  /// the reactive stage.
  void activate(std::size_t systemId);
  /// The configuration the next centre of a chain has to take, given the configuration of the centre before it.
  std::int8_t nextChirality(std::int8_t previous);
  /// One reaction scan of a system; returns the accepted reactions (logged to the bond log).
  std::vector<ReactionEvent> scan(std::size_t systemId);
  /// Applies rule 'rule' to 'groups' with the resolved configurations, records and logs the event.
  ReactionEvent react(std::size_t systemId, std::size_t rule, std::span<const std::size_t> groups,
                      std::span<const std::int8_t> chiralities, std::array<std::uint32_t, 2> captureAtoms,
                      double distance, std::optional<double> strain);
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
  /// Writes 'output/polymer_topology.s<i>.json' of a system (the current group states and bonds).
  void writeTopology(std::size_t systemId);

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const Polymerization& p);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, Polymerization& p);
};
