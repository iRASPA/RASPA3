module;

export module polymerization_settings;

import std;

import archive;

/**
 * \brief Settings of the 'Polymerization' driver: the stochastic reactive molecular dynamics with user-defined
 * reaction rules.
 *
 * Read from the 'Polymerization' object of the system block. The driver runs the spatial-decomposition MD engine
 * and, every 'scanEvery' steps of the production (curing) stage, scans the reactive groups of the components
 * (declared with 'ReactiveGroups' in the molecule files) for the reactions of the 'Reactions' list (see
 * ReactionRule and ReactionScheme): a bimolecular rule for the pairs of groups in its reactant states whose
 * capture atoms are within its capture radius, a unimolecular rule for the groups in its reactant state. The
 * rules carry their own probabilities, capture radii and orientation temperatures; the settings here are the
 * global ones.
 *
 * With an 'OrientationTemperature' T* (Kelvin) in a rule the acceptance is weighted by the geometry of the
 * approach: p = p0 exp(-dU / k_B T*), where dU is the angular strain the new bonds would carry (the bends and
 * torsions through them, evaluated with the product templates in the current geometry, relative to their minima;
 * see PolymerLinkTable::angularStrain). A monomer that approaches in the plane of its own double bond, or along
 * the sigma framework of the radical, is thereby disfavoured in proportion to the strain the force field would
 * have to relax; T* -> infinity recovers the distance-only criterion.
 *
 * A new bond is ramped in over 'relaxationSteps' MD steps (see PolymerLinkTable). The reactions stop when the
 * conversion reaches 'targetConversion' (the remaining production cycles relax the network).
 *
 * Radical trapping is measured, not imposed: for every group in an active state (a living chain end) the waiting
 * time between its successive reactions, dt = t(next reaction) - t(last reaction), is collected in a log-spaced
 * histogram per reaction, together with the open intervals from the last reaction of every still-active group
 * to the current time (its age) ('ComputeRadicalWaitingTimes' of the system block, see
 * PropertyRadicalWaitingTimes).
 *
 * The 'Tacticity' object controls the configuration of the stereocentres the rules create ('Stereo' of a rule,
 * see ReactionRule): the relative configuration of two consecutive centres of a chain is a meso dyad (isotactic
 * placement) or a racemo dyad (syndiotactic placement). With a tacticity selected, every group named as 'Next'
 * by a centre carries the configuration its own centre has to take; a candidate whose approach geometry gives
 * the other configuration is skipped ('StereoSelective', the default) or accepted and driven to the required
 * configuration by the restraint. Every created centre is held in its configuration by a harmonic restraint on
 * the normalised chiral volume (see PolymerLinkTable), since the united-atom CH bead has no hydrogen that would
 * prevent the inversion of the centre.
 */
export enum class Tacticity : std::uint8_t
{
  None = 0,          ///< No control (and no restraints); the dyad statistics are still reported.
  Isotactic = 1,     ///< All dyads meso.
  Syndiotactic = 2,  ///< All dyads racemo.
  Atactic = 3        ///< Bernoullian: a dyad is meso with probability 'mesoProbability'.
};

export struct PolymerizationSettings
{
  std::uint64_t versionNumber{1};

  std::size_t scanEvery{10};  ///< MD steps between two reaction scans ('ScanEvery').
  std::size_t startCycle{0};  ///< Production cycle at which the reactions start ('StartCycle').
  double targetConversion{1.0};  ///< Conversion at which the reactions stop ('TargetConversion').
  std::size_t relaxationSteps{1000};  ///< Force ramp of a new bond in MD steps ('RelaxationSteps').
  std::size_t maximumReactionsPerScan{std::numeric_limits<std::size_t>::max()};  ///< ('MaximumReactionsPerScan').

  bool writeBondLog{true};        ///< Write every reaction to 'output/polymerization_bonds.s<i>.txt' ('WriteBondLog').
  bool writeTopology{true};       ///< Write the final bond list to 'output/polymer_topology.s<i>.json' ('WriteTopology').
  /// Write the configuration at every reaction as a LAMMPS data file (with the polymer bonds in its topology) to
  /// 'polymerization_events/s<i>/', with a record of the reactions per frame in
  /// 'output/polymerization_events.s<i>.jsonl' ('WriteEventFrames').
  bool writeEventFrames{false};
  std::size_t eventFrameEvery{1};  ///< Write the data file of every n-th reaction frame only ('EventFrameEvery').
  std::string templateFile{};     ///< Optional molecule file whose type-based bonded tables are added ('TemplateFile').

  Tacticity tacticity{Tacticity::None};  ///< ('Tacticity' / 'Type').
  double mesoProbability{0.5};           ///< Probability of a meso dyad for 'Atactic' ('MesoProbability').
  bool stereoSelective{true};            ///< Skip candidates that approach from the wrong face ('StereoSelective').
  double chiralForceConstant{40000.0};   ///< [K] of the restraint on the normalised chiral volume ('ForceConstant').
  double chiralReferenceVolume{0.6};     ///< [-] magnitude of the normalised chiral volume held ('ReferenceVolume').

  static std::string tacticityName(Tacticity tacticity)
  {
    switch (tacticity)
    {
      case Tacticity::None:
        return "none";
      case Tacticity::Isotactic:
        return "isotactic";
      case Tacticity::Syndiotactic:
        return "syndiotactic";
      case Tacticity::Atactic:
        return "atactic";
    }
    return "none";
  }

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PolymerizationSettings &s)
  {
    archive << s.versionNumber;
    archive << s.scanEvery;
    archive << s.startCycle;
    archive << s.targetConversion;
    archive << s.relaxationSteps;
    archive << s.maximumReactionsPerScan;
    archive << s.writeBondLog;
    archive << s.writeTopology;
    archive << s.writeEventFrames;
    archive << s.eventFrameEvery;
    archive << s.templateFile;
    archive << s.tacticity;
    archive << s.mesoProbability;
    archive << s.stereoSelective;
    archive << s.chiralForceConstant;
    archive << s.chiralReferenceVolume;
    return archive;
  }
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PolymerizationSettings &s)
  {
    std::uint64_t versionNumber;
    archive >> versionNumber;
    if (versionNumber > s.versionNumber)
    {
      throw std::runtime_error(std::format("Invalid version reading 'PolymerizationSettings' at line {} in file {}\n",
                                           std::source_location::current().line(),
                                           std::source_location::current().file_name()));
    }
    archive >> s.scanEvery;
    archive >> s.startCycle;
    archive >> s.targetConversion;
    archive >> s.relaxationSteps;
    archive >> s.maximumReactionsPerScan;
    archive >> s.writeBondLog;
    archive >> s.writeTopology;
    archive >> s.writeEventFrames;
    archive >> s.eventFrameEvery;
    archive >> s.templateFile;
    archive >> s.tacticity;
    archive >> s.mesoProbability;
    archive >> s.stereoSelective;
    archive >> s.chiralForceConstant;
    archive >> s.chiralReferenceVolume;
    return archive;
  }
};
