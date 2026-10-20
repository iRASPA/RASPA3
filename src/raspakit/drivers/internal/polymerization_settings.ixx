module;

export module polymerization_settings;

import std;

import archive;

/**
 * \brief Settings of the 'Polymerization' driver: the stochastic free-radical chain-growth model.
 *
 * Read from the 'Polymerization' object of the system block. The driver runs the spatial-decomposition
 * MD engine and, every 'scanEvery' steps of the production (curing) stage, scans the reactive sites of the
 * components (declared with 'PolymerizationSites' in the molecule files):
 *
 *  - initiation:  a monomer site becomes a radical ('InitialRadicals' sites at the start of the reactive
 *                 stage, plus continuous initiation with probability 'initiationProbability' per monomer
 *                 site and scan);
 *  - propagation: a radical tail within 'propagationCaptureRadius' of a monomer head adds to it with
 *                 probability 'propagationProbability' (or 1 - exp(-k dt_scan) from a rate constant);
 *                 the monomer becomes the new radical;
 *  - termination: two radical tails within 'terminationCaptureRadius' combine with probability
 *                 'terminationProbability'.
 *
 * A new bond is ramped in over 'relaxationSteps' MD steps (see PolymerLinkTable). The reactions stop when
 * the conversion reaches 'targetConversion' (the remaining production cycles relax the network).
 *
 * The 'Tacticity' object controls the configuration of the backbone stereocentres (the tails of the reacted
 * sites, see PolymerizationSite): the relative configuration of two consecutive centres of a chain is a meso
 * dyad (isotactic placement) or a racemo dyad (syndiotactic placement). With a tacticity selected, every
 * radical carries the configuration its tail has to take when it adds to the next monomer; a candidate
 * whose approach geometry gives the other configuration is skipped ('StereoSelective', the default) or
 * accepted and driven to the required configuration by the restraint. Every reacted centre is held in its
 * configuration by a harmonic restraint on the normalised chiral volume (see PolymerLinkTable), since the
 * united-atom CH bead has no hydrogen that would prevent the inversion of the centre.
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
  std::uint64_t versionNumber{2};

  std::size_t scanEvery{10};  ///< MD steps between two reaction scans ('ScanEvery').
  std::size_t startCycle{0};  ///< Production cycle at which the reactions start ('StartCycle').
  double targetConversion{1.0};  ///< Conversion at which the reactions stop ('TargetConversion').
  std::size_t relaxationSteps{1000};  ///< Force ramp of a new bond in MD steps ('RelaxationSteps').
  std::size_t maximumReactionsPerScan{std::numeric_limits<std::size_t>::max()};  ///< ('MaximumReactionsPerScan').

  std::size_t initialRadicals{0};      ///< Sites initiated at the start of the reactive stage ('InitialRadicals').
  double initiationProbability{0.0};   ///< Per monomer site and scan ('Initiation' / 'Probability').
  std::optional<double> initiationRateConstant{};  ///< [1/ps], alternative to the probability ('RateConstant').

  double propagationCaptureRadius{4.0};  ///< [Angstrom] ('Propagation' / 'CaptureRadius').
  double propagationProbability{1.0};    ///< Per candidate pair and scan ('Propagation' / 'Probability').
  std::optional<double> propagationRateConstant{};  ///< [1/ps] ('Propagation' / 'RateConstant').
  bool allowIntramolecular{true};  ///< Whether a radical may add to a site of its own molecule (primary cycles).

  bool terminationEnabled{true};         ///< ('Termination' / 'Enabled').
  double terminationCaptureRadius{4.0};  ///< [Angstrom] ('Termination' / 'CaptureRadius').
  double terminationProbability{1.0};    ///< Per candidate pair and scan ('Termination' / 'Probability').
  std::optional<double> terminationRateConstant{};  ///< [1/ps] ('Termination' / 'RateConstant').

  bool writeBondLog{true};        ///< Write every reaction to 'output/polymerization_bonds.s<i>.txt' ('WriteBondLog').
  bool writeTopology{true};       ///< Write the final bond list to 'output/polymer_topology.s<i>.json' ('WriteTopology').
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

  /// The acceptance probability of a reaction per scan: the probability, or 1 - exp(-k dt) from the rate constant.
  static double probabilityPerScan(double probability, std::optional<double> rateConstant, double scanTime)
  {
    if (rateConstant.has_value()) return 1.0 - std::exp(-rateConstant.value() * scanTime);
    return probability;
  }

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PolymerizationSettings &s)
  {
    archive << s.versionNumber;
    archive << s.scanEvery;
    archive << s.startCycle;
    archive << s.targetConversion;
    archive << s.relaxationSteps;
    archive << s.maximumReactionsPerScan;
    archive << s.initialRadicals;
    archive << s.initiationProbability;
    archive << s.initiationRateConstant;
    archive << s.propagationCaptureRadius;
    archive << s.propagationProbability;
    archive << s.propagationRateConstant;
    archive << s.allowIntramolecular;
    archive << s.terminationEnabled;
    archive << s.terminationCaptureRadius;
    archive << s.terminationProbability;
    archive << s.terminationRateConstant;
    archive << s.writeBondLog;
    archive << s.writeTopology;
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
    archive >> s.initialRadicals;
    archive >> s.initiationProbability;
    archive >> s.initiationRateConstant;
    archive >> s.propagationCaptureRadius;
    archive >> s.propagationProbability;
    archive >> s.propagationRateConstant;
    archive >> s.allowIntramolecular;
    archive >> s.terminationEnabled;
    archive >> s.terminationCaptureRadius;
    archive >> s.terminationProbability;
    archive >> s.terminationRateConstant;
    archive >> s.writeBondLog;
    archive >> s.writeTopology;
    archive >> s.templateFile;
    if (versionNumber >= 2)
    {
      archive >> s.tacticity;
      archive >> s.mesoProbability;
      archive >> s.stereoSelective;
      archive >> s.chiralForceConstant;
      archive >> s.chiralReferenceVolume;
    }
    return archive;
  }
};
