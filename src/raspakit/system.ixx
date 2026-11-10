module;

export module system;

import std;

import ewald_structure_factor;

import archive;
import double3;
import double3x3;
import randomnumbers;
import threadpool;

import atom;
import atom_dynamics;
import molecule;
import framework;
import component;
import simulationbox;
import forcefield;
import averages;
import running_energy;
import energy_status;
import property_loading;
import sample_movies;
import property_enthalpy;
import property_partial_molar_properties;
import property_lambda_probability_histogram;
import property_simulationbox;
import property_energy;
import property_pressure;
import property_elastic_constants_fluctuation;
import property_conventional_rdf;
import property_rdf;
import property_density_grid;
import property_temperature;
import property_energy_histogram;
import property_number_of_molecules_histogram;
import property_molecule_properties;
import property_molecule_shape;
import property_molecule_backbone;
import property_msd;
import property_vacf;
import property_end_to_end_acf;
import property_number_of_molecules_evolution;
import property_volume_evolution;
import property_conserved_energy_evolution;
import units;
import bond_potential;
import move_statistics;
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_statistics;
import mc_moves_cputime;
import reaction;
import reactions;
import cross_links;
import polymer_links;
import property_radical_waiting_times;
import transition_matrix;
import equation_of_states;
import thermostat;
import thermobarostat;
import json;
import interpolation_energy_grid;
import mc_cell_list;
import spatial_decomposition_settings;
import cbmc_grow_context;
import write_lammps_data;
import minimization_cell_layout;

/**
 * \brief Represents the central system for simulations.
 *
 * The System struct holds all objects required for simulations, such as atom lists,
 * frameworks, and simulation boxes. It is passed to various simulation methods,
 * like MonteCarlo, to conduct simulations. The system can contain frameworks (zeolite/MOF)
 * or a simulation box, and it manages various simulation parameters and states.
 */
export struct System
{
  /**
   * \brief Default constructor for the System class.
   *
   * Initializes a System object with default values.
   */
  System() {};

  /**
   * \brief Constructs a System object programmatically.
   *
   * Initializes a System with the provided parameters, including system ID,
   * simulation box, temperature, pressure, force field, framework components,
   * and other simulation-related configurations.
   *
   * \param id The identifier for the system.
   * \param box The optional simulation box for the system.
   * \param T The temperature of the system.
   * \param P The optional pressure of the system.
   * \param forcefield The force field used in the simulation.
   * \param frameworkComponent The list of framework components in the system.
   * \param components The list of components in the system.
   * \param initialNumberOfMolecules The initial number of molecules per component.
   * \param numberOfBlocks The number of blocks in the simulation.
   * \param systemProbabilities The move probabilities for the Monte Carlo simulation.
   * \param sampleMoviesEvery Interval in which movies are written to PDB.
   */
  System(ForceField forcefield, std::optional<SimulationBox> box, bool hasExternalField, double T,
         std::optional<double> P, double heliumVoidFraction, std::optional<Framework> framework,
         std::vector<Component> components, std::vector<std::vector<double3>> initialPositions,
         std::vector<std::size_t> initialNumberOfMolecules, std::size_t numberOfBlocks,
         const MCMoveProbabilities& systemProbabilities = MCMoveProbabilities());

  std::uint64_t versionNumber{3};

  double temperature{300.0};
  double pressure{1e4};
  double input_pressure{1e4};
  double3 pressureTensorDiagonal{1e4, 1e4, 1e4};
  double3 input_pressureTensorDiagonal{1e4, 1e4, 1e4};
  double beta{1.0 / (Units::KB * 300.0)};
  CellMinimizationType cellMinimizationType{CellMinimizationType::Fixed};
  MonoclinicAngleType monoclinicAngleType{MonoclinicAngleType::Beta};

  double heliumVoidFraction{0.29};

  std::size_t numberOfFrameworks{0};
  std::size_t numberOfFrameworkAtoms{0};
  std::size_t numberOfRigidFrameworkAtoms{0};

  std::optional<Framework> framework;
  std::vector<Component> components;

  EquationOfState equationOfState;
  MolecularDynamicsEnsemble molecularDynamicsEnsemble{MolecularDynamicsEnsemble::NVE};
  std::optional<Thermostat> thermostat;
  std::optional<Thermobarostat> thermobarostat;

  LoadingData loadings;

  std::vector<std::size_t> swappableComponents{};
  std::vector<std::size_t> initialNumberOfMolecules{};

  // total # of molecules per component (include fractional molecules)
  std::vector<std::size_t> numberOfMoleculesPerComponent{};

  // # integer molecules
  std::vector<std::size_t> numberOfIntegerMoleculesPerComponent{};

  // # fractional molecules
  std::vector<std::size_t> numberOfFractionalMoleculesPerComponent{};

  // # fractional molecules for CFCMC grand-canonical or CFCMC Widom
  std::vector<std::size_t> numberOfGCFractionalMoleculesPerComponent_CFCMC{};

  // # fractional molecules for pair-CFCMC grand-canonical
  std::vector<std::size_t> numberOfPairGCFractionalMoleculesPerComponent_CFCMC{};

  // # fractional molecules for ion-pair insertion/deletion CFCMC (PairSwapCFCMC)
  std::vector<std::size_t> numberOfPairSwapFractionalMoleculesPerComponent_CFCMC{};

  // # fractional molecules for ion-pair insertion/deletion CB/CFCMC (PairSwapCBCFCMC)
  std::vector<std::size_t> numberOfPairSwapCBFractionalMoleculesPerComponent_CFCMC{};

  // # fractional molecules for group insertion/deletion CFCMC (GroupSwapCFCMC); a component may hold
  // more than one slot when it occurs multiple times in a group (e.g. two Cl of CaCl2)
  std::vector<std::size_t> numberOfGroupSwapFractionalMoleculesPerComponent_CFCMC{};

  // # fractional molecules for group insertion/deletion CB/CFCMC (GroupSwapCBCFCMC)
  std::vector<std::size_t> numberOfGroupSwapCBFractionalMoleculesPerComponent_CFCMC{};

  // # fractional molecules for GibbsSwapCFCMC / GibbsSwapCBCFCMC
  std::vector<std::size_t> numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC{};

  // # fractional molecules for GibbsConventionalCFCMC
  std::vector<std::size_t> numberOfGibbsFractionalMoleculesPerComponent_CFCMC{};

  // # parallel-reaction fractional molecules (ReactionConventionalCFCMC)
  std::vector<std::size_t> numberOfParallelReactionFractionalMoleculesPerComponent_CFCMC{};

  // # serial-reaction fractional molecules (ReactionCFCMC)
  std::vector<std::size_t> numberOfSerialReactionFractionalMoleculesPerComponent_CFCMC{};

  // # reactant fractional molecules for all reactions using CFCMC
  std::vector<std::vector<std::size_t>> numberOfReactionFractionalMoleculesPerComponent_CFCMC{};

  std::vector<double> idealGasEnergiesPerComponent{};

  ForceField forceField;
  bool hasExternalField{true};

  std::vector<std::vector<std::size_t>> numberOfPseudoAtoms;
  std::vector<std::size_t> totalNumberOfPseudoAtoms;

  // Brick-CFCMC-style aggregated van der Waals tail-correction accounting.
  // effectiveNumberOfPseudoAtomsVDW[type] = sum over all molecule atoms of that type of scalingVDW
  // (integer atoms contribute 1, fractional atoms contribute their lambda-scaled scalingVDW).
  std::vector<double> effectiveNumberOfPseudoAtomsVDW;
  // Per dU/dlambda group (0-based index = groupId - 1), the number of molecule atoms of each type
  // that belong to that group; used to reconstruct the tail dU/dlambda contribution.
  std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups> fractionalPseudoAtomCountsPerGroup;

  std::size_t translationalCenterOfMassConstraint{};
  std::size_t translationalDegreesOfFreedom{};
  std::size_t rotationalDegreesOfFreedom{};

  std::optional<double> frameworkMass() const;

  double timeStep{0.0005};

  SimulationBox simulationBox;
  // Reference cell h_0 of the constant-stress enthalpy sampled by the anisotropic volume move: the box the
  // system was constructed with. Only the deviatoric part of the imposed pressure tensor refers to it
  // (see 'MC_Moves::anisotropicPressureWork'); a hydrostatic pressure tensor does not depend on it.
  SimulationBox referenceSimulationBox;
  bool containsTheFractionalMolecule{true};

  // A contiguous list of adsorbate atoms per component for easy and fast looping
  // The atoms-order is defined as increasing per component and molecule.
  // Because the number of atoms is fixed per component it is easy to access the n-th molecule
  std::vector<Atom> atomData;
  // Cold per-atom molecular-dynamics fields (velocity, gradient), stored parallel to atomData with
  // identical length and indexing. Split out of Atom to keep Atom a single 64-byte cache line so the
  // energy inner loops touch only the bytes they need.
  std::vector<AtomDynamics> atomDynamics;
  // A persistent name for every atom of 'atomData', parallel to it, for clients that address atoms across the
  // insertions and deletions of a run. A tag never changes for the life of its atom: an insertion takes the lowest
  // free tags, or new ones above the highest in use; a deletion frees the tags; an exchange within the system and
  // a replica exchange move the tags with their atoms; an exchange with another system leaves them with the slot;
  // a new framework or a restart file without tags numbers the atoms afresh in storage order. The framework
  // atoms, never deleted, keep their tags for the whole run. Without deletions the tags are the storage indices.
  std::vector<std::uint32_t> atomTags;
  // The storage index of the atom holding each tag, 'noAtom' for a free tag; exactly as long as the highest tag
  // in use plus one, so it is bounded by the largest atom count the system has held.
  std::vector<std::uint32_t> atomIndexOfTag;
  static constexpr std::uint32_t noAtom = std::numeric_limits<std::uint32_t>::max();
  std::vector<Molecule> moleculeData;
  // Rigid-body state of every rigid group of every semi-flexible molecule, in molecule order
  // (each molecule contributes 'component.numberOfRigidGroups()' consecutive entries). Empty when
  // no semi-flexible molecules are present. Maintained during molecular dynamics; derived from the
  // atom positions at the molecular-dynamics/hybrid-Monte-Carlo boundary via 'initializeGroupData'.
  std::vector<GroupState> groupData;
  // Rigid-body state of every Rigid framework group (separate from molecule groupData).
  std::vector<GroupState> frameworkGroupData;
  std::vector<double> electricPotential;
  std::vector<double3> electricField;
  std::vector<double3> electricFieldNew;

  /// The conserved (extended-system) energy after the last MD step.
  double conservedEnergy{};
  /// The reference of the drift report: the conserved energy at the start of the stage, shifted by the work
  /// of the polymer force ramps (see 'accountConservedEnergy'), so that the drift is the integration error only.
  double referenceEnergy{};
  /// Sum over the steps of the stage of the relative drift |E - E_ref| / E_ref ("Average drift" per step).
  double accumulatedDrift{};
  /// The work of the polymer force ramps in the stage: the change of the conserved energy over the steps in
  /// which a ramp is active, plus the jump of the potential at each topology change (the model's heat of
  /// reaction, removed by the thermostat). Excluded from the drift.
  double polymerRampWork{};
  /// The part of 'polymerRampWork' that is the jump of the potential at the topology changes themselves (the
  /// terms a reaction changes without a ramp: the non-bonded parameters of the retyped atoms).
  double polymerReactionJumps{};
  /// Whether a polymer ramp was active in the previous step (the step that ends a ramp still changes lambda).
  bool polymerRampWasActive{false};
  RunningEnergy rigidEnergies;
  RunningEnergy runningEnergies;

  double3x3 currentExcessPressureTensor;
  EnergyStatus currentEnergyStatus;

  std::size_t numberOfHybridMCSteps{10};
  /// The hybrid MC move integrates its trajectory with the spatial-decomposition force engine (the MD engine:
  /// cell lists, particle-mesh Ewald, threads) when the engine covers the system ('HybridMCUseForceEngine'); the
  /// acceptance uses the exact energies either way. The settings are those of the MD engine ('NumberOfThreads',
  /// 'VerletSkin', 'PPPMMeshSpacing', ...).
  bool hybridMCUseForceEngine{true};
  SpatialDecompositionSettings forceEngineSettings{};
  /// The engine kept between hybrid MC moves (type-erased: the engine module depends on this one). Not archived;
  /// shared between copies of a system, which the move detects and rebuilds.
  std::shared_ptr<void> hybridMCForceEngine{};

  /// The device the multiple-try CBMC moves (SwapMultiTrialCBMC, ReinsertionMultiTrialCBMC,
  /// GibbsSwapMultiTrialCBMC) grow their trial molecules on ('CBMCDevice' in the input; CPU: the host grower). A
  /// device grows rigid molecules in single precision as a proposal and the host re-evaluates the selected one in
  /// double (see cbmc_multi_trial_device); components it does not cover are grown on the host. Not archived.
  PairDevice cbmcDevice{PairDevice::CPU};

  std::vector<std::complex<double>> eik_xy{};
  std::vector<std::complex<double>> eik_x{};
  std::vector<std::complex<double>> eik_y{};
  std::vector<std::complex<double>> eik_z{};
  EwaldStructureFactor storedEik{};
  EwaldStructureFactor fixedFrameworkStoredEik{};
  EwaldStructureFactor trialEik{};
  double CoulombicFourierEnergySingleIon{0.0};
  /// Net charge of the system as the Ewald net-charge correction sees it: the framework charge plus the
  /// Coulomb-scaled charge sum(scalingCoulomb * charge) of the molecule atoms (a fractional molecule counts with
  /// its scaled charge). Maintained incrementally by the molecule insertions and deletions and rebuilt from the
  /// atoms by 'computeTailCorrectionCounts' (every accepted scaling change calls it).
  double netCharge{0.0};
  double netChargeFramework{0.0};
  double netChargeAdsorbates{0.0};
  std::vector<double> netChargePerComponent;

  MCMoveProbabilities mc_moves_probabilities;
  MCMoveStatistics mc_moves_statistics;
  MCMoveCpuTime mc_moves_cputime;

  Reactions reactions;
  TransitionMatrix tmmc;

  /// Cross-links: inter-molecular bonds between the reactive sites of the components, with their
  /// bond types. Owned and mutated by the system (topology moves, renumbering on molecule deletion).
  CrossLinkTable crossLinks;

  /// Polymer bonds formed by the chain-growth polymerization of the 'Polymerization' driver (irreversible
  /// inter-molecular bonds with the full bonded topology of the product; see PolymerLinkTable).
  PolymerLinkTable polymerLinks;

  // property measurements
  PropertyEnergy averageEnergies;
  PropertyLoading averageLoadings;
  PropertyEnthalpy averageEnthalpiesOfAdsorption;
  PropertyPartialMolarProperties averagePartialMolarProperties;
  PropertyTemperature averageTemperature;
  PropertyTemperature averageTranslationalTemperature;
  PropertyTemperature averageRotationalTemperature;
  PropertyPressure averagePressure;
  PropertySimulationBox averageSimulationBox;
  std::optional<PropertyElasticConstantsFluctuation> propertyElasticConstantsFluctuation;
  std::size_t elasticConstantsSampleEvery{100};

  std::optional<SampleMovie> samplePDBMovie;

  std::optional<PropertyConventionalRadialDistributionFunction> propertyConventionalRadialDistributionFunction;
  std::optional<PropertyRadialDistributionFunction> propertyRadialDistributionFunction;
  std::optional<PropertyDensityGrid> propertyDensityGrid;
  std::optional<PropertyEnergyHistogram> averageEnergyHistogram;
  std::optional<PropertyNumberOfMoleculesHistogram> averageNumberOfMoleculesHistogram;
  std::optional<PropertyMoleculeProperties> propertyMoleculeProperties;
  std::optional<PropertyMoleculeShape> propertyMoleculeShape;
  std::optional<PropertyMoleculeBackbone> propertyMoleculeBackbone;
  std::optional<PropertyMeanSquaredDisplacement> propertyMSD;
  std::optional<PropertyVelocityAutoCorrelationFunction> propertyVACF;
  std::optional<PropertyEndToEndAutoCorrelationFunction> propertyEndToEndACF;
  /// Waiting times of the living radicals of a polymerization ('ComputeRadicalWaitingTimes'); the samples are
  /// added by the Polymerization driver at the reactions.
  std::optional<PropertyRadicalWaitingTimes> propertyRadicalWaitingTimes;
  std::optional<WriteLammpsData> writeLammpsData;

  std::optional<PropertyNumberOfMoleculesEvolution> propertyNumberOfMoleculesEvolution;
  std::optional<PropertyVolumeEvolution> propertyVolumeEvolution;
  std::optional<PropertyConservedEnergyEvolution> propertyConservedEnergyEvolution;

  // Breakthrough settings
  std::size_t columnNumberOfGridPoints{100};
  double columnTotalPressure{1e5};
  double columnPressureGradient{0.0};
  double columnVoidFraction{0.4};
  double columnParticleDensity{1000};
  double columnEntranceVelocity{0.1};
  double columnLength{0.3};
  double columnTimeStep{0.0005};
  std::size_t columnNumberOfTimeSteps{0};
  bool columnAutoNumberOfTimeSteps{true};

  std::vector<std::optional<InterpolationEnergyGrid>> interpolationGrids;
  std::optional<InterpolationEnergyGrid> externalFieldInterpolationGrid;

  // Fractional molecules per component are stored in Move::Types enum order, then integer molecules.
  // Regions: SwapCFCMC (GC) -> SwapCBCFCMC (pair) -> PairSwapCFCMC -> PairSwapCBCFCMC -> GroupSwapCFCMC ->
  //          GroupSwapCBCFCMC -> GibbsSwapCFCMC -> ReactionConventionalCFCMC -> ReactionCFCMC ->
  //          GibbsConventionalCFCMC
  [[nodiscard]] std::size_t indexOfGCFractionalMoleculesPerComponent_CFCMC(
      std::size_t selectedComponent) const noexcept;
  [[nodiscard]] std::size_t indexOfPairGCFractionalMoleculesPerComponent_CFCMC(
      std::size_t selectedComponent) const noexcept;
  [[nodiscard]] std::size_t indexOfPairSwapFractionalMoleculesPerComponent_CFCMC(
      std::size_t selectedComponent) const noexcept;
  [[nodiscard]] std::size_t indexOfPairSwapCBFractionalMoleculesPerComponent_CFCMC(
      std::size_t selectedComponent) const noexcept;
  [[nodiscard]] std::size_t indexOfGroupSwapFractionalMoleculesPerComponent_CFCMC(
      std::size_t selectedComponent) const noexcept;
  [[nodiscard]] std::size_t indexOfGroupSwapCBFractionalMoleculesPerComponent_CFCMC(
      std::size_t selectedComponent) const noexcept;
  [[nodiscard]] std::size_t indexOfGibbsSwapFractionalMoleculesPerComponent_CFCMC(
      std::size_t selectedComponent) const noexcept;
  [[nodiscard]] std::size_t indexOfParallelReactionFractionalMoleculesPerComponent_CFCMC(
      std::size_t selectedComponent) const noexcept;
  [[nodiscard]] std::size_t indexOfSerialReactionFractionalMoleculesPerComponent_CFCMC(
      std::size_t selectedComponent) const noexcept;
  /// The serial-reaction fractional molecules of a component that are present (the active side of each
  /// serial reaction); the Gibbs-conventional region starts after them.
  [[nodiscard]] std::size_t numberOfPresentSerialReactionFractionalMolecules(
      std::size_t selectedComponent) const noexcept;
  [[nodiscard]] std::size_t indexOfGibbsConventionalFractionalMoleculesPerComponent_CFCMC(
      std::size_t selectedComponent) const noexcept;
  [[nodiscard]] std::size_t indexOfGibbsFractionalMoleculesPerComponent_CFCMC(
      std::size_t selectedComponent) const noexcept
  {
    return indexOfGibbsConventionalFractionalMoleculesPerComponent_CFCMC(selectedComponent);
  }
  [[nodiscard]] std::size_t indexOfFractionalMoleculeForMove(Move::Types move, std::size_t selectedComponent,
                                                             std::size_t subIndex = 0) const noexcept;
  [[nodiscard]] std::size_t parallelReactionFractionalMoleculeIndex(std::size_t reactionId, std::size_t componentId,
                                                                    bool isProduct,
                                                                    std::size_t localIndex) const noexcept;
  [[nodiscard]] std::size_t serialReactionFractionalMoleculeIndex(std::size_t reactionId, std::size_t componentId,
                                                                  std::size_t localIndex) const noexcept;

  void addComponent(const Component&& component) noexcept(false);

  void createFrameworks();
  void createInitialMolecules(const std::vector<std::vector<double3>>& initialPositions);

  /**
   * \brief Retargets this system in place to a different framework and simulation box.
   *
   * Reinitializes all framework-dependent state (framework atoms and per-atom dynamics/field buffers,
   * framework atom counts, net charges, flexible-framework degrees of freedom, pseudo-atom counts, Ewald
   * parameters, and the rigid reciprocal-space precompute) so an existing System can be reused for a new
   * framework and cell without going through the full constructor. This is the single-object equivalent of
   * the constructor's framework setup. Any guest-molecule atoms (the storage after the framework prefix) are
   * preserved; because their positions are expressed in the previous cell, callers that also change the box
   * (e.g. reduction to a primitive cell for phonon dispersion) should apply this to a framework-only system
   * or a copy so the retained conventional-cell system is not disturbed.
   *
   * \param newFramework New framework whose atoms and intramolecular potentials replace the current ones.
   * \param newSimulationBox Simulation box associated with \p newFramework.
   */
  void rebuildForFramework(const Framework& newFramework, const SimulationBox& newSimulationBox);

  void checkCartesianPositions();

  /**
   * \brief Starts the conserved-energy bookkeeping of an MD stage from 'runningEnergies'.
   *
   * Sets the drift reference and the current conserved energy, and clears the ramp-work sum.
   */
  void resetConservedEnergyBookkeeping();

  /**
   * \brief Conserved-energy bookkeeping after an MD step ('runningEnergies' holds the energies of the step).
   *
   * The change of the conserved energy over a step in which a polymer force ramp is active (or that ends one)
   * is the work of the time-dependent Hamiltonian, not integration error: it is added to 'polymerRampWork' and
   * to the drift reference, so that 'conservedEnergy - referenceEnergy' stays the integration drift across
   * the reactions. The integration error of the ramp steps themselves is attributed to the ramp work (a few
   * K per ramp against the 1e4-1e5 K a ramp injects). The relative drift of the step is accumulated for the
   * "Average drift".
   */
  void accountConservedEnergy();

  /**
   * \brief Accounts a jump of the conserved energy at a topology change (a reaction) as polymer ramp work.
   *
   * 'runningEnergies' holds the energies of the new topology; the jump against the previous value of
   * 'conservedEnergy' shifts the drift reference.
   */
  void accountTopologyChange();

  void precomputeTotalRigidEnergy() noexcept;

  /// Whether molecular dynamics conserves a total momentum that can be removed from the velocities (with its three
  /// degrees of freedom): more than one mobile entity, no lab-fixed framework atoms and no external field.
  [[nodiscard]] bool conservesCenterOfMassMomentum() const noexcept;
  /// Sets the translational centre-of-mass constraint of the system and its thermostat accordingly.
  void updateCenterOfMassConstraint() noexcept;
  void precomputeTotalGradients() noexcept;
  /// Framework-molecule and molecule-molecule VDW tail corrections (energy and dU/dlambda) of the current
  /// configuration, evaluated from aggregated per-type counts (O(N)). The gradient routines do not include them.
  RunningEnergy computeTailCorrectionEnergies() const noexcept;
  RunningEnergy computeTotalEnergies() noexcept;
  RunningEnergy computePolarizationEnergy() noexcept;
  RunningEnergy computeTotalGradients() noexcept;
  void computeTotalElectrostaticPotential() noexcept;

  /**
   * \brief Change of the cross-link energy when molecule 'selectedMolecule' of 'selectedComponent'
   *        takes the positions 'newAtoms' instead of its current positions 'oldAtoms'.
   *
   * Returns zero when the system has no cross-links or the molecule carries none, so coordinate
   * moves can add it unconditionally to their energy difference before the acceptance test.
   */
  RunningEnergy crossLinkEnergyDifference(std::size_t selectedComponent, std::size_t selectedMolecule,
                                          std::span<const Atom> newAtoms, std::span<const Atom> oldAtoms) const;

  /// Cross-link energy of the current configuration (all links), in the given box.
  RunningEnergy computeCrossLinkEnergy(const SimulationBox& box, std::span<const Atom> moleculeAtoms) const;

  /// The atoms of the molecule that currently carry a cross-link (sorted, unique); empty when unlinked.
  std::vector<std::size_t> crossLinkedSiteAtoms(std::size_t selectedComponent, std::size_t selectedMolecule) const;

  /// The cross-links of the molecule as tethers for a CBMC regrowth (see CrossLinkTether); empty when
  /// unlinked. The tethers freeze the partner positions, so they are valid for one move only.
  std::vector<CrossLinkTether> crossLinkTethers(std::size_t selectedComponent, std::size_t selectedMolecule) const;

  /**
   * \brief The placed set of a CBMC regrowth of a cross-linked molecule: its linked sites plus the
   *        given fixed atoms (sorted, unique). The linked sites stay in place, so the links are kept
   *        and only the rest of the molecule is regrown ('CBMC::FirstBeadScheme::AlreadyPlaced').
   *
   * std::nullopt when nothing would be regrown (every atom is a linked site or fixed).
   */
  std::optional<std::vector<std::size_t>> crossLinkRegrowthPlacedSet(std::size_t selectedComponent,
                                                                     std::size_t selectedMolecule,
                                                                     std::span<const std::size_t> fixedAtoms = {}) const;
  void computeTotalElectricField() noexcept;

  std::size_t randomFramework(RandomNumber& random)
  {
    return std::size_t(random.uniform() * static_cast<double>(numberOfFrameworks));
  }
  std::size_t randomComponent(RandomNumber& random)
  {
    return std::size_t(random.uniform() * static_cast<double>(components.size()));
  }
  std::size_t numerOfAdsorbateComponents() { return components.size(); }
  /// The number of components of type Adsorbate (the cations are not counted).
  [[nodiscard]] std::size_t adsorbateComponentCount() const
  {
    std::size_t count = 0uz;
    for (const Component& component : components)
    {
      if (component.type == Component::Type::Adsorbate) ++count;
    }
    return count;
  }

  /// The change of the intramolecular Coulomb energy when molecules are rescaled in place (same positions, new
  /// Coulomb scaling: a CFCMC lambda change, a fractional molecule made whole or a whole one made fractional).
  /// For a flexible molecule the non-excluded intramolecular pairs carry the scaled exclusion kernel
  /// lambda_A lambda_B C q_A q_B K(r) (Potentials::intraMolecularCoulomb), which compensates the scaled pair the
  /// Ewald Fourier sum counts; the Fourier difference of such a move changes by the opposite amount, so a move
  /// that adds only that difference drifts in both terms. 'newAtoms' and 'oldAtoms' hold the same slots in the
  /// same order (several molecules are split by their molecule id); only a molecule whose positions are unchanged
  /// is rescaled in place and counted (a slot that receives another conformation, as in the Gibbs transfers, has
  /// its internal energies accounted by the move itself); zero for rigid molecules.
  [[nodiscard]] RunningEnergy intraCoulombRescaleDifference(std::span<const Atom> newAtoms,
                                                            std::span<const Atom> oldAtoms) const
  {
    RunningEnergy difference{};
    if (!forceField.useCharge || newAtoms.size() != oldAtoms.size()) return difference;
    std::size_t begin = 0;
    while (begin < newAtoms.size())
    {
      std::size_t end = begin + 1;
      while (end < newAtoms.size() && newAtoms[end].moleculeId == newAtoms[begin].moleculeId) ++end;
      const std::size_t componentId = static_cast<std::size_t>(newAtoms[begin].componentId);
      bool inPlace = true;
      for (std::size_t a = begin; a < end && inPlace; ++a)
      {
        const double3 d = newAtoms[a].position - oldAtoms[a].position;
        inPlace = d.x == 0.0 && d.y == 0.0 && d.z == 0.0;
      }
      if (inPlace && componentId < components.size())
      {
        const auto& potentials = components[componentId].intraMolecularPotentials;
        difference = difference +
                     potentials.computeInternalIntraCoulombEnergies(forceField, simulationBox,
                                                                    newAtoms.subspan(begin, end - begin)) -
                     potentials.computeInternalIntraCoulombEnergies(forceField, simulationBox,
                                                                    oldAtoms.subspan(begin, end - begin));
      }
      begin = end;
    }
    return difference;
  }

  /// The net charge the Monte Carlo moves hand to the Ewald energy difference for its net-charge correction:
  /// the whole system's, or under 'omitInterInteractions' the framework's alone (the correction is then the
  /// cross term of the moved charge with the framework, as in the full evaluations).
  [[nodiscard]] double ewaldMoveNetCharge() const
  {
    return forceField.omitInterInteractions ? netChargeFramework : netCharge;
  }

  // --- transition-matrix Monte Carlo (TMMC) ------------------------------------------------------
  // The macrostate is the molecule count of one component, the adsorbate. A cation beside it is a
  // component too, picked for its own moves, but the matrix is not of its count: its moves must neither
  // see the bias (a cation count outside the adsorbate's window would be refused as out of bound, and
  // biased by the adsorbate's ln Pi otherwise) nor write into the collection matrix.

  /// The component whose molecule count is the TMMC macrostate: the first adsorbate (component 0 when
  /// there is none).
  [[nodiscard]] std::size_t tmmcComponent() const
  {
    for (std::size_t componentId = 0uz; componentId < components.size(); ++componentId)
    {
      if (components[componentId].type == Component::Type::Adsorbate) return componentId;
    }
    return 0uz;
  }
  /// Whether the transition matrix biases and records the moves of 'componentId': TMMC is on and the
  /// component is the one whose count is the macrostate.
  [[nodiscard]] bool tmmcTracks(std::size_t componentId) const { return tmmc.doTMMC && componentId == tmmcComponent(); }
  /// The transition-matrix bias of a move of 'componentId' from 'oldN' to 'newN' molecules; 1 when the
  /// matrix does not track that component.
  [[nodiscard]] double tmmcBiasFactor(std::size_t componentId, std::size_t newN, std::size_t oldN)
  {
    return tmmcTracks(componentId) ? tmmc.biasFactor(newN, oldN) : 1.0;
  }
  /// Whether the integer swaps of 'componentId' (Swap, SwapCBMC, SwapMultiTrialCBMC, the pair and group
  /// swaps, and every Gibbs transfer: N -> N +/- 1 with the lambda of the chain unchanged) are steps of the
  /// transition matrix. On a plain N chain they are. On the (N, lambda) chain of SwapCFCMC/SwapCBCFCMC the
  /// only steps are the lambda hops (N, k) -> (N, k +/- 1) and the endpoint crossings (N, last) <-> (N+1, 0);
  /// an integer swap jumps (N, k) -> (N +/- 1, k), which the tridiagonal matrix cannot hold, so such a move
  /// is refused and not recorded (the input reader refuses the mix up front). The Gibbs CFCMC transfers are
  /// such jumps too: the direct ones turn the fractional molecule whole at its bin and make a molecule of the
  /// other box fractional at that same bin, and the conventional ones drive the Gibbs lambda, not the chain's.
  [[nodiscard]] bool tmmcTracksIntegerSwap(std::size_t componentId) const
  {
    return tmmcTracks(componentId) && !tmmc.lambdaChain();
  }
  /// The transition-matrix bias of an integer swap of 'componentId' from 'oldN' to 'newN' molecules: 1 when
  /// the matrix does not track that component, 0 (refused) on an (N, lambda) chain.
  [[nodiscard]] double tmmcIntegerSwapBiasFactor(std::size_t componentId, std::size_t newN, std::size_t oldN)
  {
    if (!tmmcTracks(componentId)) return 1.0;
    if (tmmc.lambdaChain()) return 0.0;
    return tmmc.biasFactor(newN, oldN);
  }
  /// The change of the molecule count of 'componentId' by one group insertion of the group of
  /// 'centralComponent' (the central molecule and its satellites; a group deletion removes the same).
  [[nodiscard]] std::size_t groupMoleculeCountChange(std::size_t centralComponent, std::size_t componentId) const
  {
    std::size_t change = (centralComponent == componentId) ? 1uz : 0uz;
    for (std::size_t satelliteComponent : components[centralComponent].groupComponentIds)
    {
      if (satelliteComponent == componentId) ++change;
    }
    return change;
  }
  /// Whether a group swap of 'centralComponent' fits the one-dimensional nearest-neighbour transition
  /// matrix of the tracked component 'trackedComponent': the group changes the tracked count by exactly
  /// its central molecule (a tracked central molecule without tracked satellites), or not at all (an
  /// untracked central molecule without tracked satellites). Any other group moves the macrostate by
  /// more than one, or by one on a move the matrix does not see (and the matrix would record a N -> N
  /// +/- m transition as N -> N +/- 1).
  [[nodiscard]] bool groupSwapFitsTMMC(std::size_t centralComponent, std::size_t trackedComponent) const
  {
    return groupMoleculeCountChange(centralComponent, trackedComponent) ==
           (centralComponent == trackedComponent ? 1uz : 0uz);
  }
  /// The same for the tracked component of this system; true when TMMC is off.
  [[nodiscard]] bool groupSwapFitsTMMC(std::size_t centralComponent) const
  {
    return !tmmc.doTMMC || groupSwapFitsTMMC(centralComponent, tmmcComponent());
  }
  std::size_t randomMoleculeOfComponent(RandomNumber& random, std::size_t selectedComponent);
  std::size_t randomIntegerMoleculeOfComponent(RandomNumber& random, std::size_t selectedComponent);

  std::size_t globalIndexOfComponentAndMolecule(std::size_t selectedComponent, std::size_t selectedMolecule);

  std::size_t indexOfFirstMolecule(std::size_t selectedComponent);
  std::vector<Atom>::iterator iteratorForMolecule(std::size_t selectedComponent, std::size_t selectedMolecule);
  std::vector<double3>::iterator iteratorForElectricField(std::size_t selectedComponent, std::size_t selectedMolecule);
  std::vector<Molecule>::iterator indexForMolecule(std::size_t selectedComponent, std::size_t selectedMolecule);
  std::size_t moleculeIndexOfComponent(std::size_t selectedComponent, std::size_t selectedMolecule);
  std::span<Atom> spanOfMolecule(std::size_t selectedComponent, std::size_t selectedMolecule);
  const std::span<const Atom> spanOfMolecule(std::size_t selectedComponent, std::size_t selectedMolecule) const;
  const std::span<const Atom> spanOfIntegerAtomsOfComponent(std::size_t selectedComponent) const;
  std::span<const Atom> spanOfFrameworkAtoms() const;
  std::span<Atom> spanOfFrameworkAtoms();
  std::span<const Atom> spanOfRigidFrameworkAtoms() const;
  std::span<const Atom> spanOfFlexibleAtoms() const;
  std::span<const Atom> spanOfMoleculeAtoms() const;
  std::span<Atom> spanOfMoleculeAtoms();

  /**
   * \brief The CBMC growth context of this system: its environment (force field, box, framework,
   * grids, all molecule atoms as background, beta) with the cut-offs of 'mode'.
   *
   * The default 'Growth' mode makes the dual cut-off decision (inner cut-off when the force field
   * enables the scheme, full cut-offs otherwise); the CBMC entry points correct their results to the
   * full cut-offs themselves whenever the context grew at the inner cut-off. The inter-molecular
   * background is queried through 'cellList()' (built on demand) and the framework through
   * 'frameworkCellList()'.
   */
  CBMC::GrowContext makeGrowContext(CBMC::CutOffMode mode = CBMC::CutOffMode::Growth) const;

  /**
   * \brief The CBMC growth context of this system grown against 'background' instead of its own
   * molecule atoms (e.g. the molecule atoms with a pair removed, or with already grown group
   * members appended).
   *
   * Equivalent to 'makeGrowContext(mode).withMoleculeAtoms(background)' but does not touch
   * 'cellList()': the molecule cell list describes 'spanOfMoleculeAtoms()' only, so a context over
   * another background falls back to the brute-force loop anyway and building (or rebuilding) the
   * list for it would be wasted work. The framework cell list is kept since the framework is the same.
   */
  CBMC::GrowContext makeGrowContext(std::span<const Atom> background,
                                    CBMC::CutOffMode mode = CBMC::CutOffMode::Growth) const;

  /**
   * \brief Persistent cell list over 'spanOfMoleculeAtoms()' for the Monte Carlo energy differences.
   *
   * A cache, hence mutable and not archived: 'cellList()' returns it current (rebuilt when the atom
   * array, box or cut-off changed, or after 'invalidateCellList'). Accepted single-molecule moves apply
   * their new positions with 'cellListAtomsMoved'; moves that change positions in any other way are
   * followed by 'invalidateCellList' (see 'MC_Moves::maintainsCellList'), as are molecule insertions
   * and deletions ('updateMoleculeAtomInformation'). The list is disabled (and the queries use their
   * brute-force loops) when the box holds fewer than four cells of the cut-off in every direction.
   */
  mutable MCCellList moleculeCellList;
  [[nodiscard]] const MCCellList &cellList() const;
  /// The cell list rebuilt from the current positions (for full recomputations that must not trust the
  /// incremental maintenance); O(N).
  [[nodiscard]] const MCCellList &rebuiltCellList() const;
  [[nodiscard]] double cellListCutOff() const;
  void invalidateCellList() const { moleculeCellList.invalidate(); }
  /// 'movedAtoms' must be a sub-span of 'spanOfMoleculeAtoms()' whose positions were just updated.
  void cellListAtomsMoved(std::span<const Atom> movedAtoms) const;
  /// Rebuilds a reference list and compares it with the maintained one (consistency check).
  [[nodiscard]] bool verifyCellList() const;

  /**
   * \brief Cell list over 'spanOfFrameworkAtoms()' for the framework-molecule energies of trial positions
   * (the Monte Carlo moves and the CBMC grows without interpolation grids).
   *
   * A cache like 'moleculeCellList': 'frameworkCellList()' returns it current for the framework atoms, the
   * box and the framework cut-off ('frameworkCellListCutOff'). The framework atoms never move in a Monte Carlo
   * simulation with a rigid framework, so the list is built once (and again after a box change); for a
   * framework with mobile atoms it is kept disabled (the queries then gather by a pass over the framework
   * atoms, which reads the current positions).
   */
  mutable MCCellList frameworkMCCellList;
  [[nodiscard]] const MCCellList& frameworkCellList() const;
  [[nodiscard]] double frameworkCellListCutOff() const;
  void invalidateFrameworkCellList() const { frameworkMCCellList.invalidate(); }

  /**
   * \brief Solute tempering (replica exchange with solute scaling, REST2).
   *
   * The Hamiltonian of one component (the 'solute') is scaled: its intramolecular energy and the
   * pair interactions among its atoms by lambda, its interactions with everything else by
   * sqrt(lambda), realized through the pair parameters of its pseudo-atom types, which therefore
   * must not be shared with any other component or the framework, and its partial charges (scaled
   * by sqrt(lambda)). lambda = T / T_effective: the solute samples as if at a higher temperature
   * while the solvent stays at the system temperature, and the exchange acceptance between replicas
   * depends on the solute-involving energies only (see MC_Moves::ParallelTemperingLogAcceptance).
   */
  std::optional<std::size_t> soluteTemperingComponent{};
  double soluteTemperingLambda{1.0};

  /// The pseudo-atom types used by 'componentId' (one flag per pseudo-atom of the force field).
  [[nodiscard]] std::vector<bool> pseudoAtomTypesOfComponent(std::size_t componentId) const;

  /**
   * \brief Scales the Hamiltonian of 'componentId' by 'lambda' (see 'soluteTemperingComponent').
   *
   * Applies the scaling to the force field, the component definition and the charges of the present
   * molecules, and recomputes the derived state (Ewald structure factors, running energies). Throws
   * when the component shares a pseudo-atom type with another component or the framework. Cumulative:
   * the stored lambda is multiplied, so this is normally called once on a freshly constructed system.
   */
  void scaleSoluteHamiltonian(std::size_t componentId, double lambda);

  /// Multiplies the partial charges of the solute molecules present in the system by 'factor' and
  /// updates the per-component net charges (used after a configuration exchange between replicas
  /// with different lambda: the atoms travel, the Hamiltonian stays).
  void rescaleSoluteCharges(double factor);

  /**
   * \brief Takes over the interpolation grids of 'reference' (a replica of the same system), scaled to the
   * Hamiltonian of this system.
   *
   * The replicas of a parallel-tempering run share the framework and external-field grids, which are computed
   * once. The framework Van der Waals grid of a pseudo-atom type tabulates the framework interaction with the
   * force field it was built with; under solute tempering the solute's framework interaction is scaled by
   * sqrt(lambda), so the solute-type grids of a replica with a different lambda are the reference grids times
   * sqrt(lambda / lambda_reference). The Coulomb grid is per unit charge and the solute charges carry their own
   * scaling, so it is shared as is. The running energies are not recomputed.
   */
  void copyInterpolationGridsFrom(const System& reference);

  /**
   * \brief Recomputes everything derived from the configuration after the atoms/molecules were
   * replaced wholesale (configuration exchange between replicas, restart): Ewald parameters and
   * structure factors, running energies, loadings, per-molecule bookkeeping, pseudo-atom and
   * tail-correction counts, net charge.
   */
  void rebuildConfigurationDerivedState();

  /**
   * \brief The growth context of an isolated molecule: no framework, no interpolation grids, no
   * external field, and no background molecules, evaluated at the full cut-offs and always grown with
   * configurational bias.
   *
   * For the ideal-gas grows (the conformation reservoirs and the recoil-growth reference
   * conformations): the external energies are identically zero, only the intramolecular terms act.
   * Configurational bias is forced because these conformations seed recoil growth's own reference
   * and Widom sampling, which must not depend on the production chain scheme.
   */
  CBMC::GrowContext makeIdealGasGrowContext() const;

  std::span<const AtomDynamics> spanOfFrameworkDynamics() const;
  std::span<AtomDynamics> spanOfFrameworkDynamics();
  std::span<const AtomDynamics> spanOfMoleculeDynamics() const;
  std::span<AtomDynamics> spanOfMoleculeDynamics();

  /// Rebuilds 'groupData' from the current molecule atom positions: one rigid-body state per rigid
  /// fragment of every semi-flexible molecule, recovered exactly via 'Component::deriveFragmentState'.
  void initializeGroupData();
  void initializeFrameworkGroupData();
  std::span<GroupState> spanOfGroupData();
  std::span<const GroupState> spanOfGroupData() const;
  std::span<GroupState> spanOfFrameworkGroupData();
  std::span<const GroupState> spanOfFrameworkGroupData() const;
  std::span<double> spanOfMoleculeElectrostaticPotential();
  std::span<double3> spanOfMoleculeElectricField();
  std::span<const double3> spanOfMoleculeElectricField() const;
  std::span<double3> spanOfMoleculeElectricFieldNew();
  std::span<double3> spanElectricFieldNew(std::size_t selectedComponent, std::size_t selectedMolecule);
  const std::span<const double3> spanElectricFieldNew(std::size_t selectedComponent,
                                                      std::size_t selectedMolecule) const;
  std::span<double3> spanElectricFieldOld(std::size_t selectedComponent, std::size_t selectedMolecule);
  const std::span<const double3> spanElectricFieldOld(std::size_t selectedComponent,
                                                      std::size_t selectedMolecule) const;

  std::size_t numberOfMolecules() const
  {
    return std::accumulate(numberOfMoleculesPerComponent.begin(), numberOfMoleculesPerComponent.end(), std::size_t(0),
                           [](const std::size_t& acc, const std::size_t& b) { return acc + b; });
  }

  std::size_t numberOfIntegerMolecules() const
  {
    return std::accumulate(numberOfIntegerMoleculesPerComponent.begin(), numberOfIntegerMoleculesPerComponent.end(),
                           std::size_t(0), [](const std::size_t& acc, const std::size_t& b) { return acc + b; });
  }

  // The biased-to-Boltzmann reweighting factor of a configuration removes the bias of every
  // fractional molecule simultaneously. The biasing weight functions W_i(lambda_i) add in the
  // exponent, so the per-configuration weight is exp(-sum_i W_i) = prod_i exp(-W_i), i.e. the
  // product of the individual lambda weights (each weight() already equals exp(-W_i)). This
  // reduces to a single term when only one fractional molecule is present; inactive lambda
  // coordinates contribute exp(0) = 1 and drop out of the product.
  // The ion-pair and group CFCMC lambda coordinates live on the driving component of the pair/group
  // (see componentDrivesPairSwapLambda / componentDrivesGroupSwapLambda) and are included under the
  // same gating as in sampleProperties(), so that every biased coordinate is unbiased exactly once.
  // Improving the accuracy of computing chemical potentials in CFCMC simulations
  // A. Rahbari, R. Hens, D. Dubbeldam, and T.J.H Vlugt
  // Mol. Phys.  117(23-24), 3493-3508, 2019
  double weight() const
  {
    double w = 1.0;
    for (std::size_t componentId = 0; componentId < components.size(); ++componentId)
    {
      const Component& component = components[componentId];
      if (!component.hasFractionalMolecule)
      {
        continue;
      }

      w *= component.lambdaGC.weight() * component.lambdaGibbs.weight();

      if (componentDrivesPairSwapLambda(componentId, Move::Types::PairSwapCFCMC))
      {
        w *= component.lambdaPairSwap.weight();
      }
      if (componentDrivesPairSwapLambda(componentId, Move::Types::PairSwapCBCFCMC))
      {
        w *= component.lambdaPairSwapCB.weight();
      }
      if (componentDrivesGroupSwapLambda(componentId, Move::Types::GroupSwapCFCMC))
      {
        w *= component.lambdaGroupSwap.weight();
      }
      if (componentDrivesGroupSwapLambda(componentId, Move::Types::GroupSwapCBCFCMC))
      {
        w *= component.lambdaGroupSwapCB.weight();
      }
    }

    if (usesReactionConventionalCFCMC())
    {
      for (const Reaction& reaction : reactions.list)
      {
        w *= activeReactionLambdaHistogram(reaction).weight();
      }
    }

    return w;
  }

  /// Wang-Landau updates the GC lambda bias only while this box holds a real fractional molecule.
  /// `containsTheFractionalMolecule` alone is not enough: it defaults to true on system 0 for Gibbs
  /// ownership, including conventional GCMC that has no lambda coordinate.
  bool lambdaWangLandauIsActive(std::size_t componentId) const
  {
    return containsTheFractionalMolecule && components[componentId].hasFractionalMolecule;
  }

  void removeRedundantMoves();
  void rescaleMoveProbabilities();
  void determineSwappableComponents();
  void determineFractionalComponents();
  void createReactionFractionalMolecules();
  void createParallelReactionFractionalMolecules();
  void createSerialReactionFractionalMolecules();
  void precomputeReactionFractionalLayout() noexcept;
  void syncReactionFractionalMoleculeIndices() noexcept;
  void incrementReactionFractionalMoleculeIds(std::size_t componentId) noexcept;
  void initializeReactionLambdaHistograms(std::size_t numberOfBlocks, std::size_t numberOfLambdaBins);
  [[nodiscard]] double reactionDUdlambda(const Reaction& reaction) const noexcept;

  /// Instantaneous dU/dlambda of the lambda coordinate driving dUdlambda-group 'groupId'. The
  /// VDW/charge/Ewald parts are maintained incrementally in 'runningEnergies'; the polarization part is
  /// evaluated on demand from the incrementally-maintained stored electric field.
  [[nodiscard]] double currentDUdlambda(double lambda, std::size_t groupId) const noexcept;
  void syncReactionLambdaBin(Reaction& reaction) noexcept;
  void syncReactionLambdaBins() noexcept;
  [[nodiscard]] bool usesReactionConventionalCFCMC() const noexcept;
  [[nodiscard]] bool usesSerialReactionCFCMC() const noexcept;
  [[nodiscard]] bool usesParallelReactionCFCMC() const noexcept;
  [[nodiscard]] bool usesGibbsConventionalCFCMC() const noexcept;
  void initializeGibbsConventionalFractionalMolecules() noexcept;
  void initializeFixedLambdaFractionalMolecules();
  void initializeGibbsSwapFractionalMoleculeGroupIds() noexcept;
  [[nodiscard]] bool hasReactionFractionalMolecules() const noexcept;
  [[nodiscard]] PropertyLambdaProbabilityHistogram& activeReactionLambdaHistogram(Reaction& reaction) noexcept;
  [[nodiscard]] const PropertyLambdaProbabilityHistogram& activeReactionLambdaHistogram(
      const Reaction& reaction) const noexcept;
  void reactionLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase phase) noexcept;
  void reactionLambdaSampleOccupancy() noexcept;
  void reactionLambdaClearBookkeeping() noexcept;
  void reactionLambdaFinalize() noexcept;
  [[nodiscard]] bool componentDrivesPairSwapLambda(std::size_t componentId, Move::Types move) const noexcept;
  [[nodiscard]] bool componentDrivesGroupSwapLambda(std::size_t componentId, Move::Types move) const noexcept;
  /// The component whose group-swap lambda histogram drives the first group fractional sub-slot of
  /// 'componentId' for the given group CFCMC move (the driver itself or the first component listing
  /// 'componentId' among its GroupComponents). Returns std::nullopt when there is none. A component may
  /// belong to several groups (a counter-ion shared by two salts): its group sub-slots are laid out
  /// driver by driver in component order, each driver's central slot first and then one per occurrence
  /// among that driver's satellites; see groupSwapFractionalSubSlotRange and groupSwapLambdaDriverOfSubSlot.
  [[nodiscard]] std::optional<std::size_t> groupSwapLambdaDriver(std::size_t componentId,
                                                                 Move::Types move) const noexcept;
  /// The group fractional sub-slots of 'componentId' that belong to the group of 'driverComponentId' for
  /// the given group CFCMC move: {first sub-slot, count} (count 0 when the component is not in that group).
  [[nodiscard]] std::pair<std::size_t, std::size_t> groupSwapFractionalSubSlotRange(std::size_t componentId,
                                                                                    std::size_t driverComponentId,
                                                                                    Move::Types move) const noexcept;
  /// The driver whose group owns the group fractional sub-slot 'subIndex' of 'componentId'.
  [[nodiscard]] std::optional<std::size_t> groupSwapLambdaDriverOfSubSlot(std::size_t componentId, std::size_t subIndex,
                                                                          Move::Types move) const noexcept;
  void assignDUdlambdaGroups();
  [[nodiscard]] std::uint8_t fractionalSlotDUdlambdaGroupId(std::size_t componentId,
                                                            std::size_t slotIndex) const noexcept;
  void pairSwapLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase phase) noexcept;
  void pairSwapLambdaSampleOccupancy() noexcept;
  void pairSwapLambdaClearBookkeeping() noexcept;
  [[nodiscard]] double pairSwapLambdaMinBias() const noexcept;
  void pairSwapLambdaNormalize(double minBias) noexcept;
  void pairSwapLambdaWriteBiasingFiles(std::size_t systemId);
  [[nodiscard]] double reactionLambdaMinBias() const noexcept;
  void reactionLambdaNormalize(double minBias) noexcept;
  void reactionLambdaWriteBiasingFiles(std::size_t systemId);
  void reactionLambdaSampleProductionHistograms(std::size_t blockIndex, double weight) noexcept;
  void rescaleMolarFractions();
  void computeComponentFluidProperties();
  void computeNumberOfPseudoAtoms();
  /// Rebuilds the scaling-dependent aggregates from the current molecule atoms: effectiveNumberOfPseudoAtomsVDW and
  /// fractionalPseudoAtomCountsPerGroup (tail corrections) and the Coulomb-scaled net charges (Ewald net-charge
  /// correction, see 'computeNetCharge'). Called after every accepted change of the atom scaling.
  void computeTailCorrectionCounts();
  /// Rebuilds netChargeAdsorbates, netChargePerComponent and netCharge as the Coulomb-scaled charge of the current
  /// molecule atoms plus the framework charge.
  void computeNetCharge();
  /// Coulomb-scaled charge sum(scalingCoulomb * charge) of a set of atoms.
  [[nodiscard]] static double scaledNetCharge(std::span<const Atom> atoms);
  /// Incrementally adds a single atom's contribution to the aggregated tail-correction counts.
  void addAtomToTailCorrectionCounts(const Atom& atom);
  /// Incrementally removes a single atom's contribution from the aggregated tail-correction counts.
  void removeAtomFromTailCorrectionCounts(const Atom& atom);
  void optimizeMCMoves();

  std::string writeOutputHeader() const;
  std::string writeNumberOfPseudoAtoms() const;
  std::string writePreInitializationStatusReport(std::size_t currentCycle, std::size_t numberOfProductionCycles) const;
  std::string writeInitializationStatusReport(std::size_t currentCycle, std::size_t numberOfProductionCycles) const;
  std::string writeEquilibrationStatusReportMC(std::size_t currentCycle, std::size_t numberOfProductionCycles) const;
  /// One line describing the Monte Carlo cell list (grid, neighbourhood size, rebuild count), or that it is disabled.
  std::string writeCellListStatus() const;
  std::string writeEquilibrationStatusReportMD(std::size_t currentCycle, std::size_t numberOfProductionCycles) const;
  /// Production status report; 'progress' = (current cycle, number of cycles) enables the rate/ETA line.
  std::string writeProductionStatusReportMC(
      const std::string& statusLine, std::optional<std::pair<std::size_t, std::size_t>> progress = std::nullopt) const;
  std::string writeProductionStatusReportMD(std::size_t currentCycle, std::size_t numberOfProductionCycles) const;
  /// Writes the files of the properties sampled by the Monte Carlo drivers (radial distribution functions, density
  /// grid, energy and number-of-molecules histograms, molecule properties, end-to-end ACF, shape and backbone).
  void writeSampledPropertiesMC(std::size_t systemId, std::size_t currentCycle);
  /// Writes the time-evolution properties (number of molecules, volume, conserved energy) at 'absoluteCurrentCycle';
  /// 'lastWrite' (with the number of cycles completed) writes every collected sample.
  void writePropertiesEvolution(std::size_t systemId, std::size_t absoluteCurrentCycle, bool lastWrite = false);
  std::string writeSystemStatus() const;
  std::string writeComponentStatus() const;
  std::string writeMCMoveStatistics() const;

  /// Whether the periodic status reports include the compact move-statistics table ('PrintMoveStatistics').
  bool printMoveStatistics{true};

  /// Cycle and wall-clock time of the previous periodic status report (for the rate/ETA line).
  mutable std::optional<std::pair<std::size_t, std::chrono::steady_clock::time_point>> previousStatusReport;

  /**
   * \brief Compact move-statistics table for the periodic status reports.
   *
   * One line per Monte Carlo move (per component and for the system moves) with the attempts and
   * acceptance since the previous report and cumulative, the current maximum change, and the share of
   * the move CPU time. Moves with at least 1000 attempts and no acceptance are flagged. Marks the
   * statistics as reported (start of the next window). Empty when there are no moves.
   */
  std::string writeMoveStatisticsSummary() const;

  /// Elapsed time since the previous report, cycle rate, and estimated time to the end of the stage.
  std::string writeProgressLine(std::size_t currentCycle, std::size_t numberOfCycles) const;

  nlohmann::json jsonSystemStatus() const;
  nlohmann::json jsonComponentStatus() const;
  nlohmann::json jsonMCMoveStatistics() const;

  void insertMolecule(std::size_t selectedComponent, const Molecule& molecule, std::vector<Atom> atoms);
  void insertMoleculePolarization(std::size_t selectedComponent, const Molecule& molecule, std::vector<Atom> atoms,
                                  std::span<double3> electricField);
  void insertFractionalMolecule(std::size_t selectedComponent, const Molecule& molecule, std::vector<Atom> atoms,
                                std::size_t moleculeIndex);
  void insertFractionalMoleculeAtIndex(std::size_t selectedComponent, std::size_t moleculeIndex,
                                       const Molecule& molecule, std::vector<Atom> atoms);
  void insertReactionFractionalMolecule(std::size_t selectedComponent, std::size_t moleculeIndex,
                                        const Molecule& molecule, std::vector<Atom> atoms, bool isReactant,
                                        double lambda, std::uint8_t dUdlambdaGroupId);
  void insertSerialReactionFractionalMolecule(std::size_t selectedComponent, std::size_t moleculeIndex,
                                              const Molecule& molecule, std::vector<Atom> atoms, double lambda,
                                              std::uint8_t dUdlambdaGroupId);
  void deleteMolecule(std::size_t selectedComponent, std::size_t selectedMolecule, const std::span<Atom> atoms);
  void deleteFractionalMolecule(std::size_t selectedComponent, std::size_t selectedMolecule,
                                const std::span<Atom> atoms);

  /**
   * \brief Exchanges two molecules of one component, slot for slot.
   *
   * The atoms with their dynamics and stored electric field, and the molecule record, change places; what
   * names a slot stays behind: the molecule id of the atoms and the atom index of the record. The cross-link
   * sites follow their molecules. The counts, the charges and the tail-correction bookkeeping are unchanged by
   * construction. The cell list is kept current.
   */
  void swapMolecules(std::size_t selectedComponent, std::size_t moleculeA, std::size_t moleculeB);

  /**
   * \brief Makes fractional molecule 'moleculeIndex' of the component a whole molecule in place: it moves to the
   * end of the component with its atoms, dynamics, field, tags and record (the molecules after it shift down by
   * one), and leaves the fractional count for the integer count.
   *
   * The scaling of its atoms is the caller's, and so are the net charge and the tail-correction counts that depend
   * on it ('computeTailCorrectionCounts'). The cell list is invalidated.
   */
  void makeFractionalMoleculeWhole(std::size_t selectedComponent, std::size_t moleculeIndex);

  /**
   * \brief Makes whole molecule 'moleculeIndex' of the component the fractional molecule of slot 'slotIndex' in
   * place: it moves to the slot with its atoms, dynamics, field, tags and record (the molecules from the slot up
   * to its old position shift up by one), and leaves the integer count for the fractional count. The molecule must
   * not carry cross-links. The scaling and the counts are the caller's, as in 'makeFractionalMoleculeWhole'.
   */
  void makeWholeMoleculeFractional(std::size_t selectedComponent, std::size_t moleculeIndex, std::size_t slotIndex);

  /// Moves molecule 'moleculeIndex' of the component to 'newIndex' of the component, the molecules in between
  /// shifting by one; the atoms with their dynamics, field and tags, the record and the cross-link sites follow.
  void moveMolecule(std::size_t selectedComponent, std::size_t moleculeIndex, std::size_t newIndex);

  /**
   * \brief The exchange of 'swapMolecules' between a molecule of this system and one of 'other', of the same
   * component index in both (the boxes of a Gibbs ensemble).
   *
   * The cell list of each system is refreshed for its slot; a caller that goes on editing the exchanged atoms (a
   * slot records the group id and scaling of its atom) applies the span again with 'cellListAtomsMoved'. Calling
   * it a second time with the same arguments undoes the exchange. A molecule carrying a cross-link can not change
   * system (the moves guarantee that): the call throws.
   */
  void exchangeMolecule(System& other, std::size_t selectedComponent, std::size_t molecule, std::size_t otherMolecule);

  /**
   * \brief Exchanges the configuration with 'other' (replica exchange).
   *
   * Everything that belongs to the positions changes systems: the adsorbate atoms with their dynamics, electric
   * fields and potentials, the molecule and rigid-group records, the molecule counts, the adsorbate charges and
   * degrees of freedom, the cross-link topology, the fractional-molecule state with its lambda bins, and the box
   * when there is no framework. The thermodynamic state, the force field, the learned biases, the move controls
   * and the statistics stay. The cell lists are invalidated; both systems need 'rebuildConfigurationDerivedState'
   * afterwards.
   */
  void exchangeConfiguration(System& other);

  /// Throws when the polymer link table is activated: its group bookkeeping is laid out over the molecules at
  /// activation, so no molecule may be inserted, deleted, moved or exchanged afterwards ('operation' names the
  /// attempt). Called by every primitive above that changes the molecule layout.
  void checkPolymerLayoutIsFixed(std::string_view operation) const;

  /// Hands the lowest free tags, or new ones, to 'count' atoms just inserted at 'offset' of 'atomData'.
  void insertAtomTags(std::size_t offset, std::size_t count);
  /// Frees the tags of the 'count' atoms at 'offset' of 'atomData' (before their erase) and erases them.
  void eraseAtomTags(std::size_t offset, std::size_t count);
  /// Numbers the atoms afresh in storage order.
  void resetAtomTags();
  /// Rebuilds 'atomIndexOfTag' from 'atomTags' (a restart); throws when a tag is the free-tag marker.
  void rebuildAtomIndexOfTag();
  /// Throws when the tags are not distinct with 'atomIndexOfTag' their exact inverse.
  void checkAtomTags() const;

  void updateMoleculeAtomInformation();
  void checkMoleculeIds();

  std::vector<Atom> randomConfiguration(RandomNumber& random, std::size_t selectedComponent,
                                        const std::span<const Atom> atoms);

  /**
   * \brief Generates an equilibrated ideal-gas molecule placed randomly in the simulation box.
   *
   * Produces a trial molecule at a uniformly random center-of-mass position in the box with a uniformly
   * random overall orientation. Rigid components use their fixed reference geometry. Flexible components
   * draw an internal conformation from the ideal-gas (intra-molecular) Boltzmann distribution
   * exp(-beta * U_intra) via equilibratedIdealGasConformation, so that the conventional (non-CBMC)
   * swap/CFCMC/Gibbs acceptance rules, which do not include the intra-molecular energy, sample the correct
   * distribution.
   *
   * \param random A random number generator instance.
   * \param selectedComponent The component to generate a molecule for.
   * \return A pair containing the equilibrated molecule and its corresponding atoms.
   */
  std::pair<Molecule, std::vector<Atom>> equilibratedIdealGasMoleculeRandomInBox(RandomNumber& random,
                                                                                 std::size_t selectedComponent);

  /**
   * \brief Samples an ideal-gas (isolated) internal conformation of a flexible molecule.
   *
   * Runs a short Markov chain of CBMC reinsertion moves on an isolated scratch molecule (no framework, no
   * other molecules, no external field), so the growth samples only the intra-molecular potential. The
   * accepted conformation is therefore an unbiased draw from the ideal-gas Boltzmann distribution
   * exp(-beta * U_intra). The result is returned centered on its mass center. Rigid components (or those
   * with fewer than two atoms) return the reference geometry unchanged.
   *
   * \param random A random number generator instance.
   * \param selectedComponent The component to sample an internal conformation for.
   * \return The sampled, mass-centered internal conformation.
   */
  std::vector<Atom> equilibratedIdealGasConformation(RandomNumber& random, std::size_t selectedComponent);

  /**
   * \brief Fills each flexible/semi-flexible component's ideal-gas conformation reservoir.
   *
   * For every component with more than one atom that is not fully rigid, draws a fixed number of
   * independent conformations from the ideal-gas Boltzmann distribution exp(-beta * U_intra) (via
   * equilibratedIdealGasConformation) and stores them in Component::conformationReservoir. The CBMC
   * internal Monte-Carlo seeds each grow from a random member, which decorrelates successive grows and
   * (together with the ring conformer-hopping move) samples every ring conformer. The reservoir seeds
   * the internal MC only; it does not enter any Rosenbluth weight. Built once at the system temperature
   * with a fixed local RNG so it does not perturb the simulation's random stream; a temperature change
   * (e.g. parallel tempering) would require rebuilding it. Because drivers that change the temperature
   * call this to rebuild, it first re-prepares the growth plans for the current 'beta'
   * (prepareGrowthPlans).
   */
  void buildConformationReservoirs();

  /**
   * \brief Prepares every component's cached CBMC growth plans for the system's current 'beta'
   * (Component::prepareGrowthPlans): fills the frozen base-coupling constants of the flexible attach
   * steps. Must run before the first grow, and again after any change of the temperature.
   */
  void prepareGrowthPlans();

  /**
   * \brief Fills each flexible component's recoil-growth openness reference conformations.
   *
   * When recoil growth is enabled, draws a few independent ideal-gas conformations per flexible
   * multi-atom component (via equilibratedIdealGasConformation, forced through the CBMC path) and
   * stores them in Component::recoilReferenceConformations. The recoil-growth open/closed test measures
   * every trial direction against the same growth step evaluated in these conformations (per-step
   * maximum energy), so molecules whose correctly grown chains carry unavoidable intramolecular strain
   * (crowded united-atom beads) are not rejected as 'closed'. Built once with a fixed local RNG so the
   * reference is a fixed constant of the run (grow and retrace share it; detailed balance holds).
   */
  void buildRecoilReferenceConformations();

  /**
   * \brief Works out the blocking pockets of the components that asked for them.
   *
   * A component with 'BlockingPockets' : 'auto' has its spheres measured from the framework here, for a probe
   * standing for the molecule, and the request is cleared once they are. This has to happen after the
   * framework is in place and before a single molecule is put anywhere, which is the whole of the reason it
   * is done by the system rather than by the reader: the spheres are a property of the framework, and the
   * framework does not exist while the components are being read.
   *
   * \throws std::runtime_error If a component asking for them has no framework or no atoms.
   */
  void computeAutomaticBlockingPockets();

  bool insideBlockedPockets(const Component& component, std::span<const Atom> molecule_atoms) const;

  /**
   * \brief Samples the block averages of the thermodynamic properties of the current state.
   *
   * \param kineticEnergiesFromRunningEnergies When true the temperatures are formed from the kinetic energies
   *        stored in 'runningEnergies' (set by the molecular-dynamics step) instead of being recomputed from the
   *        velocities in 'moleculeData' / the atom dynamics. The GPU-resident integrator keeps the velocities on
   *        the device and only downloads them when the host needs them (every 'PrintEvery' cycles or for a
   *        per-step property), so the host velocities are stale in between; the kinetic energies of the step are
   *        reduced on the device and are current every cycle.
   */
  void sampleProperties(std::size_t systemId, std::size_t currentBlock, std::size_t currentCycle,
                        bool kineticEnergiesFromRunningEnergies = false);
  void samplePropertiesEvolution(std::size_t absoluteCurrentCycle);

  void writeCPUTimeStatistics(std::ostream& stream) const;

  /**
   * Molecular (center-of-mass) excess pressure: evaluates framework + intermolecular + Ewald site
   * forces (no intramolecular) into a scratch gradient buffer, then applies the atomic-to-molecular
   * virial correction. Does not mutate live AtomDynamics gradients used by MD.
   */
  [[nodiscard]] std::pair<EnergyStatus, double3x3> computeMolecularPressure() noexcept;

  /// True when the force-based RDF is enabled and should sample on this cycle.
  [[nodiscard]] bool forceBasedRDFSampleDue(std::size_t currentCycle) const;

  /**
   * Accumulate the Borgis force-based RDF from site gradients already stored in atomDynamics
   * (e.g. after an MD integrator force evaluation). Does not recompute forces.
   */
  void sampleForceBasedRDFFromCurrentGradients(std::size_t currentCycle, std::size_t currentBlock);

  /**
   * Evaluate full site gradients (framework + inter + Ewald + intramolecular) and sample the
   * force-based RDF. Use this on the MC path so flexible molecules get bonded forces in Borgis.
   */
  void sampleForceBasedRDFWithFullGradients(std::size_t currentCycle, std::size_t currentBlock);

  /**
   * Rigidly translate every molecule so its center of mass scales with the box, keeping internal geometry
   * fixed. Returns the scaled molecules and the scaled *molecule* atoms only (framework atoms excluded),
   * matching the molecule-atom spans expected by the inter-molecular and Ewald energy routines.
   */
  std::pair<std::vector<Molecule>, std::vector<Atom>> scaledCenterOfMassPositions(double scale) const;

  std::pair<std::vector<Molecule>, std::vector<Atom>> scaledCenterOfMassPositions(const SimulationBox& oldBox,
                                                                                  const SimulationBox& newBox) const;

  void writeComponentFittingStatus(std::ostream& stream, const std::vector<std::pair<double, double>>& rawData) const;

  void createExternalFieldInterpolationGrid(std::ostream& stream, std::size_t systemId);
  void createFrameworkInterpolationGrids(std::ostream& stream);

  void setThermostat(const std::optional<Thermostat>& thermo);
  void setThermobarostat(const std::optional<Thermobarostat>& barostat);

  void setSamplePDBMovie(const std::optional<SampleMovie>& movie);
  void updateSamplePDBMovie(std::size_t systemId, std::size_t currentCycle);

  void setNumberOfMoleculesHistogram(const std::optional<PropertyNumberOfMoleculesHistogram>& hist);
  void setAverageEnergyHistogram(const std::optional<PropertyEnergyHistogram>& hist);
  void setPropertyDensityGrid(const std::optional<PropertyDensityGrid>& grid);

  void setPropertyNumberOfMoleculesEvolution(std::optional<PropertyNumberOfMoleculesEvolution> property);
  void setPropertyVolumeEvolution(std::optional<PropertyVolumeEvolution> property);
  void setPropertyConservedEnergyEvolution(std::optional<PropertyConservedEnergyEvolution> property);

  void setPropertyConventionalRDF(const std::optional<PropertyConventionalRadialDistributionFunction>& rdf);
  void setPropertyRDF(const std::optional<PropertyRadialDistributionFunction>& rdf);
  void setPropertyMSD(const std::optional<PropertyMeanSquaredDisplacement>& msd);
  void setPropertyVACF(const std::optional<PropertyVelocityAutoCorrelationFunction>& vacf);
  void setPropertyRadicalWaitingTimes(const std::optional<PropertyRadicalWaitingTimes>& property);

  /// Creates the intra-molecular analyses requested by the components (their 'moleculePropertiesSettings',
  /// 'moleculeShapeSettings', 'moleculeBackboneSettings' and 'endToEndACFSettings'): the molecule-properties
  /// histograms, the shape descriptors, the backbone statistics and the end-to-end autocorrelation function.
  /// Call once the components, the molecules and the time step are known; throws when a component asks for an
  /// analysis it cannot support (e.g. the autocorrelation function without end-to-end atoms).
  /// 'timeAxisInPicoseconds' selects the lag axis of the autocorrelation function: the MD time step for the
  /// molecular-dynamics drivers, cycles for the Monte Carlo drivers (where 'timeStep' is meaningless).
  void initializeMoleculeProperties(std::size_t numberOfBlocks, bool timeAxisInPicoseconds);

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const System& s);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, System& s);

  void writeRestartFile(std::size_t systemId);

  std::string repr() const;
};
