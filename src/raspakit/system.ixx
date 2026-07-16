module;

export module system;

import std;

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
import property_msd;
import property_vacf;
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
import transition_matrix;
import equation_of_states;
import thermostat;
import thermobarostat;
import json;
import interpolation_energy_grid;
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

  std::uint64_t versionNumber{5};

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
  bool containsTheFractionalMolecule{true};

  // A contiguous list of adsorbate atoms per component for easy and fast looping
  // The atoms-order is defined as increasing per component and molecule.
  // Because the number of atoms is fixed per component it is easy to access the n-th molecule
  std::vector<Atom> atomData;
  // Cold per-atom molecular-dynamics fields (velocity, gradient), stored parallel to atomData with
  // identical length and indexing. Split out of Atom to keep Atom a single 64-byte cache line so the
  // energy inner loops touch only the bytes they need.
  std::vector<AtomDynamics> atomDynamics;
  std::vector<Molecule> moleculeData;
  std::vector<double> electricPotential;
  std::vector<double3> electricField;
  std::vector<double3> electricFieldNew;

  double conservedEnergy{};
  double referenceEnergy{};
  double accumulatedDrift{};
  RunningEnergy rigidEnergies;
  RunningEnergy runningEnergies;

  double3x3 currentExcessPressureTensor;
  EnergyStatus currentEnergyStatus;

  std::size_t numberOfHybridMCSteps{10};

  std::vector<std::complex<double>> eik_xy{};
  std::vector<std::complex<double>> eik_x{};
  std::vector<std::complex<double>> eik_y{};
  std::vector<std::complex<double>> eik_z{};
  std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> storedEik{};
  std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> fixedFrameworkStoredEik{};
  std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> trialEik{};
  double CoulombicFourierEnergySingleIon{0.0};
  double netCharge{0.0};
  double netChargeFramework{0.0};
  double netChargeAdsorbates{0.0};
  std::vector<double> netChargePerComponent;

  MCMoveProbabilities mc_moves_probabilities;
  MCMoveStatistics mc_moves_statistics;
  MCMoveCpuTime mc_moves_cputime;

  Reactions reactions;
  TransitionMatrix tmmc;

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
  std::optional<PropertyMeanSquaredDisplacement> propertyMSD;
  std::optional<PropertyVelocityAutoCorrelationFunction> propertyVACF;
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
  // Regions: SwapCFCMC (GC) -> SwapCBCFCMC (pair) -> PairSwapCFCMC -> PairSwapCBCFCMC -> GibbsSwapCFCMC ->
  //          ReactionConventionalCFCMC -> ReactionCFCMC -> GibbsConventionalCFCMC
  [[nodiscard]] std::size_t indexOfGCFractionalMoleculesPerComponent_CFCMC(
      std::size_t selectedComponent) const noexcept;
  [[nodiscard]] std::size_t indexOfPairGCFractionalMoleculesPerComponent_CFCMC(
      std::size_t selectedComponent) const noexcept;
  [[nodiscard]] std::size_t indexOfPairSwapFractionalMoleculesPerComponent_CFCMC(
      std::size_t selectedComponent) const noexcept;
  [[nodiscard]] std::size_t indexOfPairSwapCBFractionalMoleculesPerComponent_CFCMC(
      std::size_t selectedComponent) const noexcept;
  [[nodiscard]] std::size_t indexOfGibbsSwapFractionalMoleculesPerComponent_CFCMC(
      std::size_t selectedComponent) const noexcept;
  [[nodiscard]] std::size_t indexOfParallelReactionFractionalMoleculesPerComponent_CFCMC(
      std::size_t selectedComponent) const noexcept;
  [[nodiscard]] std::size_t indexOfSerialReactionFractionalMoleculesPerComponent_CFCMC(
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

  void precomputeTotalRigidEnergy() noexcept;
  void precomputeTotalGradients() noexcept;
  RunningEnergy computeTotalEnergies() noexcept;
  RunningEnergy computePolarizationEnergy() noexcept;
  RunningEnergy computeTotalGradients() noexcept;
  void computeTotalElectrostaticPotential() noexcept;
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
  std::span<const AtomDynamics> spanOfFrameworkDynamics() const;
  std::span<AtomDynamics> spanOfFrameworkDynamics();
  std::span<const AtomDynamics> spanOfMoleculeDynamics() const;
  std::span<AtomDynamics> spanOfMoleculeDynamics();
  std::span<double> spanOfMoleculeElectrostaticPotential();
  std::span<double3> spanOfMoleculeElectricField();
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
  // Improving the accuracy of computing chemical potentials in CFCMC simulations
  // A. Rahbari, R. Hens, D. Dubbeldam, and T.J.H Vlugt
  // Mol. Phys.  117(23-24), 3493-3508, 2019
  double weight() const
  {
    double w = std::transform_reduce(
        components.begin(), components.end(), 1.0, [](const double& acc, const double& b) { return acc * b; },
        [](const Component& component) { return component.lambdaGC.weight() * component.lambdaGibbs.weight(); });

    if (usesReactionConventionalCFCMC())
    {
      for (const Reaction& reaction : reactions.list)
      {
        w *= activeReactionLambdaHistogram(reaction).weight();
      }
    }

    return w;
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
  void syncReactionLambdaBin(Reaction& reaction) noexcept;
  void syncReactionLambdaBins() noexcept;
  [[nodiscard]] bool usesReactionConventionalCFCMC() const noexcept;
  [[nodiscard]] bool usesSerialReactionCFCMC() const noexcept;
  [[nodiscard]] bool usesParallelReactionCFCMC() const noexcept;
  [[nodiscard]] bool usesGibbsConventionalCFCMC() const noexcept;
  void initializeGibbsConventionalFractionalMolecules() noexcept;
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
  void reactionLambdaSampleProductionHistograms(std::size_t blockIndex, double weight) noexcept;
  void rescaleMolarFractions();
  void computeComponentFluidProperties();
  void computeNumberOfPseudoAtoms();
  /// Rebuilds effectiveNumberOfPseudoAtomsVDW and fractionalPseudoAtomCountsPerGroup from the current molecule atoms.
  void computeTailCorrectionCounts();
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
  std::string writeEquilibrationStatusReportMD(std::size_t currentCycle, std::size_t numberOfProductionCycles) const;
  std::string writeProductionStatusReportMC(const std::string& statusLine) const;
  std::string writeProductionStatusReportMD(std::size_t currentCycle, std::size_t numberOfProductionCycles) const;
  std::string writeSystemStatus() const;
  std::string writeComponentStatus() const;
  std::string writeMCMoveStatistics() const;

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

  bool insideBlockedPockets(const Component& component, std::span<const Atom> molecule_atoms) const;

  void sampleProperties(std::size_t systemId, std::size_t currentBlock, std::size_t currentCycle);
  void samplePropertiesEvolution(std::size_t absoluteCurrentCycle);

  void writeCPUTimeStatistics(std::ostream& stream) const;

  /**
   * Molecular (center-of-mass) excess pressure: fills site gradients from framework + intermolecular +
   * Ewald interactions (no intramolecular), then applies the atomic-to-molecular virial correction.
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

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const System& s);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, System& s);

  void writeRestartFile(std::size_t systemId);

  std::string repr() const;
};
