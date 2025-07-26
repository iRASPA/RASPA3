module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <chrono>
#include <complex>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <ostream>
#include <span>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#endif

export module system;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;
import double3;
import double3x3;
import randomnumbers;
import threadpool;

import atom;
import molecule;
import framework;
import component;
import framework;
import simulationbox;
import forcefield;
import averages;
import running_energy;
import energy_status;
import energy_factor;
import gradient_factor;
import loadings;
import sample_movies;
import enthalpy_of_adsorption;
import property_lambda_probability_histogram;
import property_simulationbox;
import property_energy;
import property_pressure;
import property_loading;
import property_enthalpy;
import property_conventional_rdf;
import property_rdf;
import property_density_grid;
import property_temperature;
import property_energy_histogram;
import property_number_of_molecules_histogram;
import property_msd;
import property_vacf;
import multi_site_isotherm;
import pressure_range;
import units;
import bond_potential;
import isotherm;
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
import json;
import interpolation_energy_grid;

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
  System() = default;

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
  System(size_t id, ForceField forcefield, std::optional<SimulationBox> box, double T, std::optional<double> P,
         double heliumVoidFraction, std::optional<Framework> framework, std::vector<Component> components,
         std::vector<size_t> initialNumberOfMolecules, size_t numberOfBlocks,
         const MCMoveProbabilities &systemProbabilities = MCMoveProbabilities(),
         std::optional<size_t> sampleMoviesEvery = std::nullopt);

  System(size_t id, double T, std::optional<double> P, double heliumVoidFraction, std::optional<Framework> framework,
         std::vector<Component> components);

  uint64_t versionNumber{1};

  size_t systemId{};

  double temperature{300.0};
  double pressure{1e4};
  double input_pressure{1e4};
  double beta{1.0 / (Units::KB * 300.0)};

  double heliumVoidFraction{0.29};

  size_t numberOfFrameworks{0};
  size_t numberOfFrameworkAtoms{0};
  size_t numberOfRigidFrameworkAtoms{0};

  std::optional<Framework> framework;
  std::vector<Component> components;

  EquationOfState equationOfState;
  std::optional<Thermostat> thermostat;

  Loadings loadings;

  std::vector<size_t> swappableComponents{};
  std::vector<size_t> initialNumberOfMolecules{};

  // total # of molecules per component (include fractional molecules)
  std::vector<size_t> numberOfMoleculesPerComponent{};

  // # integer molecules
  std::vector<size_t> numberOfIntegerMoleculesPerComponent{};

  // # fractional molecules
  std::vector<size_t> numberOfFractionalMoleculesPerComponent{};

  // # fractional molecules for CFCMC grand-canonical or CFCMC Widom
  std::vector<size_t> numberOfGCFractionalMoleculesPerComponent_CFCMC{};

  // # fractional molecules for pair-CFCMC grand-canonical
  std::vector<size_t> numberOfPairGCFractionalMoleculesPerComponent_CFCMC{};

  // # fractional molecules for CFCMC-Gibbs
  std::vector<size_t> numberOfGibbsFractionalMoleculesPerComponent_CFCMC{};

  // # reactant fractional molecules for all reactions using CFCMC
  std::vector<std::vector<size_t>> numberOfReactionFractionalMoleculesPerComponent_CFCMC{};

  std::vector<double> idealGasEnergiesPerComponent{};

  ForceField forceField;
  bool hasExternalField;

  std::vector<std::vector<size_t>> numberOfPseudoAtoms;
  std::vector<size_t> totalNumberOfPseudoAtoms;

  // vector of pair of lamba and pseudoatom-type
  std::vector<std::pair<double, size_t>> numberOfIntegerPseudoAtoms;
  std::vector<std::pair<double, size_t>> numberOfFractionalPseudoAtoms;

  size_t translationalCenterOfMassConstraint{};
  size_t translationalDegreesOfFreedom{};
  size_t rotationalDegreesOfFreedom{};

  std::optional<double> frameworkMass() const;

  double timeStep{0.0005};

  SimulationBox simulationBox;

  // A contiguous list of adsorbate atoms per component for easy and fast looping
  // The atoms-order is defined as increasing per component and molecule.
  // Because the number of atoms is fixed per component it is easy to access the n-th molecule
  std::vector<Atom> atomPositions;
  std::vector<Molecule> moleculePositions;
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

  size_t numberOfTrialDirections{10};
  size_t numberOfHybridMCSteps{10};

  std::vector<std::complex<double>> eik_xy{};
  std::vector<std::complex<double>> eik_x{};
  std::vector<std::complex<double>> eik_y{};
  std::vector<std::complex<double>> eik_z{};
  std::vector<std::pair<std::complex<double>, std::complex<double>>> storedEik{};
  std::vector<std::pair<std::complex<double>, std::complex<double>>> fixedFrameworkStoredEik{};
  std::vector<std::pair<std::complex<double>, std::complex<double>>> totalEik{};
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

  // Breakthrough settings
  size_t columnNumberOfGridPoints{100};
  double columnTotalPressure{1e5};
  double columnPressureGradient{0.0};
  double columnVoidFraction{0.4};
  double columnParticleDensity{1000};
  double columnEntranceVelocity{0.1};
  double columnLength{0.3};
  double columnTimeStep{0.0005};
  size_t columnNumberOfTimeSteps{0};
  bool columnAutoNumberOfTimeSteps{true};
  MultiSiteIsotherm::PredictionMethod mixturePredictionMethod{MultiSiteIsotherm::PredictionMethod::IAST};
  PressureRange pressure_range;
  size_t numberOfCarrierGases{0};
  size_t carrierGasComponent{0};
  size_t maxIsothermTerms{0};

  bool containsTheFractionalMolecule{true};

  // property measurements
  PropertyEnergy averageEnergies;
  PropertyLoading averageLoadings;
  PropertyEnthalpy averageEnthalpiesOfAdsorption;
  PropertyTemperature averageTemperature;
  PropertyTemperature averageTranslationalTemperature;
  PropertyTemperature averageRotationalTemperature;
  PropertyPressure averagePressure;
  PropertySimulationBox averageSimulationBox;
  std::optional<SampleMovie> samplePDBMovie;
  std::optional<PropertyConventionalRadialDistributionFunction> propertyConventionalRadialDistributionFunction;
  std::optional<PropertyRadialDistributionFunction> propertyRadialDistributionFunction;
  std::optional<PropertyDensityGrid> propertyDensityGrid;
  std::optional<PropertyEnergyHistogram> averageEnergyHistogram;
  std::optional<PropertyNumberOfMoleculesHistogram> averageNumberOfMoleculesHistogram;
  std::optional<PropertyMeanSquaredDisplacement> propertyMSD;
  std::optional<PropertyVelocityAutoCorrelationFunction> propertyVACF;

  std::vector<std::optional<InterpolationEnergyGrid>> interpolationGrids;

  /// The fractional molecule for grand-canonical is stored first
  inline size_t indexOfGCFractionalMoleculesPerComponent_CFCMC([[maybe_unused]] size_t selectedComponent) { return 0; }

  /// The fractional molecule for grand-canonical pair-insertion is stored second
  inline size_t indexOfPairGCFractionalMoleculesPerComponent_CFCMC(size_t selectedComponent)
  {
    return numberOfGCFractionalMoleculesPerComponent_CFCMC[selectedComponent];
  }

  /// The fractional molecule for Gibbs is stored third
  inline size_t indexOfGibbsFractionalMoleculesPerComponent_CFCMC(size_t selectedComponent)
  {
    return numberOfGCFractionalMoleculesPerComponent_CFCMC[selectedComponent] +
           numberOfPairGCFractionalMoleculesPerComponent_CFCMC[selectedComponent];
  }

  void addComponent(const Component &&component) noexcept(false);

  void createFrameworks();
  void createInitialMolecules();

  void checkCartesianPositions();

  void precomputeTotalRigidEnergy() noexcept;
  void precomputeTotalGradients() noexcept;
  RunningEnergy computeTotalEnergies() noexcept;
  RunningEnergy computePolarizationEnergy() noexcept;
  RunningEnergy computeTotalGradients() noexcept;
  void computeTotalElectrostaticPotential() noexcept;
  void computeTotalElectricField() noexcept;

  size_t randomFramework(RandomNumber &random)
  {
    return size_t(random.uniform() * static_cast<double>(numberOfFrameworks));
  }
  size_t randomComponent(RandomNumber &random)
  {
    return size_t(random.uniform() * static_cast<double>(components.size()));
  }
  size_t numerOfAdsorbateComponents() { return components.size(); }
  size_t randomMoleculeOfComponent(RandomNumber &random, size_t selectedComponent);
  size_t randomIntegerMoleculeOfComponent(RandomNumber &random, size_t selectedComponent);

  size_t indexOfFirstMolecule(size_t selectedComponent);
  std::vector<Atom>::iterator iteratorForMolecule(size_t selectedComponent, size_t selectedMolecule);
  std::vector<double3>::iterator iteratorForElectricField(size_t selectedComponent, size_t selectedMolecule);
  std::vector<Molecule>::iterator indexForMolecule(size_t selectedComponent, size_t selectedMolecule);
  size_t moleculeIndexOfComponent(size_t selectedComponent, size_t selectedMolecule);
  std::span<Atom> spanOfMolecule(size_t selectedComponent, size_t selectedMolecule);
  const std::span<const Atom> spanOfMolecule(size_t selectedComponent, size_t selectedMolecule) const;
  std::span<const Atom> spanOfFrameworkAtoms() const;
  std::span<Atom> spanOfFrameworkAtoms();
  std::span<const Atom> spanOfRigidFrameworkAtoms() const;
  std::span<const Atom> spanOfFlexibleAtoms() const;
  std::span<const Atom> spanOfMoleculeAtoms() const;
  std::span<Atom> spanOfMoleculeAtoms();
  std::span<double> spanOfMoleculeElectrostaticPotential();
  std::span<double3> spanOfMoleculeElectricField();
  std::span<double3> spanOfMoleculeElectricFieldNew();
  std::span<double3> spanElectricFieldNew(size_t selectedComponent, size_t selectedMolecule);
  const std::span<const double3> spanElectricFieldNew(size_t selectedComponent, size_t selectedMolecule) const;
  std::span<double3> spanElectricFieldOld(size_t selectedComponent, size_t selectedMolecule);
  const std::span<const double3> spanElectricFieldOld(size_t selectedComponent, size_t selectedMolecule) const;

  size_t numberOfMolecules() const
  {
    return std::accumulate(numberOfMoleculesPerComponent.begin(), numberOfMoleculesPerComponent.end(), size_t(0),
                           [](const size_t &acc, const size_t &b) { return acc + b; });
  }

  size_t numberOfIntegerMolecules() const
  {
    return std::accumulate(numberOfIntegerMoleculesPerComponent.begin(), numberOfIntegerMoleculesPerComponent.end(),
                           size_t(0), [](const size_t &acc, const size_t &b) { return acc + b; });
  }

  // The system weight is the sum of the weights of all the components
  // Improving the accuracy of computing chemical potentials in CFCMC simulations
  // A. Rahbari, R. Hens, D. Dubbeldam, and T.J.H Vlugt
  // Mol. Phys.  117(23-24), 3493-3508, 2019
  double weight() const
  {
    return std::transform_reduce(
        components.begin(), components.end(), 0.0, [](const double &acc, const double &b) { return acc + b; },
        [](const Component &component) { return component.lambdaGC.weight() + component.lambdaGibbs.weight(); });
  }

  void removeRedundantMoves();
  void rescaleMoveProbabilities();
  void determineSwappableComponents();
  void determineFractionalComponents();
  void rescaleMolarFractions();
  void computeComponentFluidProperties();
  void computeNumberOfPseudoAtoms();
  void optimizeMCMoves();

  std::string writeOutputHeader() const;
  std::string writeNumberOfPseudoAtoms() const;
  std::string writeInitializationStatusReport(size_t currentCycle, size_t numberOfCycles) const;
  std::string writeEquilibrationStatusReportMC(size_t currentCycle, size_t numberOfCycles) const;
  std::string writeEquilibrationStatusReportMD(size_t currentCycle, size_t numberOfCycles) const;
  std::string writeProductionStatusReportMC(size_t currentCycle, size_t numberOfCycles) const;
  std::string writeProductionStatusReportMD(size_t currentCycle, size_t numberOfCycles) const;
  std::string writeSystemStatus() const;
  std::string writeComponentStatus() const;
  std::string writeMCMoveStatistics() const;

  nlohmann::json jsonSystemStatus() const;
  nlohmann::json jsonComponentStatus() const;
  nlohmann::json jsonMCMoveStatistics() const;

  void insertMolecule(size_t selectedComponent, const Molecule &molecule, std::vector<Atom> atoms);
  void insertMoleculePolarization(size_t selectedComponent, const Molecule &molecule, std::vector<Atom> atoms,
                                  std::span<double3> electricField);
  void insertFractionalMolecule(size_t selectedComponent, const Molecule &molecule, std::vector<Atom> atoms,
                                size_t moleculeId);
  void deleteMolecule(size_t selectedComponent, size_t selectedMolecule, const std::span<Atom> atoms);
  void checkMoleculeIds();

  std::vector<Atom> randomConfiguration(RandomNumber &random, size_t selectedComponent,
                                        const std::span<const Atom> atoms);

  bool insideBlockedPockets(const Component &component, std::span<const Atom> molecule_atoms) const;

  void sampleProperties(size_t currentBlock, size_t currentCycle);

  void writeCPUTimeStatistics(std::ostream &stream) const;

  [[nodiscard]] std::pair<EnergyStatus, double3x3> computeMolecularPressure() noexcept;

  std::pair<std::vector<Molecule>, std::vector<Atom>> scaledCenterOfMassPositions(double scale) const;

  void writeComponentFittingStatus(std::ostream &stream, const std::vector<std::pair<double, double>> &rawData) const;

  void createInterpolationGrids(std::ostream &stream);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const System &s);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, System &s);

  void writeRestartFile();
  void readRestartFile();

  std::string repr() const;
};
