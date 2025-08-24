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
import write_lammps_data;

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
  System(std::size_t id, ForceField forcefield, std::optional<SimulationBox> box, double T, std::optional<double> P,
         double heliumVoidFraction, std::optional<Framework> framework, std::vector<Component> components,
         std::vector<std::vector<double3>> initialPositions, std::vector<std::size_t> initialNumberOfMolecules,
         std::size_t numberOfBlocks, const MCMoveProbabilities &systemProbabilities = MCMoveProbabilities(),
         std::optional<std::size_t> sampleMoviesEvery = std::nullopt);

  System(std::size_t id, double T, std::optional<double> P, double heliumVoidFraction,
         std::optional<Framework> framework, std::vector<Component> components);

  std::uint64_t versionNumber{1};

  std::size_t systemId{};

  double temperature{300.0};
  double pressure{1e4};
  double input_pressure{1e4};
  double beta{1.0 / (Units::KB * 300.0)};

  double heliumVoidFraction{0.29};

  std::size_t numberOfFrameworks{0};
  std::size_t numberOfFrameworkAtoms{0};
  std::size_t numberOfRigidFrameworkAtoms{0};

  std::optional<Framework> framework;
  std::vector<Component> components;

  EquationOfState equationOfState;
  std::optional<Thermostat> thermostat;

  Loadings loadings;

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

  // # fractional molecules for CFCMC-Gibbs
  std::vector<std::size_t> numberOfGibbsFractionalMoleculesPerComponent_CFCMC{};

  // # reactant fractional molecules for all reactions using CFCMC
  std::vector<std::vector<std::size_t>> numberOfReactionFractionalMoleculesPerComponent_CFCMC{};

  std::vector<double> idealGasEnergiesPerComponent{};

  ForceField forceField;
  bool hasExternalField;

  std::vector<std::vector<std::size_t>> numberOfPseudoAtoms;
  std::vector<std::size_t> totalNumberOfPseudoAtoms;

  // vector of pair of lamba and pseudoatom-type
  std::vector<std::pair<double, std::size_t>> numberOfIntegerPseudoAtoms;
  std::vector<std::pair<double, std::size_t>> numberOfFractionalPseudoAtoms;

  std::size_t translationalCenterOfMassConstraint{};
  std::size_t translationalDegreesOfFreedom{};
  std::size_t rotationalDegreesOfFreedom{};

  std::optional<double> frameworkMass() const;

  double timeStep{0.0005};

  SimulationBox simulationBox;

  // A contiguous list of adsorbate atoms per component for easy and fast looping
  // The atoms-order is defined as increasing per component and molecule.
  // Because the number of atoms is fixed per component it is easy to access the n-th molecule
  std::vector<Atom> atomData;
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

  std::size_t numberOfTrialDirections{10};
  std::size_t numberOfHybridMCSteps{10};

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
  MultiSiteIsotherm::PredictionMethod mixturePredictionMethod{MultiSiteIsotherm::PredictionMethod::IAST};
  PressureRange pressure_range;
  std::size_t numberOfCarrierGases{0};
  std::size_t carrierGasComponent{0};
  std::size_t maxIsothermTerms{0};

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
  std::optional<WriteLammpsData> writeLammpsData;

  std::vector<std::optional<InterpolationEnergyGrid>> interpolationGrids;

  /// The fractional molecule for grand-canonical is stored first
  inline std::size_t indexOfGCFractionalMoleculesPerComponent_CFCMC([[maybe_unused]] std::size_t selectedComponent)
  {
    return 0;
  }

  /// The fractional molecule for grand-canonical pair-insertion is stored second
  inline std::size_t indexOfPairGCFractionalMoleculesPerComponent_CFCMC(std::size_t selectedComponent)
  {
    return numberOfGCFractionalMoleculesPerComponent_CFCMC[selectedComponent];
  }

  /// The fractional molecule for Gibbs is stored third
  inline std::size_t indexOfGibbsFractionalMoleculesPerComponent_CFCMC(std::size_t selectedComponent)
  {
    return numberOfGCFractionalMoleculesPerComponent_CFCMC[selectedComponent] +
           numberOfPairGCFractionalMoleculesPerComponent_CFCMC[selectedComponent];
  }

  void addComponent(const Component &&component) noexcept(false);

  void createFrameworks();
  void createInitialMolecules(const std::vector<std::vector<double3>> &initialPositions);

  void checkCartesianPositions();

  void precomputeTotalRigidEnergy() noexcept;
  void precomputeTotalGradients() noexcept;
  RunningEnergy computeTotalEnergies() noexcept;
  RunningEnergy computePolarizationEnergy() noexcept;
  RunningEnergy computeTotalGradients() noexcept;
  void computeTotalElectrostaticPotential() noexcept;
  void computeTotalElectricField() noexcept;

  std::size_t randomFramework(RandomNumber &random)
  {
    return std::size_t(random.uniform() * static_cast<double>(numberOfFrameworks));
  }
  std::size_t randomComponent(RandomNumber &random)
  {
    return std::size_t(random.uniform() * static_cast<double>(components.size()));
  }
  std::size_t numerOfAdsorbateComponents() { return components.size(); }
  std::size_t randomMoleculeOfComponent(RandomNumber &random, std::size_t selectedComponent);
  std::size_t randomIntegerMoleculeOfComponent(RandomNumber &random, std::size_t selectedComponent);

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
                           [](const std::size_t &acc, const std::size_t &b) { return acc + b; });
  }

  std::size_t numberOfIntegerMolecules() const
  {
    return std::accumulate(numberOfIntegerMoleculesPerComponent.begin(), numberOfIntegerMoleculesPerComponent.end(),
                           std::size_t(0), [](const std::size_t &acc, const std::size_t &b) { return acc + b; });
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
  std::string writeInitializationStatusReport(std::size_t currentCycle, std::size_t numberOfCycles) const;
  std::string writeEquilibrationStatusReportMC(std::size_t currentCycle, std::size_t numberOfCycles) const;
  std::string writeEquilibrationStatusReportMD(std::size_t currentCycle, std::size_t numberOfCycles) const;
  std::string writeProductionStatusReportMC(const std::string &statusLine) const;
  std::string writeProductionStatusReportMD(std::size_t currentCycle, std::size_t numberOfCycles) const;
  std::string writeSystemStatus() const;
  std::string writeComponentStatus() const;
  std::string writeMCMoveStatistics() const;

  nlohmann::json jsonSystemStatus() const;
  nlohmann::json jsonComponentStatus() const;
  nlohmann::json jsonMCMoveStatistics() const;

  void insertMolecule(std::size_t selectedComponent, const Molecule &molecule, std::vector<Atom> atoms);
  void insertMoleculePolarization(std::size_t selectedComponent, const Molecule &molecule, std::vector<Atom> atoms,
                                  std::span<double3> electricField);
  void insertFractionalMolecule(std::size_t selectedComponent, const Molecule &molecule, std::vector<Atom> atoms,
                                std::size_t moleculeId);
  void deleteMolecule(std::size_t selectedComponent, std::size_t selectedMolecule, const std::span<Atom> atoms);
  void updateMoleculeAtomInformation();
  void checkMoleculeIds();

  std::vector<Atom> randomConfiguration(RandomNumber &random, std::size_t selectedComponent,
                                        const std::span<const Atom> atoms);

  bool insideBlockedPockets(const Component &component, std::span<const Atom> molecule_atoms) const;

  void sampleProperties(std::size_t currentBlock, std::size_t currentCycle);

  void writeCPUTimeStatistics(std::ostream &stream) const;

  [[nodiscard]] std::pair<EnergyStatus, double3x3> computeMolecularPressure() noexcept;

  std::pair<std::vector<Molecule>, std::vector<Atom>> scaledCenterOfMassPositions(double scale) const;

  void writeComponentFittingStatus(std::ostream &stream, const std::vector<std::pair<double, double>> &rawData) const;

  void createInterpolationGrids(std::ostream &stream);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const System &s);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, System &s);

  void writeRestartFile();

  std::string repr() const;
};
