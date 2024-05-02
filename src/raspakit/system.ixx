module;

#ifdef USE_LEGACY_HEADERS
#include <complex>
#include <vector>
#include <array>
#include <map>
#include <utility>
#include <tuple>
#include <optional>
#include <span>
#include <fstream>
#include <ostream>
#include <iostream>
#include <numeric>
#include <chrono>
#include <algorithm>
#include <type_traits>
#endif

export module system;

#ifndef USE_LEGACY_HEADERS
import <complex>;
import <vector>;
import <array>;
import <map>;
import <utility>;
import <tuple>;
import <optional>;
import <span>;
import <fstream>;
import <ostream>;
import <iostream>;
import <numeric>;
import <chrono>;
import <algorithm>;
import <type_traits>;
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
import force_factor;
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
import multi_site_isotherm;
import pressure_range;
import units;
import bond_potential;
import isotherm;
import move_statistics;
import mc_moves_probabilities_system;
import mc_moves_statistics_system;
import mc_moves_cputime;
import mc_moves_count;
import reaction;
import reactions;
import transition_matrix;
import equation_of_states;

export struct System
{
  System() = default;

  // construct System programmatically
/*! \brief Brief description.
 *         Brief description continued.
 *
 *  Detailed description starts here.
 */
  System(size_t id, std::optional<SimulationBox> box, double T, std::optional<double> P, ForceField forcefield, 
         std::vector<Framework> frameworkComponents, std::vector<Component> components, 
         std::vector<size_t> initialNumberOfMolecules, size_t numberOfBlocks, 
         const MCMoveProbabilitiesSystem &systemProbabilities = MCMoveProbabilitiesSystem());

  uint64_t versionNumber{ 1 };

  size_t systemId{};

  double temperature{ 300.0 };
  double pressure{ 1e4 };
  double input_pressure{ 1e4 };
  double beta{ 1.0 / (Units::KB * 300.0) };

  double HeliumVoidFraction{ 0.29 };

  size_t numberOfFrameworks{ 0 };
  size_t numberOfFrameworkAtoms{ 0 };
  size_t numberOfRigidFrameworkAtoms{ 0 };

  std::vector<Framework> frameworkComponents;
  std::vector<Component> components;

  EquationOfState equationOfState;

  Loadings loadings;
  PropertyLoading averageLoadings;

  PropertyEnthalpy averageEnthalpiesOfAdsorption;

  std::vector<size_t> swapableComponents{};
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

  std::optional<double> frameworkMass{};

  double timeStep{ 0.0005 };

  SimulationBox simulationBox;
  PropertySimulationBox averageSimulationBox;

  // A contiguous list of adsorbate atoms per component for easy and fast looping
  // The atoms-order is defined as increasing per component and molecule.
  // Because the number of atoms is fixed per component it is easy to access the n-th molecule
  std::vector<Atom> atomPositions;
  std::vector<Molecule> moleculePositions;

  RunningEnergy runningEnergies;
  PropertyEnergy averageEnergies;

  double3x3 currentExcessPressureTensor;
  EnergyStatus currentEnergyStatus;
  PropertyPressure averagePressure;

  size_t numberOfTrialDirections{ 10 };

  std::vector<std::complex<double>> eik_xy{};
  std::vector<std::complex<double>> eik_x{};
  std::vector<std::complex<double>> eik_y{};
  std::vector<std::complex<double>> eik_z{};
  std::vector<std::pair<std::complex<double>, std::complex<double>>> storedEik{};
  std::vector<std::pair<std::complex<double>, std::complex<double>>> fixedFrameworkStoredEik{};
  std::vector<std::pair<std::complex<double>, std::complex<double>>> totalEik{};
  double CoulombicFourierEnergySingleIon{ 0.0 };
  std::vector<double> netCharge;

  MCMoveProbabilitiesSystem mc_moves_probabilities;
  MCMoveStatisticsSystem mc_moves_statistics;
  MCMoveCpuTime mc_moves_cputime;
  MCMoveCount mc_moves_count;

  Reactions reactions;
  TransitionMatrix tmmc;

  // Breakthrough settings
  size_t columnNumberOfGridPoints{ 100 };
  double columnTotalPressure{ 1e5 };
  double columnPressureGradient{ 0.0 };
  double columnVoidFraction{ 0.4 };
  double columnParticleDensity{ 1000 };
  double columnEntranceVelocity{ 0.1 };
  double columnLength{ 0.3 };
  double columnTimeStep{ 0.0005 };
  size_t columnNumberOfTimeSteps{ 0 };
  bool columnAutoNumberOfTimeSteps{ true };
  MultiSiteIsotherm::PredictionMethod mixturePredictionMethod{ MultiSiteIsotherm::PredictionMethod::IAST };
  PressureRange pressure_range;
  size_t numberOfCarrierGases{ 0 };
  size_t carrierGasComponent{ 0 };
  size_t maxIsothermTerms{ 0 };

  bool containsTheFractionalMolecule{ true };

  // property measurements
  std::optional<PropertyConventionalRadialDistributionFunction> propertyConventionalRadialDistributionFunction;
  std::optional<PropertyRadialDistributionFunction> propertyRadialDistributionFunction;
  std::optional<PropertyDensityGrid> propertyDensityGrid;

  /// The fractional molecule for grand-canonical is stored first
  inline size_t indexOfGCFractionalMoleculesPerComponent_CFCMC([[maybe_unused]] size_t selectedComponent) { return 0;}

  /// The fractional molecule for grand-canonical pair-insertion is stored second
  inline size_t indexOfPairGCFractionalMoleculesPerComponent_CFCMC(size_t selectedComponent) 
  { return numberOfGCFractionalMoleculesPerComponent_CFCMC[selectedComponent]; }

  /// The fractional molecule for Gibbs is stored third 
  inline size_t indexOfGibbsFractionalMoleculesPerComponent_CFCMC(size_t selectedComponent) 
  { return numberOfGCFractionalMoleculesPerComponent_CFCMC[selectedComponent] +
           numberOfPairGCFractionalMoleculesPerComponent_CFCMC[selectedComponent]; }

  void addComponent(const Component&& component) noexcept(false);

  void createFrameworks();
  void createInitialMolecules(RandomNumber &random);
  void determineSimulationBox();

  void integrate();
  void computeCenterOfMassAndQuaternionGradients();

  void precomputeTotalRigidEnergy() noexcept;
  void recomputeTotalEnergies() noexcept;
  RunningEnergy computeTotalEnergies() noexcept;
  RunningEnergy computeTotalGradients() noexcept;

  size_t randomFramework(RandomNumber &random) { return size_t(random.uniform() * static_cast<double>(numberOfFrameworks)); }
  size_t randomComponent(RandomNumber &random) { return size_t(random.uniform() * static_cast<double>(components.size())); }
  size_t numerOfAdsorbateComponents() { return components.size(); }
  size_t randomMoleculeOfComponent(RandomNumber &random, size_t selectedComponent);
  size_t randomIntegerMoleculeOfComponent(RandomNumber &random, size_t selectedComponent);

  size_t indexOfFirstMolecule(size_t selectedComponent);
  std::vector<Atom>::iterator iteratorForMolecule(size_t selectedComponent, size_t selectedMolecule);
  std::vector<Molecule>::iterator indexForMolecule(size_t selectedComponent, size_t selectedMolecule);
  std::span<Atom> spanOfMolecule(size_t selectedComponent, size_t selectedMolecule); 
  const std::span<const Atom> spanOfMolecule(size_t selectedComponent, size_t selectedMolecule) const; 
  std::span<const Atom> spanOfFrameworkAtoms() const;
  std::span<Atom> spanOfFrameworkAtoms();
  std::span<const Atom> spanOfRigidFrameworkAtoms() const;
  std::span<const Atom> spanOfFlexibleAtoms() const;
  std::span<const Atom> spanOfMoleculeAtoms() const;
  std::span<Atom> spanOfMoleculeAtoms();

  size_t numberOfMolecules() const {
      return std::reduce(numberOfMoleculesPerComponent.begin(), numberOfMoleculesPerComponent.end(), size_t(0),
          [](const size_t& acc, const size_t& b) { return acc + b; });
  }

  size_t numberOfIntegerMolecules() const {
    return std::reduce(numberOfIntegerMoleculesPerComponent.begin(), numberOfIntegerMoleculesPerComponent.end(), size_t(0),
      [](const size_t& acc, const size_t& b) { return acc + b; });
  }

  // The system weight is the sum of the weights of all the components
  // Improving the accuracy of computing chemical potentials in CFCMC simulations 
  // A. Rahbari, R. Hens, D. Dubbeldam, and T.J.H Vlugt
  // Mol. Phys.  117(23-24), 3493-3508, 2019
  double weight() const 
  {
    return std::transform_reduce(components.begin(), components.end(), 0.0,
             [](const double& acc, const double& b) { return acc + b; },
             [](const Component& component) { return component.lambdaGC.weight() + component.lambdaGibbs.weight();});
  }

  

  void removeRedundantMoves();
  void rescaleMoveProbabilities();
  void determineSwapableComponents();
  void determineFractionalComponents();
  void rescaleMolarFractions();
  void computeComponentFluidProperties();
  void computeFrameworkDensity();
  void computeNumberOfPseudoAtoms();
  void optimizeMCMoves();

  std::string writeOutputHeader() const;
  std::string writeInitializationStatusReport(size_t currentCycle, size_t numberOfCycles) const;
  std::string writeEquilibrationStatusReport(size_t currentCycle, size_t numberOfCycles) const;
  std::string writeProductionStatusReport(size_t currentCycle, size_t numberOfCycles) const;
  std::string writeSystemStatus() const;
  std::string writeComponentStatus() const;

  std::string writeMCMoveStatistics() const;

  std::vector<Component> nonFrameworkComponents()
  {
    std::vector<Component> comps{};
    std::copy_if(components.begin(), components.end(), std::back_inserter(comps), [](const Component &c){
        return c.type != Component::Type::Framework;} 
      );
    return comps;
  }

  void insertMolecule(size_t selectedComponent, const Molecule &molecule, std::vector<Atom> atoms);
  void insertFractionalMolecule(size_t selectedComponent, const Molecule &molecule, std::vector<Atom> atoms, size_t moleculeId);
  void deleteMolecule(size_t selectedComponent, size_t selectedMolecule, const std::span<Atom> atoms);
  bool checkMoleculeIds();
  
  std::vector<Atom> randomConfiguration(RandomNumber &random, size_t selectedComponent, const std::span<const Atom> atoms);

  void sampleProperties(size_t currentBlock, size_t currentCycle);
  
  void writeCPUTimeStatistics(std::ostream &stream) const;

  [[nodiscard]] std::pair<EnergyStatus, double3x3> computeMolecularPressure() noexcept;

  void clearMoveStatistics();

  std::vector<Atom> scaledCenterOfMassPositions(double scale) const;

  std::vector<Atom> equilibratedMoleculeRandomInBox(RandomNumber &random, size_t selectedComponent, double scaling, size_t moleculeId) const;

  void writeComponentFittingStatus(std::ostream &stream, const std::vector<std::pair<double, double>> &rawData) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const System &s);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, System &s);    

  std::string repr() const;
};
