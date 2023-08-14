export module system;

import <complex>;
import <vector>;
import <tuple>;
import <optional>;
import <span>;
import <fstream>;
import <ostream>;
import <iostream>;
import <numeric>;
import <chrono>;
import <algorithm>;

import double3;
import double3x3;

import randomnumbers;
import threadpool;


import atom;
import component;
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
import cbmc;
import property_simulationbox;
import property_energy;
import property_pressure;
import property_loading;
import property_enthalpy;
import property_rdf;
import multi_site_isotherm;
import pressure_range;
import units;
import bond_potential;
import isotherm;
import move_statistics;
import mc_moves_probabilities_system;
import mc_moves_cputime;
import mc_moves_count;
import reaction;
import reactions;
import transition_matrix;

export struct System
{
    System() = delete;
    System(size_t id, double T, double P, ForceField forcefield, std::vector<Component> components, std::vector<size_t> initialNumberOfMolecules, size_t numberOfBlocks);
    System(size_t s, ForceField forcefield, std::vector<Component> components, std::vector<size_t> initialNumberOfMolecules, size_t numberOfBlocks);
    ~System() = default;

    System(const System &s) = delete;
    System(const System&& s) noexcept;
    System(System&& s) noexcept;

   
    enum class FluidState : int
    {
      Unknown = 0,
      SuperCriticalFluid = 1,
      Vapor = 2,
      Liquid = 3,
      VaporLiquid = 4
    };
    FluidState fluidState = FluidState::Unknown;
    enum class EquationOfState : int
    {
      PengRobinson = 0,
      PengRobinsonGasem = 1,
      SoaveRedlichKwong = 2
    };
    EquationOfState equationOfState = EquationOfState::PengRobinson;

    enum class MultiComponentMixingRules : int
    {
      VanDerWaals = 0
    };
    MultiComponentMixingRules multiComponentMixingRules = MultiComponentMixingRules::VanDerWaals;

    size_t systemId{};

    double temperature{ 300.0 };
    double pressure{ 1e4 };
    double input_pressure{ 1e4 };
    double beta{ 1.0 / (Units::KB * 300.0) };

    double HeliumVoidFraction{ 0.29 };

    size_t numberOfFrameworks{ 0 };
    size_t numberOfFrameworkAtoms{ 0 };
    size_t numberOfRigidFrameworkAtoms{ 0 };

    std::vector<Component> components;

    Loadings loadings;
    PropertyLoading averageLoadings;

    PropertyEnthalpy averageEnthalpiesOfAdsorption;

    std::vector<size_t> swapableComponents{};
    std::vector<size_t> initialNumberOfMolecules{};

    std::vector<size_t> numberOfMoleculesPerComponent{};                                            // total # of molecules per component (include fractional molecules)
                                                                                                    
    std::vector<size_t> numberOfIntegerMoleculesPerComponent{};                                     // # integer molecules
    std::vector<size_t> numberOfFractionalMoleculesPerComponent{};                                  // # fractional molecules
                                                                                                    
    std::vector<size_t> numberOfGCFractionalMoleculesPerComponent_CFCMC{};                          // # fractional molecules for CFCMC grand-canonical or CFCMC Widom
    std::vector<size_t> numberOfPairGCFractionalMoleculesPerComponent_CFCMC{};                      // # fractional molecules for pair-CFCMC grand-canonical
    std::vector<size_t> numberOfGibbsFractionalMoleculesPerComponent_CFCMC{};                       // # fractional molecules for CFCMC-Gibbs
    std::vector<std::vector<size_t>> numberOfReactionFractionalMoleculesPerComponent_CFCMC{};       // # reactant fractional molecules for all reactions using CFCMC

    /// The fractional molecule for grand-canonical is stored first
    inline size_t indexOfGCFractionalMoleculesPerComponent_CFCMC([[maybe_unused]] size_t selectedComponent) { return 0;}

    /// The fractional molecule for grand-canonical pair-insertion is stored second
    inline size_t indexOfPairGCFractionalMoleculesPerComponent_CFCMC(size_t selectedComponent) 
    { return numberOfGCFractionalMoleculesPerComponent_CFCMC[selectedComponent]; }
  
    /// The fractional molecule for Gibbs is stored third 
    inline size_t indexOfGibbsFractionalMoleculesPerComponent_CFCMC(size_t selectedComponent) 
    { return numberOfGCFractionalMoleculesPerComponent_CFCMC[selectedComponent] +
             numberOfPairGCFractionalMoleculesPerComponent_CFCMC[selectedComponent]; }
   

    std::vector<double> idealGasEnergiesPerComponent{};

    ForceField forceField;

    std::vector<std::vector<size_t>> numberOfPseudoAtoms;
    std::vector<size_t> totalNumberOfPseudoAtoms;

    std::optional<double> frameworkMass{ std::nullopt };

    double timeStep{ 0.0005 };

    SimulationBox simulationBox;
    PropertySimulationBox averageSimulationBox;

    // A contiguous list of adsorbate atoms per component for easy and fast looping
    // The atoms-order is defined as increasing per component and molecule.
    // Because the number of atoms is fixed per component it is easy to access the n-th molecule
    std::vector<Atom> atomPositions;

    RunningEnergy runningEnergies;
    RunningEnergy rigidEnergies;
    PropertyEnergy averageEnergies;

    double3x3 currentExcessPressureTensor;
    EnergyStatus currentEnergyStatus;
    PropertyPressure averagePressure;

    size_t numberOfTrialDirections{ 10 };
    double minimumRosenbluthFactor{ 1e-150 };

    std::vector<std::complex<double>> eik_xy;
    std::vector<std::complex<double>> eik_x;
    std::vector<std::complex<double>> eik_y;
    std::vector<std::complex<double>> eik_z;
    std::vector<std::pair<std::complex<double>, std::complex<double>>> storedEik;
    std::vector<std::pair<std::complex<double>, std::complex<double>>> fixedFrameworkStoredEik;
    std::vector<std::pair<std::complex<double>, std::complex<double>>> totalEik;
    double CoulombicFourierEnergySingleIon{ 0.0 };
    std::vector<int> netCharge;

    bool noCharges{ false };
    bool omitEwaldFourier{ false };

    MCMoveProbabilitiesSystem mc_moves_probabilities;
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


    void addComponent(const Component&& component) noexcept(false);
    void initializeComponents();
    void createInitialMolecules();

    void MD_Loop();

    void precomputeTotalRigidEnergy() noexcept;
    void recomputeTotalEnergies() noexcept;
    RunningEnergy computeTotalEnergies() noexcept;
    void computeTotalGradients() noexcept;

    ForceFactor computeInterMolecularGradient() noexcept;
    ForceFactor computeFrameworkMoleculeGradient() noexcept;

    void computeFrameworkMoleculeEnergy(const SimulationBox &box, std::span<const Atom> frameworkAtomPositions, std::span<const Atom> moleculeAtomPositions, RunningEnergy &energyStatus) noexcept;
    void computeFrameworkMoleculeTailEnergy(std::span<const Atom> frameworkAtomPositions, std::span<const Atom> moleculeAtomPositions, RunningEnergy &energyStatus) noexcept;
    void computeInterMolecularEnergy(const SimulationBox &box, std::span<const Atom> moleculeAtomPositions, RunningEnergy &energyStatus) noexcept;
    void computeInterMolecularTailEnergy(std::span<const Atom> moleculeAtomPositions, RunningEnergy &energyStatus) noexcept;
    void computeTailCorrectionVDWEnergy(RunningEnergy &energyStatus) noexcept;
    void computeEwaldFourierEnergy(const SimulationBox &box, RunningEnergy &energyStatus);
    void computeEwaldFourierEnergy(const SimulationBox &box, std::span<const Atom> moleculeAtomPositions, RunningEnergy& energyStatus);
    void computeEwaldFourierRigidEnergy(const SimulationBox& box, RunningEnergy& energyStatus);
    ForceFactor computeEwaldFourierGradient();


    [[nodiscard]] std::optional<RunningEnergy> computeFrameworkMoleculeEnergyDifference(std::span<const Atom> newatoms, std::span<const Atom> oldatoms) const noexcept;
    [[nodiscard]] RunningEnergy computeFrameworkMoleculeTailEnergyDifference(std::span<const Atom> newatoms, std::span<const Atom> oldatoms) const noexcept;
    [[nodiscard]] std::optional<RunningEnergy> computeInterMolecularEnergyDifference(std::span<const Atom> newatoms, std::span<const Atom> oldatoms) const noexcept;
    [[nodiscard]] RunningEnergy computeInterMolecularTailEnergyDifference(std::span<const Atom> newatoms, std::span<const Atom> oldatoms) const noexcept;

    [[nodiscard]] std::optional<RunningEnergy> computeFrameworkMoleculeEnergy(double cutOffVDW, double cutOffCoulomb, std::span<Atom> atoms, std::make_signed_t<std::size_t> skip = -1) const noexcept;

    template <ThreadPool::ThreadingType T>
    [[nodiscard]] std::optional<RunningEnergy> computeFrameworkMoleculeEnergy(double cutOffVDW, double cutOffCoulomb, std::span<Atom> atoms, std::make_signed_t<std::size_t> skip = -1) const noexcept;


    [[nodiscard]] std::optional<RunningEnergy> computeInterMolecularEnergy(double cutOffVDW, double cutOffCoulomb, std::span<Atom> atoms, std::make_signed_t<std::size_t> skip = -1) const noexcept;
    [[nodiscard]] std::optional<EnergyStatus> computeInterMolecularSpanEnergy(std::span<const Atom>::iterator startIterator, std::span<const Atom>::iterator endIterator, std::span<Atom> atoms, std::make_signed_t<std::size_t> skip) const noexcept;

    [[nodiscard]] EnergyStatus computeTailCorrectionVDWOldEnergy() const noexcept;
    [[nodiscard]] EnergyStatus computeTailCorrectionVDWAddEnergy(size_t selectedComponent) const noexcept;
    [[nodiscard]] EnergyStatus computeTailCorrectionVDWRemoveEnergy(size_t selectedComponent) const noexcept;

    [[nodiscard]] const std::vector<std::pair<Atom, RunningEnergy>> computeExternalNonOverlappingEnergies(double cutOffVDW, double cutOffCoulomb, std::vector<Atom>& trialPositions) const noexcept;
    [[nodiscard]] const std::vector<std::pair<std::vector<Atom>,RunningEnergy>> computeExternalNonOverlappingEnergies(double cutOffVDW, double cutOffCoulomb, std::vector<std::vector<Atom>>& trialPositions, std::make_signed_t<std::size_t> skip) const noexcept;
    const std::optional<RunningEnergy> computeExternalNonOverlappingEnergyDualCutOff(double cutOffVDW, double cutOffCoulomb, std::vector<Atom>& trialPositionSet) const noexcept;

    [[nodiscard]] std::vector<double> computeInterMolecularInteractionsEnergy() const noexcept;
    [[nodiscard]] std::vector<double> computeInterMolecularInteractionsForce() noexcept;

    size_t randomFramework() { return size_t(RandomNumber::Uniform() * static_cast<double>(numberOfFrameworks)); }
    size_t randomComponent() { return size_t(RandomNumber::Uniform() * static_cast<double>((components.size() - numberOfFrameworks)) + static_cast<double>(numberOfFrameworks)); }
    size_t numerOfAdsorbateComponents() { return components.size() - numberOfFrameworks; }
    size_t randomMoleculeOfComponent(size_t selectedComponent);
    size_t randomIntegerMoleculeOfComponent(size_t selectedComponent);

    size_t indexOfFirstMolecule(size_t selectedComponent);
    std::vector<Atom>::const_iterator iteratorForMolecule(size_t selectedComponent, size_t selectedMolecule);
    std::span<Atom> spanOfMolecule(size_t selectedComponent, size_t selectedMolecule); 
    const std::span<const Atom> spanOfMolecule(size_t selectedComponent, size_t selectedMolecule) const; 
    std::span<const Atom> spanOfFrameworkAtoms() const;
    std::span<const Atom> spanOfRigidFrameworkAtoms() const;
    std::span<const Atom> spanOfFlexibleAtoms() const;
    std::span<const Atom> spanOfMoleculeAtoms() const;
    std::span<Atom> spanOfMoleculeAtoms();
    std::span<const Component> spanOfAdsorbateComponents() const {return std::span(components.cbegin() + 
          static_cast<std::vector<Component>::difference_type>(numberOfFrameworks), components.size() - numberOfFrameworks);}

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

    bool containsTheFractionalMolecule{ true };
    

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
    std::string writeOutputHeaderHardware() const;
    std::string writeInitializationStatusReport(size_t currentCycle, size_t numberOfCycles) const;
    std::string writeEquilibrationStatusReport(size_t currentCycle, size_t numberOfCycles) const;
    std::string writeProductionStatusReport(size_t currentCycle, size_t numberOfCycles) const;
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

    void insertMolecule(size_t selectedComponent, std::vector<Atom> atoms);
    void insertFractionalMolecule(size_t selectedComponent, std::vector<Atom> atoms, size_t moleculeId);
    void deleteMolecule(size_t selectedComponent, size_t selectedMolecule, const std::span<Atom>& molecule);
    bool checkMoleculeIds();
    
    //SampleMovie sampleMovie;
    
    std::vector<Atom> randomConfiguration(size_t selectedComponent, const std::span<const Atom> atoms) const;

    [[nodiscard]] std::optional<ChainData> growMoleculeSwapInsertion(double cutOffVDW, double cutOffCoulomb, size_t selectedComponent, size_t selectedMolecule, double scaling, std::vector<Atom> atoms) const noexcept;
    [[nodiscard]] std::optional<FirstBeadData> growMultipleFirstBeadSwapInsertion(double cutOffVDW, double cutOffCoulomb, const Atom& atom) const noexcept;
    [[nodiscard]] std::optional<ChainData> growChain(double cutOffVDW, double cutOffCoulomb, size_t startingBead, std::vector<Atom> atoms) const noexcept;

    [[nodiscard]] ChainData retraceMoleculeSwapDeletion(double cutOffVDW, double cutOffCoulomb, size_t selectedComponent, size_t selectedMolecule, std::span<Atom> atoms, double scaling, double storedR) const noexcept;
    [[nodiscard]] FirstBeadData retraceMultipleFirstBeadSwapDeletion(double cutOffVDW, double cutOffCoulomb, const Atom& atom, double scaling, double storedR) const noexcept;
    [[nodiscard]] ChainData retraceChain(double cutOffVDW, double cutOffCoulomb, size_t startingBead, double scaling, std::span<Atom> molecule) const noexcept;
    [[nodiscard]] ChainData retraceChainReinsertion(double cutOffVDW, double cutOffCoulomb, size_t startingBead, std::span<Atom> molecule) const noexcept;

    [[nodiscard]] std::optional<ChainData> growMoleculeReinsertion(double cutOffVDW, double cutOffCoulomb, size_t selectedComponent, size_t selectedMolecule, std::span<Atom> molecule) const noexcept;
    [[nodiscard]] ChainData retraceMoleculeReinsertion(double cutOffVDW, double cutOffCoulomb, size_t selectedComponent, size_t selectedMolecule, std::span<Atom> atoms, double storedR) const noexcept;
    [[nodiscard]] std::optional<FirstBeadData> growMultipleFirstBeadReinsertion(double cutOffVDW, double cutOffCoulomb, const Atom& atom) const noexcept;
    [[nodiscard]] FirstBeadData retraceMultipleFirstBeadReinsertion(double cutOffVDW, double cutOffCoulomb, const Atom& atom, double storedR) const noexcept;

    size_t selectTrialPosition(std::vector <double> BoltzmannFactors) const noexcept;

    RunningEnergy energyDifferenceEwaldFourier(std::vector<std::pair<std::complex<double>, std::complex<double>>> &storedWavevectors,
                                               std::span<const Atom> newatoms, std::span<const Atom> oldatoms);
    void registerEwaldFourierEnergySingleIon(double3 position, double charge);
    void acceptEwaldMove();

    void sampleProperties(size_t currentBlock);
    
    void writeCPUTimeStatistics(std::ostream &stream) const;

    [[nodiscard]] std::pair<EnergyStatus, double3x3> computeFrameworkMoleculeEnergyStrainDerivative() noexcept;
    [[nodiscard]] std::pair<EnergyStatus, double3x3> computeInterMolecularEnergyStrainDerivative() noexcept;
    [[nodiscard]] std::pair<EnergyStatus, double3x3> computeEwaldFourierEnergyStrainDerivative() noexcept;
    [[nodiscard]] std::pair<EnergyStatus, double3x3> computeMolecularPressure() noexcept;

    void clearMoveStatistics();

    std::vector<Atom> scaledCenterOfMassPositions(double scale) const;

    std::vector<Atom> equilibratedMoleculeRandomInBox(size_t selectedComponent, std::span<Atom> molecule, double scaling, size_t moleculeId) const;

    void writeComponentFittingStatus(std::ostream &stream, const std::vector<std::pair<double, double>> &rawData) const;

    
};
