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

import double3;
import double3x3;
import atom;
import component;
import simulationbox;
import forcefield;
import averages;
import running_energy;
import energy_status;
import loadings;
import sample_movies;
import enthalpy_of_adsorption;
import lambda;
import cbmc;
import randomnumbers;
import property_simulationbox;
import property_energy;
import property_pressure;
import property_loading;
import property_enthalpy;
import threadpool;

export struct System
{
    System() = delete;
    System(size_t s, ForceField forcefield, std::vector<Component> components, std::vector<size_t> initialNumberOfMolecules, size_t numberOfBlocks);
    ~System() = default;

    System(const System &s) = delete;
    System(System&& s) noexcept;

    void addComponent(const Component&& component) noexcept(false);
    void createInitialMolecules();

    void MD_Loop();

    void computeTotalEnergies() noexcept;
    void computeTotalGradients() noexcept;

    void computeFrameworkMoleculeEnergy(const SimulationBox &box, std::span<const Atom> frameworkAtomPositions, std::span<const Atom> moleculeAtomPositions, RunningEnergy &energyStatus) noexcept;
    void computeInterMolecularEnergy(const SimulationBox &box, std::span<const Atom> moleculeAtomPositions, RunningEnergy &energyStatus) noexcept;
    void computeTailCorrectionVDWEnergy(RunningEnergy &energyStatus) noexcept;
    void computeEwaldFourierEnergy(const SimulationBox &box, RunningEnergy &energyStatus);

    [[nodiscard]] std::optional<RunningEnergy> computeFrameworkMoleculeEnergyDifference(std::span<const Atom> newatoms, std::span<const Atom> oldatoms) const noexcept;
    [[nodiscard]] std::optional<RunningEnergy> computeInterMolecularEnergyDifference(std::span<const Atom> newatoms, std::span<const Atom> oldatoms) const noexcept;


    [[nodiscard]] std::optional<RunningEnergy> computeFrameworkMoleculeEnergy(std::span<Atom> atoms, std::make_signed_t<std::size_t> skip = -1) const noexcept;
    template <ThreadPool::ThreadingType T>
    [[nodiscard]] std::optional<RunningEnergy> computeFrameworkMoleculeEnergy(std::span<Atom> atoms, std::make_signed_t<std::size_t> skip = -1) const noexcept;

    RunningEnergy computeFrameworkSpanMoleculeEnergy(std::span<const Atom>::iterator startIterator, std::span<const Atom>::iterator endIterator, std::span<Atom> atoms, std::make_signed_t<std::size_t> skip, std::atomic_flag & cancel) const;

    [[nodiscard]] std::optional<RunningEnergy> computeInterMolecularEnergy(std::span<Atom> atoms, std::make_signed_t<std::size_t> skip = -1) const noexcept;
    [[nodiscard]] std::optional<EnergyStatus> computeFrameworkMoleculeEnergy(std::span<Atom> atoms, std::make_signed_t<std::size_t> skip, std::span<const Atom>::iterator startIterator, std::span<const Atom>::iterator endIterator) const noexcept;
    [[nodiscard]] std::optional<EnergyStatus> computeInterMolecularSpanEnergy(std::span<const Atom>::iterator startIterator, std::span<const Atom>::iterator endIterator, std::span<Atom> atoms, std::make_signed_t<std::size_t> skip) const noexcept;

    [[nodiscard]] EnergyStatus computeTailCorrectionVDWOldEnergy() const noexcept;
    [[nodiscard]] EnergyStatus computeTailCorrectionVDWAddEnergy(size_t selectedComponent) const noexcept;
    [[nodiscard]] EnergyStatus computeTailCorrectionVDWRemoveEnergy(size_t selectedComponent) const noexcept;

    [[nodiscard]] const std::vector<std::pair<Atom, RunningEnergy>> computeExternalNonOverlappingEnergies(std::vector<Atom>& trialPositions) const noexcept;
    [[nodiscard]] const std::vector<std::pair<std::vector<Atom>,RunningEnergy>> computeExternalNonOverlappingEnergies(std::vector<std::vector<Atom>>& trialPositions, std::make_signed_t<std::size_t> skip) const noexcept;

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
    std::span<const Atom> spanOfMoleculeAtoms() const;
    std::span<Atom> spanOfMoleculeAtoms();

    size_t numberOfMolecules() const {
        return std::reduce(numberOfMoleculesPerComponent.begin(), numberOfMoleculesPerComponent.end(), size_t(0),
            [](const size_t& acc, const size_t& b) { return acc + b; });
    }

    double weight() const 
    {
      return std::transform_reduce(components.begin(), components.end(), 0.0,
               [](const double& acc, const double& b) { return acc + b; },
               [](const Component& component) { return component.lambda.weight();});
    }

    double HeliumVoidFraction{ 0.29 };
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

    void removeRedundantMoves();
    void rescaleMoveProbabilities();
    void determineSwapableComponents();
    void determineFractionalComponents();
    void rescaleMolarFractions();
    void computeComponentFluidProperties();
    void computeFrameworkDensity();
    void computeNumberOfPseudoAtoms();

    void createOutputFile();
    void closeOutputFile();
    void writeToOutputFile(const std::string &s);
    void writeOutputHeader();
    void writeOutputHeaderHardware();
    void writeInitializationStatusReport(size_t currentCycle, size_t numberOfCycles);
    void writeEquilibrationStatusReport(size_t currentCycle, size_t numberOfCycles);
    void writeProductionStatusReport(size_t currentCycle, size_t numberOfCycles);
    void writeComponentStatus();

    void writeMCMoveStatistics();
    std::string writeEnergyAveragesStatistics() const;
    std::string writeEnthalpyOfAdsorption() const;
    std::string writePressureAveragesStatistics() const;

    size_t systemId{};

    size_t numberOfFrameworks{ 0 };
    size_t numberOfFrameworkAtoms{ 0 };

    std::vector<Component> components;

    Loadings loadings;
    PropertyLoading averageLoadings;

    PropertyEnthalpy averageEnthalpiesOfAdsorption;

    std::vector<size_t> swapableComponents{};
    std::vector<size_t> numberOfMoleculesPerComponent{};                    // includes all molecules
    std::vector<size_t> numberOfIntegerMoleculesPerComponent{};             // integer molecules
    std::vector<size_t> numberOfFractionalMoleculesPerComponent{};          // fractional molecules
    std::vector<size_t> numberOfReactionMoleculesPerComponent{};            // reaction molecules
    std::vector<size_t> numberOfReactionFractionalMoleculesPerComponent{};  // reaction fractional molecules
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
    std::vector<double3> atomVelocities;
    std::vector<double3> atomForces;

    void insertMolecule(size_t selectedComponent, std::vector<Atom> atoms);
    void insertFractionalMolecule(size_t selectedComponent, std::vector<Atom> atoms);
    void deleteMolecule(size_t selectedComponent, size_t selectedMolecule, const std::span<Atom>& molecule);
    bool checkMoleculeIds();
    
    RunningEnergy runningEnergies;
    PropertyEnergy averageEnergies;

    double3x3 currentExcessPressureTensor;
    EnergyStatus currentEnergyStatus;
    PropertyPressure averagePressure;

    SampleMovie sampleMovie;
    
    [[nodiscard]] std::optional<ChainData> growMoleculeSwapInsertion(size_t selectedComponent, size_t selectedMolecule, double scaling) const noexcept;
    [[nodiscard]] std::optional<FirstBeadData> growMultipleFirstBeadSwapInsertion(const Atom& atom) const noexcept;
    [[nodiscard]] std::optional<ChainData> growChain(size_t startingBead, std::vector<Atom> atoms) const noexcept;

    [[nodiscard]] ChainData retraceMoleculeSwapDeletion(size_t selectedComponent, size_t selectedMolecule, std::span<Atom> atoms, double scaling, double storedR) const noexcept;
    [[nodiscard]] FirstBeadData retraceMultipleFirstBeadSwapDeletion(const Atom& atom, double scaling, double storedR) const noexcept;
    [[nodiscard]] ChainData retraceChain(size_t startingBead, double scaling, std::span<Atom> molecule) const noexcept;
    [[nodiscard]] ChainData retraceChainReinsertion(size_t startingBead, std::span<Atom> molecule) const noexcept;

    [[nodiscard]] std::optional<ChainData> growMoleculeReinsertion(size_t selectedComponent, size_t selectedMolecule, std::span<Atom> molecule) const noexcept;
    [[nodiscard]] ChainData retraceMoleculeReinsertion(size_t selectedComponent, size_t selectedMolecule, std::span<Atom> atoms, double storedR) const noexcept;
    [[nodiscard]] std::optional<FirstBeadData> growMultipleFirstBeadReinsertion(const Atom& atom) const noexcept;
    [[nodiscard]] FirstBeadData retraceMultipleFirstBeadReinsertion(const Atom& atom, double storedR) const noexcept;

    size_t selectTrialPosition(std::vector <double> BoltzmannFactors) const noexcept;

    size_t numberOfTrialDirections{ 8 };
    double minimumRosenbluthFactor{ 1e-150 };

    std::size_t kx_max_unsigned{8};
    std::size_t ky_max_unsigned{8};
    std::size_t kz_max_unsigned{8};
    size_t numberOfWavevectors;
    std::make_signed_t<std::size_t> kx_max{8};
    std::make_signed_t<std::size_t> ky_max{8};
    std::make_signed_t<std::size_t> kz_max{8};
    std::vector<std::complex<double>> eik_xy;
    std::vector<std::complex<double>> eik_x;
    std::vector<std::complex<double>> eik_y;
    std::vector<std::complex<double>> eik_z;
    std::vector<std::complex<double>> storedEik;
    std::vector<std::complex<double>> fixedFrameworkStoredEik;
    std::vector<std::complex<double>> totalEik;
    double CoulombicFourierEnergySingleIon{ 0.0 };
    std::vector<int> netCharge;

    bool noCharges{ false };
    bool omitEwaldFourier { false};

    RunningEnergy energyDifferenceEwaldFourier(std::vector<std::complex<double>> &storedWavevectors, 
                                               std::span<const Atom> newatoms, std::span<const Atom> oldatoms);
    void registerEwaldFourierEnergySingleIon(double3 position, double charge);
    void acceptEwaldMove();

    void sampleProperties(size_t currentBlock);
    std::chrono::duration<double> cpuTime_Sampling{ 0.0 };
    std::chrono::duration<double> cpuTime_Pressure{ 0.0 };
    const std::string writeCPUTimeStatistics() const;

    [[nodiscard]] std::pair<EnergyStatus, double3x3> computeFrameworkMoleculeEnergyStrainDerivative() noexcept;
    [[nodiscard]] std::pair<EnergyStatus, double3x3> computeInterMolecularEnergyStrainDerivative() noexcept;
    [[nodiscard]] std::pair<EnergyStatus, double3x3> computeEwaldFourierEnergyStrainDerivative() noexcept;
    [[nodiscard]] std::pair<EnergyStatus, double3x3> computeMolecularPressure() noexcept;

    double probabilityVolumeMove{ 0.0 };
    double probabilityGibbsVolumeMove{ 0.0 };
    double probabilityGibbsSwapMove_CBMC{ 0.0 };
    double probabilityGibbsSwapMove_CFCMC{ 0.0 };
    double probabilityGibbsSwapMove_CFCMC_CBMC{ 0.0 };

    void clearMoveStatistics();
    void clearTimingStatistics();

    MoveStatistics<double> statistics_VolumeMove{.maxChange = 0.5};
    MoveStatistics<double> statistics_GibbsVolumeMove{};
    MoveStatistics<double3> statistics_GibbsSwapMove_CBMC{};
    MoveStatistics<double3> statistics_GibbsSwapMove_CFCMC{};
    MoveStatistics<double3> statistics_GibbsSwapMove_CFCMC_CBMC{};

    std::chrono::duration<double> cpuTime_VolumeMove{ 0.0 };
    std::chrono::duration<double> cpuTime_GibbsVolumeMove{ 0.0 };
    std::chrono::duration<double> cpuTime_GibbsSwapMove_CBMC{ 0.0 };
    std::chrono::duration<double> cpuTime_GibbsSwapLambdaMove_CFCMC{ 0.0 };
    std::chrono::duration<double> cpuTime_GibbsSwapLambdaMove_CFCMC_CBMC{ 0.0 };

    std::chrono::duration<double> cpuTime_VolumeMove_NonEwald{ 0.0 };
    std::chrono::duration<double> cpuTime_GibbsVolumeMove_NonEwald{ 0.0 };
    std::chrono::duration<double> cpuTime_GibbsSwapMove_CBMC_NonEwald{ 0.0 };
    std::chrono::duration<double> cpuTime_GibbsSwapLambdaMove_CFCMC_NonEwald{ 0.0 };
    std::chrono::duration<double> cpuTime_GibbsSwapLambdaMove_CFCMC_CBMC_NonEwald{ 0.0 };

    std::chrono::duration<double> cpuTime_VolumeMove_Ewald{ 0.0 };
    std::chrono::duration<double> cpuTime_GibbsVolumeMove_Ewald{ 0.0 };
    std::chrono::duration<double> cpuTime_GibbsSwapMove_CBMC_Ewald{ 0.0 };
    std::chrono::duration<double> cpuTime_GibbsSwapLambdaMove_CFCMC_Ewald{ 0.0 };
    std::chrono::duration<double> cpuTime_GibbsSwapLambdaMove_CFCMC_CBMC_Ewald{ 0.0 };

    std::vector<Atom> scaledCenterOfMassPositions(double scale) const;

    std::ofstream outputFile{};
};
