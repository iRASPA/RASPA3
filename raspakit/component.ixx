export module component;

import <format>;
import <tuple>;
import <vector>;
import <string>;
import <chrono>;
import <cstdint>;
import <fstream>;
import <ostream>;
import <vector>;
import <array>;
import <map>;
import <optional>;
import <span>;

import int3;
import double3;
import averages;
import atom;
import forcefield;
import lambda;
import simulationbox;
import property_widom;
import isotherm;
import multi_site_isotherm;
import bond_potential;

export template<typename T>
 struct MoveStatistics
{
    T counts{};
    T constructed{};
    T accepted{};
    T maxChange{};

    void clear()
    {
        counts = T();
        constructed = T();
        accepted = T();
    }
};

export struct Component
{
    enum class Type : int
    {
        Framework = 0,
        Adsorbate = 1,
        Cation = 2
    };

    

    Component(Component::Type type, size_t currentComponent, const std::string &componentName, 
              std::optional<const std::string> fileName, size_t numberOfBlocks) noexcept(false);

    void readComponent(const ForceField& forceField, const std::string& fileName);
    void readFramework(const ForceField& forceField, const std::string& fileName);

    void printStatus(std::ostream &stream, const ForceField& forceField) const;
    void printBreakthroughStatus(std::ostream &stream) const;

    void normalizeMoveProbabilties();

    const std::string writeMCMoveStatistics() const;
    
    const std::string writeMCMoveCPUTimeStatistics() const;
    
    std::vector<double3> randomlyRotatedPositionsAroundStartingBead() const;
    std::vector<Atom> newAtoms(double scaling, size_t moleculeId) const;
    std::vector<Atom> copiedAtoms(std::span<Atom> molecule) const;
    
    Type type;

    SimulationBox simulationBox{30.0, 30.0 , 30.0};
    size_t spaceGroupHallNumber{ 1 };
    int3 numberOfUnitCells{ 1, 1, 1};

    size_t componentId{ 0 };
    std::string name{};
    std::optional<std::string> filenameData{};
    std::string filename{};

    bool rigid { true };

    double criticalTemperature{ 0.0 };
    double criticalPressure{ 0.0 };
    double acentricFactor{ 0.0 };
    double molFraction{ 1.0 };
    bool swapable{ false };
    double partialPressure{ 0.0 };

    double mass{ 0.0 };
    bool computeFugacityCoefficient{ true };
    double partialFugacity{ 0.0 };
    double fugacityCoefficient{ 1.0 };
    double amountOfExcessMolecules { 0.0 };
    double bulkFluidDensity{ 0.0 };
    double compressibility{ 0.0 };

    std::optional<double> idealGasRosenbluthWeight{};
    std::optional<double> idealGasEnergy{};

    double netCharge{ 0.0 };
    size_t startingBead{ 0 };
    std::vector<Atom> definedAtoms{};
    std::vector<Atom> atoms{};

    size_t initialNumberOfMolecules{ 0 };

    Lambda lambda;
    bool hasFractionalMolecule{ false };

    std::vector<size_t> chiralCenters{};
    std::vector<BondPotential> bonds{};
    std::vector<std::pair<size_t, size_t>> bondDipoles{};
    std::vector<std::tuple<size_t, size_t, size_t>> bends{};
    std::vector<std::pair<size_t, size_t>>  UreyBradley{};
    std::vector<std::tuple<size_t, size_t, size_t, size_t>> inversionBends{};
    std::vector<std::tuple<size_t, size_t, size_t, size_t>> Torsion{};
    std::vector<std::tuple<size_t, size_t, size_t, size_t>> ImproperTorsions{};
    std::vector<std::tuple<size_t, size_t, size_t>> bondBonds{};
    std::vector<std::tuple<size_t, size_t, size_t>> stretchBends{};
    std::vector<std::tuple<size_t, size_t, size_t, size_t>> bendBends{};
    std::vector<std::tuple<size_t, size_t, size_t, size_t>> stretchTorsions{};
    std::vector<std::tuple<size_t, size_t, size_t, size_t>> bendTorsions{};
    std::vector<std::pair<size_t, size_t>> intraVDW{};
    std::vector<std::pair<size_t, size_t>> intraCoulomb{};
    std::vector<std::pair<size_t, size_t>> excludedIntraCoulomb{};
    std::vector<std::pair<size_t, std::vector<size_t>>> configMoves{};
    
    double probabilityTranslationMove{ 0.0 };
    double probabilityRandomTranslationMove{ 0.0 };
    double probabilityRotationMove{ 0.0 };
    double probabilityRandomRotationMove{ 0.0 };
    double probabilityVolumeMove{ 0.0 };
    double probabilityReinsertionMove_CBMC{ 0.0 };
    double probabilityIdentityChangeMove_CBMC{ 0.0 };
    double probabilitySwapMove_CBMC{ 0.0 };
    double probabilitySwapMove_CFCMC{ 0.0 };
    double probabilitySwapMove_CFCMC_CBMC{ 0.0 };
    double probabilityGibbsVolumeMove{ 0.0 };
    double probabilityGibbsSwapMove_CBMC{ 0.0 };
    double probabilityGibbsSwapMove_CFCMC{ 0.0 };
    double probabilityGibbsSwapMove_CFCMC_CBMC{ 0.0 };
    double probabilityWidomMove{ 0.0 };
    double probabilityWidomMove_CFCMC{ 0.0 };
    double probabilityWidomMove_CFCMC_CBMC{ 0.0 };

    double accumulatedProbabilityTranslationMove{ 0.0 };
    double accumulatedProbabilityRandomTranslationMove{ 0.0 };
    double accumulatedProbabilityRotationMove{ 0.0 };
    double accumulatedProbabilityRandomRotationMove{ 0.0 };
    double accumulatedProbabilityVolumeMove{ 0.0 };
    double accumulatedProbabilityReinsertionMove_CBMC{ 0.0 };
    double accumulatedProbabilityIdentityChangeMove_CBMC{ 0.0 };
    double accumulatedProbabilitySwapMove_CBMC{ 0.0 };
    double accumulatedProbabilitySwapMove_CFCMC{ 0.0 };
    double accumulatedProbabilitySwapMove_CFCMC_CBMC{ 0.0 };
    double accumulatedProbabilityGibbsVolumeMove{ 0.0 };
    double accumulatedProbabilityGibbsSwapMove_CBMC{ 0.0 };
    double accumulatedProbabilityGibbsSwapMove_CFCMC{ 0.0 };
    double accumulatedProbabilityGibbsSwapMove_CFCMC_CBMC{ 0.0 };
    double accumulatedProbabilityWidomMove{ 0.0 };
    double accumulatedProbabilityWidomMove_CFCMC{ 0.0 };
    double accumulatedProbabilityWidomMove_CFCMC_CBMC{ 0.0 };
    
    MoveStatistics<double3> statistics_TranslationMove{.maxChange = double3(1.0,1.0,1.0)};
    MoveStatistics<double3> statistics_RandomTranslationMove{};
    MoveStatistics<double3> statistics_RotationMove{ .maxChange = double3(1.0,1.0,1.0) };
    MoveStatistics<double3> statistics_RandomRotationMove{};
    MoveStatistics<double> statistics_ReinsertionMove_CBMC{};
    MoveStatistics<double> statistics_IdentityChangeMove_CBMC{};
    MoveStatistics<double> statistics_SwapInsertionMove_CBMC{};
    MoveStatistics<double> statistics_SwapDeletionMove_CBMC{};
    MoveStatistics<double3> statistics_SwapMove_CFCMC{};
    MoveStatistics<double3> statistics_SwapMove_CFCMC_CBMC{};
    MoveStatistics<double3> statistics_WidomMove_CBMC{};
    MoveStatistics<double3> statistics_WidomMove_CFCMC{};
    MoveStatistics<double3> statistics_WidomMove_CFCMC_CBMC{};
    void clearMoveStatistics();
    
    std::chrono::duration<double> cpuTime_TranslationMove{ 0.0 };
    std::chrono::duration<double> cpuTime_RandomTranslationMove{ 0.0 };
    std::chrono::duration<double> cpuTime_RotationMove{ 0.0 };
    std::chrono::duration<double> cpuTime_RandomRotationMove{ 0.0 };
    std::chrono::duration<double> cpuTime_ReinsertionMove_CBMC{ 0.0 };
    std::chrono::duration<double> cpuTime_IdentityChangeMove_CBMC{ 0.0 };
    std::chrono::duration<double> cpuTime_SwapInsertionMove_CBMC{ 0.0 };
    std::chrono::duration<double> cpuTime_SwapDeletionMove_CBMC{ 0.0 };
    std::chrono::duration<double> cpuTime_SwapMove_CFCMC{ 0.0 };
    std::chrono::duration<double> cpuTime_SwapMove_CFCMC_CBMC{ 0.0 };
    std::chrono::duration<double> cpuTime_WidomMove_CBMC{ 0.0 };
    std::chrono::duration<double> cpuTime_WidomMove_CFCMC{ 0.0 };
    std::chrono::duration<double> cpuTime_WidomMove_CFCMC_CBMC{ 0.0 };

    std::chrono::duration<double> cpuTime_TranslationMove_NonEwald{ 0.0 };
    std::chrono::duration<double> cpuTime_RandomTranslationMove_NonEwald{ 0.0 };
    std::chrono::duration<double> cpuTime_RotationMove_NonEwald{ 0.0 };
    std::chrono::duration<double> cpuTime_RandomRotationMove_NonEwald{ 0.0 };
    std::chrono::duration<double> cpuTime_ReinsertionGrowMove_CBMC_NonEwald{ 0.0 };
    std::chrono::duration<double> cpuTime_ReinsertionRetraceMove_CBMC_NonEwald{ 0.0 };
    std::chrono::duration<double> cpuTime_IdentityChangeMove_CBMC_NonEwald{ 0.0 };
    std::chrono::duration<double> cpuTime_SwapInsertionMove_CBMC_NonEwald{ 0.0 };
    std::chrono::duration<double> cpuTime_SwapDeletionMove_CBMC_NonEwald{ 0.0 };
    std::chrono::duration<double> cpuTime_SwapMove_CFCMC_NonEwald{ 0.0 };
    std::chrono::duration<double> cpuTime_SwapInsertionMove_CFCMC_CBMC_NonEwald{ 0.0 };
    std::chrono::duration<double> cpuTime_SwapInsertionGrowMove_CFCMC_CBMC_NonEwald{ 0.0 };
    std::chrono::duration<double> cpuTime_SwapDeletionMove_CFCMC_CBMC_NonEwald{ 0.0 };
    std::chrono::duration<double> cpuTime_SwapDeletionRetraceMove_CFCMC_CBMC_NonEwald{ 0.0 };
    std::chrono::duration<double> cpuTime_SwapLambdaMove_CFCMC_CBMC_NonEwald{ 0.0 };
    std::chrono::duration<double> cpuTime_WidomMove_CBMC_NonEwald{ 0.0 };
    std::chrono::duration<double> cpuTime_WidomMove_CFCMC_NonEwald{ 0.0 };
    std::chrono::duration<double> cpuTime_WidomMove_CFCMC_CBMC_NonEwald{ 0.0 };
 
    std::chrono::duration<double> cpuTime_TranslationMove_Ewald{ 0.0 };
    std::chrono::duration<double> cpuTime_RandomTranslationMove_Ewald{ 0.0 };
    std::chrono::duration<double> cpuTime_RotationMove_Ewald{ 0.0 };
    std::chrono::duration<double> cpuTime_RandomRotationMove_Ewald{ 0.0 };
    std::chrono::duration<double> cpuTime_ReinsertionMove_CBMC_Ewald{ 0.0 };
    std::chrono::duration<double> cpuTime_IdentityChangeMove_CBMC_Ewald{ 0.0 };
    std::chrono::duration<double> cpuTime_SwapInsertionMove_CBMC_Ewald{ 0.0 };
    std::chrono::duration<double> cpuTime_SwapDeletionMove_CBMC_Ewald{ 0.0 };
    std::chrono::duration<double> cpuTime_SwapMove_CFCMC_Ewald{ 0.0 };
    std::chrono::duration<double> cpuTime_SwapInsertionMove_CFCMC_CBMC_Ewald{ 0.0 };
    std::chrono::duration<double> cpuTime_SwapInsertionGrowMove_CFCMC_CBMC_Ewald{ 0.0 };
    std::chrono::duration<double> cpuTime_SwapDeletionMove_CFCMC_CBMC_Ewald{ 0.0 };
    std::chrono::duration<double> cpuTime_SwapDeletionRetraceMove_CFCMC_CBMC_Ewald{ 0.0 };
    std::chrono::duration<double> cpuTime_SwapLambdaMove_CFCMC_CBMC_Ewald{ 0.0 };
    std::chrono::duration<double> cpuTime_WidomMove_CBMC_Ewald{ 0.0 };
    std::chrono::duration<double> cpuTime_WidomMove_CFCMC_Ewald{ 0.0 };
    std::chrono::duration<double> cpuTime_WidomMove_CFCMC_CBMC_Ewald{ 0.0 };
    void clearTimingStatistics();

    PropertyWidom averageRosenbluthWeights;

    MultiSiteIsotherm isotherm{};      // isotherm information
    double massTransferCoefficient{ 0.0 };    // breakthrough masstransfer coefficient [1/s]
    double axialDispersionCoefficient{ 0.0 }; // breakthrough axial disperion coefficient [m^2/s]
    bool isCarrierGas{ false };        // whether or not this is the carrier-gas

    size_t columnPressure{ 0 };
    size_t columnLoading{ 1 };
    size_t columnError{ 2 };

    enum class PressureScale
    {
      Log = 0,
      Normal = 1
    };

    PressureScale pressureScale{ PressureScale::Log };
};
