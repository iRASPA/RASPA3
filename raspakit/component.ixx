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
import <optional>;
import <span>;

import double3;
import averages;
import atom;
import forcefield;
import lambda;
import simulationbox;
import property_widom;

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

	enum class BondType: int
	{
		Harmonic = 0,
		NumberOfBondTypes = 1
	};

	Component(Component::Type type, size_t currentComponent, const ForceField &forceField, const std::string& fileName, size_t numberOfBlocks) noexcept(false);

    void readComponent(const ForceField& forceField, const std::string& fileName);
    void readFramework(const ForceField& forceField, const std::string& fileName);

	std::string printStatus(const ForceField& forceField) const;

	void normalizeMoveProbabilties();

	const std::string writeMCMoveStatistics() const;
	
	const std::string writeMCMoveCPUTimeStatistics() const;
	
	std::vector<double3> randomlyRotatedPositionsAroundStartingBead() const;
	std::vector<Atom> newAtoms(double scaling, size_t moleculeId) const;
	std::vector<Atom> copiedAtoms(std::span<Atom> molecule) const;
	
	Type type;
    std::optional<SimulationBox> simulationBox{ std::nullopt };

	size_t componentId{ 0 };
	std::string name{};
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

	size_t startingBead{ 0 };
	std::vector<Atom> atoms{};
	std::vector<std::pair<size_t, size_t>> bonds{ {0,1}, {1,2} };

	size_t initialNumberOfMolecules{ 0 };

	Lambda lambda;
    bool hasFractionalMolecule{ false };
	
	double probabilityTranslationMove{ 0.0 };
	double probabilityRotationMove{ 0.0 };
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
	double accumulatedProbabilityRotationMove{ 0.0 };
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
	MoveStatistics<double3> statistics_RotationMove{ .maxChange = double3(1.0,1.0,1.0) };
	MoveStatistics<double> statistics_VolumeMove{};
	MoveStatistics<double> statistics_ReinsertionMove_CBMC{};
	MoveStatistics<double> statistics_IdentityChangeMove_CBMC{};
	MoveStatistics<double> statistics_SwapInsertionMove_CBMC{};
	MoveStatistics<double> statistics_SwapDeletionMove_CBMC{};
	MoveStatistics<double3> statistics_SwapMove_CFCMC{};
	MoveStatistics<double3> statistics_SwapMove_CFCMC_CBMC{};
	MoveStatistics<double> statistics_GibbsVolumeMove{};
	MoveStatistics<double3> statistics_GibbsSwapMove_CBMC{};
	MoveStatistics<double3> statistics_GibbsSwapMove_CFCMC{};
	MoveStatistics<double3> statistics_GibbsSwapMove_CFCMC_CBMC{};
	MoveStatistics<double3> statistics_WidomMove_CBMC{};
	MoveStatistics<double3> statistics_WidomMove_CFCMC{};
	MoveStatistics<double3> statistics_WidomMove_CFCMC_CBMC{};
    void clearMoveStatistics();
	

	std::chrono::duration<double> cpuTime_TranslationMove{ 0.0 };
	std::chrono::duration<double> cpuTime_RotationMove{ 0.0 };
	std::chrono::duration<double> cpuTime_VolumeMove{ 0.0 };
	std::chrono::duration<double> cpuTime_ReinsertionMove_CBMC{ 0.0 };
	std::chrono::duration<double> cpuTime_IdentityChangeMove_CBMC{};
	std::chrono::duration<double> cpuTime_SwapInsertionMove_CBMC{ 0.0 };
	std::chrono::duration<double> cpuTime_SwapDeletionMove_CBMC{ 0.0 };
	std::chrono::duration<double> cpuTime_SwapMove_CFCMC{ 0.0 };
	std::chrono::duration<double> cpuTime_SwapMove_CFCMC_CBMC{ 0.0 };
	std::chrono::duration<double> cpuTime_GibbsVolumeMove{ 0.0 };
	std::chrono::duration<double> cpuTime_GibbsSwapMove_CBMC{ 0.0 };
	std::chrono::duration<double> cpuTime_GibbsSwapLambdaMove_CFCMC{ 0.0 };
	std::chrono::duration<double> cpuTime_GibbsSwapLambdaMove_CFCMC_CBMC{ 0.0 };
	std::chrono::duration<double> cpuTime_WidomMove{ 0.0 };
	std::chrono::duration<double> cpuTime_WidomMove_CFCMC{ 0.0 };
	std::chrono::duration<double> cpuTime_WidomMove_CFCMC_CBMC{ 0.0 };

    PropertyWidom averageRosenbluthWeights;
};
