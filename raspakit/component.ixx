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
import property_lambda_probability_histogram;
import simulationbox;
import property_widom;
import isotherm;
import multi_site_isotherm;
import bond_potential;
import move_statistics;
import mc_moves_probabilities_particles;

export struct Component
{
  enum class Type : int
  {
    Framework = 0,
    Adsorbate = 1,
    Cation = 2
  };

    Component(size_t componentId, std::string componentName, double mass, SimulationBox simulationBox, double T_c, double P_c, double w, std::vector<Atom> definedAtoms, size_t numberOfBlocks) noexcept(false);
    Component(size_t componentId, std::string componentName, double mass, SimulationBox simulationBox, size_t spaceGroupHallNumber, std::vector<Atom> definedAtoms, int3 numberOfUnitCells, size_t numberOfBlocks) noexcept(false);

    Component(Component::Type type, size_t currentComponent, const std::string &componentName, 
              std::optional<const std::string> fileName, size_t numberOfBlocks) noexcept(false);

    Type type;

    SimulationBox simulationBox;
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

    PropertyLambdaProbabilityHistogram lambdaGC;
    PropertyLambdaProbabilityHistogram lambdaGibbs;
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
    
    MCMoveProbabilitiesParticles mc_moves_probabilities;

    PropertyWidom averageRosenbluthWeights;

    MultiSiteIsotherm isotherm{};      // isotherm information
    double massTransferCoefficient{ 0.0 };    // breakthrough masstransfer coefficient [1/s]
    double axialDispersionCoefficient{ 0.0 }; // breakthrough axial disperion coefficient [m^2/s]
    bool isCarrierGas{ false };        // whether or not this is the carrier-gas

    size_t columnPressure{ 0 };
    size_t columnLoading{ 1 };
    size_t columnError{ 2 };

    double lnPartitionFunction { 0 };

    enum class PressureScale
    {
      Log = 0,
      Normal = 1
    };

    PressureScale pressureScale{ PressureScale::Log };

    void readComponent(const ForceField& forceField, const std::string& fileName);
    void readFramework(const ForceField& forceField, const std::string& fileName);

    std::string printStatus(const ForceField& forceField) const;
    void printBreakthroughStatus(std::ostream& stream) const;

    std::vector<double3> randomlyRotatedPositionsAroundStartingBead() const;
    std::vector<Atom> newAtoms(double scaling, size_t moleculeId) const;
    std::vector<Atom> copyAtoms(std::span<Atom> molecule, double scaling, size_t moleculeId) const;
    std::vector<Atom> copyAtomsRandomlyRotatedAt(double3 position, std::span<Atom> molecule, double scaling, size_t moleculeId) const;
    std::vector<Atom> copiedAtoms(std::span<Atom> molecule) const;
};
