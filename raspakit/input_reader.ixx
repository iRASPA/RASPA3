export module input_reader;

import system;
import atom;
import component;
import simulationbox;
import forcefield;
import double3;
import loadings;
import enthalpy_of_adsorption;
import energy_status;
import averages;
import threadpool;

import <string>;
import <unordered_set>;
import <map>;
import <vector>;
import <complex>;
import <locale>;
import <algorithm>;
import <cctype>;
import <optional>;
import <fstream>;

struct Hash {
    size_t operator()([[maybe_unused]] const std::string& str) const {
        return 0;
    }
};

struct InsensitiveCompare {
    bool operator() (const std::string& str1, const std::string& str2) const {
        return str1.size() == str2.size() && std::equal(str1.begin(), str1.end(), str2.begin(), [](auto a, auto b) {return std::tolower(a) == std::tolower(b); });
    }
};

export struct InputReader
{
    enum class SimulationType : int
    {
        MonteCarlo = 0,
        MolecularDynamics = 1,
        Minimization = 2,
        Test = 3,
        Breakthrough = 4,
        IAST = 5,
        Fitting = 6
    };
    InputReader();
    ~InputReader() {};

    const std::string pseudoAtomsFileName{"pseudo_atoms.def"};
    const std::string forceFieldMixingRulesFileName{"force_field_mixing_rules.def"};
    const std::string forceFieldOverwriteFileName{"force_field.def"};
    const std::string simulationSettingsFileName{"simulation.input"};

    void requireExistingSystem(const std::string& keyword, size_t lineNumber);
    void requireExistingSystemAndComponent(const std::string& keyword, size_t lineNumber);

    SimulationType simulationType{ SimulationType::MonteCarlo };

    size_t numberOfCycles{ 0 };
    size_t numberOfInitializationCycles{ 0 };
    size_t numberOfEquilibrationCycles{ 0 };
    std::size_t printEvery{ 0 };

    size_t numberOfBlocks{ 5 };
    size_t numberOfThreads{ 1 };
    ThreadPool::ThreadingType threadingType{ThreadPool::ThreadingType::Serial};

    ForceField forceField;
    std::vector<System> systems{};

    size_t carrierGasComponent{ 0 };
    std::string displayName{"Column"};
    double temperature{ -1.0 };

    //size_t printEvery{ 1000 };
    size_t writeEvery{ 100 };
};
