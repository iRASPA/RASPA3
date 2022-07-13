module;

module input_reader;

import system;
import atom;
import component;
import simulationbox;
import forcefield;
import double3;
import units;
import print;
import sample_movies;
import threadpool;

import <filesystem>;
import <fstream>;
import <cstdlib>;
import <iostream>;
import <sstream>;
import <exception>;
import <numbers>;
import <vector>;
import <complex>;
import <ios>;

bool caseInSensStringCompare(const std::string& str1, const std::string& str2)
{
    return str1.size() == str2.size() && std::equal(str1.begin(), str1.end(), str2.begin(), [](auto a, auto b) {return std::tolower(a) == std::tolower(b); });
}

std::string trim(const std::string& s)
{
    auto start = s.begin();
    while (start != s.end() && std::isspace(*start)) {
        start++;
    }

    auto end = s.end();
    do {
        end--;
    } while (std::distance(start, end) > 0 && std::isspace(*end));

    return std::string(start, end + 1);
}

template<class T>
T parse(const std::string& arguments)
{
    T value;

    std::string str;
    std::istringstream ss(arguments);

    ss >> value;

    return value;
}

template<typename T>
std::vector<T> parseListOfSystemValues(const std::string& arguments, const std::string& keyword, size_t lineNumber)
{
    std::vector<T> list{};

    std::string str;
    std::istringstream ss(arguments);

    while (ss >> str)
    {
        if (trim(str).rfind("//", 0) == 0)
        {
            if (list.empty())
            {
                throw std::runtime_error(std::print("No values could be read for keyword '{}' at line: {}", keyword, lineNumber));
            }
            return list;
        }
        T value;
        std::istringstream s(str);
        if (s >> value)
        {
            list.push_back(value);
        }
        else
        {
            if (list.empty()) 
            {
                throw std::runtime_error(std::print("No values could be read for keyword '{}' at line: {}", keyword, lineNumber));
            }
            return list;
        }

    };

    if (list.empty()) 
    {
        throw std::runtime_error( std::print("No values could be read for keyword '{}' at line: {}", keyword, lineNumber));
    }
    return list;
}

int parseBoolean(const std::string& arguments, const std::string& keyword, size_t lineNumber)
{
    bool value{};

    std::istringstream ss(arguments);

    if (ss >> std::boolalpha >> value)
    {
        return value;
    };

    std::string str;
    std::istringstream ss2(arguments);
    if (ss2 >> str)
    {
        if (caseInSensStringCompare(str, "yes")) return true;
        if (caseInSensStringCompare(str, "no")) return false;
    };

    throw std::runtime_error(std::print("Numbers could not be read for keyword '{}' at line: {}", keyword, lineNumber));
}

size_t parseInteger(const std::string& arguments, const std::string& keyword, size_t lineNumber)
{
    size_t value{};

    std::string str;
    std::istringstream ss(arguments);

    if (ss >> value)
    {
        return value;
    };

    throw std::runtime_error(std::print("Numbers could not be read for keyword '{}' at line: {}", keyword, lineNumber));
}

double parseDouble(const std::string& arguments, const std::string& keyword, size_t lineNumber)
{
    double value{};

    std::string str;
    std::istringstream ss(arguments);

    if (ss >> value)
    {
        return value;
    };

    throw std::runtime_error(std::print("Numbers could not be read for keyword '{}' at line: {}", keyword, lineNumber));
}

double3 parseDouble3(const std::string& arguments, const std::string& keyword, size_t lineNumber)
{
    double3 value{};

    std::string str;
    std::istringstream ss(arguments);

    if (ss >> value.x >> value.y >> value.z)
    {
        return value;
    };

    throw std::runtime_error(std::print("Numbers could not be read for keyword '{}' at line: {}", keyword, lineNumber));
}

std::string parseString(const std::string& arguments, const std::string& keyword, size_t lineNumber)
{
    std::string value{};

    std::string str;
    std::istringstream ss(arguments);

    if (ss >> value)
    {
        return value;
    };

    throw std::runtime_error(std::print("Numbers could not be read for keyword '{}' at line: {}", keyword, lineNumber));
}


InputReader::InputReader(const std::string pseudoAtomsFileName, const std::string forceFieldMixingRulesFileName,
    const std::string forceFieldOverwriteFileName, const std::string simulationSettingsFileName):
    forceField(pseudoAtomsFileName, forceFieldMixingRulesFileName, forceFieldOverwriteFileName)
{
    const char* env_p = std::getenv("RASPA_DIR");

    std::filesystem::path pathfile = std::filesystem::path(simulationSettingsFileName);
    if (!std::filesystem::exists(pathfile)) pathfile = std::filesystem::path(env_p) / simulationSettingsFileName;

    if (!std::filesystem::exists(pathfile)) throw std::runtime_error("simulation.input' not found");

    std::ifstream fileInput{ pathfile };
    if (!fileInput) throw std::runtime_error("File 'simulation.input' exists, but error opening file");

    std::string line{};
    std::string keyword{};
    std::string arguments{};

    std::size_t numberOfSystems{ 0 };
    //std::size_t numberOfFrameworks{ 0 };
    std::size_t numberOfComponents{ 0 };
    size_t lineNumber{ 0 };
   
    while (std::getline(fileInput, line))
    {
        lineNumber += 1;
        if (!line.empty())
        {
            std::istringstream iss(line);

            iss >> keyword;
            keyword = trim(keyword);
            std::getline(iss, arguments);           

            if (caseInSensStringCompare(keyword, std::string("NumberOfCycles")))
            {
                numberOfCycles = parse<size_t>(arguments);
            }

            if (caseInSensStringCompare(keyword, std::string("NumberOfInitializationCycles")))
            {
                numberOfInitializationCycles = parse<size_t>(arguments);
            }

            if (caseInSensStringCompare(keyword, std::string("NumberOfEquilibrationCycles")))
            {
                numberOfEquilibrationCycles = parse<size_t>(arguments);
            }

            if (caseInSensStringCompare(keyword, std::string("PrintEvery")))
            {
                printEvery = parse<size_t>(arguments);
            }

            if (caseInSensStringCompare(keyword, std::string("NumberOfBlocks")))
            {
                numberOfBlocks = parse<size_t>(arguments);
            }

            if (caseInSensStringCompare(keyword, std::string("NumberOfThreads")))
            {
                numberOfThreads = parse<size_t>(arguments);
            }

            if (caseInSensStringCompare(keyword, std::string("ThreadingType")))
            {
                std::string str;
                std::istringstream ss2(arguments);
                if (ss2 >> str)
                {
                    if (caseInSensStringCompare(str, "Serial")) threadingType = ThreadPool::ThreadingType::Serial;
                    if (caseInSensStringCompare(str, "ThreadPool")) threadingType = ThreadPool::ThreadingType::ThreadPool;
                    if (caseInSensStringCompare(str, "OpenMP")) threadingType = ThreadPool::ThreadingType::OpenMP;
                    if (caseInSensStringCompare(str, "GPU-Offload")) threadingType = ThreadPool::ThreadingType::GPU_Offload;
                };
            }

            if (caseInSensStringCompare(keyword, std::string("Box")))
            {
                numberOfSystems += 1;
                systems.emplace_back(numberOfSystems - 1, forceField, std::vector<Component>{}, std::vector<size_t>{}, numberOfBlocks);
            }

            if (caseInSensStringCompare(keyword, "BoxLengths"))
            {
                requireExistingSystem(keyword, lineNumber);
                double3 value = parseDouble3(arguments, keyword, lineNumber);
                systems[numberOfSystems - 1].simulationBox.setBoxLengths(value);
            }

            if (caseInSensStringCompare(keyword, "BoxAngles"))
            {
                requireExistingSystem(keyword, lineNumber);
                double3 value = parseDouble3(arguments, keyword, lineNumber);
                systems[numberOfSystems - 1].simulationBox.setBoxAngles((std::numbers::pi / 180.0) * value);
            }

            if (caseInSensStringCompare(keyword, std::string("Framework")))
            {
                numberOfSystems += 1;
                systems.emplace_back(numberOfSystems - 1, forceField, std::vector<Component>{}, std::vector<size_t>{}, numberOfBlocks);
                //numberOfFrameworks = 0;
            }

            if (caseInSensStringCompare(keyword, std::string("FrameworkName")))
            {
                //numberOfFrameworks += 1;
                numberOfComponents += 1;

                std::istringstream ss(arguments);
                std::string frameworkName;
                ss >> frameworkName;
                
                systems[numberOfSystems - 1].addComponent(Component(Component::Type::Framework, numberOfComponents - 1, systems[numberOfSystems - 1].forceField, frameworkName, numberOfBlocks));
            }

            if (caseInSensStringCompare(keyword, "ExternalTemperature"))
            {
                requireExistingSystem(keyword, lineNumber);
                double value = parseDouble(arguments, keyword, lineNumber);
                systems[numberOfSystems - 1].simulationBox.temperature = value;
                systems[numberOfSystems - 1].simulationBox.Beta = 1.0/(0.8314464919 * value); 
            }

            if (caseInSensStringCompare(keyword, "ExternalPressure"))
            {
                requireExistingSystem(keyword, lineNumber);
                double value = parseDouble(arguments, keyword, lineNumber);
                systems[numberOfSystems - 1].simulationBox.pressure = value / Units::PressureConversionFactor;
            }

            if (caseInSensStringCompare(keyword, "ChargeMethod"))
            {
                requireExistingSystem(keyword, lineNumber);
                std::string value = parseString(arguments, keyword, lineNumber);
                if (caseInSensStringCompare(value, "None")) systems[numberOfSystems - 1].noCharges = true;
                if (caseInSensStringCompare(value, "Ewald")) systems[numberOfSystems - 1].noCharges = false;
            }

            if (caseInSensStringCompare(keyword, "OmitEwaldFourier"))
            {
                requireExistingSystem(keyword, lineNumber);
                bool value = parseBoolean(arguments, keyword, lineNumber);
                systems[numberOfSystems - 1].omitEwaldFourier = value;
            }

            if (caseInSensStringCompare(keyword, "probabilityVolumeMove"))
            {
                requireExistingSystem(keyword, lineNumber);
                double value = parseDouble(arguments, keyword, lineNumber);
                systems[numberOfSystems - 1].probabilityVolumeMove = value;
            }

            if (caseInSensStringCompare(keyword, "Movies"))
            {
                requireExistingSystem(keyword, lineNumber);
                bool value = parseBoolean(arguments, keyword, lineNumber);
                systems[numberOfSystems - 1].sampleMovie.sample = value;
            }

            if (caseInSensStringCompare(keyword, "WriteMoviesEvery"))
            {
                requireExistingSystem(keyword, lineNumber);
                size_t value = parseInteger(arguments, keyword, lineNumber);
                systems[numberOfSystems - 1].sampleMovie.writeEvery = value;
            }

            if (caseInSensStringCompare(keyword, std::string("Component")))
            {
                requireExistingSystem(keyword, lineNumber);
                
                std::istringstream ss(arguments);
                std::string c, moleculeNameKeyword,remainder;
                ss >> c >> moleculeNameKeyword;
                std::getline(ss, remainder);

                std::vector<std::string> values = parseListOfSystemValues<std::string>(remainder, keyword, lineNumber);
                values.resize(systems.size(), values.front());
                numberOfComponents += 1;
                for (size_t i = 0; i < systems.size(); ++i)
                {
                    systems[i].addComponent(Component(Component::Type::Adsorbate, static_cast<size_t>(numberOfComponents  - 1), systems[numberOfSystems - 1].forceField, values[i], numberOfBlocks));
                }
            }

            
            if (caseInSensStringCompare(keyword, "CreateNumberOfMolecules"))
            {
                requireExistingSystemAndComponent(keyword, lineNumber);

                std::vector<size_t> values = parseListOfSystemValues<size_t>(arguments, keyword, lineNumber);
                values.resize(systems.size(), values.front());
                for (size_t i = 0; i < systems.size(); ++i)
                {
                    systems[i].components[numberOfComponents  - 1].initialNumberOfMolecules = values[i];
                }
            }

            if (caseInSensStringCompare(keyword, "MolFraction"))
            {
                requireExistingSystem(keyword, lineNumber);

                std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
         
                values.resize(systems.size(), values.front());
                for (size_t i = 0; i < systems.size(); ++i)
                {
                    systems[i].components[numberOfComponents - 1].molFraction = values[i];
                }
            }

            if (caseInSensStringCompare(keyword, "TranslationProbability"))
            {
                requireExistingSystemAndComponent(keyword, lineNumber);
                std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
         
                values.resize(systems.size(), values.front());
                for (size_t i = 0; i < systems.size(); ++i)
                {
                    systems[i].components[numberOfComponents - 1].probabilityTranslationMove = values[i];
                }
            }

            if (caseInSensStringCompare(keyword, "RandomTranslationProbability"))
            {
                requireExistingSystemAndComponent(keyword, lineNumber);
                std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);
         
                values.resize(systems.size(), values.front());
                for (size_t i = 0; i < systems.size(); ++i)
                {
                    systems[i].components[numberOfComponents - 1].probabilityRandomTranslationMove = values[i];
                }
            }

            if (caseInSensStringCompare(keyword, "RotationProbability"))
            {
                requireExistingSystemAndComponent(keyword, lineNumber);
                std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

                values.resize(systems.size(), values.front());
                for (size_t i = 0; i < systems.size(); ++i)
                {
                    systems[i].components[numberOfComponents - 1].probabilityRotationMove = values[i];
                }
            }

            if (caseInSensStringCompare(keyword, "RandomRotationProbability"))
            {
                requireExistingSystemAndComponent(keyword, lineNumber);
                std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

                values.resize(systems.size(), values.front());
                for (size_t i = 0; i < systems.size(); ++i)
                {
                    systems[i].components[numberOfComponents - 1].probabilityRandomRotationMove = values[i];
                }
            }

            if (caseInSensStringCompare(keyword, "ReinsertionProbability"))
            {
                requireExistingSystemAndComponent(keyword, lineNumber);
                std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

                values.resize(systems.size(), values.front());
                for (size_t i = 0; i < systems.size(); ++i)
                {
                    systems[i].components[numberOfComponents - 1].probabilityReinsertionMove_CBMC = values[i];
                }
            }

            if (caseInSensStringCompare(keyword, "SwapProbability"))
            {
                requireExistingSystemAndComponent(keyword, lineNumber);
                std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

                values.resize(systems.size(), values.front());
                for (size_t i = 0; i < systems.size(); ++i)
                {
                    systems[i].components[numberOfComponents - 1].probabilitySwapMove_CBMC = values[i];
                }
            }

            if (caseInSensStringCompare(keyword, "GibbsSwapProbability"))
            {
                requireExistingSystemAndComponent(keyword, lineNumber);
                std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

                values.resize(systems.size(), values.front());
                for (size_t i = 0; i < systems.size(); ++i)
                {
                    systems[i].components[numberOfComponents - 1].probabilityGibbsSwapMove_CBMC = values[i];
                }
            }

            
            if (caseInSensStringCompare(keyword, "CFCMC_CBMC_SwapProbability"))
            {
                requireExistingSystemAndComponent(keyword, lineNumber);
                std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

                values.resize(systems.size(), values.front());
                for (size_t i = 0; i < systems.size(); ++i)
                {
                    systems[i].components[numberOfComponents - 1].probabilitySwapMove_CFCMC_CBMC = values[i];
                }
            }

            if (caseInSensStringCompare(keyword, "WidomProbability"))
            {
                requireExistingSystemAndComponent(keyword, lineNumber);
                std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

                values.resize(systems.size(), values.front());
                for (size_t i = 0; i < systems.size(); ++i)
                {
                    systems[i].components[numberOfComponents - 1].probabilityWidomMove = values[i];
                }
            }

            if (caseInSensStringCompare(keyword, "CFCMCWidomProbability"))
            {
                requireExistingSystemAndComponent(keyword, lineNumber);
                std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

                values.resize(systems.size(), values.front());
                for (size_t i = 0; i < systems.size(); ++i)
                {
                    systems[i].components[numberOfComponents - 1].probabilityWidomMove_CFCMC = values[i];
                }
            }

            if (caseInSensStringCompare(keyword, "CBCFCMCWidomProbability"))
            {
                requireExistingSystemAndComponent(keyword, lineNumber);
                std::vector<double> values = parseListOfSystemValues<double>(arguments, keyword, lineNumber);

                values.resize(systems.size(), values.front());
                for (size_t i = 0; i < systems.size(); ++i)
                {
                    systems[i].components[numberOfComponents - 1].probabilityWidomMove_CFCMC_CBMC = values[i];
                }
            }
           
        }

    }

}

void InputReader::requireExistingSystem(const std::string& keyword, size_t lineNumber)
{
    if (systems.empty()) {
        throw std::runtime_error(std::print("No system (Framework or Box) defined yet at keyword '{}' at line: {}", keyword, lineNumber));
    }
}

void InputReader::requireExistingSystemAndComponent(const std:: string &keyword, size_t lineNumber)
{
    if (systems.empty()) {
        throw std::runtime_error(
            std::print("No system (Framework or Box) defined yet at keyword '{}' at line: {}", keyword, lineNumber)
            );
    }
    if (systems[0].components.empty()) {
        throw std::runtime_error(
            std::print("No component defined yet at keyword '{}' at line: {}", keyword, lineNumber)
            );
    }
}
