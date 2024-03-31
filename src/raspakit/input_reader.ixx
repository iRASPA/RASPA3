export module input_reader;

import <string>;
import <map>;
import <unordered_set>;
import <vector>;
import <complex>;
import <locale>;
import <algorithm>;
import <cctype>;
import <optional>;
import <fstream>;
import <istream>;
import <string>;
import <format>;

import stringutils;

import int3;
import double3;
import threadpool;

import system;
import atom;
import framework;
import component;
import simulationbox;
import forcefield;
import loadings;
import enthalpy_of_adsorption;
import energy_status;
import averages;

struct InputDataSystem
{
  // sampling the conventional radial distribution function (RDF)
  bool computeConventionalRadialDistributionFunction{ false };
  size_t sampleConventionalRadialDistributionFunctionEvery{ 10 };
  size_t writeConventionalRadialDistributionFunctionEvery { 5000 };
  size_t conventionalRadialDistributionFunctionHistogramSize{ 128 };
  double conventionalRadialDistributionFunctionRange{ 12 };

  // sampling the radial distribution function (RDF)
  bool computeRadialDistributionFunction{ false };
  size_t sampleRadialDistributionFunctionEvery{ 10 };
  size_t writeRadialDistributionFunctionEvery{ 5000 };
  size_t radialDistributionFunctionHistogramSize{ 128 };
  double radialDistributionFunctionRange{ 15.0 };

  // sampling the 3D density
  bool computeDensityGrid{ false };
  size_t sampleDensityGridEvery{ 1 };
  size_t writeDensityGridEvery{ 5000 };
  int3 densityGridSize{ 128, 128, 128 };
};


export struct InputReader
{
  enum class SimulationType : int
  {
    MonteCarlo = 0,
    MonteCarloTransitionMatrix = 1,
    MolecularDynamics = 2,
    Minimization = 3,
    Test = 4,
    Breakthrough = 5,
    MixturePrediction = 6,
    Fitting = 7
  };

  InputReader(const std::string inputFile);

  std::ifstream inputStream;
  std::string inputString;

  InputReader::SimulationType simulationType{ SimulationType::MonteCarlo };
  InputReader::SimulationType scanSimulationType();

  size_t numberOfBlocks{ 5 };
  size_t scanNumberOfBlocks();
  size_t totalNumberOfSystems;
  size_t scanNumberOfSystems();
  size_t totalNumberOfComponents;
  size_t scanNumberOfComponents();

  std::vector<ForceField> forceFields;
  std::vector<ForceField> scanForceFields();
  void scanForceFieldParameters(std::vector<ForceField> &scannedForceFields);

  std::vector<std::optional<SimulationBox>> scanBoxes();
  std::vector<std::vector<Framework>> scanFrameworkComponents();

  std::vector<std::vector<Component>> scanComponents();

  std::vector<std::vector<size_t>> scanInitialNumberOfMolecules();
  std::vector<double> scanSystemTemperatures();
  std::vector<std::optional<double>> scanSystemPressures();

  void scanGeneralSettings();
  void scanSystemSettings();
  void scanComponentSettings();

  std::vector<InputDataSystem> inputDataSystem;

  void requireExistingSystem(const std::string& keyword, size_t lineNumber);
  void requireExistingSystemAndComponent(const std::string& keyword, size_t lineNumber);



  size_t numberOfCycles{ 0 };
  size_t numberOfInitializationCycles{ 0 };
  size_t numberOfEquilibrationCycles{ 0 };
  std::size_t printEvery{ 5000 };
  size_t writeBinaryRestartEvery{ 5000 };
  size_t rescaleWangLandauEvery{ 5000 };
  size_t optimizeMCMovesEvery{ 5000 };
  size_t writeEvery{ 100 };
  bool restartFromBinary{ false };

  std::optional<unsigned long long> randomSeed{ std::nullopt };

  size_t numberOfThreads{ 1 };
  ThreadPool::ThreadingType threadingType{ThreadPool::ThreadingType::Serial};

  ForceField forceField;
  std::vector<System> systems{};

  // size_t carrierGasComponent{ 0 };
  std::string displayName{"Column"};
  double temperature{ -1.0 };
};


template<class T>
T parse(const std::string& arguments, [[maybe_unused]] const std::string& keyword, [[maybe_unused]] size_t lineNumber)
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
        throw std::runtime_error(std::format("No values could be read for keyword '{}' at line: {}\n",
                                             keyword, lineNumber));
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
        throw std::runtime_error(std::format("No values could be read for keyword '{}' at line: {}\n",
                                             keyword, lineNumber));
      }
      return list;
    }

  };

  if (list.empty())
  {
    throw std::runtime_error(std::format("No values could be read for keyword '{}' at line: {}",
                                         keyword, lineNumber));
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

  throw std::runtime_error(std::format("Numbers could not be read for keyword '{}' at line: {}\n", keyword, lineNumber));
}

template<>
std::vector<bool> parseListOfSystemValues(const std::string& arguments, const std::string& keyword, size_t lineNumber)
{
  std::vector<bool> list{};

  std::string str, boolString;
  std::istringstream ss(arguments);

  while (ss >> str)
  {
    if (trim(str).rfind("//", 0) == 0)
    {
      if (list.empty())
      {
        throw std::runtime_error(std::format("No values could be read for keyword '{}' at line: {}\n",
                                             keyword, lineNumber));
      }
      return list;
    }
    bool value = false;
    std::istringstream s(str);
    std::istringstream s2(str);
    if (s >> std::boolalpha >>  value)
    {
      list.push_back(value);
    }
    else if (s2 >> boolString)
    {
      if (caseInSensStringCompare(boolString, "yes"))
      {
        list.push_back(true);
      }
      else if (caseInSensStringCompare(boolString, "no"))
      {
        list.push_back(false);
      }
      else
      {
        if (list.empty())
        {
          throw std::runtime_error(std::format("No values could be read for keyword '{}' at line: {}\n",
                                               keyword, lineNumber));
        }
        return list;
      }
    }
  };

  if (list.empty())
  {
    throw std::runtime_error( std::format("No values could be read for keyword '{}' at line: {}\n", keyword, lineNumber));
  }
  return list;
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

  throw std::runtime_error(std::format("Numbers could not be read for keyword '{}' at line: {}\n", keyword, lineNumber));
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

  throw std::runtime_error(std::format("Numbers could not be read for keyword '{}' at line: {}\n", keyword, lineNumber));
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

  throw std::runtime_error(std::format("Numbers could not be read for keyword '{}' at line: {}\n", keyword, lineNumber));
}

int3 parseInt3(const std::string& arguments, const std::string& keyword, size_t lineNumber)
{
  int3 value{};

  std::string str;
  std::istringstream ss(arguments);

  if (ss >> value.x >> value.y >> value.z)
  {
    return value;
  };

  throw std::runtime_error(std::format("Numbers could not be read for keyword '{}' at line: {}\n", keyword, lineNumber));
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

  throw std::runtime_error(std::format("Numbers could not be read for keyword '{}' at line: {}\n", keyword, lineNumber));
}



struct InsensitiveCompare 
{
  bool operator() (const std::string& a, const std::string& b) const 
  {
    return caseInSensStringCompare(a,b);
    //return strcasecmp(a.c_str(), b.c_str()) < 0;
  }
};

static std::unordered_set<std::string, std::hash<std::string>, InsensitiveCompare> generalKeywords 
{
  "RestartFromBinary",
  "RandomSeed",
  "SimulationType",
  "NumberOfBlocks"
  "NumberOfCycles",
  "NumberOfInitializationCycles",
  "NumberOfEquilibrationCycles",
  "PrintEvery",
  "WriteBinaryRestartEvery",
  "RescaleWangLandauEvery",
  "OptimizeMCMovesEvery",
  "NumberOfLambdaBins",
  "NumberOfThreads",
  "ThreadingType"
};

static std::unordered_set<std::string, std::hash<std::string>, InsensitiveCompare> systemKeywords 
{
  "ForceField",
  "Box",
  "BoxLengths",
  "BoxAngles",
  "Framework",
  "FrameworkName",
  "NumberOfUnitCells",
  "ExternalTemperature",
  "ExternalPressure",
  "ChargeMethod",
  "EwaldPrecision",
  "EwaldParameters",
  "OmitEwaldFourier",
  "ComputeConventionalRDF",
  "WriteConventionalRDFEvery",
  "ConventionalRDFistogramSize",
  "ConventionalRDFRange",
  "ComputeRDF",
  "WriteRDFEvery",
  "RDFistogramSize",
  "RDFRange",
  "ComputeDensityGrid",
  "SampleDensityGridEvery",
  "WriteDensityGridEvery",
  "DensityGridSize",
  "ExternalField",
  "VolumeMoveProbability",
  "GibbsVolumeMoveProbability",
  "UseBiasOnMacrostate",
  "TMMCMin",
  "TMMCMax",
  "Reaction",
  "Movies",
  "WriteMoviesEvery",
  "PressureStart",
  "PressureEnd",
  "NumberOfPressurePoints",
  "PressureScale",
  "ColumnVoidFraction",
  "ParticleDensity",
  "TotalPressure",
  "PressureGradient",
  "ColumnEntranceVelocity",
  "TimeStep",
  "NumberOfTimeSteps",
  "WriteEvery",
  "ColumnLength",
  "NumberOfGridPoints",
  "MixturePredictionMethod",
  "ColumnPressure",
  "ColumnLoading",
  "ColumnError"
};

static std::unordered_set<std::string, std::hash<std::string>, InsensitiveCompare> componentKeywords 
{
  "Component",
  "MoleculeName",
  "MoleculeDefinition",
  "CreateNumberOfMolecules",
  "MolFraction",
  "FugacityCoefficient",
  "ThermodynamicIntegration",
  "LnPartitionFunction",
  "FileName",
  "CarrierGas",
  "MassTransferCoefficient",
  "AxialDispersionCoefficient",
  "TranslationProbability",
  "RandomTranslationProbability",
  "RotationProbability",
  "RandomRotationProbability",
  "ReinsertionProbability",
  "SwapConventionalProbability",
  "SwapProbability",
  "GibbsSwapProbability",
  "GibbsCFCMCSwapProbability",
  "CFCMC_SwapProbability",
  "CFCMC_CBMC_SwapProbability",
  "WidomProbability",
  "CFCMCWidomProbability",
  "CBCFCMCWidomProbability",
  "MassTransferCoefficient",
  "AxialDispersionCoefficient",
  "NumberOfIsothermSites",
  "Langmuir",
  "Anti-Langmuir",
  "BET",
  "Henry",
  "Freundlich",
  "Sips",
  "Langmuir-Freundlich",
  "Redlich-Peterson",
  "Toth",
  "Unilan",
  "O'Brien&Myers",
  "Quadratic",
  "Temkin"
};

static std::unordered_set<std::string, std::hash<std::string>, InsensitiveCompare> allKeywords 
{
  "RestartFromBinary",
  "RandomSeed",
  "SimulationType",
  "NumberOfBlocks",
  "NumberOfCycles",
  "NumberOfInitializationCycles",
  "NumberOfEquilibrationCycles",
  "PrintEvery",
  "WriteBinaryRestartEvery",
  "RescaleWangLandauEvery",
  "OptimizeMCMovesEvery",
  "NumberOfLambdaBins",
  "NumberOfThreads",
  "ThreadingType",

  "ForceField",
  "Box",
  "BoxLengths",
  "BoxAngles",
  "Framework",
  "FrameworkName",
  "NumberOfUnitCells",
  "ExternalTemperature",
  "ExternalPressure",
  "ChargeMethod",
  "EwaldPrecision",
  "EwaldParameters",
  "OmitEwaldFourier",
  "ComputeConventionalRDF",
  "WriteConventionalRDFEvery",
  "ConventionalRDFistogramSize",
  "ConventionalRDFRange",
  "ComputeRDF",
  "WriteRDFEvery",
  "RDFistogramSize",
  "RDFRange",
  "ComputeDensityGrid",
  "SampleDensityGridEvery",
  "WriteDensityGridEvery",
  "DensityGridSize",
  "ExternalField",
  "VolumeMoveProbability",
  "GibbsVolumeMoveProbability",
  "UseBiasOnMacrostate",
  "TMMCMin",
  "TMMCMax",
  "Reaction",
  "Movies",
  "WriteMoviesEvery",
  "PressureStart",
  "PressureEnd",
  "NumberOfPressurePoints",
  "PressureScale",
  "ColumnVoidFraction",
  "ParticleDensity",
  "TotalPressure",
  "PressureGradient",
  "ColumnEntranceVelocity",
  "TimeStep",
  "NumberOfTimeSteps",
  "WriteEvery",
  "ColumnLength",
  "NumberOfGridPoints",
  "MixturePredictionMethod",
  "ColumnPressure",
  "ColumnLoading",
  "ColumnError",

  "Component",
  "MoleculeName",
  "MoleculeDefinition",
  "CreateNumberOfMolecules",
  "MolFraction",
  "FugacityCoefficient",
  "ThermodynamicIntegration",
  "LnPartitionFunction",
  "FileName",
  "CarrierGas",
  "MassTransferCoefficient",
  "AxialDispersionCoefficient",
  "TranslationProbability",
  "RandomTranslationProbability",
  "RotationProbability",
  "RandomRotationProbability",
  "ReinsertionProbability",
  "SwapConventionalProbability",
  "SwapProbability",
  "GibbsSwapProbability",
  "GibbsCFCMCSwapProbability",
  "CFCMC_SwapProbability",
  "CFCMC_CBMC_SwapProbability",
  "WidomProbability",
  "CFCMCWidomProbability",
  "CBCFCMCWidomProbability",
  "MassTransferCoefficient",
  "AxialDispersionCoefficient",
  "NumberOfIsothermSites",
  "Langmuir",
  "Anti-Langmuir",
  "BET",
  "Henry",
  "Freundlich",
  "Sips",
  "Langmuir-Freundlich",
  "Redlich-Peterson",
  "Toth",
  "Unilan",
  "O'Brien&Myers",
  "Quadratic",
  "Temkin"
};
