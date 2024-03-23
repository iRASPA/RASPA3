export module input_reader;

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

import stringutils;

import int3;
import double3;
import threadpool;

import system;
import atom;
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
  InputReader(size_t totalNumberOfSystems, size_t numberOfComponents, size_t numberOfBlocks);

  static std::tuple<size_t, size_t, size_t> readNumberOfSystemsComponentsAndBlocks();

  std::vector<ForceField> forceFields;
  std::vector<std::string> componentNames;
  void readForceFields();
  void readComponentNames();

  std::vector<InputDataSystem> inputDataSystem;

  void requireExistingSystem(const std::string& keyword, size_t lineNumber);
  void requireExistingSystemAndComponent(const std::string& keyword, size_t lineNumber);

  SimulationType simulationType{ SimulationType::MonteCarlo };

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

  size_t numberOfBlocks{ 5 };
  size_t numberOfLambdaBins{ 41 };
  size_t numberOfThreads{ 1 };
  ThreadPool::ThreadingType threadingType{ThreadPool::ThreadingType::Serial};

  ForceField forceField;
  std::vector<System> systems{};

  // size_t carrierGasComponent{ 0 };
  std::string displayName{"Column"};
  double temperature{ -1.0 };
};
