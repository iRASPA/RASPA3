module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cctype>
#include <complex>
#include <format>
#include <fstream>
#include <istream>
#include <locale>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>
#endif

export module input_reader;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <map>;
import <set>;
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
#endif

import stringutils;

import int3;
import double3;
import threadpool;
import json;

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
  bool computeConventionalRadialDistributionFunction{false};
  size_t sampleConventionalRadialDistributionFunctionEvery{10};
  size_t writeConventionalRadialDistributionFunctionEvery{5000};
  size_t conventionalRadialDistributionFunctionHistogramSize{128};
  double conventionalRadialDistributionFunctionRange{12};

  // sampling the radial distribution function (RDF)
  bool computeRadialDistributionFunction{false};
  size_t sampleRadialDistributionFunctionEvery{10};
  size_t writeRadialDistributionFunctionEvery{5000};
  size_t radialDistributionFunctionHistogramSize{128};
  double radialDistributionFunctionRange{15.0};

  // sampling the 3D density
  bool computeDensityGrid{false};
  size_t sampleDensityGridEvery{1};
  size_t writeDensityGridEvery{5000};
  int3 densityGridSize{128, 128, 128};
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
    Fitting = 7,
    ParallelTempering = 8
  };

  InputReader(const std::string inputFile);

  std::string inputString;

  InputReader::SimulationType simulationType{SimulationType::MonteCarlo};

  size_t numberOfBlocks{5};

  size_t numberOfCycles{0};
  size_t numberOfInitializationCycles{0};
  size_t numberOfEquilibrationCycles{0};
  std::size_t printEvery{5000};
  size_t writeBinaryRestartEvery{5000};
  size_t rescaleWangLandauEvery{5000};
  size_t optimizeMCMovesEvery{5000};
  size_t writeEvery{100};
  bool restartFromBinary{false};

  std::optional<unsigned long long> randomSeed{std::nullopt};

  size_t numberOfThreads{1};
  ThreadPool::ThreadingType threadingType{ThreadPool::ThreadingType::Serial};

  ForceField forceField;
  std::vector<System> systems{};

  // size_t carrierGasComponent{ 0 };
  std::string displayName{"Column"};
  double temperature{-1.0};

  void parseMolecularSimulations(const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);
  void parseFitting(const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);
  void parseMixturePrediction(const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);
  void parseBreakthrough(const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);

  void validateInput(const nlohmann::basic_json<nlohmann::raspa_map> &parsed_data);

  struct InsensitiveCompare
  {
    bool operator()(const std::string &a, const std::string &b) const 
    { 
    #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
      return _stricmp(a.c_str(), b.c_str()) < 0; 
    #else
      return strcasecmp(a.c_str(), b.c_str()) < 0; 
    #endif
    }
  };

  static const std::set<std::string, InsensitiveCompare> generalOptions;
  static const std::set<std::string, InsensitiveCompare> systemOptions;
  static const std::set<std::string, InsensitiveCompare> componentOptions;
};
