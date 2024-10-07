module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <fstream>
#include <iostream>
#include <optional>
#include <vector>
#endif

export module monte_carlo;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <iostream>;
import <fstream>;
import <chrono>;
import <optional>;
#endif

import randomnumbers;
import threadpool;

import averages;
import system;
import mc_moves;
import input_reader;
import energy_status;
import archive;
import json;

export struct MonteCarlo
{
  enum class WeightingMethod : size_t
  {
    LambdaZero = 0,
    AllLambdas = 1
  };

  enum class SimulationStage : size_t
  {
    Uninitialized = 0,
    Initialization = 1,
    Equilibration = 2,
    Production = 3
  };

  MonteCarlo();

  MonteCarlo(InputReader &reader) noexcept;
  MonteCarlo(size_t numberOfCycles, size_t numberOfInitializationCycles, size_t numberOfEquilibrationCycles,
             size_t printEvery, size_t writeBinaryRestartEvery, size_t rescaleWangLandauEvery,
             size_t optimizeMCMovesEvery, std::vector<System> &systems, RandomNumber &randomSeed,
             size_t numberOfBlocks);

  uint64_t versionNumber{1};

  size_t numberOfCycles;
  size_t numberOfSteps;
  size_t numberOfInitializationCycles;
  size_t numberOfEquilibrationCycles;
  size_t printEvery;
  size_t writeBinaryRestartEvery;
  size_t rescaleWangLandauEvery;
  size_t optimizeMCMovesEvery;

  size_t currentCycle{0};
  SimulationStage simulationStage{SimulationStage::Uninitialized};

  std::vector<System> systems;
  RandomNumber random;
  size_t fractionalMoleculeSystem{0};  // the system where the fractional molecule is located

  std::vector<std::ofstream> streams;
  std::vector<std::string> outputJsonFileNames;
  std::vector<nlohmann::json> outputJsons;

  BlockErrorEstimation estimation;

  std::chrono::duration<double> totalInitializationSimulationTime{0};
  std::chrono::duration<double> totalEquilibrationSimulationTime{0};
  std::chrono::duration<double> totalProductionSimulationTime{0};
  std::chrono::duration<double> totalSimulationTime{0};

  void createOutputFiles();
  void run();
  void cycle();
  void initialize();
  void equilibrate();
  void production();
  void output();
  System &randomSystem();

  std::string repr() const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MonteCarlo &mc);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MonteCarlo &mc);
};
