module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <fstream>
#include <iostream>
#include <optional>
#include <vector>
#endif

export module parallel_tempering;

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

export struct ParallelTempering
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

  ParallelTempering();

  ParallelTempering(InputReader &reader) noexcept;

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

  BlockErrorEstimation estimation;

  std::chrono::duration<double> totalSimulationTime{0};
  ThreadPool threadPool;

  void runSystemCycleInitialize(System &system);
  void runSystemCycleEquilibrate(System &system);
  void runSystemCycleProduction(System &system);
  void createOutputFiles();
  void run();
  void initialize();
  void equilibrate();
  void production();
  void output();
  System &randomSystem();

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const ParallelTempering &pt);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, ParallelTempering &pt);
};
