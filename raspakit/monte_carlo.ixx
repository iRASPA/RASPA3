export module monte_carlo;

import randomnumbers;
import threadpool;

import <vector>;
import <iostream>;
import <fstream>;
import <chrono>;
import <optional>;

import averages;
import system;
import mc_moves;
import input_reader;
import energy_status;
import archive;

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

  MonteCarlo(): random(std::nullopt) {};

  MonteCarlo(InputReader& reader) noexcept;

  bool operator==(MonteCarlo const &rhs) const
  {
    return (versionNumber == rhs.versionNumber) &&
           (numberOfCycles == rhs.numberOfCycles) &&
           (numberOfSteps == rhs.numberOfSteps) &&
           (numberOfInitializationCycles == rhs.numberOfInitializationCycles) &&
           (numberOfEquilibrationCycles == rhs.numberOfEquilibrationCycles) &&
           (optimizeMCMovesEvery == rhs.optimizeMCMovesEvery) &&
           (printEvery == rhs.printEvery) &&
           (currentCycle == rhs.currentCycle) &&
           (simulationStage == rhs.simulationStage) &&
           (random == rhs.random) &&
           (systems == rhs.systems) &&
           (fractionalMoleculeSystem == rhs.fractionalMoleculeSystem) &&
           (estimation == rhs.estimation) &&
           (particleMoves == rhs.particleMoves);
  }



  uint64_t versionNumber{ 1 };

  size_t numberOfCycles;
  size_t numberOfSteps;
  size_t numberOfInitializationCycles;
  size_t numberOfEquilibrationCycles;
  size_t optimizeMCMovesEvery{ 5000 };
  size_t printEvery;

  size_t currentCycle{ 0 };
  SimulationStage simulationStage{SimulationStage::Uninitialized};

  std::vector<System> systems;
  RandomNumber random;
  size_t fractionalMoleculeSystem{ 0 };   // the system where the fractional molecule is located

  std::vector<std::ofstream> streams;

  BlockErrorEstimation estimation;

  MC_Moves particleMoves;

  std::chrono::system_clock::time_point t1;
  std::chrono::system_clock::time_point t2;

  void createOutputFiles();
  void run();
  void initialize();
  void equilibrate();
  void production();
  void output();
  System& randomSystem();

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MonteCarlo &mc);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MonteCarlo &mc);
};
