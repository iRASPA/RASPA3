export module monte_carlo;

import randomnumbers;
import threadpool;

import averages;
import system;
import mc_moves;
import input_reader;
import energy_status;
import archive;

import <vector>;
import <iostream>;
import <fstream>;
import <chrono>;

export struct MonteCarlo
{
  int64_t versionNumber{ 1 };
  enum class WeightingMethod : size_t
  {
      LambdaZero = 0,
      AllLambdas = 1
  };

  MonteCarlo(InputReader& reader) noexcept;

  void run();
  void initialize();
  void equilibrate();
  void production();
  void output();

  System& randomSystem();

  size_t numberOfCycles;
  size_t numberOfSteps;
  size_t numberOfInitializationCycles;
  size_t numberOfEquilibrationCycles;
  size_t optimizeMCMovesEvery{ 5000 };
  size_t printEvery;

  std::vector<System> systems;
  size_t fractionalMoleculeSystem{ 0 };   // the system where the fractional molecule is located

  std::vector<std::ofstream> streams;

  BlockErrorEstimation estimation;

  MC_Moves particleMoves;

  std::chrono::system_clock::time_point t1;
  std::chrono::system_clock::time_point t2;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const MonteCarlo &mc);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, MonteCarlo &mc);
};
