export module monte_carlo_transition_matrix;

import randomnumbers;
import threadpool;

import averages;
import system;
import mc_moves;
import input_reader;
import energy_status;

import <vector>;
import <iostream>;
import <fstream>;
import <chrono>;

export struct MonteCarloTransitionMatrix
{
  enum class WeightingMethod : size_t
  {
    LambdaZero = 0,
    AllLambdas = 1
  };

  MonteCarloTransitionMatrix(InputReader& reader) noexcept;

  void run();
  void initialize();
  void equilibrate();
  void production();
  void output();
  void cleanup();

  System& randomSystem();

  size_t numberOfCycles;
  size_t numberOfInitializationCycles;
  size_t numberOfEquilibrationCycles;
  size_t optimizeMCMovesEvery{ 5000 };
  size_t printEvery;

  std::vector<System> systems;
  RandomNumber random;
  size_t fractionalMoleculeSystem{ 0 };   // the system where the fractional molecule is located

  std::vector<std::ofstream> streams;

  BlockErrorEstimation estimation;

  std::chrono::system_clock::time_point t1;
  std::chrono::system_clock::time_point t2;
};
