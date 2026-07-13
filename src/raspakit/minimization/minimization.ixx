module;

export module minimization;

import std;

import input_reader;
import system;
import generalized_hessian;
import minimization_dof_layout;
import minimization_evaluate_derivatives;
import minimization_options;
import minimization_baker;
import minimization_generalized_coordinates;

export struct MinimizationSystemResult
{
  bool converged{};
  std::size_t iterations{};
  double initialEnergy{};
  double finalEnergy{};
  double rmsGradient{};
  double maxGradient{};
  std::size_t negativeModes{};
  std::size_t zeroModes{};
};

/** Fixed-cell Baker eigenvector-following minimization driver. */
export struct Minimization
{
  Minimization() = default;
  explicit Minimization(InputReader &inputReader);
  Minimization(const MinimizationOptions &options, std::vector<System> systems, bool outputToFiles = false);

  void run();
  void setup();
  void output();
  void tearDown();

  MinimizationOptions options{};
  bool outputToFiles{true};
  std::vector<System> systems;
  std::vector<MinimizationSystemResult> results;
  std::vector<std::ofstream> streams;
};
