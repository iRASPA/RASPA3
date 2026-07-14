module;

export module minimization;

import std;
import randomnumbers;

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

/** Baker eigenvector-following minimization driver with optional logarithmic cell DOFs. */
export struct Minimization
{
  Minimization() = default;
  explicit Minimization(InputReader &inputReader);
  Minimization(const MinimizationOptions &options, std::vector<System> systems, bool outputToFiles = false);

  enum class SimulationStage : std::size_t
  {
    Uninitialized = 0,
    PreInitialization = 1,
    Initialization = 2,
    Equilibration = 3,
    Run = 4
  };

  void run();
  void setup();
  void performCycle();
  void preInitialize();
  void initialize();
  void equilibrate();
  void runPhase();
  void output();
  void tearDown();

  MinimizationOptions options{};
  bool outputToFiles{true};

  std::size_t numberOfPreInitializationCycles{0};
  std::size_t numberOfInitializationCycles{0};
  std::size_t numberOfEquilibrationCycles{0};

  std::size_t printEvery{5000};
  std::size_t writeBinaryRestartEvery{5000};
  std::size_t rescaleWangLandauEvery{5000};
  std::size_t optimizeMCMovesEvery{5000};

  std::size_t currentCycle{0};
  std::size_t absoluteCurrentCycle{0};

  SimulationStage simulationStage{SimulationStage::Uninitialized};

  std::optional<std::size_t> randomSeed{std::nullopt};
  RandomNumber random{std::nullopt};
  std::size_t fractionalMoleculeSystem{0};

  std::vector<System> systems;
  std::vector<MinimizationSystemResult> results;
  std::vector<std::ofstream> streams;
};
