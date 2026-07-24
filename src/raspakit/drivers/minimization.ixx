module;

export module minimization;

import std;
import randomnumbers;

import archive;
import input_reader;
import system;
import generalized_hessian;
import minimization_dof_layout;
import minimization_evaluate_derivatives;
import minimization_options;
import minimization_baker;
import minimization_generalized_coordinates;
import elastic_constants;
import normal_modes;
import phonon_dynamical_matrix;

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
  std::optional<ElasticConstantsResult> elasticConstants{};
  std::optional<NormalModesResult> normalModes{};
  std::optional<PhononDispersionResult> phononDispersion{};
  std::optional<PhononDensityOfStates> phononDensityOfStates{};

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const MinimizationSystemResult& r);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, MinimizationSystemResult& r);
};

/** Baker eigenvector-following minimization driver with optional logarithmic cell DOFs. */
export struct Minimization
{
  Minimization() = default;
  explicit Minimization(InputReader& inputReader);
  Minimization(const MinimizationOptions& options, std::vector<System> systems, bool outputToFiles = false);

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

  /**
   * \brief Writes the binary restart file. The Baker temporaries (generalized Hessian, gradient,
   *        DOF layout) are not stored; they are rebuilt from the serialized system state.
   */
  void writeBinaryRestartFile() noexcept;

  /**
   * \brief Writes the periodic binary restart file and services a pending shutdown signal.
   */
  void checkpointIfDue(std::size_t currentCycle);

  std::uint64_t versionNumber{1};  ///< Version number for serialization.

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

  /// Binary-restart resume point of the Run stage: the system currently being minimized and the
  /// next Baker iteration of that system (systems before 'runSystemIndex' are already done and
  /// their results are stored).
  std::size_t runSystemIndex{0};
  std::size_t runIteration{0};

  std::optional<std::size_t> randomSeed{std::nullopt};
  RandomNumber random{std::nullopt};
  std::size_t fractionalMoleculeSystem{0};

  std::vector<System> systems;
  std::vector<MinimizationSystemResult> results;
  std::vector<std::ofstream> streams;

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const Minimization& m);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, Minimization& m);
};
