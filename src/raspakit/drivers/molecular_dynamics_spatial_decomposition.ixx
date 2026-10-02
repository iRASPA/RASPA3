module;

export module molecular_dynamics_spatial_decomposition;

import std;

import randomnumbers;
import averages;
import system;
import input_reader;
import archive;
import json;
import double3x3;
import running_energy;
import spatial_decomposition_settings;
import force_engine;

/**
 * \brief Molecular dynamics with the multithreaded spatial-decomposition force engine.
 *
 * Selected with 'SimulationType' : 'MolecularDynamicsSpatialDecomposition'. The pre-initialization and
 * initialization stages are the ordinary serial Monte Carlo stages of the MolecularDynamics driver (exact O(N^2)
 * pair sums and direct Ewald); at the start of the equilibration stage one SpatialDecompositionForceEngine per
 * system is built and from then on every force evaluation (velocity-Verlet steps in NVE / NVT, the thermobarostat
 * steps in NPT / NPT-PR, and the stage-boundary recomputations) goes through the engine: cell-list Verlet
 * neighbour lists over `NumberOfThreads` sub-domains and a particle-mesh Ewald sum. With one thread it is the
 * serial O(N) cell-list + PPPM variant.
 *
 * The pressure is sampled every cycle from the engine's virial; with a thermobarostat the reported tensor is the
 * estimator of the barostat's coupling ('BarostatCoupling', molecular by default), so its average is the external
 * pressure. The per-component energy decomposition comes from the engine's running energies for single-component
 * systems and from the exact code every 'PrintEvery' cycles otherwise. The engine is not archived: after a binary
 * restart the neighbour lists and the mesh are rebuilt from the settings.
 */
export struct MolecularDynamicsSpatialDecomposition
{
  enum class SimulationStage : std::size_t
  {
    Uninitialized = 0,
    PreInitialization = 1,
    Initialization = 2,
    Equilibration = 3,
    Production = 4
  };

  MolecularDynamicsSpatialDecomposition();
  MolecularDynamicsSpatialDecomposition(const MolecularDynamicsSpatialDecomposition&) = delete;
  MolecularDynamicsSpatialDecomposition& operator=(const MolecularDynamicsSpatialDecomposition&) = delete;

  MolecularDynamicsSpatialDecomposition(InputReader& reader) noexcept;

  std::uint64_t versionNumber{1};

  bool outputToFiles{true};
  RandomNumber random;

  std::size_t numberOfProductionCycles{0};
  std::size_t numberOfSteps{0};
  std::size_t numberOfPreInitializationCycles{0};
  std::size_t numberOfInitializationCycles{0};
  std::size_t numberOfEquilibrationCycles{0};

  std::size_t printEvery{5000};
  std::size_t writeRestartEvery{5000};
  std::size_t writeBinaryRestartEvery{5000};
  std::size_t rescaleWangLandauEvery{5000};
  std::size_t optimizeMCMovesEvery{5000};

  std::size_t currentCycle{};
  std::size_t absoluteCurrentCycle{};
  SimulationStage simulationStage{SimulationStage::Uninitialized};

  std::vector<System> systems;
  std::size_t fractionalMoleculeSystem{0};

  SpatialDecompositionSettings engineSettings{};
  std::vector<ForceEngine> engines;  ///< One per system, built lazily.

  /// Running sum and count (per system) of the pressure tensor the barostat drives to the external pressure,
  /// over the steps since the last status report (printed with every report, then reset). Not part of the restart.
  std::vector<double3x3> barostatPressureWindowSum;
  std::vector<std::size_t> barostatPressureWindowCount;

  std::vector<std::ofstream> streams;
  std::vector<nlohmann::json> outputJsons;

  BlockErrorEstimation estimation;

  std::chrono::duration<double> totalSimulationTime{0};

  void createOutputFiles();
  void checkpointIfDue(std::size_t currentCycle);

  void run();
  void setup();
  void tearDown();

  void preInitialize();
  void initialize();
  void equilibrate();
  void production();
  void output();

  /// Builds the engines (if not yet built), evaluates the forces and reports the engine status and the
  /// comparison with the exact Ewald / pair code at the start of an MD stage.
  void startEngines(std::string_view stageName);

  /// Builds and initializes missing engines silently (after a binary restart into an MD stage) and refreshes the
  /// forces and the pressure tensor through them.
  void ensureEngines();

  /// One MD step of a system with the engine forces (velocity Verlet, or the thermobarostat step in NPT).
  RunningEnergy molecularDynamicsStep(std::size_t systemId);

  /// Recomputes forces, energies and the pressure tensor of a system through its engine.
  void recomputeGradients(std::size_t systemId);

  /// Sets 'currentExcessPressureTensor' of a system from the engine's molecular virial; with a thermobarostat the
  /// tensor is the one conjugate to the barostat coupling, and it is accumulated into the status-report window.
  void updateReportedPressure(std::size_t systemId, bool accumulate);

  /// The window-averaged barostat pressure tensor block of the status report (empty without a thermobarostat).
  std::string writeBarostatPressureWindow(std::size_t systemId);

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive,
                                            const MolecularDynamicsSpatialDecomposition& md);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, MolecularDynamicsSpatialDecomposition& md);
};
