module;

export module force_engine;

import std;

import double3x3;
import running_energy;
import system;
import spatial_decomposition_settings;
import spatial_decomposition_force_engine;

/**
 * \brief The MD force engine of one System as a value: one of the available implementations behind a single
 * interface.
 *
 * The driver holds a ForceEngine per system, in a plain vector, and never deals with the implementation it
 * wraps. The implementations are the alternatives of a std::variant (every one a movable value type); the calls
 * are forwarded with std::visit, so adding an engine (e.g. an OpenCL engine for the same spatial-decomposition
 * data layout) means adding an alternative and a constructor case, not a pointer-based class hierarchy.
 *
 * All alternatives share the data types of the interface: the Validation and Timings records and the
 * SpatialDecompositionSettings they are built from.
 */
export class ForceEngine
{
 public:
  using Validation = SpatialDecompositionForceEngine::Validation;
  using Timings = SpatialDecompositionForceEngine::Timings;

  /// Builds the engine for `settings` (currently the multithreaded CPU engine).
  explicit ForceEngine(const SpatialDecompositionSettings& settings) : implementation(std::in_place_index<0>, settings)
  {
  }

  ForceEngine(const ForceEngine&) = delete;
  ForceEngine& operator=(const ForceEngine&) = delete;
  ForceEngine(ForceEngine&&) noexcept = default;
  ForceEngine& operator=(ForceEngine&&) noexcept = default;

  /// Returns whether an engine covers the system; on false, `reason` names the unsupported feature.
  static bool supports(const System& system, std::string& reason)
  {
    return SpatialDecompositionForceEngine::supports(system, reason);
  }

  /// Builds the engine's data structures for the system's current box and atoms.
  void initialize(System& system)
  {
    std::visit([&](auto& engine) { engine.initialize(system); }, implementation);
  }

  /// Evaluates all forces and energies of the system; with `withVirial` also the molecular pressure tensor.
  RunningEnergy computeGradients(System& system, bool withVirial = true)
  {
    return std::visit([&](auto& engine) { return engine.computeGradients(system, withVirial); }, implementation);
  }

  /// Configurational molecular pressure tensor of the last computeGradients(system, true).
  const double3x3& molecularPressureTensor() const
  {
    return std::visit([](const auto& engine) -> const double3x3& { return engine.molecularPressureTensor(); },
                      implementation);
  }

  /// Whether the engine integrates on the device with the MD state resident there (see
  /// SpatialDecompositionForceEngine::usesResident).
  bool usesResident() const
  {
    return std::visit([](const auto& engine) { return engine.usesResident(); }, implementation);
  }
  /// One velocity-Verlet step on the device (usesResident() must hold).
  RunningEnergy residentVelocityVerlet(System& system)
  {
    return std::visit([&](auto& engine) { return engine.residentVelocityVerlet(system); }, implementation);
  }
  /// Refreshes the host state of the system from the device (no-op when the host copy is current).
  void downloadResidentState(System& system)
  {
    std::visit([&](auto& engine) { engine.downloadResidentState(system); }, implementation);
  }

  /// Compares the engine against Integrators::updateGradients on the system's current configuration.
  Validation validate(System& system)
  {
    return std::visit([&](auto& engine) { return engine.validate(system); }, implementation);
  }

  std::size_t numberOfThreads() const
  {
    return std::visit([](const auto& engine) { return engine.numberOfThreads(); }, implementation);
  }
  const Timings& timings() const
  {
    return std::visit([](const auto& engine) -> const Timings& { return engine.timings(); }, implementation);
  }
  bool initialized() const
  {
    return std::visit([](const auto& engine) { return engine.initialized(); }, implementation);
  }
  bool usesFastKernel() const
  {
    return std::visit([](const auto& engine) { return engine.usesFastKernel(); }, implementation);
  }

  std::string writeStatus() const
  {
    return std::visit([](const auto& engine) { return engine.writeStatus(); }, implementation);
  }
  std::string writeTimings() const
  {
    return std::visit([](const auto& engine) { return engine.writeTimings(); }, implementation);
  }

 private:
  std::variant<SpatialDecompositionForceEngine> implementation;
};
