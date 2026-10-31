module;

export module cbmc_device_environment;

import std;

import atom;
import forcefield;
import simulationbox;
import spatial_decomposition_settings;
import spatial_decomposition_device_context;

/**
 * \file environment.ixx
 * \brief The device side of the multiple-try CBMC grower: the background of a move and the batch evaluation of
 * trial molecules against it, in single precision, on an OpenCL, Metal or CUDA device.
 *
 * The environment is a proposal engine only (the "option 2" precision design): it evaluates the external
 * energies of many complete trial molecules at once in float; the grower built on it (cbmc_multi_trial_device)
 * selects among them with those energies, and the host re-evaluates the one selected trial in double, so that
 * the system's running energies stay in double and the acceptance carries the ratio of the two Boltzmann
 * factors as a correction (see MultiTrialGrowResult::logProposalCorrection).
 *
 * Per move the grower uploads the background (the framework atoms, cached while they and the box do not change,
 * and the molecule atoms), the trial atoms of the batch, and reads back three partials per trial set: the
 * Lennard-Jones energy, the real-space Ewald energy, and the number of overlapping pairs.
 */
export namespace CBMC
{
/// Mirrors the GrowParameters struct of the kernel source (4-byte members only, so the layouts agree).
struct DeviceGrowParameters
{
  float cell[9]{};
  float inverseCell[9]{};
  float cutOffFrameworkVDWSquared{0.0f};
  float cutOffMoleculeVDWSquared{0.0f};
  float cutOffChargeSquared{0.0f};
  float alpha{0.0f};
  float coulombFactor{0.0f};
  float overlapCriteria{0.0f};
  std::uint32_t useCharge{0};
  std::uint32_t orthorhombic{0};
  std::uint32_t numberOfTypes{0};
  std::uint32_t numberOfBackground{0};
  std::uint32_t atomsPerSet{0};
  std::uint32_t trialMoleculeId{0};
  std::uint32_t skipMoleculeId{0};
  std::uint32_t padding[5]{};
};
static_assert(sizeof(DeviceGrowParameters) == 144);

/// The device energies of one trial set, as doubles (converted from the float partials).
struct DeviceTrialEnergy
{
  double vdw{0.0};
  double charge{0.0};
  bool overlap{false};
  [[nodiscard]] double total() const { return vdw + charge; }
};

class DeviceGrowEnvironment
{
 public:
  /// The largest trial molecule the kernel holds in local memory (GROW_MAX_ATOMS).
  static constexpr std::size_t maximumAtomsPerSet = 32;
  /// Lanes per trial set (GROW_GROUP) and partials per set (GROW_PARTIALS) of the kernel source.
  static constexpr std::size_t groupSize = 64;
  static constexpr std::size_t partials = 3;

  /// Whether a device of the kind is available (its runtime is initialized on first use).
  static bool available(PairDevice device);
  /// Why no device of the kind is available (empty when one is).
  static std::string unavailableReason(PairDevice device);

  /// Whether the device kernel evaluates the pair potentials of the force field as the host does for the given
  /// types: plain (unswitched) 12-6 Lennard-Jones or no interaction between every pair of 'presentTypes', and
  /// Ewald or no electrostatics. Else 'reason' says why not.
  static bool supportsForceField(const ForceField& forceField, std::span<const std::uint8_t> presentTypes,
                                 std::string& reason);

  /// Creates the context and compiles the kernel; throws when no device of the kind is available.
  explicit DeviceGrowEnvironment(PairDevice device);
  ~DeviceGrowEnvironment();
  DeviceGrowEnvironment(const DeviceGrowEnvironment&) = delete;
  DeviceGrowEnvironment& operator=(const DeviceGrowEnvironment&) = delete;

  PairDevice device() const { return kind; }
  std::string deviceName() const;

  /// Fixes the pair tables, the cutoffs and the Ewald parameters from the force field, and the box. The pair
  /// table is rebuilt on the host every call and re-uploaded only when its contents changed: the environment is
  /// shared by every system of the process, so the identity of the force-field object says nothing about what
  /// the device holds (another system, or a force field rebuilt in place, may have the same address).
  void setForceField(const ForceField& forceField, const SimulationBox& box);

  /// Uploads the background of a move: the framework atoms (re-uploaded only when their packed contents changed)
  /// followed by the molecule atoms (always). Call after 'setForceField'.
  void setBackground(std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms);

  /// Evaluates the external energies of 'trialAtoms.size() / atomsPerSet' trial sets (set-major) against the
  /// background, leaving out the background atoms of molecule 'trialMoleculeId' and of 'skipMolecule'. Blocking.
  void evaluate(std::span<const Atom> trialAtoms, std::size_t atomsPerSet, std::uint32_t trialMoleculeId,
                std::optional<std::size_t> skipMolecule, std::vector<DeviceTrialEnergy>& energies);

 private:
  PairDevice kind;
  std::unique_ptr<DeviceContext> context;
  DeviceKernel kernel{};

  DeviceGrowParameters parameters{};
  DeviceBufferOwner parameterBuffer;
  DeviceBufferOwner lennardJonesBuffer;
  DeviceBufferOwner backgroundBuffer;
  DeviceBufferOwner backgroundTypeBuffer;
  DeviceBufferOwner backgroundMoleculeBuffer;
  DeviceBufferOwner trialBuffer;
  DeviceBufferOwner trialTypeBuffer;
  DeviceBufferOwner partialBuffer;

  // What the device currently holds, by content (the environment is a process-wide singleton: keys by object
  // address or span size alone would alias between systems and miss in-place changes).
  std::vector<std::array<float, 4>> uploadedLennardJones;
  std::size_t uploadedTypes{0};
  std::vector<std::array<float, 4>> uploadedFrameworkPositions;
  std::vector<std::uint32_t> uploadedFrameworkTypes;
  bool frameworkUploaded{false};
  std::size_t backgroundCapacity{0};
  std::size_t trialCapacity{0};
  std::size_t partialCapacity{0};
  std::size_t numberOfFramework{0};

  // host staging (kept alive across the asynchronous writes of one evaluation)
  std::vector<std::array<float, 4>> hostPositions;
  std::vector<std::uint32_t> hostTypes;
  std::vector<std::uint32_t> hostMolecules;
  std::vector<std::array<float, 4>> hostTrial;
  std::vector<std::uint32_t> hostTrialTypes;
  std::vector<float> hostPartials;

  void ensureBackgroundCapacity(std::size_t atoms);
};
}  // namespace CBMC
