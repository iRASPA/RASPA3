module;

export module spatial_decomposition_device_resident;

import std;

import double3;
import simd_quatd;
import simulationbox;
import system;
import spatial_decomposition_cell_list;
import spatial_decomposition_device_context;
import spatial_decomposition_device_step;

/**
 * \brief The MD state of one System resident on the device, and the velocity-Verlet integrator over it.
 *
 * With the device pair kernel, mesh and bonded terms the forces of a step are complete on the device; this class
 * keeps the integration variables there as well, so that a step moves no atom data between host and device:
 *
 *  - Positions, velocities, centers of mass, quaternions and their conjugate momenta are double-float numbers
 *    (hi + lo floats, ~48 bits of mantissa: "emulated double"), the forces are the single-precision device
 *    forces; the per-step arithmetic is that of Integrators::updateVelocities / updatePositions /
 *    noSquishFreeRotorOrderTwo / createCartesianPositions / updateCenterOfMassAndQuaternionGradients for rigid
 *    and flexible molecules (resident_kernel_source.cpp).
 *  - The first half of a step (thermostat scaling, kick, drift, free rotor, cartesian positions) also packs the
 *    slot positions for the force kernels and reduces the largest displacement since the list build and since
 *    the list compaction; the host reads those two numbers and decides whether to rebuild or compact, exactly
 *    as the host path does with the atom positions.
 *  - The second half (molecular gradients and torques, kick, kinetic energies) kicks into a second set of
 *    velocity buffers that the host accepts with acceptVelocities(); a lane-list overflow then just repeats the
 *    chain and the second half with the corrected forces.
 *  - The Nose-Hoover chain stays on the host: it needs the two kinetic-energy sums of the step, which are
 *    reduced on the device in double-float, and returns the scale factors that enter the next first half.
 *
 * The host state of the System is refreshed by download() (positions, velocities, molecule records, gradients)
 * when the driver needs it (property sampling, restart files); upload() makes the device state from the host
 * state. The slot layout follows the list builds of DeviceStep through setLayout(); the constant tables (masses,
 * molecule and component records, body-fixed reference positions) are written by upload().
 *
 * Covers rigid and flexible molecules (no semi-flexible components, no framework, no particle exchange) in an
 * NVE or NVT ensemble. Movable value type.
 */
export class DeviceResident
{
 public:
  DeviceResident() = default;
  DeviceResident(const DeviceResident&) = delete;
  DeviceResident& operator=(const DeviceResident&) = delete;
  DeviceResident(DeviceResident&&) noexcept = default;
  DeviceResident& operator=(DeviceResident&&) noexcept = default;

  /// Whether the resident integrator covers the system (else `reason`).
  static bool supports(const System& system, std::string& reason);

  /// Compiles the kernels on the context of the device step.
  void initialize(DeviceContext& context);
  bool initialized() const { return context != nullptr; }

  /// Host -> device: the constant tables and the dynamic state of the system (positions, velocities, molecule
  /// records), the current gradients into the force buffer of the step, and the slot layout of the step.
  void upload(const System& system, DeviceStep& step, const CellList& cells);
  /// The slot layout and the box translations of the atoms after a list build (the device state is unchanged).
  void setLayout(DeviceStep& step, const CellList& cells, const SimulationBox& box);

  /// Scale factors of the velocities (translational, rotational) applied at the start of the next first half.
  struct Scaling
  {
    double translational{1.0};
    double rotational{1.0};
  };

  /// First half of the step: scaling and kick with the forces in the step's force buffer, drift, free rotor,
  /// cartesian positions, slot positions and the displacement reduction (read back; complete after the wait on
  /// the returned event).
  DeviceEvent enqueueFirstHalf(const Scaling& scaling, double timeStep);
  /// Largest squared displacements since the list build and since the list compaction, after the wait.
  struct Displacement
  {
    double sinceBuild{0.0};
    double sinceCompaction{0.0};
  };
  Displacement collectDisplacement() const;
  /// Slot positions and displacements only (after a list build: the positions did not change, the layout did).
  DeviceEvent enqueuePack();

  /// Second half of the step with the forces in the step's force buffer: molecular gradients and torques, kick
  /// into the pending velocity buffers, kinetic energies (read back; complete after the wait).
  DeviceEvent enqueueSecondHalf(double timeStep);
  /// Kinetic energies of the pending velocities, after the wait.
  struct Kinetic
  {
    double translational{0.0};
    double rotational{0.0};
  };
  Kinetic collectKinetic() const;
  /// Makes the pending velocities the current ones (the second half is accepted).
  void acceptVelocities();

  /// Scales the current velocities (translational, rotational) on the device.
  void enqueueScale(const Scaling& scaling);

  /// Device -> host: positions, velocities, molecule records and the gradients of the last step into the system
  /// (blocking).
  void download(System& system, const DeviceStep& step, const CellList& cells);
  /// Device -> host: the positions only (blocking; for the binning of a list rebuild).
  void downloadPositions(System& system);

  std::size_t numberOfAtoms() const { return atoms; }
  std::size_t numberOfMolecules() const { return molecules; }

 private:
  DeviceContext* context{nullptr};
  DeviceKernel atomsA{}, moleculesA{}, pack{}, torques{}, atomsB{}, moleculesB{}, scaleAtoms{}, scaleMolecules{};
  DeviceStep::ResidentTargets targets{};

  std::size_t atoms{0}, molecules{0}, components{0}, references{0};
  std::size_t atomGroups{0}, moleculeGroups{0};

  // atoms (system order); the velocities in two sets (current and pending, see acceptVelocities)
  DeviceBufferOwner positionHi{}, positionLo{};
  std::array<DeviceBufferOwner, 2> velocityHi{}, velocityLo{};
  DeviceBufferOwner atomMass{}, atomInfo{}, slotOfAtom{}, translationHi{}, translationLo{};
  // molecules
  DeviceBufferOwner comHi{}, comLo{};
  std::array<DeviceBufferOwner, 2> moleculeVelocityHi{}, moleculeVelocityLo{};
  DeviceBufferOwner orientationHi{}, orientationLo{};
  std::array<DeviceBufferOwner, 2> momentumHi{}, momentumLo{};
  DeviceBufferOwner moleculeGradient{}, moleculeTorque{}, moleculeMass{}, moleculeInfo{};
  // components
  DeviceBufferOwner componentInertia{}, componentReferenceOffset{}, referenceHi{}, referenceLo{};
  // reductions (float2 per work-group)
  DeviceBufferOwner packPartials{}, kineticPartials{}, moleculePartials{};
  std::size_t current{0};  ///< index of the current velocity buffers (the other set is the pending one)

  std::vector<float> packHost{}, kineticHost{};
  std::vector<float> stage{};           ///< host staging of the uploads and downloads (float4 arrays)
  std::vector<std::uint32_t> stageU{};  ///< host staging of the uint arrays

  void ensureBuffers(const System& system);
  void writeFloat4(DeviceBufferOwner& buffer, std::span<const float> values);
  void readFloat4(const DeviceBufferOwner& buffer, std::size_t count, std::vector<float>& into);
};
