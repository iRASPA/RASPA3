module;

export module spatial_decomposition_force_engine;

import std;

import double3;
import double3x3;
import running_energy;
import system;
import spatial_decomposition_settings;
import spatial_decomposition_cell_list;
import spatial_decomposition_pppm;
import spatial_decomposition_pair_kernel;
import spatial_decomposition_cluster_kernel;
import spatial_decomposition_device_step;
import spatial_decomposition_device_resident;
import spatial_decomposition_worker_team;

/**
 * \brief Multithreaded spatial-decomposition evaluation of the MD forces of one System.
 *
 * A drop-in replacement for Integrators::updateGradients for molecules in a periodic box: it fills the atomic
 * gradients of `system.atomDynamics` and returns the same RunningEnergy contributions (molecule-molecule VDW and
 * real-space Coulomb, Ewald reciprocal / self / exclusion, intramolecular terms). Instead of the O(N^2) pair loop
 * and the direct k-space sum it uses
 *
 *  - a cell list with per-sub-domain Verlet lists (CellList): each worker thread owns a balanced share of the
 *    atoms and evaluates every pair of the system exactly once with Newton's third law on a compact private copy
 *    of its atoms and their ghost images (periodic shifts resolved at list build, so the kernel has no
 *    minimum-image operation); the forces on the images are collected by the owners in a reduction phase, so no
 *    two threads ever write the same force;
 *  - specialised pair kernels (Lennard-Jones + Ewald real space) when the force field is plain 12-6
 *    Lennard-Jones with Ewald (or no) electrostatics and all atoms are fully coupled: a scalar kernel over the
 *    half lists in double with the tabulated erfc (PairPrecision::Double) or the SIMD cluster kernel
 *    (ClusterPairKernel: spatially sorted clusters, a pruned dual list, closed-form erfc) in mixed precision
 *    (PairPrecision::Mixed); otherwise the generic Potentials::potentialVDW / potentialCoulomb kernels of the
 *    rest of the code;
 *  - a smooth particle-mesh Ewald sum (PPPM) for the reciprocal-space Coulomb energy with the force field's
 *    Ewald alpha, so the self and intramolecular exclusion terms of the exact Ewald code apply unchanged;
 *  - the bonded terms and the self / exclusion corrections per molecule, handed out to the threads in chunks;
 *    with the mesh this work overlaps with the FFTs, which only thread 0 drives.
 *
 * The threads form a persistent WorkerTeam; one call of computeGradients is one task in which the phases
 * (position refresh and rebuild check, optional rebuild, pairs + mesh spreading, ghost-force and mesh reduction,
 * forward FFT || bonded, influence function, backward FFT || bonded, interpolation + scatter) are separated by
 * barriers. With one thread the same phases run inline on the calling thread, which makes the serial cell-list +
 * PPPM run the reference for the parallel ones.
 *
 * Scope (checked by supports()): no framework, external field, polarization, cross-links or fractional molecules,
 * no MD-stage particle exchange, no 'OmitInterInteractions' or dual cutoff.
 */
export class SpatialDecompositionForceEngine
{
 public:
  explicit SpatialDecompositionForceEngine(const SpatialDecompositionSettings& settings);
  ~SpatialDecompositionForceEngine();
  SpatialDecompositionForceEngine(const SpatialDecompositionForceEngine&) = delete;
  SpatialDecompositionForceEngine& operator=(const SpatialDecompositionForceEngine&) = delete;
  /// Movable (the worker team and the FFTW resources are owned through handles), so the engine is a value that
  /// can live in a container or a ForceEngine variant.
  SpatialDecompositionForceEngine(SpatialDecompositionForceEngine&&) noexcept;
  SpatialDecompositionForceEngine& operator=(SpatialDecompositionForceEngine&&) noexcept;

  struct Validation
  {
    double engineEnergy{};
    double referenceEnergy{};
    double engineReciprocalEnergy{};
    double referenceReciprocalEnergy{};
    double maximumGradientDifference{};
    double rmsGradientDifference{};
    double rmsGradient{};
    double referenceSeconds{};  ///< Wall time of the exact reference evaluation.
    double engineSeconds{};     ///< Wall time of the engine evaluation (includes the first list build).
  };

  struct Timings
  {
    std::chrono::duration<double> total{};
    std::chrono::duration<double> rebuild{};
    std::chrono::duration<double> influence{};  ///< influence-function recomputation after a cell change (NPT)
    std::chrono::duration<double> pack{};       ///< device pairs: staging of the positions for the device (all threads)
    std::chrono::duration<double> pairs{};
    std::chrono::duration<double> mesh{};
    std::chrono::duration<double> bonded{};
    std::chrono::duration<double> device{};        ///< execution time of the OpenCL pair kernel on the device
    std::chrono::duration<double> devicePrune{};   ///< execution time of the OpenCL list compaction kernel
    std::chrono::duration<double> deviceMesh{};    ///< execution time of the OpenCL mesh chain (spread, FFTs, ...)
    std::chrono::duration<double> deviceBonded{};  ///< execution time of the OpenCL bonded / exclusion kernel
    std::chrono::duration<double> deviceBuild{};   ///< wall time of the device list builds (layout, upload, build)
    std::chrono::duration<double> deviceWait{};    ///< time thread 0 waited for the device results
    std::size_t steps{};
    std::size_t rebuilds{};
    std::size_t influenceUpdates{};
  };

  /// Returns whether the engine covers the system; on false, `reason` names the unsupported feature.
  static bool supports(const System& system, std::string& reason);

  /// Builds the worker team, the cell list and the mesh for the system's current box and atoms.
  void initialize(System& system);

  /// Evaluates all forces and energies of the system (the replacement of Integrators::updateGradients). With
  /// `withVirial` the configurational molecular pressure tensor is accumulated in the same pass.
  RunningEnergy computeGradients(System& system, bool withVirial = true);

  /// Configurational part of the molecular pressure tensor at the configuration of the last
  /// computeGradients(system, true): the same quantity as System::computeMolecularPressure().second (strain
  /// derivative of the potential energy, molecular center-of-mass convention, tail correction included).
  const double3x3& molecularPressureTensor() const { return pressureTensor; }

  /// Whether the MD state is resident on the device (DeviceResident): the driver then integrates with
  /// residentVelocityVerlet and refreshes the host state with downloadResidentState when it needs it.
  bool usesResident() const { return residentEnabled; }
  /// One velocity-Verlet step with the Nose-Hoover thermostat and, when the system has one, the isotropic
  /// barostat of the system (the NVT / NPT steps of the driver with the engine forces), entirely on the device:
  /// positions, velocities and molecule records stay there in double-float arithmetic; the thermostat and
  /// barostat chains run on the host from the device reductions. The first call after a host-side evaluation
  /// (computeGradients) uploads the host state. Returns the energies of the step (potential terms, kinetic
  /// energies, thermostat and barostat energies); the molecular pressure tensor is available as after
  /// computeGradients(system, true), and the simulation box of the system follows the barostat.
  RunningEnergy residentVelocityVerlet(System& system);
  /// Copies the device state (positions, velocities, molecule records, gradients) into the system; no-op when
  /// the host copy is current.
  void downloadResidentState(System& system);

  /// Compares the engine against Integrators::updateGradients on the system's current configuration; the
  /// system's gradients are left as computed by the engine.
  Validation validate(System& system);

  std::size_t numberOfThreads() const { return settings.numberOfThreads; }
  const Timings& timings() const { return timing; }
  bool initialized() const { return initializedFlag; }

  /// Runs `body(member, numberOfMembers)` once on every member of the engine's worker team (the caller is
  /// member 0) and returns when all have finished: the driver uses it to run the host-side integrator passes
  /// over disjoint molecule ranges in parallel between two force evaluations. Before the team exists the body
  /// runs inline as a team of one.
  void runOnTeam(const std::function<void(std::size_t, std::size_t)>& body);

  std::string writeStatus() const;
  std::string writeTimings() const;

  /// Whether the specialised Lennard-Jones + tabulated Ewald kernel is in use (else the generic kernels).
  bool usesFastKernel() const { return fastKernel; }
  /// Whether the pairs are evaluated on the OpenCL device.
  bool usesDevice() const { return deviceKernel; }
  /// Whether the particle-mesh Ewald sum runs on the OpenCL device.
  bool usesDeviceMesh() const { return deviceMesh; }
  /// Whether the bonded terms and the self / exclusion corrections run on the OpenCL device.
  bool usesDeviceBonded() const { return deviceBonded; }

 private:
  SpatialDecompositionSettings settings;
  std::unique_ptr<WorkerTeam> team{};
  CellList cellList{};
  PPPM pppm{};
  bool useMesh{false};
  bool initializedFlag{false};

  // shared force accumulators in sorted order; every slot is written by the owner of the atom only
  std::vector<double> fx{}, fy{}, fz{};
  // per-thread private force buffers over the compact local atoms (owned atoms, then ghost images, padded to the
  // clusters of the cluster kernel); the image part is collected by the owners of the atoms after the pair phase
  std::vector<std::vector<LocalForce>> localForce{};

  // specialised pair kernels: per pair-type Lennard-Jones table and the tabulated Ewald real-space term, evaluated
  // by the scalar kernel (double) or on clusters (ClusterPairKernel, mixed precision; the double instantiation
  // on request for validation)
  bool fastKernel{false};
  bool fastCoulomb{false};
  std::size_t numberOfPseudoAtomTypes{0};
  std::vector<LennardJonesPair> lennardJones{};
  EwaldRealSpaceTable ewaldTable{};
  ClusterPairKernel<double> clusterKernelDouble{};
  ClusterPairKernel<float> clusterKernelMixed{};
  // the pairs on a device (settings.pairDevice != CPU): the device builds its own list from the cell binning, the
  // per-domain Verlet lists and ghost images are not built, and the device work overlaps with the mesh and bonded
  // phases of the threads
  bool deviceKernel{false};
  DeviceStep devicePairs{};
  // with the device pairs: the mesh (settings.deviceMesh) and the per-molecule terms (settings.deviceBonded, when
  // the device bonded kernels cover the system's intramolecular potentials) run on the device as well; the device results of
  // the mesh that computeGradients needs after the step
  bool deviceMesh{false};
  bool deviceBonded{false};
  std::string deviceBondedFallback{};  ///< why the bonded work stayed on the host (status line)
  double deviceSingleIonSum{0.0};
  double3x3 deviceReciprocalStrain{};
  double3x3 deviceSingleIonStrain{};
  // the resident integrator (settings.resident, with the complete step on the device)
  bool residentEnabled{false};
  bool residentValid{false};  ///< the device holds the authoritative state (else the host does)
  bool residentHostCurrent{true};  ///< the host copy equals the device state
  std::string residentFallback{};  ///< why the integration stayed on the host (status line)
  DeviceResident resident{};
  DeviceResident::Scaling residentPendingScale{};  ///< thermostat factor of the last step, not yet applied
  DeviceResident::Kinetic residentKinetic{};       ///< kinetic energies and virial of the current (scaled) velocities
  std::size_t residentSteps{0};
  /// The isotropic barostat of the resident step: the chain and the first half-kick of the cell velocity, the
  /// propagation of the cell (and the force-field parameters that follow it); returns the coupling of the step.
  DeviceResident::Coupling residentBarostatFirstHalf(System& system, double pressureVirialTrace);
  bool mixedPrecision() const { return settings.pairPrecision == PairPrecision::Mixed; }
  bool usesClusterKernel() const { return mixedPrecision() || settings.clusterKernelForDouble; }
  std::vector<std::uint8_t> rebuildRequested{};
  std::uint8_t influenceUpdateRequested{0};  ///< set by thread 0 in phase 0, read by all after the barrier
  std::vector<RunningEnergy> threadEnergies{};
  std::vector<double3x3> threadStrain{};      ///< pair + exclusion strain derivatives per thread
  std::vector<double3x3> threadCorrection{};  ///< atomic-to-molecular virial correction per thread
  bool virialRequested{false};
  double3x3 pressureTensor{};

  // bonded / exclusion work by molecule, handed out in chunks through the team's work counter so that it can be
  // done by whichever threads are free (it overlaps with the FFTs of thread 0 when the mesh is in use)
  std::size_t bondedChunkSize{1};
  std::size_t partitionedMolecules{0};
  /// Energies, strain derivatives and virial corrections per chunk: whichever thread takes a chunk, the sums are
  /// reduced in chunk order, so that the results do not depend on the scheduling (bit-reproducible runs).
  std::vector<RunningEnergy> chunkEnergies{};
  std::vector<double3x3> chunkStrain{}, chunkCorrection{};
  /// Mass-weighted center of mass of the molecule of every atom (original atom order), set in the bonded work and
  /// used for the atomic-to-molecular virial correction of the pair + mesh gradients in the scatter phase.
  std::vector<double3> atomCenterOfMass{};

  double reciprocalEnergy{0.0};
  Timings timing{};

  /// Position-independent sums over the atoms used by finishStep every step (net charge for the Bogusz
  /// correction, scaled atom count per type for the tail virial), valid for `staticSumsAtoms` atoms.
  double staticNetCharge{0.0};
  std::vector<double> staticScaledCountPerType{};
  std::size_t staticSumsAtoms{std::numeric_limits<std::size_t>::max()};
  void refreshStaticSums(const System& system);

  void prepareBondedWork(const System& system);
  void refreshCutoffs(System& system);
  RunningEnergy finishStep(const System& system, RunningEnergy total, double3x3 strain, double3x3 correction);
  void residentRebuild(System& system);
  void step(std::size_t thread, System& system);
  void pairPhase(std::size_t thread, const System& system, RunningEnergy& energy);
  template <bool Fast>
  void pairLoop(std::size_t thread, const System& system, RunningEnergy& energy);
  void clusterPairLoop(std::size_t thread, RunningEnergy& energy);
  void buildClusterLists(std::size_t thread);
  void collectGhostForces(std::size_t thread);
  void prepareKernel(const System& system);
  void bondedWork(System& system);
  void scatterPhase(std::size_t thread, System& system);
};
