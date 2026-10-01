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
 *  - a specialised Lennard-Jones + tabulated Ewald real-space kernel (EwaldRealSpaceTable, in r^2) when the force field
 *    is plain 12-6 Lennard-Jones with Ewald (or no) electrostatics and all atoms are fully coupled; otherwise the
 *    generic Potentials::potentialVDW / potentialCoulomb kernels of the rest of the code;
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
    std::chrono::duration<double> pairs{};
    std::chrono::duration<double> mesh{};
    std::chrono::duration<double> bonded{};
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

  /// Compares the engine against Integrators::updateGradients on the system's current configuration; the
  /// system's gradients are left as computed by the engine.
  Validation validate(System& system);

  std::size_t numberOfThreads() const { return settings.numberOfThreads; }
  const Timings& timings() const { return timing; }
  bool initialized() const { return initializedFlag; }

  std::string writeStatus() const;
  std::string writeTimings() const;

  /// Whether the specialised Lennard-Jones + tabulated Ewald kernel is in use (else the generic kernels).
  bool usesFastKernel() const { return fastKernel; }

 private:
  SpatialDecompositionSettings settings;
  std::unique_ptr<WorkerTeam> team{};
  CellList cellList{};
  PPPM pppm{};
  bool useMesh{false};
  bool initializedFlag{false};

  // shared force accumulators in sorted order; every slot is written by the owner of the atom only
  std::vector<double> fx{}, fy{}, fz{};
  // per-thread private force buffers over the compact local atoms (owned atoms, then ghost images); the image
  // part is collected by the owners of the atoms after the pair phase
  struct LocalForce
  {
    double x{0.0};
    double y{0.0};
    double z{0.0};
  };
  std::vector<std::vector<LocalForce>> localForce{};

  // specialised pair kernel: per pair-type Lennard-Jones table and the tabulated Ewald real-space term
  bool fastKernel{false};
  bool fastCoulomb{false};
  std::size_t numberOfPseudoAtomTypes{0};
  std::vector<LennardJonesPair> lennardJones{};
  EwaldRealSpaceTable ewaldTable{};
  std::vector<std::uint8_t> rebuildRequested{};
  std::uint8_t influenceUpdateRequested{0};  ///< set by thread 0 in phase 0, read by all after the barrier
  std::vector<RunningEnergy> threadEnergies{};
  std::vector<double3x3> threadStrain{};      ///< pair + exclusion strain derivatives per thread
  std::vector<double3x3> threadCorrection{};  ///< atomic-to-molecular virial correction per thread
  bool virialRequested{false};
  double3x3 pressureTensor{};

  // bonded / exclusion work by molecule, handed out in chunks through an atomic counter so that it can be done by
  // whichever threads are free (it overlaps with the FFTs of thread 0 when the mesh is in use)
  std::atomic<std::size_t> nextBondedChunk{0};
  std::size_t bondedChunkSize{1};
  std::size_t partitionedMolecules{0};
  /// Mass-weighted center of mass of the molecule of every atom (original atom order), set in the bonded work and
  /// used for the atomic-to-molecular virial correction of the pair + mesh gradients in the scatter phase.
  std::vector<double3> atomCenterOfMass{};

  double reciprocalEnergy{0.0};
  Timings timing{};

  void prepareBondedWork(const System& system);
  void step(std::size_t thread, System& system);
  void pairPhase(std::size_t thread, const System& system, RunningEnergy& energy);
  template <bool Fast>
  void pairLoop(std::size_t thread, const System& system, RunningEnergy& energy);
  void collectGhostForces(std::size_t thread);
  void prepareKernel(const System& system);
  void bondedWork(std::size_t thread, System& system, RunningEnergy& energy);
  void scatterPhase(std::size_t thread, System& system);
};
