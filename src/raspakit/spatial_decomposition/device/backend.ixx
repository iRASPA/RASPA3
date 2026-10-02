module;

export module spatial_decomposition_device_backend;

import std;

import double3x3;
import int3;
import simulationbox;
import spatial_decomposition_device_bonded_topology;

/**
 * \file backend.ixx
 * \brief The interface between the backend-neutral device step (DeviceStep, step.ixx) and a device API.
 *
 * DeviceStep owns everything that does not depend on the API: the slot layout, the host staging of positions,
 * the decision when to rebuild, compact or sample, the growth of the list capacities, the double reductions of
 * the partial sums, the timers and the status report. A DeviceBackend owns the API objects (context, queue,
 * program, buffers, events) and turns the step's requests into uploads and kernel launches. The kernels
 * themselves are the shared sources of spatial_decomposition_device_kernels, compiled by the backend with its
 * dialect header.
 *
 * Conventions:
 *  - All calls are made from one thread (thread 0 of the engine). Enqueue calls are asynchronous; waitBuild /
 *    waitStep complete them, finish() completes everything (used for the synchronous timing samples).
 *  - The step tracks the capacities; an allocate* call replaces the buffers of that family (contents are lost)
 *    and invalidates the mapped host pointers of the family (allocateSlots: positions and forces,
 *    allocateRelative: relative positions). The step unmaps before growing and maps again afterwards.
 *  - The mapped pointers are host-visible staging of the device buffers: writable (positions, relative
 *    positions) between unmapAll() of the previous step and the enqueue calls, readable (forces) after waitStep().
 *    A backend without unified memory implements them with pinned staging and copies in unmapAll / mapForces.
 *  - The pair, mesh and bonded chains write their per-group partial sums to device buffers; enqueueReadPartials
 *    reads them back (the pair partials into the span given, the mesh and bonded partials into backend-owned
 *    storage read by collectMesh / collectBonded) and marks the completion of the step.
 */

/// Launch geometry fixed by the pair kernel source (CLUSTER_I, CLUSTER_J, GROUP_SIZE, the partial sums per group).
export namespace DeviceKernelLayout
{
constexpr std::size_t clusterI = 8;
constexpr std::size_t clusterJ = 4;
constexpr std::size_t pairGroupSize = 32;  ///< work-items (lanes) per i-cluster: one per (a, b) of 8 x 4
constexpr std::size_t pairPartials = 11;   ///< partial sums per i-cluster: 2 energies + 9 strain components
constexpr std::uint32_t noAtom = std::numeric_limits<std::uint32_t>::max();
}  // namespace DeviceKernelLayout

/// Mirrors the Parameters struct of the pair kernel source (4-byte members only, so the layouts agree).
export struct DevicePairParameters
{
  float cell[9];
  float inverseCell[9];
  float cutOffVDWSquared{0.0f};
  float cutOffChargeSquared{0.0f};
  float alpha{0.0f};
  float alphaSquared{0.0f};
  float alphaOverSqrtPi{0.0f};
  float coulombFactor{0.0f};
  float innerCutoffSquared{0.0f};
  std::uint32_t useCharge{0};
  std::uint32_t orthorhombic{0};
  std::uint32_t numberOfTypes{0};
  std::uint32_t blocksPerCluster{0};
  std::uint32_t pairsPerLane{0};
  std::uint32_t padding[2]{};
};
static_assert(sizeof(DevicePairParameters) == 128);

/// Mirrors the BuildParameters struct of the pair kernel source.
export struct DeviceBuildParameters
{
  float cell[9];
  float inverseCell[9];
  float listCutoffSquared{0.0f};
  std::uint32_t gridX{1}, gridY{1}, gridZ{1};
  std::uint32_t blocksPerCluster{0};
  std::uint32_t orthorhombic{0};
};
static_assert(sizeof(DeviceBuildParameters) == 96);

/// The per-molecule results of a step (bonded chain), summed over all molecules.
export struct DeviceBondedResults
{
  /// The self energy is evaluated in double on the host (it depends on the charges and alpha only); the
  /// device evaluates self + exclusion without their mutual cancellation and `exclusion` is the difference.
  double self{0.0}, exclusion{0.0}, bond{0.0}, bend{0.0}, torsion{0.0}, improperTorsion{0.0};
  double intraVDW{0.0}, intraCoulomb{0.0};
  double3x3 exclusionStrain{};
  double3x3 correction{};
};

/// The host arrays of a list build the device needs (the slot layout).
export struct DeviceLayoutUpload
{
  std::span<const float> buildPositions;         ///< x, y, z, molecule bits per slot (4 floats)
  std::span<const std::uint32_t> types;          ///< pseudo-atom type per slot
  std::span<const std::uint32_t> cellSlotStart;  ///< first slot per coarse cell (cells + 1)
  std::span<const std::uint32_t> cellOfCluster;  ///< coarse cell per i-cluster
};

export class DeviceBackend
{
 public:
  virtual ~DeviceBackend() = default;

  /// Name of the device in use (status report).
  virtual std::string deviceName() const = 0;

  // ---- tables and parameters of the pair kernel ----
  /// The Lennard-Jones table (4 floats per type pair); blocking.
  virtual void setLennardJones(std::span<const float> table) = 0;
  /// The pair-kernel parameters, written on the step queue (ordered before the following launches).
  virtual void writeParameters(const DevicePairParameters& parameters, bool blocking) = 0;
  /// The list-build parameters, written on the step queue.
  virtual void writeBuildParameters(const DeviceBuildParameters& parameters, bool blocking) = 0;

  // ---- device allocations (contents lost; see the conventions above) ----
  virtual void allocateSlots(std::size_t slots) = 0;         ///< positions, build positions, types, forces, bounds
  virtual void allocateRelative(std::size_t floats) = 0;     ///< relative positions (bonded terms)
  virtual void allocateClusters(std::size_t iClusters) = 0;  ///< per i-cluster: cell, row counts, lane counts, partials
  virtual void allocateCells(std::size_t entries) = 0;       ///< cell slot starts (entries of the cells + 1 array)
  virtual void allocateBlocks(std::size_t words) = 0;        ///< outer list rows (cluster + mask words)
  virtual void allocateLanes(std::size_t words) = 0;         ///< lane lists

  // ---- host-visible staging ----
  virtual float* mapPositions(std::size_t floats) = 0;
  virtual float* mapRelative(std::size_t floats) = 0;
  virtual const float* mapForces(std::size_t floats) = 0;
  virtual void unmapAll() = 0;

  // ---- list build ----
  /// Uploads the slot layout of a build (asynchronous; the host arrays stay valid until waitBuild).
  virtual void uploadLayout(const DeviceLayoutUpload& layout) = 0;
  /// Enqueues the bounds and list-build kernels and the read-back of the row counts into `outerCount`.
  virtual void enqueueBuild(std::size_t iClusters, std::size_t jClusters, std::span<std::uint32_t> outerCount) = 0;
  /// Completes the last enqueueBuild (the row counts are then readable).
  virtual void waitBuild() = 0;

  // ---- the step ----
  virtual void enqueueCompaction(std::size_t iClusters) = 0;
  virtual void enqueuePairs(std::size_t iClusters) = 0;
  /// Read-back of the lane counts of the last compaction (asynchronous, completed by waitStep).
  virtual void enqueueReadLaneCounts(std::span<std::uint32_t> laneCount) = 0;
  /// The mesh chain (spreading, transforms, influence function, interpolation into the forces).
  virtual void enqueueMesh() = 0;
  /// The bonded chain (term and atom kernels; adds to the forces).
  virtual void enqueueBonded() = 0;
  /// Read-back of the pair partials (11 floats per i-cluster) and of the mesh / bonded partials when enabled;
  /// marks the completion of the step for waitStep.
  virtual void enqueueReadPartials(std::span<float> pairPartials) = 0;
  virtual void waitStep() = 0;
  /// Starts the device on the work enqueued so far (non-blocking).
  virtual void flush() = 0;
  /// Completes all work enqueued so far (blocking).
  virtual void finish() = 0;

  // ---- the mesh (particle-mesh Ewald) ----
  virtual void enableMesh(int3 mesh, std::size_t order, double alpha, double conversionFactor, std::size_t slots) = 0;
  virtual void setMeshSlots(std::size_t slots) = 0;
  virtual void setMeshAlpha(double alpha) = 0;
  virtual void updateMeshBox(const SimulationBox& box) = 0;
  /// Times the stages of one mesh step synchronously (the mesh and the forces are left modified).
  virtual void profileMesh() = 0;
  /// After waitStep: the reciprocal energy, its strain derivative, the single-ion sum and tensor.
  virtual void collectMesh(double& energy, double3x3& strain, double& singleIonSum,
                           double3x3& singleIonStrain) const = 0;
  virtual std::string meshStatus() const = 0;

  // ---- the per-molecule terms ----
  virtual void enableBonded(const BondedTopology& topology, double alpha, double conversionFactor,
                            bool useCharge) = 0;
  virtual void setBondedAlpha(double alpha) = 0;
  /// The slot layout of a build: (molecule << 8) | index in the molecule per slot, or BondedTopology::noAtom
  /// for a dummy slot (asynchronous upload; the span stays valid until waitBuild).
  virtual void setBondedLayout(std::span<const std::uint32_t> slotMolecule) = 0;
  /// After waitStep: the sums over all molecules.
  virtual DeviceBondedResults collectBonded() const = 0;
  virtual std::string bondedStatus() const = 0;
};
