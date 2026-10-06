module;

export module spatial_decomposition_device_step;

import std;

import double3;
import double3x3;
import int3;
import simulationbox;
import system;
import spatial_decomposition_cell_list;
import spatial_decomposition_pair_kernel;
import spatial_decomposition_settings;
import spatial_decomposition_device_context;
import spatial_decomposition_device_mesh;
import spatial_decomposition_device_bonded;
import spatial_decomposition_device_bonded_topology;

/// Launch geometry fixed by the pair kernel source (CLUSTER_I, CLUSTER_J, GROUP_SIZE, the partial sums per group).
export namespace DeviceKernelLayout
{
constexpr std::size_t clusterI = 8;
constexpr std::size_t clusterJ = 4;
constexpr std::size_t pairGroupSize = 32;  ///< work-items (lanes) per i-cluster: one per (a, b) of 8 x 4
constexpr std::size_t pairPartials = 11;   ///< partial sums per i-cluster: 2 energies + 9 strain components
constexpr std::uint32_t noAtom = std::numeric_limits<std::uint32_t>::max();
/// Bit of an atom's exclusion start marking that every pair of its molecule is excluded (ALL_EXCLUDED_BIT).
constexpr std::uint32_t allExcludedBit = 0x80000000u;
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
  // the global Lennard-Jones switching (LennardJonesSwitching; mode 0: none)
  std::uint32_t switchMode{0};
  float switchDistanceSquared{0.0f};
  float switchDistance{0.0f};
  float switchInverseWidth{0.0f};
  float switchInverseCutOff3{0.0f};
  float switchA12{1.0f};
  float switchA6{1.0f};
  std::uint32_t padding[3]{};
};
static_assert(sizeof(DevicePairParameters) == 160);

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

/// Creates the context of a device kind; throws when none is available. Defined by the backends
/// (spatial_decomposition_opencl_context, spatial_decomposition_metal_context, spatial_decomposition_cuda_context)
/// and dispatched by DeviceStep.
export std::unique_ptr<DeviceContext> createDeviceContext(PairDevice device);
/// Whether a device of the kind is available, and its name (empty when none).
export bool deviceAvailable(PairDevice device);
export std::string deviceNameOf(PairDevice device);
/// Why no device of the kind is available (empty when one is, or when the backend gives no reason).
export std::string deviceUnavailableReason(PairDevice device);

/**
 * \brief The specialised pair kernel (Lennard-Jones + Ewald real space) on a device, and the façade of the other
 * device work of a step (the particle-mesh Ewald sum and the per-molecule terms). Backend-neutral: the device API
 * is DeviceContext (device/context.ixx), the kernels are the shared sources of
 * spatial_decomposition_device_kernels.
 *
 * The device evaluates the short-range pairs of the whole system in single precision: the engine packs the
 * positions after the position refresh, enqueues the step (upload, pruning when due, kernel, read-back) from
 * thread 0 and collects the forces, energies and strain derivative before the scatter phase. The device also
 * builds its own cluster pair lists; the host only lays out the slots after the binning of the CellList (no
 * per-domain Verlet lists and no ghost images are needed when the device does the pairs):
 *
 *  - Slots: the atoms in the order of a coarse grid of cells at least the list cutoff wide, Morton-sorted within
 *    a cell and padded per cell to a multiple of 8, so that the i-clusters (8 consecutive slots) and j-clusters
 *    (4 consecutive slots) are compact and never span two cells. Dummy slots carry charge 0 and no list bits.
 *  - Outer list (device, at the rebuild): per i-cluster the j-clusters of the 27 surrounding cells with at least
 *    one pair within cutoff + Verlet skin at the positions of the binning (bounding-box prefilter), with a 32-bit
 *    mask of the pairs in range between different molecules. Rows of a fixed capacity; the host grows the rows
 *    and rebuilds when one overflows.
 *  - Inner list (device, every few steps): per work-item (a, b) of the i-cluster the j-clusters of its listed
 *    pairs within cutoff + prune skin at the current positions, so that the pair kernel evaluates listed pairs
 *    only. Reused while no atom moved more than prune skin / 2 since the compaction (checked on the host when
 *    the positions are packed); without pruning, compacted once per build within cutoff + Verlet skin. The count
 *    is exact; the host grows the lane lists and compacts again when one overflows.
 *  - Both lists are full lists (every pair in the lists of both its clusters: no atomics and no j-forces, the
 *    host halves the energies and the strain).
 *  - Positions: the wrapped positions of the atoms (position minus the box translation recorded at the binning),
 *    the device applies the minimum image per pair, so the lists stay valid while no atom moved more than
 *    skin / 2, the cell list's rebuild criterion.
 *
 * With the mesh enabled the same step also spreads the charges of the slots, transforms, applies the influence
 * function and interpolates the reciprocal gradients into the device forces; with the bonded terms enabled the
 * self / exclusion corrections, the bonded terms and the virial correction follow, on positions relative to the
 * first atom of the molecule that the host packs alongside the wrapped positions. The device forces are then the
 * complete gradients.
 *
 * Forces (gradients) come back per slot in single precision; energies and the strain tensor are reduced per
 * work-group on the device and summed on the host in double. Movable value type.
 */
export class DeviceStep
{
 public:
  static constexpr std::size_t clusterI = DeviceKernelLayout::clusterI;
  static constexpr std::size_t clusterJ = DeviceKernelLayout::clusterJ;
  static constexpr std::uint32_t noAtom = DeviceKernelLayout::noAtom;

  /// Whether a device of the given kind is available (initializes its runtime on first use).
  static bool available(PairDevice device);
  /// Name of the device of the given kind (empty when none).
  static std::string deviceName(PairDevice device);
  /// Why no device of the given kind is available (empty when one is).
  static std::string unavailableReason(PairDevice device);
  /// Whether the device bonded kernels cover the intramolecular terms of the system (else `reason`).
  static bool supportsBonded(const System& system, std::string& reason);

  DeviceStep() = default;
  ~DeviceStep();
  DeviceStep(const DeviceStep&) = delete;
  DeviceStep& operator=(const DeviceStep&) = delete;
  DeviceStep(DeviceStep&&) noexcept = default;
  DeviceStep& operator=(DeviceStep&&) noexcept = default;

  /// Creates the context of the given kind (queue, programs); throws when no device is available or a build fails.
  void initialize(PairDevice device);
  bool initialized() const { return context != nullptr; }
  /// The device context (for the resident integrator, which shares the position and force buffers).
  DeviceContext& deviceContext() { return *context; }

  /// Fixes the same-molecule pairs the list build leaves out: the excluded pairs and the scaled (1-4) pairs of
  /// every molecule (IntraMolecularExclusions::isExcludedFromPairList), as a CSR table over the atoms in system
  /// order. Every other same-molecule pair is listed and evaluated by the pair kernel at full strength. Must be
  /// called before the first build and whenever the molecules of the system change.
  void setExclusions(const System& system);

  /// Fixes the pair-type tables, the cutoffs, the Ewald parameters (`useCharge`: Ewald real space on) and the
  /// pruning skin (0 or at least the Verlet skin: no pruning, the lane lists hold the whole outer list).
  void setParameters(std::span<const LennardJonesPair> lennardJones, std::size_t types, bool useCharge,
                     double cutOffVDW, double cutOffCharge, double conversionFactor, double alpha, double verletSkin,
                     double pruneSkinValue);
  /// Sets the global Lennard-Jones switching of the pairs (after setParameters; mode 0: none).
  void setSwitching(const LennardJonesSwitching& switching);

  /// Adds the particle-mesh Ewald sum to the device step (mesh of PPPM::chooseMesh, B-spline order 3..7).
  void enableMesh(int3 meshSize, std::size_t order, double alpha, double conversionFactor);
  /// Adds the per-molecule terms to the device step (supportsBonded must hold for the system).
  void enableBonded(const System& system, double alpha, double conversionFactor, bool useCharge);
  bool meshEnabled() const { return useMesh; }
  bool bondedEnabled() const { return useBonded; }
  /// Follows a change of the Ewald alpha in the mesh and bonded parameters (setParameters covers the pairs).
  void setEwaldAlpha(double alpha);

  /// List build after CellList::bin (one thread): the slot layout on the host, then the uploads and the build
  /// kernels enqueued on the device; `parts` is the number of threads that pack positions.
  void beginBuild(const CellList& cells, const SimulationBox& box, std::size_t parts);
  /// Waits for the build, grows the rows and rebuilds on an overflow, and binds the new lists (one thread).
  void finishBuild();

  /// Writes the current wrapped positions and charges of the given sorted atoms into the host staging array (any
  /// thread, for the atoms it owns; after the position refresh of the step), and records how far they moved
  /// since the last prune. With the bonded terms enabled it also stages the positions relative to the first
  /// atom of the molecule.
  void packPositions(std::size_t part, std::span<const std::uint32_t> sortedAtoms, const CellList& cells,
                     const SimulationBox& box);

  /// Enqueues the step: position upload, list compaction when due, pair kernel, mesh and bonded chains when
  /// enabled, read-back of forces and partial sums (one thread, non-blocking).
  void enqueue(const SimulationBox& box);

  /// Resident mode (the resident integrator, DeviceResident): the slot positions and the relative positions are
  /// written by device kernels, the host staging pointers are not mapped, the forces are read from the device
  /// buffer by the integrator and the reference positions of the list compaction are kept on the device.
  void setResident(bool resident);
  bool resident() const { return residentMode; }
  /// The device buffers of the step that the resident kernels read and write (valid after beginBuild; the
  /// buffers may be reallocated by a build, so query them after every build).
  struct ResidentTargets
  {
    DeviceBuffer positions{};         ///< float4 per slot: wrapped position, charge (written by the integrator)
    DeviceBuffer forces{};            ///< float4 per slot: the complete gradients after the step
    DeviceBuffer relative{};          ///< float4 per atom (system order): position relative to the first atom
    DeviceBuffer buildPositions{};    ///< float4 per slot: the wrapped positions at the list build
    DeviceBuffer compactReference{};  ///< float4 per slot: the positions at the last list compaction
    bool relativeEnabled{false};      ///< the bonded kernels are enabled (the relative positions are needed)
  };
  ResidentTargets residentTargets() const;
  /// Slot of every sorted atom of the current layout.
  std::span<const std::uint32_t> slotsOfSorted() const { return slotOfSorted; }
  /// Squared displacement since the last compaction above which the lane lists must be compacted again.
  double pruneDisplacementThreshold() const { return pruneDisplacementSquared; }
  /// Enqueues the step in resident mode: compaction when the lane lists are new or when `compactionDue` (an atom
  /// moved more than prune skin / 2 since the last compaction, measured on the device), then the chain.
  void enqueueResident(const SimulationBox& box, bool compactionDue);

  /// Results of a device step (the strain derivatives are the sums g (x) dr; the pair energies and strain are
  /// already halved for the full lists).
  struct Results
  {
    double energyVDW{0.0};
    double energyCharge{0.0};
    double3x3 pairStrain{};
    // mesh
    double reciprocalEnergy{0.0};
    double3x3 reciprocalStrain{};
    double singleIonSum{0.0};
    double3x3 singleIonStrain{};
    // bonded
    DeviceBondedResults bonded{};
  };
  /// Waits for the step and returns its energies and strain derivatives; the forces are then readable. When a
  /// lane list overflowed, the chain is enqueued again and `afterRetry` is called (the resident integrator
  /// re-enqueues the work that follows the chain).
  Results wait(const std::function<void()>& afterRetry = {});

  /// Gradient on a sorted atom from the last completed step.
  double3 force(std::uint32_t sorted) const
  {
    const float* f = mappedForces + 4 * static_cast<std::size_t>(slotOfSorted[sorted]);
    return double3(static_cast<double>(f[0]), static_cast<double>(f[1]), static_cast<double>(f[2]));
  }

  bool pruning() const { return pruneSkin > 0.0; }
  std::size_t numberOfSlots() const { return padded; }
  std::size_t numberOfIClusters() const { return padded / clusterI; }
  std::size_t numberOfBlocks() const { return outerBlocks; }
  /// Pairs in the lane lists of the last compaction (each pair once).
  std::size_t numberOfListedPairs() const { return listedPairs; }
  std::size_t numberOfCompactions() const { return compactions; }
  /// Estimated execution time of the pair kernel over all steps: one synchronous run is timed after every list
  /// build and counted for the steps up to the next build (the device APIs' profiling events are not relied on);
  /// likewise for the compaction kernel, the mesh chain and the bonded kernel.
  std::chrono::duration<double> deviceTime() const { return kernelTime; }
  std::chrono::duration<double> pruneTime() const { return prunedTime; }
  std::chrono::duration<double> meshTime() const { return meshedTime; }
  std::chrono::duration<double> bondedTime() const { return bondedRunTime; }
  /// Wall time of the list builds (host layout, uploads and the wait for the device build).
  std::chrono::duration<double> buildTime() const { return builtTime; }

  std::string status() const;

 private:
  std::unique_ptr<DeviceContext> context{};
  PairDevice deviceKind{PairDevice::CPU};

  // kernels and device buffers of the pair step
  DeviceKernel boundsKernel{}, buildKernel{}, compactKernel{}, pairKernel{};
  DeviceBufferOwner parameterBuffer{}, buildParameterBuffer{}, lennardJonesBuffer{};
  DeviceBufferOwner positionBuffer{}, buildPositionBuffer{}, typeBuffer{}, forceBuffer{}, relativeBuffer{};
  DeviceBufferOwner clusterMinBuffer{}, clusterMaxBuffer{}, compactReferenceBuffer{};
  DeviceBufferOwner cellOfClusterBuffer{}, outerCountBuffer{}, laneCountBuffer{}, partialBuffer{};
  DeviceBufferOwner cellSlotStartBuffer{}, outerClusterBuffer{}, outerMaskBuffer{}, pairListBuffer{};
  DeviceBufferOwner slotAtomBuffer{}, exclusionStartBuffer{}, exclusionPartnerBuffer{};
  std::size_t lennardJonesCapacity{0};
  std::size_t exclusionStartCapacity{0}, exclusionPartnerCapacity{0};
  DeviceEvent buildEvent{}, stepEvent{};
  DeviceMesh mesh{};
  DeviceBonded bonded{};

  // device capacities
  std::size_t slotCapacity{0}, iClusterCapacity{0}, cellCapacity{0}, relativeCapacity{0};
  std::size_t blocksPerCluster{0};  ///< Row capacity of the outer list on the device.
  std::size_t blockWords{0};        ///< Allocated words per outer list array (iClusterCapacity rows).
  std::size_t pairsPerLane{0};      ///< Capacity of a lane list.
  std::size_t laneWords{0};         ///< Allocated words of the lane lists (iClusterCapacity * 32 lanes).

  bool useMesh{false};
  bool meshProfiled{false};  ///< the stage profile of the mesh ran (once, at the first timing sample)
  bool useBonded{false};
  BondedTopology bondedTopology{};
  std::vector<std::uint32_t> referenceOfSorted{};

  DevicePairParameters parameters{};
  DeviceBuildParameters buildParameters{};
  std::vector<float> lennardJonesTable{};  ///< 4 epsilon, sigma^6, shift per type pair.
  bool parametersChanged{true};
  double pruneSkin{0.0};
  double innerCutoff{0.0};  ///< cutoff + prune skin (or + Verlet skin without pruning) of the lane lists
  double pruneDisplacementSquared{0.0};

  // slot layout of the current list
  std::size_t padded{0};
  int3 grid{1, 1, 1};
  std::vector<std::uint32_t> cellSlotStart{};
  std::vector<std::uint32_t> cellOfCluster{};
  std::vector<std::uint32_t> slotOfSorted{};
  std::vector<std::uint32_t> sortedOfSlot{};
  std::vector<std::uint32_t> typeOfSlot{};
  std::vector<float> buildPosition{};                          ///< x, y, z, molecule bits per slot at the binning.
  std::vector<std::uint32_t> slotAtom{};                       ///< System atom index per slot (noAtom for dummies).
  /// Per atom (system order): CSR start of its excluded + scaled same-molecule partners (bit 31: every pair of
  /// its molecule is excluded, rigid molecule), and the partner atoms (system order).
  std::vector<std::uint32_t> atomExclusionStart{}, atomExclusionPartner{};
  std::vector<std::uint32_t> sortKey{}, order{}, cellCount{};  ///< Layout scratch.
  std::vector<std::uint32_t> outerCount{}, laneCount{};
  std::size_t outerBlocks{0}, listedPairs{0}, maximumRow{0};

  // host staging: positions, relative positions and forces live in host-visible (Shared) device buffers
  float* mappedPositions{nullptr};
  float* mappedRelative{nullptr};
  const float* mappedForces{nullptr};
  std::size_t mappedPositionFloats{0}, mappedRelativeFloats{0}, mappedForceFloats{0};
  std::vector<float> listPositions{}, hostPartials{}, zeroPositions{};
  std::vector<double> partDisplacement{};  ///< Per packing thread: largest squared move since the last compaction.
  bool residentMode{false};                ///< see setResident
  bool laneListsValid{false};              ///< the lane lists (and listPositions) refer to the current slot layout.
  bool compactedThisStep{false};           ///< the enqueued step compacts the lists (lane counts are read back).
  bool sampleRequested{false};             ///< Time the kernels synchronously at the next step (after a build).
  std::chrono::steady_clock::time_point buildStart{};
  std::chrono::duration<double> kernelTime{}, prunedTime{}, builtTime{}, meshedTime{}, bondedRunTime{};
  std::chrono::duration<double> sampledKernelTime{}, sampledPruneTime{}, sampledMeshTime{}, sampledBondedTime{};
  std::size_t steps{0}, compactions{0}, builds{0};

  void ensureBuffers();
  void ensureBlockBuffers();
  void ensureLaneBuffers();
  void writeParameters(bool blocking);
  void writeBuildParameters(bool blocking);
  void enqueueBuild();
  void enqueueCompaction();
  void enqueuePairs();
  void enqueueMesh();
  void enqueueBonded();
  void enqueueChain();
  void enqueueStep(const SimulationBox& box, bool compact);
  void updateParameters(const SimulationBox& box);
  /// Maps the position (and relative-position) buffers for the host packing; `discard` when the device has not
  /// written them since the host last unmapped them (the host-integrator steps), not after the resident
  /// integrator wrote them on the device.
  void mapInputs(bool discard);
  void mapForces();
  void unmapHost();
};
