module;

#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

export module spatial_decomposition_opencl_backend;

import std;

import double3x3;
import int3;
import simulationbox;
import spatial_decomposition_device_backend;
import spatial_decomposition_device_bonded_topology;
import spatial_decomposition_opencl_handles;
import spatial_decomposition_opencl_mesh;
import spatial_decomposition_opencl_bonded;

/// Whether an OpenCL device is available (initializes the shared OpenCL context on first use).
export bool openclAvailable();
/// Name of the OpenCL device in use (empty when none).
export std::string openclDeviceName();
/// Creates the OpenCL backend of the device step: the command queue, the pair-kernel program and the parameter
/// buffers; throws when no device is available or the build fails.
export std::unique_ptr<DeviceBackend> createOpenCLBackend();

/**
 * \brief The OpenCL (1.2) implementation of DeviceBackend.
 *
 * One in-order command queue carries the step; the shared kernel sources (spatial_decomposition_device_kernels)
 * are compiled with the OpenCL dialect header. The host-visible staging of the positions, the relative positions
 * and the forces are device buffers allocated with CL_MEM_ALLOC_HOST_PTR and mapped (zero-copy on a unified-memory
 * device). The mesh (OpenCLMesh) and the per-molecule terms (OpenCLBonded) share the queue and the position /
 * force buffers. The arguments of the compaction and pair kernels are bound at the next launch after a buffer of
 * theirs was (re)allocated.
 */
export class OpenCLBackend final : public DeviceBackend
{
 public:
  OpenCLBackend();
  ~OpenCLBackend() override;
  OpenCLBackend(const OpenCLBackend&) = delete;
  OpenCLBackend& operator=(const OpenCLBackend&) = delete;

  std::string deviceName() const override;

  void setLennardJones(std::span<const float> table) override;
  void writeParameters(const DevicePairParameters& parameters, bool blocking) override;
  void writeBuildParameters(const DeviceBuildParameters& parameters, bool blocking) override;

  void allocateSlots(std::size_t slots) override;
  void allocateRelative(std::size_t floats) override;
  void allocateClusters(std::size_t iClusters) override;
  void allocateCells(std::size_t entries) override;
  void allocateBlocks(std::size_t words) override;
  void allocateLanes(std::size_t words) override;

  float* mapPositions(std::size_t floats) override;
  float* mapRelative(std::size_t floats) override;
  const float* mapForces(std::size_t floats) override;
  void unmapAll() override;

  void uploadLayout(const DeviceLayoutUpload& layout) override;
  void enqueueBuild(std::size_t iClusters, std::size_t jClusters, std::span<std::uint32_t> outerCount) override;
  void waitBuild() override;

  void enqueueCompaction(std::size_t iClusters) override;
  void enqueuePairs(std::size_t iClusters) override;
  void enqueueReadLaneCounts(std::span<std::uint32_t> laneCount) override;
  void enqueueMesh() override;
  void enqueueBonded() override;
  void enqueueReadPartials(std::span<float> pairPartials) override;
  void waitStep() override;
  void flush() override;
  void finish() override;

  void enableMesh(int3 mesh, std::size_t order, double alpha, double conversionFactor, std::size_t slots) override;
  void setMeshSlots(std::size_t slots) override;
  void setMeshAlpha(double alpha) override;
  void updateMeshBox(const SimulationBox& box) override;
  void profileMesh() override;
  void collectMesh(double& energy, double3x3& strain, double& singleIonSum, double3x3& singleIonStrain) const override;
  std::string meshStatus() const override;

  void enableBonded(const BondedTopology& topology, double alpha, double conversionFactor, bool useCharge) override;
  void setBondedAlpha(double alpha) override;
  void setBondedLayout(std::span<const std::uint32_t> slotMolecule) override;
  DeviceBondedResults collectBonded() const override;
  std::string bondedStatus() const override;

 private:
  cl_context context{nullptr};
  cl_device_id device{nullptr};
  OpenCLDevice::QueueHandle queue{};
  OpenCLDevice::ProgramHandle program{};
  OpenCLDevice::KernelHandle boundsKernel{}, buildKernel{}, compactKernel{}, pairKernel{};
  OpenCLDevice::MemHandle positionBuffer{}, relativeBuffer{}, buildPositionBuffer{}, typeBuffer{}, forceBuffer{};
  OpenCLDevice::MemHandle cellSlotStartBuffer{}, cellOfClusterBuffer{}, clusterMinBuffer{}, clusterMaxBuffer{};
  OpenCLDevice::MemHandle outerClusterBuffer{}, outerMaskBuffer{}, outerCountBuffer{};
  OpenCLDevice::MemHandle pairListBuffer{}, laneCountBuffer{};
  OpenCLDevice::MemHandle lennardJonesBuffer{}, parameterBuffer{}, buildParameterBuffer{}, partialBuffer{};
  std::size_t lennardJonesCapacity{0};
  OpenCLDevice::EventHandle stepEvent{}, buildEvent{};
  bool listArgumentsDirty{true};

  OpenCLMesh mesh{};
  OpenCLBonded bonded{};
  bool meshEnabled{false};
  bool bondedEnabled{false};

  float* mappedPositions{nullptr};
  float* mappedRelative{nullptr};
  float* mappedForces{nullptr};
  std::size_t mappedPositionFloats{0}, mappedRelativeFloats{0}, mappedForceFloats{0};

  void bindListArguments();
  void unmap(cl_mem buffer, float*& pointer, std::string_view what);
  float* map(cl_mem buffer, float*& pointer, std::size_t& mappedFloats, std::size_t floats, cl_map_flags flags,
             std::string_view what);
};
