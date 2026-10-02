module;

#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

module spatial_decomposition_opencl_backend;

import std;

import double3x3;
import int3;
import simulationbox;
import opencl;
import spatial_decomposition_device_kernels;
import spatial_decomposition_device_backend;
import spatial_decomposition_device_bonded_topology;
import spatial_decomposition_opencl_handles;
import spatial_decomposition_opencl_mesh;
import spatial_decomposition_opencl_bonded;

using OpenCLDevice::check;
using OpenCLDevice::roundUp;
using DeviceKernelLayout::clusterJ;
using DeviceKernelLayout::pairGroupSize;
using DeviceKernelLayout::pairPartials;

namespace
{
constexpr std::size_t boundsGroupSize = 64;
constexpr cl_mem_flags hostReadFlags = CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_ONLY;
constexpr cl_mem_flags hostWriteFlags = CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE;
}  // namespace

bool openclAvailable()
{
  OpenCL::initialize();
  return OpenCL::clContext.has_value() && OpenCL::clDeviceId.has_value();
}

std::string openclDeviceName()
{
  if (!openclAvailable()) return {};
  char name[256] = {};
  std::size_t length = 0;
  if (clGetDeviceInfo(OpenCL::clDeviceId.value(), CL_DEVICE_NAME, sizeof(name) - 1, name, &length) != CL_SUCCESS)
  {
    return "unknown device";
  }
  return std::string(name);
}

std::unique_ptr<DeviceBackend> createOpenCLBackend() { return std::make_unique<OpenCLBackend>(); }

OpenCLBackend::OpenCLBackend()
{
  if (!openclAvailable())
  {
    throw std::runtime_error("[OpenCL pair kernel]: no OpenCL device is available\n");
  }
  cl_int error = CL_SUCCESS;
  context = OpenCL::clContext.value();
  device = OpenCL::clDeviceId.value();

  queue.reset(clCreateCommandQueue(context, device, 0, &error));
  check(error, "clCreateCommandQueue");

  program = OpenCLDevice::buildDeviceProgram(context, device, deviceKernelPairSource,
                                             "-cl-mad-enable -cl-no-signed-zeros -cl-fast-relaxed-math",
                                             "OpenCL pair kernel");
  boundsKernel = OpenCLDevice::createKernel(program.get(), "clusterBounds");
  buildKernel = OpenCLDevice::createKernel(program.get(), "buildList");
  compactKernel = OpenCLDevice::createKernel(program.get(), "compactList");
  pairKernel = OpenCLDevice::createKernel(program.get(), "clusterPairs");

  parameterBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_ONLY, sizeof(DevicePairParameters)));
  buildParameterBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_ONLY, sizeof(DeviceBuildParameters)));
}

OpenCLBackend::~OpenCLBackend()
{
  if (queue.get() != nullptr) unmapAll();
}

std::string OpenCLBackend::deviceName() const { return openclDeviceName(); }

// ---- tables and parameters ----

void OpenCLBackend::setLennardJones(std::span<const float> table)
{
  if (table.size() > lennardJonesCapacity)
  {
    lennardJonesCapacity = table.size();
    lennardJonesBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_ONLY, lennardJonesCapacity * sizeof(float)));
    listArgumentsDirty = true;
  }
  if (!table.empty())
  {
    OpenCL::writeBuffer(lennardJonesBuffer.get(), table.size() * sizeof(float), table.data());
  }
}

void OpenCLBackend::writeParameters(const DevicePairParameters& parameters, bool blocking)
{
  check(clEnqueueWriteBuffer(queue.get(), parameterBuffer.get(), blocking ? CL_TRUE : CL_FALSE, 0,
                             sizeof(DevicePairParameters), &parameters, 0, nullptr, nullptr),
        "clEnqueueWriteBuffer (parameters)");
}

void OpenCLBackend::writeBuildParameters(const DeviceBuildParameters& parameters, bool blocking)
{
  check(clEnqueueWriteBuffer(queue.get(), buildParameterBuffer.get(), blocking ? CL_TRUE : CL_FALSE, 0,
                             sizeof(DeviceBuildParameters), &parameters, 0, nullptr, nullptr),
        "clEnqueueWriteBuffer (build parameters)");
}

// ---- allocations ----

void OpenCLBackend::allocateSlots(std::size_t slots)
{
  unmapAll();
  positionBuffer.reset(OpenCL::createBuffer(hostReadFlags, slots * 4 * sizeof(float)));
  buildPositionBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_ONLY, slots * 4 * sizeof(float)));
  typeBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_ONLY, slots * sizeof(std::uint32_t)));
  forceBuffer.reset(OpenCL::createBuffer(hostWriteFlags, slots * 4 * sizeof(float)));
  clusterMinBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_WRITE, (slots / clusterJ) * 4 * sizeof(float)));
  clusterMaxBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_WRITE, (slots / clusterJ) * 4 * sizeof(float)));
  listArgumentsDirty = true;
}

void OpenCLBackend::allocateRelative(std::size_t floats)
{
  unmap(relativeBuffer.get(), mappedRelative, "clEnqueueUnmapMemObject (relative positions)");
  relativeBuffer.reset(OpenCL::createBuffer(hostReadFlags, floats * sizeof(float)));
}

void OpenCLBackend::allocateClusters(std::size_t iClusters)
{
  cellOfClusterBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_ONLY, iClusters * sizeof(std::uint32_t)));
  outerCountBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_WRITE, iClusters * sizeof(std::uint32_t)));
  laneCountBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_WRITE, iClusters * pairGroupSize * sizeof(std::uint32_t)));
  partialBuffer.reset(OpenCL::createBuffer(CL_MEM_WRITE_ONLY, iClusters * pairPartials * sizeof(float)));
  listArgumentsDirty = true;
}

void OpenCLBackend::allocateCells(std::size_t entries)
{
  cellSlotStartBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_ONLY, entries * sizeof(std::uint32_t)));
}

void OpenCLBackend::allocateBlocks(std::size_t words)
{
  const std::size_t bytes = words * sizeof(std::uint32_t);
  outerClusterBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_WRITE, bytes));
  outerMaskBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_WRITE, bytes));
  listArgumentsDirty = true;
}

void OpenCLBackend::allocateLanes(std::size_t words)
{
  pairListBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_WRITE, words * sizeof(std::uint32_t)));
  listArgumentsDirty = true;
}

// ---- host-visible staging ----

void OpenCLBackend::unmap(cl_mem buffer, float*& pointer, std::string_view what)
{
  if (pointer == nullptr) return;
  if (buffer != nullptr)
  {
    check(clEnqueueUnmapMemObject(queue.get(), buffer, pointer, 0, nullptr, nullptr), what);
  }
  pointer = nullptr;
}

float* OpenCLBackend::map(cl_mem buffer, float*& pointer, std::size_t& mappedFloats, std::size_t floats,
                          cl_map_flags flags, std::string_view what)
{
  if (buffer == nullptr || floats == 0) return pointer;
  if (pointer != nullptr && mappedFloats >= floats) return pointer;
  // mapped for a smaller layout: map again for the current one
  unmap(buffer, pointer, "clEnqueueUnmapMemObject");
  cl_int error = CL_SUCCESS;
  mappedFloats = floats;
  pointer = static_cast<float*>(
      clEnqueueMapBuffer(queue.get(), buffer, CL_TRUE, flags, 0, floats * sizeof(float), 0, nullptr, nullptr, &error));
  check(error, what);
  return pointer;
}

float* OpenCLBackend::mapPositions(std::size_t floats)
{
  return map(positionBuffer.get(), mappedPositions, mappedPositionFloats, floats, CL_MAP_WRITE,
             "clEnqueueMapBuffer (positions)");
}

float* OpenCLBackend::mapRelative(std::size_t floats)
{
  return map(relativeBuffer.get(), mappedRelative, mappedRelativeFloats, floats, CL_MAP_WRITE,
             "clEnqueueMapBuffer (relative positions)");
}

const float* OpenCLBackend::mapForces(std::size_t floats)
{
  return map(forceBuffer.get(), mappedForces, mappedForceFloats, floats, CL_MAP_READ, "clEnqueueMapBuffer (forces)");
}

void OpenCLBackend::unmapAll()
{
  unmap(positionBuffer.get(), mappedPositions, "clEnqueueUnmapMemObject (positions)");
  unmap(relativeBuffer.get(), mappedRelative, "clEnqueueUnmapMemObject (relative positions)");
  unmap(forceBuffer.get(), mappedForces, "clEnqueueUnmapMemObject (forces)");
}

// ---- list build ----

void OpenCLBackend::uploadLayout(const DeviceLayoutUpload& layout)
{
  auto upload = [&](cl_mem buffer, std::size_t bytes, const void* host, std::string_view what)
  { check(clEnqueueWriteBuffer(queue.get(), buffer, CL_FALSE, 0, bytes, host, 0, nullptr, nullptr), what); };
  upload(buildPositionBuffer.get(), layout.buildPositions.size_bytes(), layout.buildPositions.data(),
         "clEnqueueWriteBuffer (build positions)");
  upload(typeBuffer.get(), layout.types.size_bytes(), layout.types.data(), "clEnqueueWriteBuffer (types)");
  upload(cellSlotStartBuffer.get(), layout.cellSlotStart.size_bytes(), layout.cellSlotStart.data(),
         "clEnqueueWriteBuffer (cell slots)");
  upload(cellOfClusterBuffer.get(), layout.cellOfCluster.size_bytes(), layout.cellOfCluster.data(),
         "clEnqueueWriteBuffer (cell of cluster)");
}

void OpenCLBackend::enqueueBuild(std::size_t iClusters, std::size_t jClusters, std::span<std::uint32_t> outerCount)
{
  const cl_uint numberOfJClusters = static_cast<cl_uint>(jClusters);
  const cl_mem boundsBuffers[2] = {clusterMinBuffer.get(), clusterMaxBuffer.get()};
  cl_mem buildPositions = buildPositionBuffer.get();
  check(clSetKernelArg(boundsKernel.get(), 0, sizeof(cl_mem), &buildPositions), "clSetKernelArg (clusterBounds)");
  check(clSetKernelArg(boundsKernel.get(), 1, sizeof(cl_uint), &numberOfJClusters), "clSetKernelArg (clusterBounds)");
  check(clSetKernelArg(boundsKernel.get(), 2, sizeof(cl_mem), &boundsBuffers[0]), "clSetKernelArg (clusterBounds)");
  check(clSetKernelArg(boundsKernel.get(), 3, sizeof(cl_mem), &boundsBuffers[1]), "clSetKernelArg (clusterBounds)");
  const std::size_t boundsGlobal = roundUp(std::max<std::size_t>(jClusters, 1), boundsGroupSize);
  const std::size_t boundsLocal = boundsGroupSize;
  check(clEnqueueNDRangeKernel(queue.get(), boundsKernel.get(), 1, nullptr, &boundsGlobal, &boundsLocal, 0, nullptr,
                               nullptr),
        "clEnqueueNDRangeKernel (clusterBounds)");

  const cl_mem buffers[9] = {buildPositionBuffer.get(), cellSlotStartBuffer.get(), cellOfClusterBuffer.get(),
                             clusterMinBuffer.get(),    clusterMaxBuffer.get(),    buildParameterBuffer.get(),
                             outerClusterBuffer.get(),  outerMaskBuffer.get(),     outerCountBuffer.get()};
  OpenCLDevice::setBufferArguments(buildKernel.get(), buffers, "clSetKernelArg (buildList)");
  const std::size_t global = std::max<std::size_t>(iClusters, 1) * pairGroupSize;
  const std::size_t local = pairGroupSize;
  check(clEnqueueNDRangeKernel(queue.get(), buildKernel.get(), 1, nullptr, &global, &local, 0, nullptr, nullptr),
        "clEnqueueNDRangeKernel (buildList)");
  cl_event event = nullptr;
  check(clEnqueueReadBuffer(queue.get(), outerCountBuffer.get(), CL_FALSE, 0, outerCount.size_bytes(),
                            outerCount.data(), 0, nullptr, &event),
        "clEnqueueReadBuffer (row counts)");
  buildEvent.reset(event);
}

void OpenCLBackend::waitBuild()
{
  cl_event event = buildEvent.get();
  check(clWaitForEvents(1, &event), "clWaitForEvents (build)");
  buildEvent.reset();
}

// ---- the step ----

void OpenCLBackend::bindListArguments()
{
  if (!listArgumentsDirty) return;
  listArgumentsDirty = false;
  const cl_mem compactBuffers[7] = {positionBuffer.get(),   outerClusterBuffer.get(), outerMaskBuffer.get(),
                                    outerCountBuffer.get(), parameterBuffer.get(),    pairListBuffer.get(),
                                    laneCountBuffer.get()};
  OpenCLDevice::setBufferArguments(compactKernel.get(), compactBuffers, "clSetKernelArg (compactList)");
  const cl_mem pairBuffers[8] = {positionBuffer.get(),  typeBuffer.get(),         pairListBuffer.get(),
                                 laneCountBuffer.get(), lennardJonesBuffer.get(), parameterBuffer.get(),
                                 forceBuffer.get(),     partialBuffer.get()};
  OpenCLDevice::setBufferArguments(pairKernel.get(), pairBuffers, "clSetKernelArg (clusterPairs)");
}

void OpenCLBackend::enqueueCompaction(std::size_t iClusters)
{
  bindListArguments();
  const std::size_t global = std::max<std::size_t>(iClusters, 1) * pairGroupSize;
  const std::size_t local = pairGroupSize;
  check(clEnqueueNDRangeKernel(queue.get(), compactKernel.get(), 1, nullptr, &global, &local, 0, nullptr, nullptr),
        "clEnqueueNDRangeKernel (compactList)");
}

void OpenCLBackend::enqueuePairs(std::size_t iClusters)
{
  bindListArguments();
  const std::size_t global = std::max<std::size_t>(iClusters, 1) * pairGroupSize;
  const std::size_t local = pairGroupSize;
  check(clEnqueueNDRangeKernel(queue.get(), pairKernel.get(), 1, nullptr, &global, &local, 0, nullptr, nullptr),
        "clEnqueueNDRangeKernel (clusterPairs)");
}

void OpenCLBackend::enqueueReadLaneCounts(std::span<std::uint32_t> laneCount)
{
  check(clEnqueueReadBuffer(queue.get(), laneCountBuffer.get(), CL_FALSE, 0, laneCount.size_bytes(), laneCount.data(),
                            0, nullptr, nullptr),
        "clEnqueueReadBuffer (lane counts)");
}

void OpenCLBackend::enqueueMesh() { mesh.enqueue(queue.get(), positionBuffer.get(), forceBuffer.get(), true); }

void OpenCLBackend::enqueueBonded() { bonded.enqueue(queue.get(), relativeBuffer.get(), forceBuffer.get()); }

void OpenCLBackend::enqueueReadPartials(std::span<float> pairPartialSums)
{
  cl_event event = nullptr;
  check(clEnqueueReadBuffer(queue.get(), partialBuffer.get(), CL_FALSE, 0, pairPartialSums.size_bytes(),
                            pairPartialSums.data(), 0, nullptr, &event),
        "clEnqueueReadBuffer (partials)");
  // the in-order queue: the event of the last read covers the earlier ones
  if (meshEnabled)
  {
    clReleaseEvent(event);
    event = nullptr;
    mesh.enqueueRead(queue.get(), &event);
  }
  if (bondedEnabled)
  {
    clReleaseEvent(event);
    event = nullptr;
    bonded.enqueueRead(queue.get(), &event);
  }
  stepEvent.reset(event);
}

void OpenCLBackend::waitStep()
{
  cl_event event = stepEvent.get();
  check(clWaitForEvents(1, &event), "clWaitForEvents");
  stepEvent.reset();
}

void OpenCLBackend::flush() { check(clFlush(queue.get()), "clFlush"); }

void OpenCLBackend::finish() { check(clFinish(queue.get()), "clFinish"); }

// ---- the mesh ----

void OpenCLBackend::enableMesh(int3 meshSize, std::size_t order, double alpha, double conversionFactor,
                               std::size_t slots)
{
  if (!mesh.initialized() || mesh.interpolationOrder() != std::clamp<std::size_t>(order, 3, 7))
  {
    mesh.initialize(context, device, order);
  }
  mesh.setup(meshSize, alpha, conversionFactor);
  mesh.setSlots(slots);
  meshEnabled = true;
}

void OpenCLBackend::setMeshSlots(std::size_t slots) { mesh.setSlots(slots); }

void OpenCLBackend::setMeshAlpha(double alpha) { mesh.setAlpha(alpha); }

void OpenCLBackend::updateMeshBox(const SimulationBox& box) { mesh.updateBox(box); }

void OpenCLBackend::profileMesh() { mesh.profile(queue.get(), positionBuffer.get(), forceBuffer.get()); }

void OpenCLBackend::collectMesh(double& energy, double3x3& strain, double& singleIonSum,
                                double3x3& singleIonStrain) const
{
  mesh.collect(energy, strain, singleIonSum, singleIonStrain);
}

std::string OpenCLBackend::meshStatus() const { return mesh.status(); }

// ---- the per-molecule terms ----

void OpenCLBackend::enableBonded(const BondedTopology& topology, double alpha, double conversionFactor,
                                 bool useCharge)
{
  if (!bonded.initialized()) bonded.initialize(context, device);
  bonded.setTopology(topology);
  bonded.setParameters(alpha, conversionFactor, useCharge);
  bondedEnabled = true;
}

void OpenCLBackend::setBondedAlpha(double alpha) { bonded.setAlpha(alpha); }

void OpenCLBackend::setBondedLayout(std::span<const std::uint32_t> slotMolecule)
{
  bonded.setLayout(queue.get(), slotMolecule);
}

DeviceBondedResults OpenCLBackend::collectBonded() const { return bonded.collect(); }

std::string OpenCLBackend::bondedStatus() const { return bonded.status(); }
