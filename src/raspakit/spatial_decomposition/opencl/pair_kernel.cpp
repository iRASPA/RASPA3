module;

#define CL_TARGET_OPENCL_VERSION 120
#ifdef __APPLE__
#include <OpenCL/cl.h>
#elif _WIN32
#include <CL/cl.h>
#else
#include <CL/opencl.h>
#endif

module spatial_decomposition_opencl_pair_kernel;

import std;

import double3;
import double3x3;
import int3;
import simulationbox;
import system;
import opencl;
import spatial_decomposition_cell_list;
import spatial_decomposition_pair_kernel;
import spatial_decomposition_opencl_handles;
import spatial_decomposition_opencl_mesh;
import spatial_decomposition_opencl_bonded;

using OpenCLDevice::check;
using OpenCLDevice::roundUp;

namespace
{
constexpr std::size_t partialsPerCluster = 11;
constexpr std::size_t groupSize = 32;
constexpr std::size_t boundsGroupSize = 64;
constexpr float farAway = 1e30f;         // build position of a dummy slot: never within the list cutoff of a real atom
constexpr std::uint32_t mortonBits = 3;  // 8 x 8 x 8 sub-cells order the atoms within a cell

float bitsAsFloat(std::uint32_t bits) { return std::bit_cast<float>(bits); }

std::uint32_t mortonKey(std::uint32_t x, std::uint32_t y, std::uint32_t z)
{
  std::uint32_t key = 0;
  for (std::uint32_t b = 0; b < mortonBits; ++b)
  {
    key |= ((x >> b) & 1u) << (3 * b);
    key |= ((y >> b) & 1u) << (3 * b + 1);
    key |= ((z >> b) & 1u) << (3 * b + 2);
  }
  return key;
}

void setArguments(cl_kernel kernel, std::span<const cl_mem> buffers, std::string_view what)
{
  OpenCLDevice::setBufferArguments(kernel, buffers, what);
}
}  // namespace

bool OpenCLPairKernel::available()
{
  OpenCL::initialize();
  return OpenCL::clContext.has_value() && OpenCL::clDeviceId.has_value();
}

std::string OpenCLPairKernel::deviceName()
{
  if (!available()) return {};
  char name[256] = {};
  std::size_t length = 0;
  if (clGetDeviceInfo(OpenCL::clDeviceId.value(), CL_DEVICE_NAME, sizeof(name) - 1, name, &length) != CL_SUCCESS)
  {
    return "unknown device";
  }
  return std::string(name);
}

void OpenCLPairKernel::initialize()
{
  if (!available())
  {
    throw std::runtime_error("[OpenCL pair kernel]: no OpenCL device is available\n");
  }
  cl_int error = CL_SUCCESS;
  cl_context context = OpenCL::clContext.value();
  cl_device_id device = OpenCL::clDeviceId.value();

  queue.reset(clCreateCommandQueue(context, device, 0, &error));
  check(error, "clCreateCommandQueue");

  program =
      OpenCLDevice::buildProgram(context, device, openclPairKernelSource,
                                 "-cl-mad-enable -cl-no-signed-zeros -cl-fast-relaxed-math", "OpenCL pair kernel");
  boundsKernel = OpenCLDevice::createKernel(program.get(), "clusterBounds");
  buildKernel = OpenCLDevice::createKernel(program.get(), "buildList");
  compactKernel = OpenCLDevice::createKernel(program.get(), "compactList");
  pairKernel = OpenCLDevice::createKernel(program.get(), "clusterPairs");

  parameterBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_ONLY, sizeof(Parameters)));
  buildParameterBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_ONLY, sizeof(BuildParameters)));
  useMesh = false;
  useBonded = false;
  kernelTime = {};
  prunedTime = {};
  builtTime = {};
  meshedTime = {};
  bondedRunTime = {};
  sampledKernelTime = {};
  sampledPruneTime = {};
  sampledMeshTime = {};
  sampledBondedTime = {};
  steps = 0;
  compactions = 0;
  builds = 0;
}

void OpenCLPairKernel::enableMesh(int3 meshSize, std::size_t interpolationOrder, double alpha, double conversionFactor)
{
  if (!mesh.initialized() || mesh.interpolationOrder() != std::clamp<std::size_t>(interpolationOrder, 3, 7))
  {
    mesh.initialize(OpenCL::clContext.value(), OpenCL::clDeviceId.value(), interpolationOrder);
  }
  mesh.setup(meshSize, alpha, conversionFactor);
  mesh.setSlots(padded);
  useMesh = true;
}

void OpenCLPairKernel::enableBonded(const System& system, double alpha, double conversionFactor, bool useCharge)
{
  if (!bonded.initialized()) bonded.initialize(OpenCL::clContext.value(), OpenCL::clDeviceId.value());
  bonded.setTopology(system);
  bonded.setParameters(alpha, conversionFactor, useCharge);
  useBonded = true;
}

void OpenCLPairKernel::setEwaldAlpha(double alpha)
{
  if (useMesh) mesh.setAlpha(alpha);
  if (useBonded) bonded.setAlpha(alpha);
}

void OpenCLPairKernel::setParameters(std::span<const LennardJonesPair> lennardJones, std::size_t types, bool useCharge,
                                     double cutOffVDW, double cutOffCharge, double conversionFactor, double alpha,
                                     double verletSkin, double pruneSkinValue)
{
  lennardJonesTable.resize(3 * types * types);
  for (std::size_t k = 0; k < types * types; ++k)
  {
    const LennardJonesPair& pair = lennardJones[k];
    const double sigma2 = 1.0 / pair.inverseSigma2;
    lennardJonesTable[3 * k] = static_cast<float>(pair.epsilon4);
    lennardJonesTable[3 * k + 1] = static_cast<float>(sigma2 * sigma2 * sigma2);
    lennardJonesTable[3 * k + 2] = static_cast<float>(pair.shift);
  }
  parameters.cutOffVDWSquared = static_cast<float>(cutOffVDW * cutOffVDW);
  parameters.cutOffChargeSquared = static_cast<float>(cutOffCharge * cutOffCharge);
  parameters.alpha = static_cast<float>(alpha);
  parameters.alphaSquared = static_cast<float>(alpha * alpha);
  parameters.alphaOverSqrtPi = static_cast<float>(alpha / std::sqrt(std::numbers::pi));
  parameters.coulombFactor = static_cast<float>(conversionFactor);
  parameters.useCharge = useCharge ? 1u : 0u;
  parameters.numberOfTypes = static_cast<std::uint32_t>(types);

  const double cutoff = useCharge ? std::max(cutOffVDW, cutOffCharge) : cutOffVDW;
  pruneSkin = (pruneSkinValue > 0.0 && pruneSkinValue < verletSkin) ? pruneSkinValue : 0.0;
  // the lane lists hold the pairs within cutoff + prune skin (without pruning: the whole outer list, within
  // cutoff + Verlet skin), a hair above: a pair rounded (in float) to just beyond the threshold at the compaction
  // cannot come within the cutoff before the next compaction by more than that rounding
  innerCutoff = cutoff + (pruning() ? pruneSkin : verletSkin);
  const double inner = innerCutoff + 1e-4;
  parameters.innerCutoffSquared = static_cast<float>(inner * inner);
  pruneDisplacementSquared = 0.25 * pruneSkin * pruneSkin;
  parametersChanged = true;

  if (lennardJonesTable.size() > lennardJonesCapacity)
  {
    lennardJonesCapacity = lennardJonesTable.size();
    lennardJonesBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_ONLY, lennardJonesCapacity * sizeof(float)));
  }
  if (!lennardJonesTable.empty())
  {
    OpenCL::writeBuffer(lennardJonesBuffer.get(), lennardJonesTable.size() * sizeof(float), lennardJonesTable.data());
  }
}

void OpenCLPairKernel::beginBuild(const CellList& cells, const SimulationBox& box, std::size_t parts)
{
  buildStart = std::chrono::steady_clock::now();
  ++builds;

  // the coarse grid: cells at least the list cutoff wide (27-cell stencil), 1 cell where the box is smaller
  const double listCutoff = cells.listCutoff;
  const double3 widths = box.perpendicularWidths();
  auto cellsAlong = [&](double width)
  { return static_cast<std::int32_t>(std::clamp(std::floor(width / listCutoff), 1.0, 1024.0)); };
  grid = int3(cellsAlong(widths.x), cellsAlong(widths.y), cellsAlong(widths.z));
  const std::size_t numberOfCells =
      static_cast<std::size_t>(grid.x) * static_cast<std::size_t>(grid.y) * static_cast<std::size_t>(grid.z);
  const std::size_t numberOfAtoms = cells.numberOfAtoms;
  const double3x3& inverseCell = box.inverseCell;
  const double gx = static_cast<double>(grid.x), gy = static_cast<double>(grid.y), gz = static_cast<double>(grid.z);
  const double subdivisions = static_cast<double>(1u << mortonBits);

  // sort key per atom: coarse cell, then the Morton code of the sub-cell within it
  sortKey.resize(numberOfAtoms);
  cellCount.assign(numberOfCells + 1, 0);
  for (std::size_t k = 0; k < numberOfAtoms; ++k)
  {
    double3 s = inverseCell * double3(cells.wrappedX[k], cells.wrappedY[k], cells.wrappedZ[k]);
    s = double3(s.x - std::floor(s.x), s.y - std::floor(s.y), s.z - std::floor(s.z));
    const double fx = s.x * gx, fy = s.y * gy, fz = s.z * gz;
    const std::int32_t cx = std::clamp(static_cast<std::int32_t>(fx), 0, grid.x - 1);
    const std::int32_t cy = std::clamp(static_cast<std::int32_t>(fy), 0, grid.y - 1);
    const std::int32_t cz = std::clamp(static_cast<std::int32_t>(fz), 0, grid.z - 1);
    const std::uint32_t ux = static_cast<std::uint32_t>(
        std::clamp(static_cast<std::int32_t>((fx - static_cast<double>(cx)) * subdivisions), 0, 7));
    const std::uint32_t uy = static_cast<std::uint32_t>(
        std::clamp(static_cast<std::int32_t>((fy - static_cast<double>(cy)) * subdivisions), 0, 7));
    const std::uint32_t uz = static_cast<std::uint32_t>(
        std::clamp(static_cast<std::int32_t>((fz - static_cast<double>(cz)) * subdivisions), 0, 7));
    const std::uint32_t c = static_cast<std::uint32_t>((cz * grid.y + cy) * grid.x + cx);
    sortKey[k] = (c << (3 * mortonBits)) | mortonKey(ux, uy, uz);
    ++cellCount[c + 1];
  }
  // counting sort by cell (cellCount[c] becomes the end of cell c), then the Morton order within every cell
  for (std::size_t c = 0; c < numberOfCells; ++c) cellCount[c + 1] += cellCount[c];
  order.resize(numberOfAtoms);
  for (std::size_t k = 0; k < numberOfAtoms; ++k)
  {
    const std::uint32_t c = sortKey[k] >> (3 * mortonBits);
    order[cellCount[c]++] = static_cast<std::uint32_t>(k);
  }
  cellSlotStart.resize(numberOfCells + 1);
  cellSlotStart[0] = 0;
  for (std::size_t c = 0; c < numberOfCells; ++c)
  {
    const std::uint32_t begin = c == 0 ? 0u : cellCount[c - 1];
    const std::uint32_t end = cellCount[c];
    std::sort(order.begin() + begin, order.begin() + end,
              [&](std::uint32_t p, std::uint32_t q) { return sortKey[p] < sortKey[q]; });
    cellSlotStart[c + 1] = cellSlotStart[c] + static_cast<std::uint32_t>(roundUp(end - begin, clusterI));
  }
  padded = cellSlotStart[numberOfCells];
  const std::size_t iClusters = padded / clusterI;

  // the slot layout and the build data per slot (dummy slots: far away, molecule noAtom)
  slotOfSorted.resize(numberOfAtoms);
  sortedOfSlot.assign(padded, noAtom);
  typeOfSlot.assign(padded, 0);
  buildPosition.assign(4 * padded, farAway);
  for (std::size_t slot = 0; slot < padded; ++slot) buildPosition[4 * slot + 3] = bitsAsFloat(noAtom);
  cellOfCluster.assign(iClusters, 0);
  for (std::size_t c = 0; c < numberOfCells; ++c)
  {
    const std::uint32_t begin = c == 0 ? 0u : cellCount[c - 1];
    const std::uint32_t end = cellCount[c];
    for (std::uint32_t m = begin; m < end; ++m)
    {
      const std::uint32_t k = order[m];
      const std::size_t slot = cellSlotStart[c] + (m - begin);
      slotOfSorted[k] = static_cast<std::uint32_t>(slot);
      sortedOfSlot[slot] = k;
      typeOfSlot[slot] = cells.type[k];
      buildPosition[4 * slot] = static_cast<float>(cells.wrappedX[k]);
      buildPosition[4 * slot + 1] = static_cast<float>(cells.wrappedY[k]);
      buildPosition[4 * slot + 2] = static_cast<float>(cells.wrappedZ[k]);
      buildPosition[4 * slot + 3] = bitsAsFloat(cells.moleculeId[k]);
    }
    for (std::size_t I = cellSlotStart[c] / clusterI; I < cellSlotStart[c + 1] / clusterI; ++I)
    {
      cellOfCluster[I] = static_cast<std::uint32_t>(c);
    }
  }

  // host staging for the new layout
  hostPositions.assign(4 * padded, 0.0f);
  hostRelative.assign(useBonded ? 4 * padded : 0, 0.0f);
  hostForces.assign(4 * padded, 0.0f);
  hostPartials.assign(partialsPerCluster * iClusters, 0.0f);
  partDisplacement.assign(parts, 0.0);
  laneListsValid = false;
  outerCount.assign(iClusters, 0);
  laneCount.assign(iClusters * groupSize, 0);

  // row and lane capacities: from the density on the first build, otherwise kept (grown on an overflow in
  // finishBuild / wait)
  const double density = static_cast<double>(numberOfAtoms) / box.volume;
  auto atomsWithin = [&](double reach) { return density * (4.0 / 3.0) * std::numbers::pi * reach * reach * reach; };
  if (blocksPerCluster == 0)
  {
    // plus the extent of a cluster
    blocksPerCluster =
        roundUp(std::max<std::size_t>(64, static_cast<std::size_t>(1.25 * atomsWithin(listCutoff + 4.0) / 4.0)), 32);
  }
  if (pairsPerLane == 0)
  {
    // a lane holds the pairs of one atom with the j-cluster slot b: a quarter of the atoms within the inner cutoff
    pairsPerLane = roundUp(
        std::max<std::size_t>(32, static_cast<std::size_t>(1.25 * atomsWithin(innerCutoff + 0.5) / clusterJ)), 32);
    parameters.pairsPerLane = static_cast<std::uint32_t>(pairsPerLane);
    parametersChanged = true;
  }

  const double3x3& cell = box.cell;
  const float cellValues[9] = {static_cast<float>(cell.ax), static_cast<float>(cell.ay), static_cast<float>(cell.az),
                               static_cast<float>(cell.bx), static_cast<float>(cell.by), static_cast<float>(cell.bz),
                               static_cast<float>(cell.cx), static_cast<float>(cell.cy), static_cast<float>(cell.cz)};
  const float inverseValues[9] = {
      static_cast<float>(inverseCell.ax), static_cast<float>(inverseCell.ay), static_cast<float>(inverseCell.az),
      static_cast<float>(inverseCell.bx), static_cast<float>(inverseCell.by), static_cast<float>(inverseCell.bz),
      static_cast<float>(inverseCell.cx), static_cast<float>(inverseCell.cy), static_cast<float>(inverseCell.cz)};
  std::copy(std::begin(cellValues), std::end(cellValues), std::begin(buildParameters.cell));
  std::copy(std::begin(inverseValues), std::end(inverseValues), std::begin(buildParameters.inverseCell));
  buildParameters.listCutoffSquared = static_cast<float>(listCutoff * listCutoff);
  buildParameters.gridX = static_cast<std::uint32_t>(grid.x);
  buildParameters.gridY = static_cast<std::uint32_t>(grid.y);
  buildParameters.gridZ = static_cast<std::uint32_t>(grid.z);
  buildParameters.blocksPerCluster = static_cast<std::uint32_t>(blocksPerCluster);
  buildParameters.orthorhombic =
      (cell.ay == 0.0 && cell.az == 0.0 && cell.bx == 0.0 && cell.bz == 0.0 && cell.cx == 0.0 && cell.cy == 0.0) ? 1u
                                                                                                                 : 0u;

  // uploads (the host arrays stay untouched until finishBuild has waited for the device) and the build kernels
  ensureBuffers();
  ensureBlockBuffers();
  ensureLaneBuffers();
  auto upload = [&](cl_mem buffer, std::size_t bytes, const void* host, std::string_view what)
  { check(clEnqueueWriteBuffer(queue.get(), buffer, CL_FALSE, 0, bytes, host, 0, nullptr, nullptr), what); };
  upload(buildPositionBuffer.get(), 4 * padded * sizeof(float), buildPosition.data(),
         "clEnqueueWriteBuffer (build positions)");
  upload(typeBuffer.get(), padded * sizeof(std::uint32_t), typeOfSlot.data(), "clEnqueueWriteBuffer (types)");
  upload(cellSlotStartBuffer.get(), (numberOfCells + 1) * sizeof(std::uint32_t), cellSlotStart.data(),
         "clEnqueueWriteBuffer (cell slots)");
  upload(cellOfClusterBuffer.get(), iClusters * sizeof(std::uint32_t), cellOfCluster.data(),
         "clEnqueueWriteBuffer (cell of cluster)");
  upload(buildParameterBuffer.get(), sizeof(BuildParameters), &buildParameters,
         "clEnqueueWriteBuffer (build parameters)");
  enqueueBuild();
  if (useMesh) mesh.setSlots(padded);
  if (useBonded) bonded.setLayout(queue.get(), slotOfSorted, cells.originalToSorted, padded);
  check(clFlush(queue.get()), "clFlush");
}

void OpenCLPairKernel::enqueueBuild()
{
  const std::size_t iClusters = padded / clusterI;
  const std::size_t jClusters = padded / clusterJ;
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

  setBuildArguments();
  const std::size_t global = std::max<std::size_t>(iClusters, 1) * groupSize;
  const std::size_t local = groupSize;
  check(clEnqueueNDRangeKernel(queue.get(), buildKernel.get(), 1, nullptr, &global, &local, 0, nullptr, nullptr),
        "clEnqueueNDRangeKernel (buildList)");
  cl_event event = nullptr;
  check(clEnqueueReadBuffer(queue.get(), outerCountBuffer.get(), CL_FALSE, 0, iClusters * sizeof(std::uint32_t),
                            outerCount.data(), 0, nullptr, &event),
        "clEnqueueReadBuffer (row counts)");
  buildEvent.reset(event);
}

void OpenCLPairKernel::finishBuild()
{
  for (;;)
  {
    cl_event event = buildEvent.get();
    check(clWaitForEvents(1, &event), "clWaitForEvents (build)");
    buildEvent.reset();
    maximumRow = 0;
    outerBlocks = 0;
    for (const std::uint32_t count : outerCount)
    {
      maximumRow = std::max<std::size_t>(maximumRow, count);
      outerBlocks += count;
    }
    if (maximumRow <= blocksPerCluster) break;
    // a row overflowed: grow the rows and build again
    blocksPerCluster = roundUp(maximumRow + maximumRow / 4, 32);
    buildParameters.blocksPerCluster = static_cast<std::uint32_t>(blocksPerCluster);
    check(clEnqueueWriteBuffer(queue.get(), buildParameterBuffer.get(), CL_TRUE, 0, sizeof(BuildParameters),
                               &buildParameters, 0, nullptr, nullptr),
          "clEnqueueWriteBuffer (build parameters)");
    ensureBlockBuffers();
    enqueueBuild();
  }
  parameters.blocksPerCluster = static_cast<std::uint32_t>(blocksPerCluster);
  parametersChanged = true;
  setListArguments();
  sampleRequested = true;
  builtTime += std::chrono::steady_clock::now() - buildStart;
}

void OpenCLPairKernel::ensureBuffers()
{
  const std::size_t iClusters = padded / clusterI;
  const std::size_t numberOfCells = cellSlotStart.size();
  bool blocksAffected = false;
  if (padded > slotCapacity)
  {
    slotCapacity = roundUp(padded + padded / 4, clusterI);
    positionBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_ONLY, slotCapacity * 4 * sizeof(float)));
    relativeBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_ONLY, slotCapacity * 4 * sizeof(float)));
    buildPositionBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_ONLY, slotCapacity * 4 * sizeof(float)));
    typeBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_ONLY, slotCapacity * sizeof(std::uint32_t)));
    forceBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_WRITE, slotCapacity * 4 * sizeof(float)));
    clusterMinBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_WRITE, (slotCapacity / clusterJ) * 4 * sizeof(float)));
    clusterMaxBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_WRITE, (slotCapacity / clusterJ) * 4 * sizeof(float)));
  }
  if (iClusters > iClusterCapacity)
  {
    iClusterCapacity = iClusters + iClusters / 4 + 1;
    cellOfClusterBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_ONLY, iClusterCapacity * sizeof(std::uint32_t)));
    outerCountBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_WRITE, iClusterCapacity * sizeof(std::uint32_t)));
    laneCountBuffer.reset(
        OpenCL::createBuffer(CL_MEM_READ_WRITE, iClusterCapacity * groupSize * sizeof(std::uint32_t)));
    partialBuffer.reset(OpenCL::createBuffer(CL_MEM_WRITE_ONLY, iClusterCapacity * partialsPerCluster * sizeof(float)));
    blocksAffected = true;
  }
  if (numberOfCells > cellCapacity)
  {
    cellCapacity = numberOfCells + numberOfCells / 4;
    cellSlotStartBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_ONLY, cellCapacity * sizeof(std::uint32_t)));
  }
  if (blocksAffected)
  {
    // the rows and lane lists are allocated for iClusterCapacity i-clusters
    blockWords = 0;
    laneWords = 0;
  }
}

void OpenCLPairKernel::ensureBlockBuffers()
{
  const std::size_t required = iClusterCapacity * blocksPerCluster;
  if (required <= blockWords) return;
  blockWords = required;
  const std::size_t bytes = blockWords * sizeof(std::uint32_t);
  outerClusterBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_WRITE, bytes));
  outerMaskBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_WRITE, bytes));
}

void OpenCLPairKernel::ensureLaneBuffers()
{
  const std::size_t required = iClusterCapacity * groupSize * pairsPerLane;
  if (required <= laneWords) return;
  laneWords = required;
  pairListBuffer.reset(OpenCL::createBuffer(CL_MEM_READ_WRITE, laneWords * sizeof(std::uint32_t)));
}

void OpenCLPairKernel::setBuildArguments()
{
  const cl_mem buffers[9] = {buildPositionBuffer.get(), cellSlotStartBuffer.get(), cellOfClusterBuffer.get(),
                             clusterMinBuffer.get(),    clusterMaxBuffer.get(),    buildParameterBuffer.get(),
                             outerClusterBuffer.get(),  outerMaskBuffer.get(),     outerCountBuffer.get()};
  setArguments(buildKernel.get(), buffers, "clSetKernelArg (buildList)");
}

void OpenCLPairKernel::setListArguments()
{
  const cl_mem compactBuffers[7] = {positionBuffer.get(),   outerClusterBuffer.get(), outerMaskBuffer.get(),
                                    outerCountBuffer.get(), parameterBuffer.get(),    pairListBuffer.get(),
                                    laneCountBuffer.get()};
  setArguments(compactKernel.get(), compactBuffers, "clSetKernelArg (compactList)");
  const cl_mem pairBuffers[8] = {positionBuffer.get(),  typeBuffer.get(),         pairListBuffer.get(),
                                 laneCountBuffer.get(), lennardJonesBuffer.get(), parameterBuffer.get(),
                                 forceBuffer.get(),     partialBuffer.get()};
  setArguments(pairKernel.get(), pairBuffers, "clSetKernelArg (clusterPairs)");
}

void OpenCLPairKernel::packPositions(std::size_t part, std::span<const std::uint32_t> sortedAtoms,
                                     const CellList& cells, const SimulationBox& box)
{
  const double3x3& cell = box.cell;
  float* out = hostPositions.data();
  float* relativeOut = hostRelative.data();
  const bool track = pruning() && laneListsValid;
  const float* reference = listPositions.data();
  const std::span<const std::uint32_t> referenceAtoms =
      useBonded ? bonded.referenceAtoms() : std::span<const std::uint32_t>{};
  double maximum = 0.0;
  for (const std::uint32_t k : sortedAtoms)
  {
    const int3 w = cells.wrap[k];
    const double3 translation =
        cell * double3(static_cast<double>(w.x), static_cast<double>(w.y), static_cast<double>(w.z));
    const std::size_t offset = 4 * static_cast<std::size_t>(slotOfSorted[k]);
    float* slot = out + offset;
    slot[0] = static_cast<float>(cells.x[k] - translation.x);
    slot[1] = static_cast<float>(cells.y[k] - translation.y);
    slot[2] = static_cast<float>(cells.z[k] - translation.z);
    slot[3] = static_cast<float>(cells.charge[k]);
    if (track)
    {
      const double dx = static_cast<double>(slot[0]) - static_cast<double>(reference[offset]);
      const double dy = static_cast<double>(slot[1]) - static_cast<double>(reference[offset + 1]);
      const double dz = static_cast<double>(slot[2]) - static_cast<double>(reference[offset + 2]);
      maximum = std::max(maximum, dx * dx + dy * dy + dz * dz);
    }
    if (useBonded)
    {
      // relative to the first atom of the molecule, in double: the bonded terms see ~1e-7 A of rounding
      const std::uint32_t r = referenceAtoms[k];
      float* relative = relativeOut + offset;
      relative[0] = static_cast<float>(cells.x[k] - cells.x[r]);
      relative[1] = static_cast<float>(cells.y[k] - cells.y[r]);
      relative[2] = static_cast<float>(cells.z[k] - cells.z[r]);
      relative[3] = 0.0f;
    }
  }
  partDisplacement[part] = maximum;
}

void OpenCLPairKernel::updateParameters(const SimulationBox& box)
{
  const double3x3& cell = box.cell;
  const double3x3& inverse = box.inverseCell;
  const float cellValues[9] = {static_cast<float>(cell.ax), static_cast<float>(cell.ay), static_cast<float>(cell.az),
                               static_cast<float>(cell.bx), static_cast<float>(cell.by), static_cast<float>(cell.bz),
                               static_cast<float>(cell.cx), static_cast<float>(cell.cy), static_cast<float>(cell.cz)};
  const float inverseValues[9] = {
      static_cast<float>(inverse.ax), static_cast<float>(inverse.ay), static_cast<float>(inverse.az),
      static_cast<float>(inverse.bx), static_cast<float>(inverse.by), static_cast<float>(inverse.bz),
      static_cast<float>(inverse.cx), static_cast<float>(inverse.cy), static_cast<float>(inverse.cz)};
  const std::uint32_t orthorhombic =
      (cell.ay == 0.0 && cell.az == 0.0 && cell.bx == 0.0 && cell.bz == 0.0 && cell.cx == 0.0 && cell.cy == 0.0) ? 1u
                                                                                                                 : 0u;
  bool changed = parametersChanged || parameters.orthorhombic != orthorhombic;
  for (std::size_t k = 0; k < 9; ++k)
  {
    changed = changed || parameters.cell[k] != cellValues[k] || parameters.inverseCell[k] != inverseValues[k];
  }
  if (!changed) return;
  std::copy(std::begin(cellValues), std::end(cellValues), std::begin(parameters.cell));
  std::copy(std::begin(inverseValues), std::end(inverseValues), std::begin(parameters.inverseCell));
  parameters.orthorhombic = orthorhombic;
  parametersChanged = false;
  check(clEnqueueWriteBuffer(queue.get(), parameterBuffer.get(), CL_FALSE, 0, sizeof(Parameters), &parameters, 0,
                             nullptr, nullptr),
        "clEnqueueWriteBuffer (parameters)");
}

void OpenCLPairKernel::enqueueMesh(bool interpolate)
{
  mesh.enqueue(queue.get(), positionBuffer.get(), forceBuffer.get(), interpolate);
}

void OpenCLPairKernel::enqueueBonded()
{
  bonded.enqueue(queue.get(), positionBuffer.get(), relativeBuffer.get(), typeBuffer.get(), forceBuffer.get());
}

void OpenCLPairKernel::enqueueCompaction()
{
  const std::size_t iClusters = padded / clusterI;
  const std::size_t global = std::max<std::size_t>(iClusters, 1) * groupSize;
  const std::size_t local = groupSize;
  check(clEnqueueNDRangeKernel(queue.get(), compactKernel.get(), 1, nullptr, &global, &local, 0, nullptr, nullptr),
        "clEnqueueNDRangeKernel (compactList)");
}

// the pair kernel, the mesh and bonded chains and the read-back of the forces and partial sums
void OpenCLPairKernel::enqueueChain()
{
  const std::size_t iClusters = padded / clusterI;
  const std::size_t global = std::max<std::size_t>(iClusters, 1) * groupSize;
  const std::size_t local = groupSize;
  check(clEnqueueNDRangeKernel(queue.get(), pairKernel.get(), 1, nullptr, &global, &local, 0, nullptr, nullptr),
        "clEnqueueNDRangeKernel (clusterPairs)");
  // start the device on the pair kernel now: the host-side cost of enqueueing the mesh and bonded chains (a dozen
  // launches) then overlaps with its execution
  check(clFlush(queue.get()), "clFlush");
  if (useMesh) enqueueMesh(true);
  if (useBonded) enqueueBonded();
  check(clEnqueueReadBuffer(queue.get(), forceBuffer.get(), CL_FALSE, 0, padded * 4 * sizeof(float), hostForces.data(),
                            0, nullptr, nullptr),
        "clEnqueueReadBuffer (forces)");
  if (compactedThisStep)
  {
    check(clEnqueueReadBuffer(queue.get(), laneCountBuffer.get(), CL_FALSE, 0,
                              iClusters * groupSize * sizeof(std::uint32_t), laneCount.data(), 0, nullptr, nullptr),
          "clEnqueueReadBuffer (lane counts)");
  }
  cl_event event = nullptr;
  check(clEnqueueReadBuffer(queue.get(), partialBuffer.get(), CL_FALSE, 0,
                            iClusters * partialsPerCluster * sizeof(float), hostPartials.data(), 0, nullptr, &event),
        "clEnqueueReadBuffer (partials)");
  if (useMesh)
  {
    clReleaseEvent(event);
    event = nullptr;
    mesh.enqueueRead(queue.get(), &event);
  }
  if (useBonded)
  {
    clReleaseEvent(event);
    event = nullptr;
    bonded.enqueueRead(queue.get(), &event);
  }
  readEvent.reset(event);
  check(clFlush(queue.get()), "clFlush");
}

void OpenCLPairKernel::enqueue(const SimulationBox& box)
{
  updateParameters(box);
  check(clEnqueueWriteBuffer(queue.get(), positionBuffer.get(), CL_FALSE, 0, padded * 4 * sizeof(float),
                             hostPositions.data(), 0, nullptr, nullptr),
        "clEnqueueWriteBuffer (positions)");
  if (useBonded)
  {
    check(clEnqueueWriteBuffer(queue.get(), relativeBuffer.get(), CL_FALSE, 0, padded * 4 * sizeof(float),
                               hostRelative.data(), 0, nullptr, nullptr),
          "clEnqueueWriteBuffer (relative positions)");
  }
  if (useMesh) mesh.updateBox(box);
  const std::size_t iClusters = padded / clusterI;
  const std::size_t global = std::max<std::size_t>(iClusters, 1) * groupSize;
  const std::size_t local = groupSize;
  // the fastest of a few back-to-back synchronous runs: the first run after an idle period also pays for the
  // clock ramp-up of the device
  auto timedChain = [&](auto&& enqueueWork) -> std::chrono::duration<double>
  {
    std::chrono::duration<double> best = std::chrono::duration<double>::max();
    for (std::size_t run = 0; run < 3; ++run)
    {
      check(clFinish(queue.get()), "clFinish");
      const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
      enqueueWork();
      check(clFinish(queue.get()), "clFinish");
      best = std::min(best, std::chrono::duration<double>(std::chrono::steady_clock::now() - start));
    }
    return best;
  };

  // the mesh chain and the bonded kernel add to the forces: their timing samples run before the pair kernel,
  // which overwrites the forces of every slot
  if (sampleRequested)
  {
    if (useMesh && !meshProfiled)
    {
      meshProfiled = true;
      mesh.profile(queue.get(), positionBuffer.get(), forceBuffer.get());
    }
    if (useMesh) sampledMeshTime = timedChain([&] { enqueueMesh(true); });
    if (useBonded) sampledBondedTime = timedChain([&] { enqueueBonded(); });
  }

  // compact the lane lists when the outer list is new or (with pruning) an atom moved more than prune skin / 2
  // since the last compaction
  bool compact = !laneListsValid;
  if (pruning())
  {
    for (const double displacement : partDisplacement) compact = compact || displacement > pruneDisplacementSquared;
  }
  compactedThisStep = compact;
  if (compact)
  {
    listPositions = hostPositions;
    laneListsValid = true;
    ++compactions;
    if (sampleRequested)
    {
      sampledPruneTime = timedChain([&] { enqueueCompaction(); });
    }
    else
    {
      enqueueCompaction();
    }
    prunedTime += sampledPruneTime;
  }

  // one synchronous timed run after every list build: the execution time of the kernel per step until the next
  // build (Apple's OpenCL profiling events do not report usable device times)
  if (sampleRequested)
  {
    sampledKernelTime = timedChain(
        [&]
        {
          check(clEnqueueNDRangeKernel(queue.get(), pairKernel.get(), 1, nullptr, &global, &local, 0, nullptr, nullptr),
                "clEnqueueNDRangeKernel (clusterPairs, timing sample)");
        });
    sampleRequested = false;
  }
  if (useMesh) meshedTime += sampledMeshTime;
  if (useBonded) bondedRunTime += sampledBondedTime;
  enqueueChain();
}

OpenCLPairKernel::Results OpenCLPairKernel::wait()
{
  for (;;)
  {
    cl_event event = readEvent.get();
    check(clWaitForEvents(1, &event), "clWaitForEvents");
    readEvent.reset();
    if (!compactedThisStep) break;
    // the lane counts of the compaction: the statistics, and the overflow check (a lane beyond its capacity:
    // grow the lists, compact again and evaluate the step again)
    std::size_t maximumLane = 0;
    std::size_t total = 0;
    for (const std::uint32_t count : laneCount)
    {
      maximumLane = std::max<std::size_t>(maximumLane, count);
      total += count;
    }
    listedPairs = total / 2;
    if (maximumLane <= pairsPerLane) break;
    pairsPerLane = roundUp(maximumLane + maximumLane / 4, 32);
    parameters.pairsPerLane = static_cast<std::uint32_t>(pairsPerLane);
    check(clEnqueueWriteBuffer(queue.get(), parameterBuffer.get(), CL_TRUE, 0, sizeof(Parameters), &parameters, 0,
                               nullptr, nullptr),
          "clEnqueueWriteBuffer (parameters)");
    ensureLaneBuffers();
    setListArguments();
    enqueueCompaction();
    enqueueChain();
  }
  compactedThisStep = false;
  kernelTime += sampledKernelTime;
  ++steps;

  // every pair is counted from both clusters: half the sums
  Results results{};
  const std::size_t iClusters = padded / clusterI;
  double sums[partialsPerCluster] = {};
  for (std::size_t I = 0; I < iClusters; ++I)
  {
    const float* partial = hostPartials.data() + I * partialsPerCluster;
    for (std::size_t q = 0; q < partialsPerCluster; ++q) sums[q] += static_cast<double>(partial[q]);
  }
  results.energyVDW = 0.5 * sums[0];
  results.energyCharge = 0.5 * sums[1];
  results.pairStrain.ax = 0.5 * sums[2];
  results.pairStrain.bx = 0.5 * sums[3];
  results.pairStrain.cx = 0.5 * sums[4];
  results.pairStrain.ay = 0.5 * sums[5];
  results.pairStrain.by = 0.5 * sums[6];
  results.pairStrain.cy = 0.5 * sums[7];
  results.pairStrain.az = 0.5 * sums[8];
  results.pairStrain.bz = 0.5 * sums[9];
  results.pairStrain.cz = 0.5 * sums[10];
  if (useMesh)
  {
    mesh.collect(results.reciprocalEnergy, results.reciprocalStrain, results.singleIonSum, results.singleIonStrain);
  }
  if (useBonded) results.bonded = bonded.collect();
  return results;
}

std::string OpenCLPairKernel::status() const
{
  std::string result;
  result += std::format("    pair kernel: OpenCL cluster {} x {} (single precision) Lennard-Jones{} on {}\n", clusterI,
                        clusterJ, parameters.useCharge ? " + analytic Ewald real space" : "", deviceName());
  const std::size_t iClusters = padded / clusterI;
  result += std::format(
      "    device list: {} atoms in {} slots ({} i-clusters in {} x {} x {} cells), rows of {} "
      "blocks, {} builds on the device\n",
      slotOfSorted.size(), padded, iClusters, grid.x, grid.y, grid.z, blocksPerCluster, builds);
  result +=
      std::format("    blocks: {} in the outer list ({:.1f} per i-cluster); lane lists of {} entries: {} pairs",
                  outerBlocks, iClusters > 0 ? static_cast<double>(outerBlocks) / static_cast<double>(iClusters) : 0.0,
                  pairsPerLane, listedPairs);
  if (pruning())
  {
    result += std::format(" within cutoff + {:.3f} A ({} compactions)\n", pruneSkin, compactions);
  }
  else
  {
    result += " (no pruning: the whole outer list)\n";
  }
  if (useMesh) result += mesh.status();
  if (useBonded) result += bonded.status();
  return result;
}
