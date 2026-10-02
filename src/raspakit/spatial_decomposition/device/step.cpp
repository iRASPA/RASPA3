module;

module spatial_decomposition_device_step;

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
import spatial_decomposition_device_kernels;
import spatial_decomposition_device_mesh;
import spatial_decomposition_device_bonded;
import spatial_decomposition_device_bonded_topology;
import spatial_decomposition_opencl_context;
#ifdef RASPA_DEVICE_METAL
import spatial_decomposition_metal_context;
#endif

using DeviceKernelLayout::clusterJ;
using DeviceKernelLayout::pairGroupSize;
using DeviceKernelLayout::pairPartials;

namespace
{
constexpr std::size_t sampleInterval = 8;  // list builds between two synchronous timing samples of the kernels
constexpr std::size_t boundsGroupSize = 64;
constexpr float farAway = 1e30f;         // build position of a dummy slot: never within the list cutoff of a real atom
constexpr std::uint32_t mortonBits = 3;  // 8 x 8 x 8 sub-cells order the atoms within a cell

float bitsAsFloat(std::uint32_t bits) { return std::bit_cast<float>(bits); }

std::size_t roundUp(std::size_t value, std::size_t multiple) { return ((value + multiple - 1) / multiple) * multiple; }

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

template <typename Matrix>
void matrixToFloats(const Matrix& m, float* out)
{
  out[0] = static_cast<float>(m.ax);
  out[1] = static_cast<float>(m.ay);
  out[2] = static_cast<float>(m.az);
  out[3] = static_cast<float>(m.bx);
  out[4] = static_cast<float>(m.by);
  out[5] = static_cast<float>(m.bz);
  out[6] = static_cast<float>(m.cx);
  out[7] = static_cast<float>(m.cy);
  out[8] = static_cast<float>(m.cz);
}

std::uint32_t isOrthorhombic(const double3x3& cell)
{
  return (cell.ay == 0.0 && cell.az == 0.0 && cell.bx == 0.0 && cell.bz == 0.0 && cell.cx == 0.0 && cell.cy == 0.0)
             ? 1u
             : 0u;
}
}  // namespace

std::unique_ptr<DeviceContext> createDeviceContext(PairDevice device)
{
  switch (device)
  {
    case PairDevice::OpenCL:
      return createOpenCLContext();
    case PairDevice::Metal:
#ifdef RASPA_DEVICE_METAL
      return createMetalContext();
#else
      throw std::runtime_error("[Device pair kernel]: this build has no Metal backend\n");
#endif
    case PairDevice::CPU:
      break;
  }
  throw std::runtime_error("[Device pair kernel]: 'PairDevice' CPU has no device backend\n");
}

bool deviceAvailable(PairDevice device)
{
  switch (device)
  {
    case PairDevice::OpenCL:
      return openclAvailable();
    case PairDevice::Metal:
#ifdef RASPA_DEVICE_METAL
      return metalAvailable();
#else
      return false;
#endif
    case PairDevice::CPU:
      return false;
  }
  return false;
}

std::string deviceNameOf(PairDevice device)
{
  switch (device)
  {
    case PairDevice::OpenCL:
      return openclDeviceName();
    case PairDevice::Metal:
#ifdef RASPA_DEVICE_METAL
      return metalDeviceName();
#else
      return {};
#endif
    case PairDevice::CPU:
      return {};
  }
  return {};
}

DeviceStep::~DeviceStep()
{
  if (context) unmapHost();
}

bool DeviceStep::available(PairDevice device) { return deviceAvailable(device); }

std::string DeviceStep::deviceName(PairDevice device) { return deviceNameOf(device); }

bool DeviceStep::supportsBonded(const System& system, std::string& reason)
{
  return BondedTopology::supports(system, reason);
}

void DeviceStep::initialize(PairDevice device)
{
  if (context) unmapHost();
  // the buffers and kernels belong to the old context: release them first
  mesh = DeviceMesh{};
  bonded = DeviceBonded{};
  for (DeviceBufferOwner* buffer :
       {&parameterBuffer, &buildParameterBuffer, &lennardJonesBuffer, &positionBuffer, &buildPositionBuffer,
        &typeBuffer, &forceBuffer, &relativeBuffer, &clusterMinBuffer, &clusterMaxBuffer, &compactReferenceBuffer,
        &cellOfClusterBuffer, &outerCountBuffer, &laneCountBuffer, &partialBuffer, &cellSlotStartBuffer,
        &outerClusterBuffer, &outerMaskBuffer, &pairListBuffer})
  {
    buffer->reset();
  }
  context.reset();
  mappedPositions = nullptr;
  mappedRelative = nullptr;
  mappedForces = nullptr;
  mappedPositionFloats = mappedRelativeFloats = mappedForceFloats = 0;
  deviceKind = device;
  context = createDeviceContext(device);

  const std::string_view program = "pair";
  boundsKernel = context->compileKernel(program, deviceKernelPairSource, DeviceMath::Fast, "clusterBounds");
  buildKernel = context->compileKernel(program, deviceKernelPairSource, DeviceMath::Fast, "buildList");
  compactKernel = context->compileKernel(program, deviceKernelPairSource, DeviceMath::Fast, "compactList");
  pairKernel = context->compileKernel(program, deviceKernelPairSource, DeviceMath::Fast, "clusterPairs");
  parameterBuffer.allocate(*context, sizeof(DevicePairParameters), DeviceMemory::Device);
  buildParameterBuffer.allocate(*context, sizeof(DeviceBuildParameters), DeviceMemory::Device);
  lennardJonesCapacity = 0;

  // the device buffers are new: the capacities start over
  slotCapacity = 0;
  iClusterCapacity = 0;
  cellCapacity = 0;
  relativeCapacity = 0;
  blockWords = 0;
  laneWords = 0;
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

void DeviceStep::enableMesh(int3 meshSize, std::size_t interpolationOrder, double alpha, double conversionFactor)
{
  if (!mesh.initialized() || mesh.interpolationOrder() != std::clamp<std::size_t>(interpolationOrder, 3, 7))
  {
    mesh.initialize(*context, interpolationOrder);
  }
  mesh.setup(meshSize, alpha, conversionFactor);
  mesh.setSlots(padded);
  useMesh = true;
}

void DeviceStep::enableBonded(const System& system, double alpha, double conversionFactor, bool useCharge)
{
  bondedTopology.build(system);
  if (!bonded.initialized()) bonded.initialize(*context);
  bonded.setTopology(bondedTopology);
  bonded.setParameters(alpha, conversionFactor, useCharge);
  useBonded = true;
}

void DeviceStep::setEwaldAlpha(double alpha)
{
  if (useMesh) mesh.setAlpha(alpha);
  if (useBonded) bonded.setAlpha(alpha);
}

void DeviceStep::writeParameters(bool blocking)
{
  context->write(parameterBuffer.get(), 0, sizeof(DevicePairParameters), &parameters, blocking);
}

void DeviceStep::writeBuildParameters(bool blocking)
{
  context->write(buildParameterBuffer.get(), 0, sizeof(DeviceBuildParameters), &buildParameters, blocking);
}

void DeviceStep::setParameters(std::span<const LennardJonesPair> lennardJones, std::size_t types, bool useCharge,
                               double cutOffVDW, double cutOffCharge, double conversionFactor, double alpha,
                               double verletSkin, double pruneSkinValue)
{
  // float4 per type pair: 4 epsilon, sigma^6, shift (one 16-byte load per pair in the kernel)
  lennardJonesTable.resize(4 * types * types);
  for (std::size_t k = 0; k < types * types; ++k)
  {
    const LennardJonesPair& pair = lennardJones[k];
    const double sigma2 = 1.0 / pair.inverseSigma2;
    lennardJonesTable[4 * k] = static_cast<float>(pair.epsilon4);
    lennardJonesTable[4 * k + 1] = static_cast<float>(sigma2 * sigma2 * sigma2);
    lennardJonesTable[4 * k + 2] = static_cast<float>(pair.shift);
    lennardJonesTable[4 * k + 3] = 0.0f;
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
    lennardJonesBuffer.allocate(*context, lennardJonesCapacity * sizeof(float), DeviceMemory::Device);
  }
  if (!lennardJonesTable.empty())
  {
    context->write(lennardJonesBuffer.get(), 0, lennardJonesTable.size() * sizeof(float), lennardJonesTable.data(),
                   true);
  }
}

void DeviceStep::beginBuild(const CellList& cells, const SimulationBox& box, std::size_t parts)
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

  // host staging for the new layout (mapped device buffers; dummy slots stay at 0)
  partDisplacement.assign(parts, 0.0);
  laneListsValid = false;
  outerCount.assign(iClusters, 0);
  laneCount.assign(iClusters * pairGroupSize, 0);
  hostPartials.assign(pairPartials * iClusters, 0.0f);

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
  matrixToFloats(cell, buildParameters.cell);
  matrixToFloats(inverseCell, buildParameters.inverseCell);
  buildParameters.listCutoffSquared = static_cast<float>(listCutoff * listCutoff);
  buildParameters.gridX = static_cast<std::uint32_t>(grid.x);
  buildParameters.gridY = static_cast<std::uint32_t>(grid.y);
  buildParameters.gridZ = static_cast<std::uint32_t>(grid.z);
  buildParameters.blocksPerCluster = static_cast<std::uint32_t>(blocksPerCluster);
  buildParameters.orthorhombic = isOrthorhombic(cell);

  // uploads (the host arrays stay untouched until finishBuild has waited for the device) and the build kernels
  ensureBuffers();
  ensureBlockBuffers();
  ensureLaneBuffers();
  if (residentMode)
  {
    // the dummy slots must read as 0 (the integrator's pack kernel writes the real slots only)
    zeroPositions.assign(4 * padded, 0.0f);
    context->write(positionBuffer.get(), 0, 4 * padded * sizeof(float), zeroPositions.data(), false);
  }
  else
  {
    mapInputs();
    std::fill(mappedPositions, mappedPositions + 4 * padded, 0.0f);
  }
  context->write(buildPositionBuffer.get(), 0, 4 * padded * sizeof(float), buildPosition.data(), false);
  context->write(typeBuffer.get(), 0, padded * sizeof(std::uint32_t), typeOfSlot.data(), false);
  context->write(cellSlotStartBuffer.get(), 0, (numberOfCells + 1) * sizeof(std::uint32_t), cellSlotStart.data(),
                 false);
  context->write(cellOfClusterBuffer.get(), 0, iClusters * sizeof(std::uint32_t), cellOfCluster.data(), false);
  writeBuildParameters(false);
  enqueueBuild();
  if (useMesh) mesh.setSlots(padded);
  if (useBonded)
  {
    bondedTopology.layout(slotOfSorted, cells.originalToSorted, padded, slotMolecule, referenceOfSorted);
    bonded.setLayout(slotMolecule);
  }
  context->flush();
}

void DeviceStep::enqueueBuild()
{
  const std::size_t iClusters = padded / clusterI;
  const std::size_t jClusters = padded / clusterJ;
  const std::uint32_t numberOfJClusters = static_cast<std::uint32_t>(jClusters);
  {
    const DeviceArg arguments[] = {DeviceArg::of(buildPositionBuffer.get()), DeviceArg::value(numberOfJClusters),
                                   DeviceArg::of(clusterMinBuffer.get()), DeviceArg::of(clusterMaxBuffer.get())};
    context->launch(boundsKernel, arguments, roundUp(std::max<std::size_t>(jClusters, 1), boundsGroupSize) / boundsGroupSize,
                    boundsGroupSize);
  }
  {
    const DeviceArg arguments[] = {
        DeviceArg::of(buildPositionBuffer.get()),  DeviceArg::of(cellSlotStartBuffer.get()),
        DeviceArg::of(cellOfClusterBuffer.get()),  DeviceArg::of(clusterMinBuffer.get()),
        DeviceArg::of(clusterMaxBuffer.get()),     DeviceArg::of(buildParameterBuffer.get()),
        DeviceArg::of(outerClusterBuffer.get()),   DeviceArg::of(outerMaskBuffer.get()),
        DeviceArg::of(outerCountBuffer.get())};
    context->launch(buildKernel, arguments, std::max<std::size_t>(iClusters, 1), pairGroupSize);
  }
  context->read(outerCountBuffer.get(), 0, outerCount.size() * sizeof(std::uint32_t), outerCount.data());
  buildEvent = context->mark();
}

void DeviceStep::finishBuild()
{
  for (;;)
  {
    context->wait(buildEvent);
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
    writeBuildParameters(true);
    ensureBlockBuffers();
    enqueueBuild();
  }
  parameters.blocksPerCluster = static_cast<std::uint32_t>(blocksPerCluster);
  parametersChanged = true;
  // the synchronous timing samples cost a few chain executions each: the first build and then every few builds
  sampleRequested = (builds == 1) || (builds % sampleInterval == 0);
  builtTime += std::chrono::steady_clock::now() - buildStart;
}

void DeviceStep::ensureBuffers()
{
  const std::size_t iClusters = padded / clusterI;
  const std::size_t cellEntries = cellSlotStart.size();
  const std::size_t numberOfAtoms = slotOfSorted.size();
  bool blocksAffected = false;
  if (padded > slotCapacity)
  {
    // the host pointers into the slot buffers end with them
    unmapHost();
    slotCapacity = roundUp(padded + padded / 4, clusterI);
    positionBuffer.allocate(*context, slotCapacity * 4 * sizeof(float), DeviceMemory::Shared);
    buildPositionBuffer.allocate(*context, slotCapacity * 4 * sizeof(float), DeviceMemory::Device);
    typeBuffer.allocate(*context, slotCapacity * sizeof(std::uint32_t), DeviceMemory::Device);
    forceBuffer.allocate(*context, slotCapacity * 4 * sizeof(float), DeviceMemory::Shared);
    clusterMinBuffer.allocate(*context, (slotCapacity / clusterJ) * 4 * sizeof(float), DeviceMemory::Device);
    clusterMaxBuffer.allocate(*context, (slotCapacity / clusterJ) * 4 * sizeof(float), DeviceMemory::Device);
    compactReferenceBuffer.allocate(*context, slotCapacity * 4 * sizeof(float), DeviceMemory::Device);
  }
  if (useBonded)
  {
    const std::size_t needed = 4 * numberOfAtoms;
    if (needed > relativeCapacity)
    {
      if (mappedRelative != nullptr) context->unmap(relativeBuffer.get());
      mappedRelative = nullptr;
      mappedRelativeFloats = 0;
      relativeCapacity = needed + needed / 4;
      relativeBuffer.allocate(*context, relativeCapacity * sizeof(float), DeviceMemory::Shared);
    }
  }
  if (iClusters > iClusterCapacity)
  {
    iClusterCapacity = iClusters + iClusters / 4 + 1;
    cellOfClusterBuffer.allocate(*context, iClusterCapacity * sizeof(std::uint32_t), DeviceMemory::Device);
    outerCountBuffer.allocate(*context, iClusterCapacity * sizeof(std::uint32_t), DeviceMemory::Device);
    laneCountBuffer.allocate(*context, iClusterCapacity * pairGroupSize * sizeof(std::uint32_t),
                             DeviceMemory::Device);
    partialBuffer.allocate(*context, iClusterCapacity * pairPartials * sizeof(float), DeviceMemory::Device);
    blocksAffected = true;
  }
  if (cellEntries > cellCapacity)
  {
    cellCapacity = cellEntries + cellEntries / 4;
    cellSlotStartBuffer.allocate(*context, cellCapacity * sizeof(std::uint32_t), DeviceMemory::Device);
  }
  if (blocksAffected)
  {
    // the rows and lane lists are allocated for iClusterCapacity i-clusters
    blockWords = 0;
    laneWords = 0;
  }
}

void DeviceStep::ensureBlockBuffers()
{
  const std::size_t required = iClusterCapacity * blocksPerCluster;
  if (required <= blockWords) return;
  blockWords = required;
  outerClusterBuffer.allocate(*context, blockWords * sizeof(std::uint32_t), DeviceMemory::Device);
  outerMaskBuffer.allocate(*context, blockWords * sizeof(std::uint32_t), DeviceMemory::Device);
}

void DeviceStep::ensureLaneBuffers()
{
  const std::size_t required = iClusterCapacity * pairGroupSize * pairsPerLane;
  if (required <= laneWords) return;
  laneWords = required;
  pairListBuffer.allocate(*context, laneWords * sizeof(std::uint32_t), DeviceMemory::Device);
}

void DeviceStep::enqueueCompaction()
{
  const DeviceArg arguments[] = {DeviceArg::of(positionBuffer.get()),   DeviceArg::of(outerClusterBuffer.get()),
                                 DeviceArg::of(outerMaskBuffer.get()),  DeviceArg::of(outerCountBuffer.get()),
                                 DeviceArg::of(parameterBuffer.get()),  DeviceArg::of(pairListBuffer.get()),
                                 DeviceArg::of(laneCountBuffer.get())};
  context->launch(compactKernel, arguments, std::max<std::size_t>(padded / clusterI, 1), pairGroupSize);
}

void DeviceStep::enqueuePairs()
{
  const DeviceArg arguments[] = {DeviceArg::of(positionBuffer.get()),  DeviceArg::of(typeBuffer.get()),
                                 DeviceArg::of(pairListBuffer.get()),  DeviceArg::of(laneCountBuffer.get()),
                                 DeviceArg::of(lennardJonesBuffer.get()), DeviceArg::of(parameterBuffer.get()),
                                 DeviceArg::of(forceBuffer.get()),     DeviceArg::of(partialBuffer.get())};
  context->launch(pairKernel, arguments, std::max<std::size_t>(padded / clusterI, 1), pairGroupSize);
}

void DeviceStep::enqueueMesh() { mesh.enqueue(positionBuffer.get(), forceBuffer.get(), true); }

void DeviceStep::enqueueBonded() { bonded.enqueue(relativeBuffer.get(), forceBuffer.get()); }

void DeviceStep::packPositions(std::size_t part, std::span<const std::uint32_t> sortedAtoms, const CellList& cells,
                               const SimulationBox& box)
{
  const double3x3& cell = box.cell;
  float* out = mappedPositions;
  float* relativeOut = mappedRelative;
  const bool track = pruning() && laneListsValid;
  const float* reference = listPositions.data();
  const std::span<const std::uint32_t> referenceAtoms =
      useBonded ? std::span<const std::uint32_t>(referenceOfSorted) : std::span<const std::uint32_t>{};
  double maximum = 0.0;
  for (const std::uint32_t k : sortedAtoms)
  {
    const int3 w = cells.wrap[k];
    const double tx = cell.ax * static_cast<double>(w.x) + cell.bx * static_cast<double>(w.y) +
                      cell.cx * static_cast<double>(w.z);
    const double ty = cell.ay * static_cast<double>(w.x) + cell.by * static_cast<double>(w.y) +
                      cell.cy * static_cast<double>(w.z);
    const double tz = cell.az * static_cast<double>(w.x) + cell.bz * static_cast<double>(w.y) +
                      cell.cz * static_cast<double>(w.z);
    const std::size_t offset = 4 * static_cast<std::size_t>(slotOfSorted[k]);
    float* slot = out + offset;
    slot[0] = static_cast<float>(cells.x[k] - tx);
    slot[1] = static_cast<float>(cells.y[k] - ty);
    slot[2] = static_cast<float>(cells.z[k] - tz);
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
      // relative to the first atom of the molecule, in double: the bonded terms see ~1e-7 A of rounding; stored
      // per atom in the system order (the atoms of a molecule consecutive) with the charge
      const std::uint32_t r = referenceAtoms[k];
      float* relative = relativeOut + 4 * static_cast<std::size_t>(cells.sortedToOriginal[k]);
      relative[0] = static_cast<float>(cells.x[k] - cells.x[r]);
      relative[1] = static_cast<float>(cells.y[k] - cells.y[r]);
      relative[2] = static_cast<float>(cells.z[k] - cells.z[r]);
      relative[3] = slot[3];
    }
  }
  partDisplacement[part] = maximum;
}

void DeviceStep::updateParameters(const SimulationBox& box)
{
  float cellValues[9];
  float inverseValues[9];
  matrixToFloats(box.cell, cellValues);
  matrixToFloats(box.inverseCell, inverseValues);
  const std::uint32_t orthorhombic = isOrthorhombic(box.cell);
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
  writeParameters(false);
}

// the pair kernel, the mesh and bonded chains and the read-back of the partial sums
void DeviceStep::enqueueChain()
{
  enqueuePairs();
  // start the device on the pair kernel now: the host-side cost of enqueueing the mesh and bonded chains (a dozen
  // launches) then overlaps with its execution
  context->flush();
  if (compactedThisStep)
  {
    context->read(laneCountBuffer.get(), 0, laneCount.size() * sizeof(std::uint32_t), laneCount.data());
  }
  if (useMesh) enqueueMesh();
  if (useBonded) enqueueBonded();
  context->read(partialBuffer.get(), 0, hostPartials.size() * sizeof(float), hostPartials.data());
  if (useMesh) mesh.enqueueRead();
  if (useBonded) bonded.enqueueRead();
  stepEvent = context->mark();
  context->flush();
}

void DeviceStep::enqueue(const SimulationBox& box)
{
  // compact the lane lists when the outer list is new or (with pruning) an atom moved more than prune skin / 2
  // since the last compaction
  bool compact = !laneListsValid;
  if (pruning())
  {
    for (const double displacement : partDisplacement) compact = compact || displacement > pruneDisplacementSquared;
  }
  if (compact) listPositions.assign(mappedPositions, mappedPositions + 4 * padded);
  // the mapped host pointers must not be live while the device uses the buffers
  unmapHost();
  enqueueStep(box, compact);
}

void DeviceStep::setResident(bool resident)
{
  if (resident == residentMode) return;
  residentMode = resident;
  // the reference positions of the compaction live on the other side after a switch: compact at the next step
  laneListsValid = false;
  if (resident)
  {
    unmapHost();
  }
  else if (padded > 0)
  {
    mapForces();
    mapInputs();
  }
}

DeviceStep::ResidentTargets DeviceStep::residentTargets() const
{
  ResidentTargets targets{};
  targets.positions = positionBuffer.get();
  targets.forces = forceBuffer.get();
  targets.relative = relativeBuffer.get();
  targets.buildPositions = buildPositionBuffer.get();
  targets.compactReference = compactReferenceBuffer.get();
  targets.relativeEnabled = useBonded;
  return targets;
}

void DeviceStep::enqueueResident(const SimulationBox& box, bool compactionDue)
{
  const bool compact = !laneListsValid || (pruning() && compactionDue);
  if (compact)
  {
    // the reference positions of the compaction stay on the device
    context->copy(positionBuffer.get(), 0, compactReferenceBuffer.get(), 0, 4 * padded * sizeof(float));
  }
  enqueueStep(box, compact);
}

void DeviceStep::enqueueStep(const SimulationBox& box, bool compact)
{
  updateParameters(box);
  if (useMesh) mesh.updateBox(box);

  compactedThisStep = compact;
  if (compact)
  {
    laneListsValid = true;
    ++compactions;
  }

  // the fastest of a few back-to-back synchronous runs: the first run after an idle period also pays for the
  // clock ramp-up of the device
  auto timedChain = [&](auto&& enqueueWork) -> std::chrono::duration<double>
  {
    std::chrono::duration<double> best = std::chrono::duration<double>::max();
    for (std::size_t run = 0; run < 3; ++run)
    {
      context->finish();
      const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
      enqueueWork();
      context->finish();
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
      mesh.profile(positionBuffer.get(), forceBuffer.get());
    }
    if (useMesh) sampledMeshTime = timedChain([&] { enqueueMesh(); });
    if (useBonded) sampledBondedTime = timedChain([&] { enqueueBonded(); });
  }

  if (compact)
  {
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
  // build
  if (sampleRequested)
  {
    sampledKernelTime = timedChain([&] { enqueuePairs(); });
    sampleRequested = false;
  }
  if (useMesh) meshedTime += sampledMeshTime;
  if (useBonded) bondedRunTime += sampledBondedTime;
  enqueueChain();
}

DeviceStep::Results DeviceStep::wait(const std::function<void()>& afterRetry)
{
  for (;;)
  {
    context->wait(stepEvent);
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
    writeParameters(true);
    ensureLaneBuffers();
    enqueueCompaction();
    enqueueChain();
    if (afterRetry) afterRetry();
  }
  compactedThisStep = false;
  kernelTime += sampledKernelTime;
  ++steps;
  if (!residentMode)
  {
    mapForces();
    mapInputs();
  }

  // every pair is counted from both clusters: half the sums
  Results results{};
  const std::size_t iClusters = padded / clusterI;
  double sums[pairPartials] = {};
  for (std::size_t I = 0; I < iClusters; ++I)
  {
    const float* partial = hostPartials.data() + I * pairPartials;
    for (std::size_t q = 0; q < pairPartials; ++q) sums[q] += static_cast<double>(partial[q]);
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

std::string DeviceStep::status() const
{
  std::string result;
  result += std::format("    pair kernel: {} cluster {} x {} (single precision) Lennard-Jones{} on {}\n",
                        pairDeviceName(deviceKind), clusterI, clusterJ,
                        parameters.useCharge ? " + analytic Ewald real space" : "", context->deviceName());
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
    result += " (no pruning: compacted once per build)\n";
  }
  if (useMesh) result += mesh.status();
  if (useBonded) result += bonded.status();
  return result;
}

void DeviceStep::unmapHost()
{
  if (mappedPositions != nullptr) context->unmap(positionBuffer.get());
  if (mappedRelative != nullptr) context->unmap(relativeBuffer.get());
  if (mappedForces != nullptr) context->unmap(forceBuffer.get());
  mappedPositions = nullptr;
  mappedRelative = nullptr;
  mappedForces = nullptr;
  mappedPositionFloats = mappedRelativeFloats = mappedForceFloats = 0;
}

void DeviceStep::mapInputs()
{
  if (padded > 0 && (mappedPositions == nullptr || mappedPositionFloats < 4 * padded))
  {
    if (mappedPositions != nullptr) context->unmap(positionBuffer.get());
    mappedPositionFloats = 4 * padded;
    mappedPositions = static_cast<float*>(context->map(positionBuffer.get(), 4 * padded * sizeof(float), true));
  }
  if (useBonded && !slotOfSorted.empty())
  {
    const std::size_t floats = 4 * slotOfSorted.size();
    if (mappedRelative == nullptr || mappedRelativeFloats < floats)
    {
      if (mappedRelative != nullptr) context->unmap(relativeBuffer.get());
      mappedRelativeFloats = floats;
      mappedRelative = static_cast<float*>(context->map(relativeBuffer.get(), floats * sizeof(float), true));
    }
  }
}

void DeviceStep::mapForces()
{
  if (padded > 0 && (mappedForces == nullptr || mappedForceFloats < 4 * padded))
  {
    if (mappedForces != nullptr) context->unmap(forceBuffer.get());
    mappedForceFloats = 4 * padded;
    mappedForces = static_cast<const float*>(context->map(forceBuffer.get(), 4 * padded * sizeof(float), false));
  }
}
