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
import spatial_decomposition_device_backend;
import spatial_decomposition_device_bonded_topology;
import spatial_decomposition_opencl_backend;

using DeviceKernelLayout::pairGroupSize;
using DeviceKernelLayout::pairPartials;

namespace
{
constexpr std::size_t sampleInterval = 8;  // list builds between two synchronous timing samples of the kernels
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

DeviceStep::~DeviceStep()
{
  if (backend) unmapHost();
}

bool DeviceStep::available(PairDevice device)
{
  switch (device)
  {
    case PairDevice::OpenCL:
      return openclAvailable();
    case PairDevice::CPU:
      return false;
  }
  return false;
}

std::string DeviceStep::deviceName(PairDevice device)
{
  switch (device)
  {
    case PairDevice::OpenCL:
      return openclDeviceName();
    case PairDevice::CPU:
      return {};
  }
  return {};
}

bool DeviceStep::supportsBonded(const System& system, std::string& reason)
{
  return BondedTopology::supports(system, reason);
}

void DeviceStep::initialize(PairDevice device)
{
  if (backend) unmapHost();
  backend.reset();
  mappedPositions = nullptr;
  mappedRelative = nullptr;
  mappedForces = nullptr;
  deviceKind = device;
  switch (device)
  {
    case PairDevice::OpenCL:
      backend = createOpenCLBackend();
      break;
    case PairDevice::CPU:
      throw std::runtime_error("[Device pair kernel]: 'PairDevice' CPU has no device backend\n");
  }
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
  backend->enableMesh(meshSize, interpolationOrder, alpha, conversionFactor, padded);
  useMesh = true;
}

void DeviceStep::enableBonded(const System& system, double alpha, double conversionFactor, bool useCharge)
{
  bondedTopology.build(system);
  backend->enableBonded(bondedTopology, alpha, conversionFactor, useCharge);
  useBonded = true;
}

void DeviceStep::setEwaldAlpha(double alpha)
{
  if (useMesh) backend->setMeshAlpha(alpha);
  if (useBonded) backend->setBondedAlpha(alpha);
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

  backend->setLennardJones(lennardJonesTable);
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
  mapInputs();
  std::fill(mappedPositions, mappedPositions + 4 * padded, 0.0f);
  backend->uploadLayout(DeviceLayoutUpload{
      .buildPositions = std::span<const float>(buildPosition.data(), 4 * padded),
      .types = std::span<const std::uint32_t>(typeOfSlot.data(), padded),
      .cellSlotStart = std::span<const std::uint32_t>(cellSlotStart.data(), numberOfCells + 1),
      .cellOfCluster = std::span<const std::uint32_t>(cellOfCluster.data(), iClusters),
  });
  backend->writeBuildParameters(buildParameters, false);
  enqueueBuild();
  if (useMesh) backend->setMeshSlots(padded);
  if (useBonded)
  {
    bondedTopology.layout(slotOfSorted, cells.originalToSorted, padded, slotMolecule, referenceOfSorted);
    backend->setBondedLayout(slotMolecule);
  }
  backend->flush();
}

void DeviceStep::enqueueBuild() { backend->enqueueBuild(padded / clusterI, padded / clusterJ, outerCount); }

void DeviceStep::finishBuild()
{
  for (;;)
  {
    backend->waitBuild();
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
    backend->writeBuildParameters(buildParameters, true);
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
    // the backend unmaps the host pointers of the slot buffers before it replaces them
    mappedPositions = nullptr;
    mappedRelative = nullptr;
    mappedForces = nullptr;
    slotCapacity = roundUp(padded + padded / 4, clusterI);
    backend->allocateSlots(slotCapacity);
  }
  if (useBonded)
  {
    const std::size_t needed = 4 * numberOfAtoms;
    if (needed > relativeCapacity)
    {
      mappedRelative = nullptr;
      relativeCapacity = needed + needed / 4;
      backend->allocateRelative(relativeCapacity);
    }
  }
  if (iClusters > iClusterCapacity)
  {
    iClusterCapacity = iClusters + iClusters / 4 + 1;
    backend->allocateClusters(iClusterCapacity);
    blocksAffected = true;
  }
  if (cellEntries > cellCapacity)
  {
    cellCapacity = cellEntries + cellEntries / 4;
    backend->allocateCells(cellCapacity);
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
  backend->allocateBlocks(blockWords);
}

void DeviceStep::ensureLaneBuffers()
{
  const std::size_t required = iClusterCapacity * pairGroupSize * pairsPerLane;
  if (required <= laneWords) return;
  laneWords = required;
  backend->allocateLanes(laneWords);
}

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
  backend->writeParameters(parameters, false);
}

// the pair kernel, the mesh and bonded chains and the read-back of the partial sums
void DeviceStep::enqueueChain()
{
  const std::size_t iClusters = padded / clusterI;
  backend->enqueuePairs(iClusters);
  // start the device on the pair kernel now: the host-side cost of enqueueing the mesh and bonded chains (a dozen
  // launches) then overlaps with its execution
  backend->flush();
  if (compactedThisStep) backend->enqueueReadLaneCounts(laneCount);
  if (useMesh) backend->enqueueMesh();
  if (useBonded) backend->enqueueBonded();
  backend->enqueueReadPartials(hostPartials);
  backend->flush();
}

void DeviceStep::enqueue(const SimulationBox& box)
{
  updateParameters(box);
  if (useMesh) backend->updateMeshBox(box);

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
    listPositions.assign(mappedPositions, mappedPositions + 4 * padded);
    laneListsValid = true;
    ++compactions;
  }
  // the mapped host pointers must not be live while the device uses the buffers
  unmapHost();

  const std::size_t iClusters = padded / clusterI;
  // the fastest of a few back-to-back synchronous runs: the first run after an idle period also pays for the
  // clock ramp-up of the device
  auto timedChain = [&](auto&& enqueueWork) -> std::chrono::duration<double>
  {
    std::chrono::duration<double> best = std::chrono::duration<double>::max();
    for (std::size_t run = 0; run < 3; ++run)
    {
      backend->finish();
      const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
      enqueueWork();
      backend->finish();
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
      backend->profileMesh();
    }
    if (useMesh) sampledMeshTime = timedChain([&] { backend->enqueueMesh(); });
    if (useBonded) sampledBondedTime = timedChain([&] { backend->enqueueBonded(); });
  }

  if (compact)
  {
    if (sampleRequested)
    {
      sampledPruneTime = timedChain([&] { backend->enqueueCompaction(iClusters); });
    }
    else
    {
      backend->enqueueCompaction(iClusters);
    }
    prunedTime += sampledPruneTime;
  }

  // one synchronous timed run after every list build: the execution time of the kernel per step until the next
  // build
  if (sampleRequested)
  {
    sampledKernelTime = timedChain([&] { backend->enqueuePairs(iClusters); });
    sampleRequested = false;
  }
  if (useMesh) meshedTime += sampledMeshTime;
  if (useBonded) bondedRunTime += sampledBondedTime;
  enqueueChain();
}

DeviceStep::Results DeviceStep::wait()
{
  for (;;)
  {
    backend->waitStep();
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
    backend->writeParameters(parameters, true);
    ensureLaneBuffers();
    backend->enqueueCompaction(padded / clusterI);
    enqueueChain();
  }
  compactedThisStep = false;
  kernelTime += sampledKernelTime;
  ++steps;
  mapForces();
  mapInputs();

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
    backend->collectMesh(results.reciprocalEnergy, results.reciprocalStrain, results.singleIonSum,
                         results.singleIonStrain);
  }
  if (useBonded) results.bonded = backend->collectBonded();
  return results;
}

std::string DeviceStep::status() const
{
  std::string result;
  result += std::format("    pair kernel: {} cluster {} x {} (single precision) Lennard-Jones{} on {}\n",
                        pairDeviceName(deviceKind), clusterI, clusterJ,
                        parameters.useCharge ? " + analytic Ewald real space" : "", backend->deviceName());
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
  if (useMesh) result += backend->meshStatus();
  if (useBonded) result += backend->bondedStatus();
  return result;
}

void DeviceStep::unmapHost()
{
  backend->unmapAll();
  mappedPositions = nullptr;
  mappedRelative = nullptr;
  mappedForces = nullptr;
}

void DeviceStep::mapInputs()
{
  if (padded > 0) mappedPositions = backend->mapPositions(4 * padded);
  if (useBonded && !slotOfSorted.empty()) mappedRelative = backend->mapRelative(4 * slotOfSorted.size());
}

void DeviceStep::mapForces()
{
  if (padded > 0) mappedForces = backend->mapForces(4 * padded);
}
