module;

module spatial_decomposition_cell_list;

import std;

import int3;
import double3;
import double3x3;
import atom;
import component;
import intra_molecular_exclusions;
import simulationbox;
import spatial_decomposition_domain_decomposition;

namespace
{
constexpr std::uint32_t noIndex = std::numeric_limits<std::uint32_t>::max();

// periodic shift (nx, ny, nz) <-> packed code; the shifts are small integers but not limited to -1..1 because the
// atom positions are not wrapped into the box
std::uint32_t packShift(int nx, int ny, int nz)
{
  return static_cast<std::uint32_t>(nx + 128) | (static_cast<std::uint32_t>(ny + 128) << 8) |
         (static_cast<std::uint32_t>(nz + 128) << 16);
}

double3 unpackShift(std::uint32_t code)
{
  return double3(static_cast<double>(static_cast<int>(code & 0xFFu) - 128),
                 static_cast<double>(static_cast<int>((code >> 8) & 0xFFu) - 128),
                 static_cast<double>(static_cast<int>((code >> 16) & 0xFFu) - 128));
}
}  // namespace

// Neighbour cells whose pairs a cell evaluates, as a CSR table. The stencil spans +-stencilRange cells per axis
// (the cells are at least listCutoff / subdivision wide, so that many cells cover the list cutoff). Every unordered
// cell pair must be evaluated by exactly one of the two cells: a cell takes the neighbours reached with an offset
// in the upper half of the stencil (dz > 0, or dz == 0 and dy > 0, or dz == dy == 0 and dx > 0), which mirrors to
// the lower half for the neighbour. When the stencil wraps around a small axis a neighbour can be reached with
// both an upper and a lower offset; such a pair is assigned to the cell with the lower index instead (the same
// decision is reached from both sides). The minimum-image distance then selects the interacting image.
void CellList::buildStencil(const SimulationBox& box)
{
  const int3 cells = numberOfCells;
  const double3 widths = box.perpendicularWidths();
  auto range = [&](double width, int n) -> int
  {
    const double cellWidth = width / static_cast<double>(n);
    return std::max(1, static_cast<int>(std::ceil(listCutoff / cellWidth - 1e-9)));
  };
  stencilRange = int3(range(widths.x, cells.x), range(widths.y, cells.y), range(widths.z, cells.z));

  const std::size_t numberOfCellsTotal =
      static_cast<std::size_t>(cells.x) * static_cast<std::size_t>(cells.y) * static_cast<std::size_t>(cells.z);
  const std::size_t stencilSize = static_cast<std::size_t>(2 * stencilRange.x + 1) *
                                  static_cast<std::size_t>(2 * stencilRange.y + 1) *
                                  static_cast<std::size_t>(2 * stencilRange.z + 1);
  stencilStart.resize(numberOfCellsTotal + 1);
  stencilList.clear();
  stencilList.reserve(numberOfCellsTotal * (stencilSize / 2 + 1));

  // per candidate neighbour: the wrap of the offset (box translations that place the neighbour cell next to this
  // cell) and bit 0 = reached with an upper offset, bit 1 = reached with a lower offset
  struct Candidate
  {
    std::uint32_t neighbour;
    std::int8_t wrapX, wrapY, wrapZ;
    std::uint8_t flag;
    bool operator<(const Candidate& other) const { return neighbour < other.neighbour; }
  };
  std::vector<Candidate> scratch(stencilSize);
  auto wrapOf = [](int c, int d, int n) -> int
  {
    // c + d = nWrapped + w n with 0 <= nWrapped < n
    const int sum = c + d;
    return sum >= 0 ? sum / n : -((-sum + n - 1) / n);
  };
  for (int cz = 0; cz < cells.z; ++cz)
  {
    for (int cy = 0; cy < cells.y; ++cy)
    {
      for (int cx = 0; cx < cells.x; ++cx)
      {
        const std::uint32_t cellIndex = static_cast<std::uint32_t>((cz * cells.y + cy) * cells.x + cx);
        stencilStart[cellIndex] = static_cast<std::uint32_t>(stencilList.size());
        std::size_t count = 0;
        for (int dz = -stencilRange.z; dz <= stencilRange.z; ++dz)
        {
          const int wz = wrapOf(cz, dz, cells.z);
          const int nz = cz + dz - wz * cells.z;
          for (int dy = -stencilRange.y; dy <= stencilRange.y; ++dy)
          {
            const int wy = wrapOf(cy, dy, cells.y);
            const int ny = cy + dy - wy * cells.y;
            for (int dx = -stencilRange.x; dx <= stencilRange.x; ++dx)
            {
              const int wx = wrapOf(cx, dx, cells.x);
              const int nx = cx + dx - wx * cells.x;
              const std::uint32_t neighbour = static_cast<std::uint32_t>((nz * cells.y + ny) * cells.x + nx);
              if (neighbour == cellIndex) continue;
              const bool upper = (dz > 0) || (dz == 0 && dy > 0) || (dz == 0 && dy == 0 && dx > 0);
              scratch[count++] = {neighbour, static_cast<std::int8_t>(wx), static_cast<std::int8_t>(wy),
                                  static_cast<std::int8_t>(wz), static_cast<std::uint8_t>(upper ? 1 : 2)};
            }
          }
        }
        std::sort(scratch.begin(), scratch.begin() + static_cast<std::ptrdiff_t>(count));
        std::size_t k = 0;
        while (k < count)
        {
          const Candidate first = scratch[k];
          std::uint8_t flags = 0;
          std::size_t reached = 0;
          while (k < count && scratch[k].neighbour == first.neighbour)
          {
            flags |= scratch[k].flag;
            ++reached;
            ++k;
          }
          const bool evaluateHere = (flags == 1) || (flags == 3 && first.neighbour > cellIndex);
          if (evaluateHere)
          {
            // reached through more than one offset: the image is not fixed by the cell pair (minimum image)
            stencilList.push_back({first.neighbour, first.wrapX, first.wrapY, first.wrapZ,
                                   static_cast<std::uint8_t>(reached > 1 ? 1 : 0)});
          }
        }
      }
    }
  }
  stencilStart[numberOfCellsTotal] = static_cast<std::uint32_t>(stencilList.size());
}

int3 CellList::cellGridFor(const SimulationBox& box) const { return cellGridFor(box, subdivision); }

int3 CellList::cellGridFor(const SimulationBox& box, int cellsPerCutoff) const
{
  const double3 widths = box.perpendicularWidths();
  auto count = [&](double width) -> std::int32_t
  {
    const double n = std::floor(width * static_cast<double>(cellsPerCutoff) / listCutoff);
    return static_cast<std::int32_t>(std::clamp(n, 1.0, 4096.0));
  };
  return int3(count(widths.x), count(widths.y), count(widths.z));
}

double3 CellList::wrappedFractional(const SimulationBox& box, const double3& position)
{
  double3 s = box.inverseCell * position;
  s.x -= std::floor(s.x);
  s.y -= std::floor(s.y);
  s.z -= std::floor(s.z);
  // guard against rounding to exactly 1
  s.x = std::min(s.x, std::nextafter(1.0, 0.0));
  s.y = std::min(s.y, std::nextafter(1.0, 0.0));
  s.z = std::min(s.z, std::nextafter(1.0, 0.0));
  return s;
}

std::uint32_t CellList::cellIndexOf(const SimulationBox& box, const double3& position) const
{
  const double3 s = wrappedFractional(box, position);
  const int cx = std::clamp(static_cast<int>(s.x * static_cast<double>(numberOfCells.x)), 0, numberOfCells.x - 1);
  const int cy = std::clamp(static_cast<int>(s.y * static_cast<double>(numberOfCells.y)), 0, numberOfCells.y - 1);
  const int cz = std::clamp(static_cast<int>(s.z * static_cast<double>(numberOfCells.z)), 0, numberOfCells.z - 1);
  return static_cast<std::uint32_t>((cz * numberOfCells.y + cy) * numberOfCells.x + cx);
}

double3 CellList::shiftVector(const SimulationBox& box, std::uint32_t code) { return box.cell * unpackShift(code); }

void CellList::setup(const SimulationBox& box, double cutoffValue, double skinValue, std::size_t numberOfThreads,
                     std::optional<int3> requested)
{
  cutoff = cutoffValue;
  skin = std::max(0.0, skinValue);
  listCutoff = cutoff + skin;
  numberOfDomains = std::max<std::size_t>(1, numberOfThreads);
  requestedGrid = requested;

  const double3 widths = box.perpendicularWidths();
  const double halfWidth = 0.5 * std::min({widths.x, widths.y, widths.z});
  if (listCutoff > halfWidth)
  {
    throw std::runtime_error(
        std::format("[Spatial decomposition]: cutoff + skin ({:.3f} + {:.3f} = {:.3f} A) exceeds half the smallest "
                    "perpendicular width of the box ({:.3f} A); the neighbour lists need cutoff + skin <= half the "
                    "box; use a smaller 'VerletSkin', set an explicit 'CutOffCoulomb' / 'CutOffVDW' smaller than "
                    "half the box (the automatic Coulomb cutoff is exactly half the box), or use a larger box\n",
                    cutoff, skin, listCutoff, halfWidth));
  }

  // Cells of a third of the list cutoff (7^3 stencil): the stencil encloses the cutoff sphere three times tighter
  // than cells of the full list cutoff (27 cells), which makes the list build cheaper, while the cells still hold
  // enough atoms to keep the per-cell overhead small. The sub-domains are cut on the atom positions and do not
  // depend on the cell grid.
  subdivision = 3;
  numberOfCells = cellGridFor(box, subdivision);
  const int3 grid = DomainDecomposition::chooseGrid(numberOfDomains, numberOfCells, widths, requestedGrid);
  decomposition.build(grid, {});
  domains.assign(numberOfDomains, DomainLists{});
  buildStencil(box);
  numberOfBuilds = 0;
}

void CellList::updateCellGrid(const SimulationBox& box)
{
  const int3 cells = cellGridFor(box);
  if (cells.x == numberOfCells.x && cells.y == numberOfCells.y && cells.z == numberOfCells.z) return;
  numberOfCells = cells;
  buildStencil(box);
}

void CellList::bin(const SimulationBox& box, std::span<const Atom> atoms, std::span<const Component> componentList)
{
  components = componentList;
  numberOfAtoms = atoms.size();
  const std::size_t numberOfCellsTotal = static_cast<std::size_t>(numberOfCells.x) *
                                         static_cast<std::size_t>(numberOfCells.y) *
                                         static_cast<std::size_t>(numberOfCells.z);

  sortedToOriginal.resize(numberOfAtoms);
  originalToSorted.resize(numberOfAtoms);
  x.resize(numberOfAtoms);
  y.resize(numberOfAtoms);
  z.resize(numberOfAtoms);
  refX.resize(numberOfAtoms);
  refY.resize(numberOfAtoms);
  refZ.resize(numberOfAtoms);
  charge.resize(numberOfAtoms);
  scalingVDW.resize(numberOfAtoms);
  scalingCoulomb.resize(numberOfAtoms);
  type.resize(numberOfAtoms);
  moleculeId.resize(numberOfAtoms);
  atomInMolecule.resize(numberOfAtoms);
  componentOfAtom.resize(numberOfAtoms);
  cellOfAtom.resize(numberOfAtoms);
  ownerOfAtom.resize(numberOfAtoms);
  wrappedX.resize(numberOfAtoms);
  wrappedY.resize(numberOfAtoms);
  wrappedZ.resize(numberOfAtoms);
  wrap.resize(numberOfAtoms);
  cellStart.assign(numberOfCellsTotal + 1, 0);

  // the index of every atom within its molecule (the atoms of a molecule are consecutive)
  std::vector<std::uint32_t> inMolecule(numberOfAtoms);
  for (std::size_t start = 0; start < numberOfAtoms;)
  {
    std::size_t end = start + 1;
    while (end < numberOfAtoms && atoms[end].componentId == atoms[start].componentId &&
           atoms[end].moleculeId == atoms[start].moleculeId)
    {
      ++end;
    }
    if (!components.empty())
    {
      const std::size_t componentId = static_cast<std::size_t>(atoms[start].componentId);
      if (componentId >= components.size() ||
          components[componentId].intraMolecularPotentials.exclusions.numberOfAtoms != end - start)
      {
        throw std::runtime_error(std::format(
            "[Spatial decomposition]: molecule {} of component {} has {} atoms in the system but {} in the component\n",
            atoms[start].moleculeId, componentId, end - start,
            componentId < components.size() ? components[componentId].intraMolecularPotentials.exclusions.numberOfAtoms : 0));
      }
    }
    for (std::size_t k = start; k < end; ++k) inMolecule[k] = static_cast<std::uint32_t>(k - start);
    start = end;
  }

  // wrapped fractional positions: cell binning and the balanced sub-domain cuts
  std::vector<double3> fractional(numberOfAtoms);
  std::vector<std::uint32_t> cellOfOriginal(numberOfAtoms);
  for (std::size_t i = 0; i < numberOfAtoms; ++i)
  {
    const double3 s = wrappedFractional(box, atoms[i].position);
    fractional[i] = s;
    const int cx = std::clamp(static_cast<int>(s.x * static_cast<double>(numberOfCells.x)), 0, numberOfCells.x - 1);
    const int cy = std::clamp(static_cast<int>(s.y * static_cast<double>(numberOfCells.y)), 0, numberOfCells.y - 1);
    const int cz = std::clamp(static_cast<int>(s.z * static_cast<double>(numberOfCells.z)), 0, numberOfCells.z - 1);
    const std::uint32_t cell = static_cast<std::uint32_t>((cz * numberOfCells.y + cy) * numberOfCells.x + cx);
    cellOfOriginal[i] = cell;
    ++cellStart[cell + 1];
  }
  decomposition.build(decomposition.grid, fractional);

  // counting sort by cell
  for (std::size_t cell = 0; cell < numberOfCellsTotal; ++cell)
  {
    cellStart[cell + 1] += cellStart[cell];
  }
  std::vector<std::uint32_t> fill(cellStart.begin(), cellStart.end() - 1);
  for (std::size_t i = 0; i < numberOfAtoms; ++i)
  {
    const std::uint32_t sorted = fill[cellOfOriginal[i]]++;
    sortedToOriginal[sorted] = static_cast<std::uint32_t>(i);
    originalToSorted[i] = sorted;
  }

  for (DomainLists& domain : domains)
  {
    domain.ownedAtoms.clear();
  }
  for (std::size_t sorted = 0; sorted < numberOfAtoms; ++sorted)
  {
    const std::uint32_t original = sortedToOriginal[sorted];
    const Atom& atom = atoms[original];
    x[sorted] = atom.position.x;
    y[sorted] = atom.position.y;
    z[sorted] = atom.position.z;
    refX[sorted] = atom.position.x;
    refY[sorted] = atom.position.y;
    refZ[sorted] = atom.position.z;
    charge[sorted] = atom.charge;
    scalingVDW[sorted] = atom.scalingVDW;
    scalingCoulomb[sorted] = atom.scalingCoulomb;
    type[sorted] = atom.type;
    moleculeId[sorted] = atom.moleculeId;
    atomInMolecule[sorted] = inMolecule[original];
    componentOfAtom[sorted] = atom.componentId;
    cellOfAtom[sorted] = cellOfOriginal[original];
    // wrapped position and the integer translation removed: position = wrapped + cell * wrap
    const double3 s = box.inverseCell * atom.position;
    const int3 w(static_cast<std::int32_t>(std::floor(s.x)), static_cast<std::int32_t>(std::floor(s.y)),
                 static_cast<std::int32_t>(std::floor(s.z)));
    const double3 translation =
        box.cell * double3(static_cast<double>(w.x), static_cast<double>(w.y), static_cast<double>(w.z));
    wrap[sorted] = w;
    wrappedX[sorted] = atom.position.x - translation.x;
    wrappedY[sorted] = atom.position.y - translation.y;
    wrappedZ[sorted] = atom.position.z - translation.z;
    const std::uint32_t owner = decomposition.ownerOf(fractional[original]);
    ownerOfAtom[sorted] = owner;
    domains[owner].ownedAtoms.push_back(static_cast<std::uint32_t>(sorted));
  }

  cellAtBuild = box.cell;
  ++numberOfBuilds;
}

void CellList::buildLists(std::size_t domainIndex, const SimulationBox& box)
{
  DomainLists& domain = domains[domainIndex];
  const std::size_t owned = domain.ownedAtoms.size();
  const std::uint32_t me = static_cast<std::uint32_t>(domainIndex);
  domain.neighbourStart.resize(owned + 1);
  domain.neighbourList.clear();
  domain.imageAtom.clear();
  domain.imageShift.clear();
  domain.imageOwner.clear();
  domain.imageNext.clear();
  domain.maximumDisplacementSquared = 0.0;

  domain.localIndexOfAtom.assign(numberOfAtoms, noIndex);
  for (std::size_t k = 0; k < owned; ++k) domain.localIndexOfAtom[domain.ownedAtoms[k]] = static_cast<std::uint32_t>(k);
  domain.imageHead.assign(numberOfAtoms, noIndex);

  const double listCutoffSquared = listCutoff * listCutoff;
  const std::uint32_t ownedCount = static_cast<std::uint32_t>(owned);
  const double3x3& cell = box.cell;
  const double* wx = wrappedX.data();
  const double* wy = wrappedY.data();
  const double* wz = wrappedZ.data();

  // appends the image (j shifted by the box translation n) of an accepted pair
  auto add = [&](std::uint32_t j, int nx, int ny, int nz)
  {
    if (nx == 0 && ny == 0 && nz == 0 && ownerOfAtom[j] == me)
    {
      domain.neighbourList.push_back(domain.localIndexOfAtom[j]);
      return;
    }
    // the image as a ghost: reuse its slot when an earlier owned atom already found it
    const std::uint32_t code = packShift(nx, ny, nz);
    std::uint32_t slot = domain.imageHead[j];
    while (slot != noIndex && domain.imageShift[slot] != code) slot = domain.imageNext[slot];
    if (slot == noIndex)
    {
      slot = static_cast<std::uint32_t>(domain.imageAtom.size());
      domain.imageAtom.push_back(j);
      domain.imageShift.push_back(code);
      domain.imageOwner.push_back(ownerOfAtom[j]);
      domain.imageNext.push_back(domain.imageHead[j]);
      domain.imageHead[j] = slot;
    }
    domain.neighbourList.push_back(ownedCount + slot);
  };

  // box translations of the stencil entries of the current cell (the owned atoms are in cell order)
  std::uint32_t currentCell = noIndex;
  std::vector<double3> translation;

  for (std::size_t k = 0; k < owned; ++k)
  {
    const std::uint32_t i = domain.ownedAtoms[k];
    const std::uint32_t cellI = cellOfAtom[i];
    const double xi = wx[i];
    const double yi = wy[i];
    const double zi = wz[i];
    const int3 wrapI = wrap[i];
    const std::uint32_t moleculeI = moleculeId[i];
    const std::uint32_t inMoleculeI = atomInMolecule[i];
    const IntraMolecularExclusions* exclusions =
        components.empty() ? nullptr : &components[componentOfAtom[i]].intraMolecularPotentials.exclusions;
    // a same-molecule candidate is left out when it is excluded or scaled (or always, without components)
    auto skipped = [&](std::uint32_t j)
    {
      return moleculeId[j] == moleculeI &&
             (exclusions == nullptr || exclusions->isExcludedFromPairList(inMoleculeI, atomInMolecule[j]));
    };
    domain.neighbourStart[k] = static_cast<std::uint32_t>(domain.neighbourList.size());

    if (cellI != currentCell)
    {
      currentCell = cellI;
      translation.clear();
      for (std::uint32_t n = stencilStart[cellI]; n < stencilStart[cellI + 1]; ++n)
      {
        const StencilEntry& entry = stencilList[n];
        translation.push_back(cell * double3(static_cast<double>(entry.wrapX), static_cast<double>(entry.wrapY),
                                             static_cast<double>(entry.wrapZ)));
      }
    }

    // same cell: pairs stored once (j > i), no box translation between the wrapped positions
    for (std::uint32_t j = i + 1; j < cellStart[cellI + 1]; ++j)
    {
      const double dx = xi - wx[j];
      const double dy = yi - wy[j];
      const double dz = zi - wz[j];
      if (dx * dx + dy * dy + dz * dz >= listCutoffSquared || skipped(j)) continue;
      add(j, wrapI.x - wrap[j].x, wrapI.y - wrap[j].y, wrapI.z - wrap[j].z);
    }

    // the neighbour cells this cell evaluates (half stencil): every pair of the system is built exactly once
    for (std::uint32_t n = stencilStart[cellI], e = 0; n < stencilStart[cellI + 1]; ++n, ++e)
    {
      const StencilEntry& entry = stencilList[n];
      const std::uint32_t other = entry.cell;
      const std::uint32_t begin = cellStart[other];
      const std::uint32_t end = cellStart[other + 1];
      if (entry.ambiguous == 0) [[likely]]
      {
        // the neighbour cell's atoms translated next to this cell: image = wrapped_j + cell * w
        const double3 t = translation[e];
        const double xs = xi - t.x;
        const double ys = yi - t.y;
        const double zs = zi - t.z;
        for (std::uint32_t j = begin; j < end; ++j)
        {
          const double dx = xs - wx[j];
          const double dy = ys - wy[j];
          const double dz = zs - wz[j];
          if (dx * dx + dy * dy + dz * dz >= listCutoffSquared || skipped(j)) continue;
          // unwrapped: r_i - (r_j + cell n) with n = w + wrap_i - wrap_j
          add(j, entry.wrapX + wrapI.x - wrap[j].x, entry.wrapY + wrapI.y - wrap[j].y,
              entry.wrapZ + wrapI.z - wrap[j].z);
        }
      }
      else
      {
        // the neighbour cell is reached through more than one offset (few cells along an axis): minimum image
        const double3 ri(x[i], y[i], z[i]);
        for (std::uint32_t j = begin; j < end; ++j)
        {
          if (skipped(j)) continue;
          const double3 sv = box.inverseCell * (ri - double3(x[j], y[j], z[j]));
          const double3 nv(std::round(sv.x), std::round(sv.y), std::round(sv.z));
          const double3 dr = cell * (sv - nv);
          if (double3::dot(dr, dr) >= listCutoffSquared) continue;
          add(j, static_cast<int>(nv.x), static_cast<int>(nv.y), static_cast<int>(nv.z));
        }
      }
    }
  }
  domain.neighbourStart[owned] = static_cast<std::uint32_t>(domain.neighbourList.size());

  domain.imagesByOwner.assign(numberOfDomains, {});
  for (std::uint32_t slot = 0; slot < domain.imageAtom.size(); ++slot)
  {
    domain.imagesByOwner[domain.imageOwner[slot]].push_back(slot);
  }

  // compact per-domain atom data: owned atoms, then the images
  const std::size_t local = domain.numberOfLocalAtoms();
  domain.positions.resize(local);
  domain.localType.resize(local);
  domain.localScalingVDW.resize(local);
  domain.localScalingCoulomb.resize(local);
  for (std::size_t l = 0; l < local; ++l)
  {
    const std::uint32_t atom = l < owned ? domain.ownedAtoms[l] : domain.imageAtom[l - owned];
    domain.positions[l].charge = charge[atom];
    domain.localType[l] = type[atom];
    domain.localScalingVDW[l] = scalingVDW[atom];
    domain.localScalingCoulomb[l] = scalingCoulomb[atom];
  }
  gatherPositions(domainIndex, box);
}

void CellList::gatherPositions(std::size_t domainIndex, const SimulationBox& box)
{
  DomainLists& domain = domains[domainIndex];
  const std::size_t owned = domain.ownedAtoms.size();
  LocalAtom* positions = domain.positions.data();
  for (std::size_t k = 0; k < owned; ++k)
  {
    const std::uint32_t i = domain.ownedAtoms[k];
    positions[k].x = x[i];
    positions[k].y = y[i];
    positions[k].z = z[i];
  }
  const std::size_t images = domain.imageAtom.size();
  const double3x3& cell = box.cell;
  for (std::size_t slot = 0; slot < images; ++slot)
  {
    const std::uint32_t j = domain.imageAtom[slot];
    const double3 shift = cell * unpackShift(domain.imageShift[slot]);
    positions[owned + slot].x = x[j] + shift.x;
    positions[owned + slot].y = y[j] + shift.y;
    positions[owned + slot].z = z[j] + shift.z;
  }
}

bool CellList::refreshPositionsAndCheck(std::size_t domainIndex, std::span<const Atom> atoms)
{
  DomainLists& domain = domains[domainIndex];
  double maximumSquared = 0.0;
  for (const std::uint32_t i : domain.ownedAtoms)
  {
    const double3 p = atoms[sortedToOriginal[i]].position;
    x[i] = p.x;
    y[i] = p.y;
    z[i] = p.z;
    const double dx = p.x - refX[i];
    const double dy = p.y - refY[i];
    const double dz = p.z - refZ[i];
    maximumSquared = std::max(maximumSquared, dx * dx + dy * dy + dz * dz);
  }
  domain.maximumDisplacementSquared = maximumSquared;
  return maximumSquared > 0.25 * skin * skin;
}

bool CellList::boxChanged(const SimulationBox& box) const
{
  const double3x3& a = box.cell;
  const double3x3& b = cellAtBuild;
  return a.ax != b.ax || a.ay != b.ay || a.az != b.az || a.bx != b.bx || a.by != b.by || a.bz != b.bz || a.cx != b.cx ||
         a.cy != b.cy || a.cz != b.cz;
}

std::size_t CellList::totalPairs() const
{
  std::size_t total = 0;
  for (const DomainLists& domain : domains) total += domain.neighbourList.size();
  return total;
}

std::size_t CellList::totalImages() const
{
  std::size_t total = 0;
  for (const DomainLists& domain : domains) total += domain.imageAtom.size();
  return total;
}

std::string CellList::status() const
{
  std::string result;
  result += std::format("    cutoff {:.4f} A, skin {:.4f} A, list cutoff {:.4f} A\n", cutoff, skin, listCutoff);
  result += std::format("    cell grid {} x {} x {} ({} cells, {} per list cutoff, stencil {} x {} x {})\n    {}\n",
                        numberOfCells.x, numberOfCells.y, numberOfCells.z,
                        static_cast<std::size_t>(numberOfCells.x) * static_cast<std::size_t>(numberOfCells.y) *
                            static_cast<std::size_t>(numberOfCells.z),
                        subdivision, 2 * stencilRange.x + 1, 2 * stencilRange.y + 1, 2 * stencilRange.z + 1,
                        decomposition.status());
  for (std::size_t d = 0; d < domains.size(); ++d)
  {
    result += std::format("    domain {:>3}: {:>8} atoms, {:>10} pairs, {:>8} ghost images\n", d,
                          domains[d].ownedAtoms.size(), domains[d].neighbourList.size(), domains[d].imageAtom.size());
  }
  return result;
}
