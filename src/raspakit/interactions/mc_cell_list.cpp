module;

module mc_cell_list;

import std;

import int3;
import double3;
import double3x3;
import atom;
import simulationbox;

namespace
{
// Upper bound on the number of cells (a huge, dilute box would otherwise allocate a grid far larger than
// the atom count). Cells larger than the cut-off remain correct, only less selective.
constexpr std::size_t maximumNumberOfCells = 32768;
// A direction with fewer cells than this is collapsed to a single cell: with three cells the 3-wide
// stencil already visits every cell, so nothing is gained.
constexpr int minimumCellsPerDirection = 4;
}  // namespace

int3 MCCellList::gridFor(const SimulationBox& simulationBox, double cutOff)
{
  if (!(cutOff > 0.0) || simulationBox.volume <= 0.0) return int3(1, 1, 1);

  const double3 widths = simulationBox.perpendicularWidths();
  int3 n(static_cast<int>(std::floor(widths.x / cutOff)), static_cast<int>(std::floor(widths.y / cutOff)),
         static_cast<int>(std::floor(widths.z / cutOff)));

  for (int i = 0; i != 3; ++i)
  {
    if (n[i] < minimumCellsPerDirection) n[i] = 1;
  }

  // Scale an oversized grid down uniformly (keeps the cells >= cut-off).
  while (static_cast<std::size_t>(n.x) * static_cast<std::size_t>(n.y) * static_cast<std::size_t>(n.z) >
         maximumNumberOfCells)
  {
    for (int i = 0; i != 3; ++i)
    {
      if (n[i] > 1) n[i] = std::max(1, n[i] / 2 >= minimumCellsPerDirection ? n[i] / 2 : 1);
    }
    if (n.x == 1 && n.y == 1 && n.z == 1) break;
  }

  return n;
}

bool MCCellList::wouldBeEnabled(const SimulationBox& simulationBox, double cutOff)
{
  const int3 n = gridFor(simulationBox, cutOff);
  return n.x > 1 || n.y > 1 || n.z > 1;
}

void MCCellList::build(std::span<const Atom> atoms, const SimulationBox& simulationBox, double cutOffDistance)
{
  cutOff = cutOffDistance;
  cellMatrix = simulationBox.cell;
  inverseCellMatrix = simulationBox.inverseCell;
  numberOfAtoms = atoms.size();
  numberOfCells = gridFor(simulationBox, cutOff);
  totalNumberOfCells = static_cast<std::size_t>(numberOfCells.x) * static_cast<std::size_t>(numberOfCells.y) *
                       static_cast<std::size_t>(numberOfCells.z);
  enabled = totalNumberOfCells > 1;
  ++numberOfBuilds;

  if (!enabled)
  {
    counts.clear();
    slots.clear();
    cellOfAtom.clear();
    slotOfAtom.clear();
    capacity = 0;
    valid = true;
    return;
  }

  // First pass: occupancy per cell, which sets the bucket capacity with a margin for density fluctuations.
  cellOfAtom.resize(atoms.size());
  counts.assign(totalNumberOfCells, 0u);
  std::uint32_t maximumCount = 0;
  for (std::size_t i = 0; i != atoms.size(); ++i)
  {
    const std::size_t cell = cellIndexOf(atoms[i].position);
    cellOfAtom[i] = static_cast<std::uint32_t>(cell);
    maximumCount = std::max(maximumCount, ++counts[cell]);
  }
  capacity = std::max<std::size_t>(
      16, static_cast<std::size_t>(maximumCount) + static_cast<std::size_t>(maximumCount) / 2 + 8);

  // Second pass: fill the buckets.
  slots.resize(totalNumberOfCells * capacity);
  slotOfAtom.resize(atoms.size());
  std::fill(counts.begin(), counts.end(), 0u);
  for (std::size_t i = 0; i != atoms.size(); ++i)
  {
    const std::size_t cell = cellOfAtom[i];
    const std::uint32_t slot = counts[cell]++;
    slotOfAtom[i] = slot;
    slots[cell * capacity + slot] = Record{atoms[i].position, static_cast<std::uint32_t>(i), atoms[i].moleculeId};
  }

  valid = true;
}

bool MCCellList::isCurrent(std::span<const Atom> atoms, const SimulationBox& simulationBox, double cutOffDistance) const
{
  return valid && atoms.size() == numberOfAtoms && cutOffDistance == cutOff && simulationBox.cell == cellMatrix;
}

void MCCellList::removeAtom(std::size_t atomIndex)
{
  const std::size_t cell = cellOfAtom[atomIndex];
  const std::uint32_t slot = slotOfAtom[atomIndex];
  const std::uint32_t last = --counts[cell];
  if (slot != last)
  {
    // move the last record of the bucket into the freed slot and fix its back-reference
    Record& moved = slots[cell * capacity + last];
    slots[cell * capacity + slot] = moved;
    slotOfAtom[moved.atomIndex] = slot;
  }
}

bool MCCellList::insertAtom(const Atom& atom, std::size_t atomIndex)
{
  const std::size_t cell = cellIndexOf(atom.position);
  const std::uint32_t slot = counts[cell];
  if (slot >= capacity) return false;
  counts[cell] = slot + 1;
  slots[cell * capacity + slot] = Record{atom.position, static_cast<std::uint32_t>(atomIndex), atom.moleculeId};
  cellOfAtom[atomIndex] = static_cast<std::uint32_t>(cell);
  slotOfAtom[atomIndex] = slot;
  return true;
}

void MCCellList::updateAtoms(std::span<const Atom> atoms, std::size_t first, std::size_t count)
{
  if (!valid || !enabled) return;
  if (atoms.size() != numberOfAtoms || first + count > numberOfAtoms)
  {
    valid = false;
    return;
  }

  for (std::size_t i = first; i != first + count; ++i)
  {
    const Atom& atom = atoms[i];
    const std::size_t newCell = cellIndexOf(atom.position);
    const std::size_t oldCell = cellOfAtom[i];
    Record& record = slots[oldCell * capacity + slotOfAtom[i]];
    if (newCell == oldCell)
    {
      record.setPosition(atom.position);
      record.moleculeId = atom.moleculeId;
      continue;
    }
    removeAtom(i);
    if (!insertAtom(atom, i))
    {
      // bucket overflow: the list is rebuilt (with a larger capacity) by the owner on the next query
      valid = false;
      return;
    }
  }
  numberOfAtomUpdates += count;
}

bool MCCellList::verify(std::span<const Atom> atoms, const SimulationBox& simulationBox) const
{
  if (!valid) return false;
  if (!enabled) return true;
  if (atoms.size() != numberOfAtoms) return false;
  if (!(simulationBox.cell == cellMatrix)) return false;

  std::size_t total = 0;
  for (std::size_t cell = 0; cell != totalNumberOfCells; ++cell)
  {
    if (counts[cell] > capacity) return false;
    total += counts[cell];
    for (std::uint32_t k = 0; k != counts[cell]; ++k)
    {
      const Record& record = slots[cell * capacity + k];
      if (record.atomIndex >= atoms.size()) return false;
      if (cellOfAtom[record.atomIndex] != cell || slotOfAtom[record.atomIndex] != k) return false;
      const Atom& atom = atoms[record.atomIndex];
      if (!(record.x == atom.position.x && record.y == atom.position.y && record.z == atom.position.z)) return false;
      if (record.moleculeId != atom.moleculeId) return false;
      if (cellIndexOf(atom.position) != cell) return false;
    }
  }
  return total == atoms.size();
}

double MCCellList::averageNeighbourhoodSize() const
{
  if (!valid || !enabled || totalNumberOfCells == 0) return static_cast<double>(numberOfAtoms);
  double sum = 0.0;
  for (int ix = 0; ix != numberOfCells.x; ++ix)
  {
    for (int iy = 0; iy != numberOfCells.y; ++iy)
    {
      for (int iz = 0; iz != numberOfCells.z; ++iz)
      {
        const int dxMin = numberOfCells.x > 1 ? -1 : 0, dxMax = numberOfCells.x > 1 ? 1 : 0;
        const int dyMin = numberOfCells.y > 1 ? -1 : 0, dyMax = numberOfCells.y > 1 ? 1 : 0;
        const int dzMin = numberOfCells.z > 1 ? -1 : 0, dzMax = numberOfCells.z > 1 ? 1 : 0;
        for (int dx = dxMin; dx <= dxMax; ++dx)
          for (int dy = dyMin; dy <= dyMax; ++dy)
            for (int dz = dzMin; dz <= dzMax; ++dz)
              sum += static_cast<double>(counts[cellIndex(
                  wrap(ix + dx, numberOfCells.x), wrap(iy + dy, numberOfCells.y), wrap(iz + dz, numberOfCells.z))]);
      }
    }
  }
  return sum / static_cast<double>(totalNumberOfCells);
}
