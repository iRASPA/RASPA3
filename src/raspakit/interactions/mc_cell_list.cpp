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
// the atom count). Cells larger than cutOff / subdivision remain correct, only less selective.
constexpr std::size_t maximumNumberOfCells = 32768;

// Whether a query of radius 'cutOff' along an axis of 'n' cells over the perpendicular width 'width' leaves
// cells unvisited (the stencil is narrower than the axis).
bool axisReduces(int n, double width, double cutOff)
{
  if (n <= 1) return false;
  const int range = static_cast<int>(std::ceil(cutOff * static_cast<double>(n) / width));
  return 2 * range + 1 < n;
}
}  // namespace

int3 MCCellList::gridForSubdivision(const SimulationBox& simulationBox, double cutOff, int subdivision)
{
  if (!(cutOff > 0.0) || simulationBox.volume <= 0.0 || subdivision < 1) return int3(1, 1, 1);

  const double3 widths = simulationBox.perpendicularWidths();
  const double cellWidth = cutOff / static_cast<double>(subdivision);
  int3 n(static_cast<int>(std::floor(widths.x / cellWidth)), static_cast<int>(std::floor(widths.y / cellWidth)),
         static_cast<int>(std::floor(widths.z / cellWidth)));

  // Scale an oversized grid down uniformly (keeps the cells >= cellWidth).
  while (static_cast<std::size_t>(std::max(n.x, 1)) * static_cast<std::size_t>(std::max(n.y, 1)) *
             static_cast<std::size_t>(std::max(n.z, 1)) >
         maximumNumberOfCells)
  {
    n.x = std::max(1, n.x / 2);
    n.y = std::max(1, n.y / 2);
    n.z = std::max(1, n.z / 2);
  }

  // A direction along which a cut-off query visits every cell gives no reduction: collapse it.
  for (int i = 0; i != 3; ++i)
  {
    if (!axisReduces(n[i], widths[static_cast<std::size_t>(i)], cutOff)) n[i] = 1;
  }

  return n;
}

std::pair<int3, int> MCCellList::gridFor(const SimulationBox& simulationBox, double cutOff, std::size_t numberOfAtoms)
{
  auto cellsOf = [](const int3& n)
  { return static_cast<double>(n.x) * static_cast<double>(n.y) * static_cast<double>(n.z); };
  auto enabled = [](const int3& n) { return n.x > 1 || n.y > 1 || n.z > 1; };

  std::optional<std::pair<int3, int>> coarsest{};
  for (int subdivision = maximumSubdivision; subdivision >= 1; --subdivision)
  {
    const int3 n = gridForSubdivision(simulationBox, cutOff, subdivision);
    if (!enabled(n)) continue;
    if (static_cast<double>(numberOfAtoms) >= minimumAtomsPerCell * cellsOf(n)) return {n, subdivision};
    coarsest = {n, subdivision};
  }
  return coarsest.value_or(std::pair<int3, int>{int3(1, 1, 1), maximumSubdivision});
}

bool MCCellList::wouldBeEnabled(const SimulationBox& simulationBox, double cutOff)
{
  const int3 n = gridForSubdivision(simulationBox, cutOff, maximumSubdivision);
  return n.x > 1 || n.y > 1 || n.z > 1;
}

void MCCellList::build(std::span<const Atom> atoms, const SimulationBox& simulationBox, double cutOffDistance)
{
  cutOff = cutOffDistance;
  cellMatrix = simulationBox.cell;
  inverseCellMatrix = simulationBox.inverseCell;
  perpendicularWidths = simulationBox.perpendicularWidths();
  orthorhombic = cellMatrix.ay == 0.0 && cellMatrix.az == 0.0 && cellMatrix.bx == 0.0 && cellMatrix.bz == 0.0 &&
                 cellMatrix.cx == 0.0 && cellMatrix.cy == 0.0;
  numberOfAtoms = atoms.size();
  std::tie(numberOfCells, subdivision) = gridFor(simulationBox, cutOff, numberOfAtoms);
  totalNumberOfCells = static_cast<std::size_t>(numberOfCells.x) * static_cast<std::size_t>(numberOfCells.y) *
                       static_cast<std::size_t>(numberOfCells.z);
  enabled = totalNumberOfCells > 1;
  ++numberOfBuilds;

  if (!enabled)
  {
    counts.clear();
    resizeSlots(0);
    cellOfAtom.clear();
    slotOfAtom.clear();
    capacity = 0;
    numberOfPartiallyCoupledAtoms = 0;
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
  capacity =
      std::max<std::size_t>(8, static_cast<std::size_t>(maximumCount) + static_cast<std::size_t>(maximumCount) / 2 + 4);

  // Second pass: fill the buckets.
  resizeSlots(totalNumberOfCells * capacity);
  slotOfAtom.resize(atoms.size());
  std::fill(counts.begin(), counts.end(), 0u);
  numberOfPartiallyCoupledAtoms = 0;
  for (std::size_t i = 0; i != atoms.size(); ++i)
  {
    const std::size_t cell = cellOfAtom[i];
    const std::uint32_t slot = counts[cell]++;
    slotOfAtom[i] = slot;
    writeSlot(cell * capacity + slot, atoms[i], i);
    numberOfPartiallyCoupledAtoms += slotPartiallyCoupled(cell * capacity + slot);
  }

  valid = true;
}

void MCCellList::resizeSlots(std::size_t numberOfSlots)
{
  slotX.resize(numberOfSlots);
  slotY.resize(numberOfSlots);
  slotZ.resize(numberOfSlots);
  slotCharge.resize(numberOfSlots);
  slotScalingVDW.resize(numberOfSlots);
  slotScalingCoulomb.resize(numberOfSlots);
  slotMoleculeIdValue.resize(numberOfSlots);
  slotType.resize(numberOfSlots);
  slotMoleculeId.resize(numberOfSlots);
  slotAtomIndex.resize(numberOfSlots);
  slotGroupId.resize(numberOfSlots);
}

void MCCellList::writeSlot(std::size_t slot, const Atom& atom, std::size_t atomIndex)
{
  const double3 wrapped = wrappedPosition(atom.position);
  slotX[slot] = wrapped.x;
  slotY[slot] = wrapped.y;
  slotZ[slot] = wrapped.z;
  slotCharge[slot] = atom.charge;
  slotScalingVDW[slot] = atom.scalingVDW;
  slotScalingCoulomb[slot] = atom.scalingCoulomb;
  slotMoleculeIdValue[slot] = static_cast<double>(atom.moleculeId);
  slotType[slot] = static_cast<std::uint32_t>(atom.type);
  slotMoleculeId[slot] = atom.moleculeId;
  slotAtomIndex[slot] = static_cast<std::uint32_t>(atomIndex);
  slotGroupId[slot] = static_cast<std::uint8_t>(atom.groupId);
}

bool MCCellList::slotPartiallyCoupled(std::size_t slot) const
{
  return slotScalingVDW[slot] != 1.0 || slotScalingCoulomb[slot] != 1.0 || slotGroupId[slot] != 0;
}

void MCCellList::copySlot(std::size_t to, std::size_t from)
{
  slotX[to] = slotX[from];
  slotY[to] = slotY[from];
  slotZ[to] = slotZ[from];
  slotCharge[to] = slotCharge[from];
  slotScalingVDW[to] = slotScalingVDW[from];
  slotScalingCoulomb[to] = slotScalingCoulomb[from];
  slotMoleculeIdValue[to] = slotMoleculeIdValue[from];
  slotType[to] = slotType[from];
  slotMoleculeId[to] = slotMoleculeId[from];
  slotAtomIndex[to] = slotAtomIndex[from];
  slotGroupId[to] = slotGroupId[from];
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
  numberOfPartiallyCoupledAtoms -= slotPartiallyCoupled(cell * capacity + slot);
  if (slot != last)
  {
    // move the last record of the bucket into the freed slot and fix its back-reference
    copySlot(cell * capacity + slot, cell * capacity + last);
    slotOfAtom[slotAtomIndex[cell * capacity + slot]] = slot;
  }
}

bool MCCellList::insertAtom(const Atom& atom, std::size_t atomIndex)
{
  const std::size_t cell = cellIndexOf(atom.position);
  const std::uint32_t slot = counts[cell];
  if (slot >= capacity) return false;
  counts[cell] = slot + 1;
  writeSlot(cell * capacity + slot, atom, atomIndex);
  numberOfPartiallyCoupledAtoms += slotPartiallyCoupled(cell * capacity + slot);
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
    if (newCell == oldCell)
    {
      // the atom stayed in its cell: refresh the slot in place (the position, and whatever else a move may
      // have changed, e.g. the scaling of a fractional molecule)
      const std::size_t slot = oldCell * capacity + slotOfAtom[i];
      numberOfPartiallyCoupledAtoms -= slotPartiallyCoupled(slot);
      writeSlot(slot, atom, i);
      numberOfPartiallyCoupledAtoms += slotPartiallyCoupled(slot);
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
  std::size_t partiallyCoupled = 0;
  for (std::size_t cell = 0; cell != totalNumberOfCells; ++cell)
  {
    if (counts[cell] > capacity) return false;
    total += counts[cell];
    for (std::uint32_t k = 0; k != counts[cell]; ++k)
    {
      const std::size_t slot = cell * capacity + k;
      partiallyCoupled += slotPartiallyCoupled(slot);
      const std::uint32_t atomIndex = slotAtomIndex[slot];
      if (atomIndex >= atoms.size()) return false;
      if (cellOfAtom[atomIndex] != cell || slotOfAtom[atomIndex] != k) return false;
      const Atom& atom = atoms[atomIndex];
      const double3 wrapped = wrappedPosition(atom.position);
      if (!(slotX[slot] == wrapped.x && slotY[slot] == wrapped.y && slotZ[slot] == wrapped.z)) return false;
      if (slotMoleculeId[slot] != atom.moleculeId) return false;
      if (slotMoleculeIdValue[slot] != static_cast<double>(atom.moleculeId)) return false;
      if (slotCharge[slot] != atom.charge || slotScalingVDW[slot] != atom.scalingVDW ||
          slotScalingCoulomb[slot] != atom.scalingCoulomb)
      {
        return false;
      }
      if (slotType[slot] != static_cast<std::uint32_t>(atom.type) ||
          slotGroupId[slot] != static_cast<std::uint8_t>(atom.groupId))
      {
        return false;
      }
      if (cellIndexOf(atom.position) != cell) return false;
    }
  }
  return total == atoms.size() && partiallyCoupled == numberOfPartiallyCoupledAtoms;
}

double MCCellList::averageNeighbourhoodSize() const
{
  if (!valid || !enabled || totalNumberOfCells == 0) return static_cast<double>(numberOfAtoms);

  int ranges[3];
  bool whole[3];
  for (std::size_t axis = 0; axis < 3; ++axis)
  {
    const int n = numberOfCells[static_cast<int>(axis)];
    ranges[axis] = stencilRange(axis, cutOff);
    whole[axis] = (n <= 1) || (2 * ranges[axis] + 1 >= n);
  }
  const int beginX = whole[0] ? 0 : -ranges[0], endX = whole[0] ? numberOfCells.x - 1 : ranges[0];
  const int beginY = whole[1] ? 0 : -ranges[1], endY = whole[1] ? numberOfCells.y - 1 : ranges[1];
  const int beginZ = whole[2] ? 0 : -ranges[2], endZ = whole[2] ? numberOfCells.z - 1 : ranges[2];

  double sum = 0.0;
  for (int ix = 0; ix != numberOfCells.x; ++ix)
  {
    for (int iy = 0; iy != numberOfCells.y; ++iy)
    {
      for (int iz = 0; iz != numberOfCells.z; ++iz)
      {
        for (int dx = beginX; dx <= endX; ++dx)
        {
          const int jx = whole[0] ? dx : wrap(ix + dx, numberOfCells.x);
          for (int dy = beginY; dy <= endY; ++dy)
          {
            const int jy = whole[1] ? dy : wrap(iy + dy, numberOfCells.y);
            for (int dz = beginZ; dz <= endZ; ++dz)
            {
              const int jz = whole[2] ? dz : wrap(iz + dz, numberOfCells.z);
              sum += static_cast<double>(counts[cellIndex(jx, jy, jz)]);
            }
          }
        }
      }
    }
  }
  return sum / static_cast<double>(totalNumberOfCells);
}
