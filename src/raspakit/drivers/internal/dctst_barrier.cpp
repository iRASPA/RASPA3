module;

module dctst_barrier;

import std;

import int3;
import double3;
import double3x3;
import simulationbox;
import randomnumbers;

DcTSTBarrier::DcTSTBarrier(const SimulationBox& simulationBox, double3 fractionalPosition, double3 cartesianNormal,
                           double latticeDistance)
    : point(simulationBox.cell * fractionalPosition),
      normal(double3::normalize(cartesianNormal)),
      lambda(latticeDistance)
{
  // two unit vectors spanning the plane: start from the Cartesian axis least aligned with the normal
  double3 seed = (std::abs(normal.x) <= std::abs(normal.y) && std::abs(normal.x) <= std::abs(normal.z))
                     ? double3(1.0, 0.0, 0.0)
                 : (std::abs(normal.y) <= std::abs(normal.z)) ? double3(0.0, 1.0, 0.0)
                                                              : double3(0.0, 0.0, 1.0);
  inPlaneU = double3::normalize(double3::cross(normal, seed));
  inPlaneV = double3::cross(normal, inPlaneU);
}

double DcTSTBarrier::projection(const SimulationBox& simulationBox, double3 r) const
{
  return double3::dot(simulationBox.applyPeriodicBoundaryConditions(r - point), normal);
}

double DcTSTBarrier::displacement(const SimulationBox& simulationBox, double3 rOld, double3 rNew) const
{
  return double3::dot(simulationBox.applyPeriodicBoundaryConditions(rNew - rOld), normal);
}

std::int64_t DcTSTBarrier::cageIndex(double q) const { return static_cast<std::int64_t>(std::floor(q / lambda)); }

std::size_t DcTSTBarrier::shellOfCage(std::int64_t cage)
{
  return static_cast<std::size_t>(cage >= 0 ? cage : -1 - cage);
}

double DcTSTBarrier::reducedCoordinate(double q) const
{
  double reduced = q - lambda * std::floor(q / lambda + 0.5);
  if (reduced >= 0.5 * lambda) reduced -= lambda;
  if (reduced < -0.5 * lambda) reduced += lambda;
  return reduced;
}

std::size_t DcTSTBarrier::bin(double q, std::size_t numberOfBins) const
{
  // bins of width lambda / numberOfBins, bin 0 centered on the barrier (q = 0)
  const double width = lambda / static_cast<double>(numberOfBins);
  const double reduced = reducedCoordinate(q);
  std::int64_t index = static_cast<std::int64_t>(std::floor(reduced / width + 0.5));
  const std::int64_t n = static_cast<std::int64_t>(numberOfBins);
  index = ((index % n) + n) % n;
  return static_cast<std::size_t>(index);
}

double DcTSTBarrier::binCenter(std::size_t index, std::size_t numberOfBins) const
{
  const double width = lambda / static_cast<double>(numberOfBins);
  double q = static_cast<double>(index) * width;
  if (q >= 0.5 * lambda) q -= lambda;
  return q;
}

double3 DcTSTBarrier::randomPointOnPlane(RandomNumber& random, const SimulationBox& simulationBox, double offset) const
{
  // a uniformly random point of the box, shifted along the normal onto the requested plane; the
  // in-plane coordinates stay uniformly distributed over the periodic plane
  return projectOntoPlane(simulationBox, simulationBox.randomPosition(random), offset);
}

double3 DcTSTBarrier::projectOntoPlane(const SimulationBox& simulationBox, double3 r, double offset) const
{
  const double q = projection(simulationBox, r);
  const double reduced = reducedCoordinate(q);
  return r + (offset - reduced) * normal;
}

// =================================================================================================
// Cage lattice
// =================================================================================================

DcTSTCageLattice::DcTSTCageLattice(const SimulationBox& unitCell, int3 numberOfUnitCells, double3 cageCenter)
    : cell(unitCell.cell),
      inverseCell(unitCell.inverseCell),
      numberOfUnitCells(std::max(1, numberOfUnitCells.x), std::max(1, numberOfUnitCells.y),
                        std::max(1, numberOfUnitCells.z)),
      origin(cageCenter)
{
}

std::size_t DcTSTCageLattice::numberOfCages() const
{
  return static_cast<std::size_t>(numberOfUnitCells.x) * static_cast<std::size_t>(numberOfUnitCells.y) *
         static_cast<std::size_t>(numberOfUnitCells.z);
}

std::size_t DcTSTCageLattice::index(int3 n) const
{
  for (int i = 0; i < 3; ++i)
  {
    const int size = numberOfUnitCells[i];
    n[i] = ((n[i] % size) + size) % size;
  }
  return (static_cast<std::size_t>(n.x) * static_cast<std::size_t>(numberOfUnitCells.y) +
          static_cast<std::size_t>(n.y)) *
             static_cast<std::size_t>(numberOfUnitCells.z) +
         static_cast<std::size_t>(n.z);
}

std::size_t DcTSTCageLattice::cageOf(double3 r) const
{
  const double3 f = inverseCell * (r - origin);
  return index(int3(static_cast<int>(std::floor(f.x + 0.5)), static_cast<int>(std::floor(f.y + 0.5)),
                    static_cast<int>(std::floor(f.z + 0.5))));
}

std::size_t DcTSTCageLattice::neighborAlongAxis(std::size_t cage, int axis, int step) const
{
  int3 n = indices(cage);
  n[axis] += step;
  return index(n);
}

double3 DcTSTCageLattice::fractionalOffset(const SimulationBox& simulationBox, double3 r, std::size_t cage) const
{
  return inverseCell * simulationBox.applyPeriodicBoundaryConditions(r - center(cage));
}

DcTSTCageLattice::AxisAlignment DcTSTCageLattice::axisAlongDirection(double3 direction) const
{
  const double3 unit = double3::normalize(direction);
  AxisAlignment best{0, 1, -1.0};
  for (int axis = 0; axis < 3; ++axis)
  {
    const double3 vector = cell[static_cast<std::size_t>(axis)];
    const double cosine = double3::dot(vector, unit) / vector.length();
    if (std::abs(cosine) > best.cosine)
    {
      best = AxisAlignment{axis, (cosine >= 0.0) ? 1 : -1, std::abs(cosine)};
    }
  }
  return best;
}

double DcTSTCageLattice::axisLength(int axis) const { return cell[static_cast<std::size_t>(axis)].length(); }

int3 DcTSTCageLattice::indices(std::size_t cage) const
{
  const std::size_t nz = static_cast<std::size_t>(numberOfUnitCells.z);
  const std::size_t ny = static_cast<std::size_t>(numberOfUnitCells.y);
  return int3(static_cast<int>(cage / (ny * nz)), static_cast<int>((cage / nz) % ny), static_cast<int>(cage % nz));
}

double3 DcTSTCageLattice::center(std::size_t cage) const
{
  const int3 n = indices(cage);
  return origin + cell * double3(static_cast<double>(n.x), static_cast<double>(n.y), static_cast<double>(n.z));
}

std::size_t DcTSTCageLattice::neighbor(std::size_t cage, double3 shift) const { return cageOf(center(cage) + shift); }

double3 DcTSTCageLattice::randomPointInCage(RandomNumber& random, std::size_t cage) const
{
  return center(cage) + cell * double3(random.uniform() - 0.5, random.uniform() - 0.5, random.uniform() - 0.5);
}
