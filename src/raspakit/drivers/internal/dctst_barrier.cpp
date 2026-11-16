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

  // the Miller indices of the plane in the basis of the box: n . (cell m) = (cell^T n) . m for a lattice
  // translation cell m (m integer), so the normal of the lattice plane (h k l) is G = h1 b1 + h2 b2 + h3 b3
  // (cell^T G = h) and n . (cell m) = (h . m) / |G|: for coprime h the box repeats every 1 / |G| along it. The
  // given normal is snapped onto the lowest-index plane normal within ~0.8 degrees (e.g. the x axis of a nearly
  // orthogonal triclinic cell, which is not exactly normal to a lattice plane)
  const double3 h = transposedMultiply(simulationBox.cell, normal);
  const double largest = std::max({std::abs(h.x), std::abs(h.y), std::abs(h.z)});
  constexpr int maximumIndex = 12;
  constexpr double alignmentTolerance = 1e-4;  // 1 - cos of the largest accepted angle
  bool found = false;
  for (int c = 1; c <= maximumIndex && !found; ++c)
  {
    const double3 scaled = (static_cast<double>(c) / largest) * h;
    int3 rounded(static_cast<int>(std::lround(scaled.x)), static_cast<int>(std::lround(scaled.y)),
                 static_cast<int>(std::lround(scaled.z)));
    const int divisor = std::gcd(std::gcd(std::abs(rounded.x), std::abs(rounded.y)), std::abs(rounded.z));
    if (divisor == 0) continue;
    rounded = int3(rounded.x / divisor, rounded.y / divisor, rounded.z / divisor);
    const double3 reciprocal = transposedMultiply(
        simulationBox.inverseCell,
        double3(static_cast<double>(rounded.x), static_cast<double>(rounded.y), static_cast<double>(rounded.z)));
    const double cosine = double3::dot(reciprocal, normal) / reciprocal.length();
    if (cosine > 1.0 - alignmentTolerance)
    {
      millerIndices = rounded;
      period = 1.0 / reciprocal.length();
      normal = double3::normalize(reciprocal);
      found = true;
    }
  }
  if (!found)
  {
    throw std::invalid_argument(std::format(
        "dcTST: the barrier normal ({:.6f}, {:.6f}, {:.6f}) is not the normal of a lattice plane (h k l) of the "
        "periodic box (|h|, |k|, |l| <= {}); its periodic images would not be planes\n",
        normal.x, normal.y, normal.z, maximumIndex));
  }
  // the in-plane vectors of the (possibly snapped) normal
  {
    const double3 seed = (std::abs(normal.x) <= std::abs(normal.y) && std::abs(normal.x) <= std::abs(normal.z))
                             ? double3(1.0, 0.0, 0.0)
                         : (std::abs(normal.y) <= std::abs(normal.z)) ? double3(0.0, 1.0, 0.0)
                                                                      : double3(0.0, 0.0, 1.0);
    inPlaneU = double3::normalize(double3::cross(normal, seed));
    inPlaneV = double3::cross(normal, inPlaneU);
  }
  const double cagesPerPeriod = period / lambda;
  if (!(lambda > 0.0) || std::abs(cagesPerPeriod - std::round(cagesPerPeriod)) > 1e-3 || std::round(cagesPerPeriod) < 1.0)
  {
    throw std::invalid_argument(
        std::format("dcTST: the periodic box repeats every {:.6f} [A] along the barrier normal (plane ({} {} {})), "
                    "which is not a whole number of 'LatticeDistance' {:.6f} [A]\n",
                    period, millerIndices.x, millerIndices.y, millerIndices.z, lambda));
  }
  // The lattice distance is that whole fraction of the period, not the number as it was typed. The reaction
  // coordinate is wrapped on the period and then reduced on lambda: with 'LatticeDistance' 12.28 for a period of
  // 49.128 (four cages) the seam of the box at +-period/2 does not fold onto a plane q = k lambda but 0.004 A
  // beside it, and with an even number of cages that plane is the barrier top: a slab of that width about q* = 0
  // had five preimages where every other point has four (a spike of 25% in P(q*)), and one fewer for a lattice
  // distance typed too large.
  lambda = period / std::round(cagesPerPeriod);
}

double DcTSTBarrier::projection([[maybe_unused]] const SimulationBox& simulationBox, double3 r) const
{
  // the scalar distance from the plane, wrapped on the repeat distance of the box along the normal (the
  // scalar s . h wrapped in one dimension, not a component-wise minimum image)
  const double q = double3::dot(r - point, normal);
  double reduced = q - period * std::floor(q / period + 0.5);
  if (reduced >= 0.5 * period) reduced -= period;
  if (reduced < -0.5 * period) reduced += period;
  return reduced;
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

bool DcTSTBarrier::inPeriodicWindow(double x, double center, double halfWidth) const
{
  return std::abs(reducedCoordinate(x - center)) <= halfWidth;
}

double DcTSTBarrier::samplePeriodicWindow(RandomNumber& random, double center, double halfWidth) const
{
  // Uniform on the circle-arc of radius halfWidth about center. Sampling the clipped Cartesian
  // interval [center - w, center + w] ∩ [−λ/2, λ/2) would drop the wrapped cap, and the reverse
  // of a jump into that cap would have proposal density zero.
  if (halfWidth >= 0.5 * lambda)
  {
    return reducedCoordinate((random.uniform() - 0.5) * lambda);
  }
  return reducedCoordinate(center + (2.0 * random.uniform() - 1.0) * halfWidth);
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
  // onto the plane of this barrier (q = offset), not onto the nearest of its periodic images q = k lambda + offset:
  // the cages of the other molecules are counted from q = 0
  const double q = projection(simulationBox, r);
  return r + (offset - q) * normal;
}

// =================================================================================================
// Cage lattice
// =================================================================================================

DcTSTCageLattice::DcTSTCageLattice(const SimulationBox& unitCell, int3 numberOfUnitCells, double3 cageCenter,
                                   int3 cagesPerUnitCell)
    : origin(cageCenter)
{
  // (a count below one was taken as one without a word)
  if (cagesPerUnitCell.x < 1 || cagesPerUnitCell.y < 1 || cagesPerUnitCell.z < 1)
  {
    throw std::runtime_error(std::format("dcTST: 'CagesPerUnitCell' must be three positive integers (got {}, {}, {})\n",
                                         cagesPerUnitCell.x, cagesPerUnitCell.y, cagesPerUnitCell.z));
  }
  const int3 divisions(std::max(1, cagesPerUnitCell.x), std::max(1, cagesPerUnitCell.y),
                       std::max(1, cagesPerUnitCell.z));
  cell = double3x3(unitCell.cell[0] / static_cast<double>(divisions.x),
                   unitCell.cell[1] / static_cast<double>(divisions.y),
                   unitCell.cell[2] / static_cast<double>(divisions.z));
  inverseCell = cell.inverse();
  this->numberOfUnitCells = int3(std::max(1, numberOfUnitCells.x) * divisions.x,
                                 std::max(1, numberOfUnitCells.y) * divisions.y,
                                 std::max(1, numberOfUnitCells.z) * divisions.z);
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
  // the reciprocal vectors b_i are the rows of the inverse cell: the cage boundaries s_i = const are
  // perpendicular to them
  const double3 unit = double3::normalize(direction);
  const double3x3 reciprocal = inverseCell.transpose();
  AxisAlignment best{0, 1, -1.0};
  for (int axis = 0; axis < 3; ++axis)
  {
    const double3 vector = reciprocal[static_cast<std::size_t>(axis)];
    const double cosine = double3::dot(vector, unit) / vector.length();
    if (std::abs(cosine) > best.cosine)
    {
      best = AxisAlignment{axis, (cosine >= 0.0) ? 1 : -1, std::abs(cosine)};
    }
  }
  return best;
}

double DcTSTCageLattice::planeSpacing(int axis) const
{
  return 1.0 / inverseCell.transpose()[static_cast<std::size_t>(axis)].length();
}

std::optional<std::string> DcTSTCageLattice::inconsistencyWith(const SimulationBox& simulationBox,
                                                               const DcTSTBarrier& barrier) const
{
  const AxisAlignment alignment = axisAlongDirection(barrier.normal);
  if (alignment.cosine < 0.999)
  {
    return std::format("the barrier normal is not perpendicular to the cage planes of any cage axis (best "
                       "|cos| = {:.4f}); for a non-orthogonal cell use the reciprocal-lattice direction",
                       alignment.cosine);
  }
  // The barrier is a lattice plane (h1 h2 h3) of the box, whose axes are those of the cages: it is a cage boundary
  // only as (h 0 0), (0 h 0) or (0 0 h). The cosine alone lets a high-index plane through ((12 1 0) in a cell with
  // b = 2a is 2.4 degrees from (1 0 0), |cos| = 0.99913), which cuts across the cages beside the column.
  const int transverse1 = barrier.millerIndices[(alignment.axis + 1) % 3];
  const int transverse2 = barrier.millerIndices[(alignment.axis + 2) % 3];
  if (transverse1 != 0 || transverse2 != 0)
  {
    return std::format("the barrier plane ({} {} {}) is not a boundary plane of the cages (it is inclined to the "
                       "cage planes of axis {}: only one Miller index may differ from zero)",
                       barrier.millerIndices.x, barrier.millerIndices.y, barrier.millerIndices.z, alignment.axis);
  }
  const double spacing = planeSpacing(alignment.axis);
  if (std::abs(spacing - barrier.lambda) > 1e-3 * barrier.lambda)
  {
    return std::format("'LatticeDistance' {:.5f} [A] is not the spacing {:.5f} [A] of the cage planes along the "
                       "normal (set 'CagesPerUnitCell' to the number of cages per unit cell along each axis)",
                       barrier.lambda, spacing);
  }
  if (numberOfUnitCells[alignment.axis] < 2)
  {
    return std::string("the cage lattice has a single cage along the normal (cages A and B coincide)");
  }
  const double offset = std::abs(barrier.reducedCoordinate(barrier.projection(simulationBox, origin)));
  if (std::abs(offset - 0.5 * barrier.lambda) > 1e-3 * barrier.lambda)
  {
    return std::format("the cage centers lie {:.4f} [A] from the barrier plane, not half a 'LatticeDistance'",
                       offset);
  }
  return std::nullopt;
}

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

std::size_t DcTSTCageLattice::neighbor(std::size_t cage, double3 shift) const
{
  // integer steps along the cage axis whose planes are perpendicular to 'shift' (cageOf(center + shift) rounds
  // the components along the other axes too, and lands on a diagonal cage in a non-orthogonal cell)
  const AxisAlignment alignment = axisAlongDirection(shift);
  const double distance = shift.length() * alignment.cosine;
  const int steps = static_cast<int>(std::lround(distance / planeSpacing(alignment.axis)));
  return neighborAlongAxis(cage, alignment.axis, alignment.sign * steps);
}

std::pair<std::size_t, std::size_t> DcTSTCageLattice::cagesAcrossNearestPlane(double3 r, double3 direction) const
{
  const AxisAlignment alignment = axisAlongDirection(direction);
  const double3 f = inverseCell * (r - origin);
  const int3 n(static_cast<int>(std::floor(f.x + 0.5)), static_cast<int>(std::floor(f.y + 0.5)),
               static_cast<int>(std::floor(f.z + 0.5)));
  const std::size_t own = index(n);
  // where the point sits in its cage along the axis, in [-1/2, 1/2): toward 'direction' the nearest boundary
  // plane is the one ahead, and the cage is the one behind it
  const double within = f[alignment.axis] - static_cast<double>(n[alignment.axis]);
  if (static_cast<double>(alignment.sign) * within >= 0.0)
  {
    return {own, neighborAlongAxis(own, alignment.axis, alignment.sign)};
  }
  return {neighborAlongAxis(own, alignment.axis, -alignment.sign), own};
}

double DcTSTCageLattice::fractionAlongDirection(double3 r, double3 direction) const
{
  const AxisAlignment alignment = axisAlongDirection(direction);
  const double f = (inverseCell * (r - origin))[static_cast<std::size_t>(alignment.axis)];
  return static_cast<double>(alignment.sign) * (f - std::floor(f + 0.5));
}

double3 DcTSTCageLattice::randomPointInCage(RandomNumber& random, std::size_t cage) const
{
  return center(cage) + cell * double3(random.uniform() - 0.5, random.uniform() - 0.5, random.uniform() - 0.5);
}
