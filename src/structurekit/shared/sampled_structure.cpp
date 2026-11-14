module;

module sampled_structure;

import std;

import double3;
import int3;
import double3x3;
import unit_cell;
import skspacegroupsetting;
import skspacegroupdatabase;
import units;

namespace
{
double maximumRadiusOf(const std::vector<double> &radii)
{
  return radii.empty() ? 0.0 : *std::ranges::max_element(radii);
}

// Splits a possibly out-of-range bin coordinate into a wrapped bin index and the periodic image it
// came from (same convention as PoreAccessibility).
std::pair<int, int> binAndImage(int coordinate, int gridExtent)
{
  int image = (coordinate >= 0) ? coordinate / gridExtent : -((-coordinate + gridExtent - 1) / gridExtent);
  return {coordinate - image * gridExtent, image};
}

int3 binOfFractional(const double3 &fractional, const int3 &gridSize)
{
  return int3(std::min(gridSize.x - 1, static_cast<int>(fractional.x * static_cast<double>(gridSize.x))),
              std::min(gridSize.y - 1, static_cast<int>(fractional.y * static_cast<double>(gridSize.y))),
              std::min(gridSize.z - 1, static_cast<int>(fractional.z * static_cast<double>(gridSize.z))));
}

// Bin of an unwrapped fractional coordinate (any real), for queries that must stay in the caller's frame.
int3 binOfUnwrappedFractional(const double3 &fractional, const int3 &gridSize)
{
  return int3(static_cast<int>(std::floor(fractional.x * static_cast<double>(gridSize.x))),
              static_cast<int>(std::floor(fractional.y * static_cast<double>(gridSize.y))),
              static_cast<int>(std::floor(fractional.z * static_cast<double>(gridSize.z))));
}

// The bins, built if they are not there or are not this structure's.
//
// Nearly every call finds them there, and many threads call at once from inside their own loops, so that way
// through takes no lock. It used to be a flag saying the bins were ready, followed by a look at what they
// were built for. That is sound for bins built once. It is not for bins built a second time: a thread that
// read the flag just before another began rebuilding went on to read the fields the builder was by then
// writing, could find them already those of the new build, and walked off into bins that were half made.
// Thousands of queries from eight threads after each change to a structure gave a few wrong answers in every
// thousand.
//
// So the flag is the number of the build, nothing while one is running, and it is read on both sides of the
// look. The builder puts it to nothing before it writes anything else, and each thing a reader looks at is
// handed over the way the number is, so a reader that saw anything of a new build cannot then read the old
// number again: it finds nothing, or a later number, and goes to the lock like any other that found the bins
// wanting. One that reads the same number twice with the fields right has the bins of a build that was
// finished before it looked and was not begun again while it looked.
//
// What comes back is the largest radius of the structure, which every query bounds its search by and which
// has been worked out here already. The queries take it from here and not from the grid, where it is one of
// the fields read on the way in and so is not to be read any other way.
double ensureNeighborGrid(const SampledStructure &structure)
{
  SampledNeighborGrid &grid = structure.neighborGrid;
  const std::size_t atomCount = structure.positions.size();

  std::atomic_ref<std::uint32_t> published{grid.published};
  std::atomic_ref<std::size_t> builtForAtoms{grid.atomCount};
  std::atomic_ref<double> builtForRadius{grid.maximumRadius};

  // The way of every query: the bins are there and were built for this many atoms. The largest radius is read
  // from the build, not worked out again: finding it is a pass over every atom, which on this way would be
  // made for every point asked about (see the header: a change of the radii at the same atom count is one the
  // structure's owner announces by putting a fresh grid in place).
  const std::uint32_t build = published.load(std::memory_order_acquire);
  if (build != 0 && builtForAtoms.load(std::memory_order_acquire) == atomCount)
  {
    const double builtRadius = builtForRadius.load(std::memory_order_acquire);
    if (published.load(std::memory_order_acquire) == build) return builtRadius;
  }

  // Building is serialized, and so is finding out under the lock that another thread has built meanwhile.
  // Queries on bins that are there stay on the way above, so the worker threads in SampledRoadmap::build do
  // not share this mutex.
  static std::mutex buildMutex;
  std::lock_guard<std::mutex> lock(buildMutex);
  const double maxRadius = maximumRadiusOf(structure.radii);
  if (published.load(std::memory_order_acquire) != 0 &&
      builtForAtoms.load(std::memory_order_acquire) == atomCount &&
      builtForRadius.load(std::memory_order_acquire) == maxRadius)
  {
    return maxRadius;
  }

  published.store(0, std::memory_order_release);
  builtForAtoms.store(atomCount, std::memory_order_release);
  builtForRadius.store(maxRadius, std::memory_order_release);

  // The widths across the cell, as the cell gives them, which is the same for a left-handed set of vectors as
  // for a right-handed one. Worked out here from a volume that carried the sign of the determinant they came
  // out negative for a left-handed cell, the least bin width with them, and the searches below, which stop
  // once a shell is further off than matters, had shells that were never any distance off.
  const UnitCell &unitCell = structure.unitCell;
  const double3x3 &cell = unitCell.cell;
  const double volume = std::fabs(unitCell.volume);
  const double3 perpendicularWidths = unitCell.perpendicularWidths();

  // Roughly four atoms per bin; bin counts proportional to the perpendicular widths so the bins are
  // approximately metrically cubic (same construction as PoreAccessibility / SKVoronoi).
  const double targetBinSize =
      std::cbrt(volume / std::max(1.0, static_cast<double>(atomCount) / 4.0));
  grid.gridSize = int3(std::max(1, static_cast<int>(perpendicularWidths.x / targetBinSize)),
                       std::max(1, static_cast<int>(perpendicularWidths.y / targetBinSize)),
                       std::max(1, static_cast<int>(perpendicularWidths.z / targetBinSize)));
  grid.minimumBinWidth = std::min({perpendicularWidths.x / static_cast<double>(grid.gridSize.x),
                                   perpendicularWidths.y / static_cast<double>(grid.gridSize.y),
                                   perpendicularWidths.z / static_cast<double>(grid.gridSize.z)});

  const int3 gridSize = grid.gridSize;
  grid.wrappedPositions.resize(structure.positions.size());
  grid.bins.assign(static_cast<std::size_t>(gridSize.x) * static_cast<std::size_t>(gridSize.y) *
                       static_cast<std::size_t>(gridSize.z),
                   {});
  for (std::size_t i = 0; i < structure.positions.size(); ++i)
  {
    double3 fractional = double3::fract(unitCell.inverseCell * structure.positions[i]);
    grid.wrappedPositions[i] = cell * fractional;
    int3 bin = binOfFractional(fractional, gridSize);
    grid.bins[static_cast<std::size_t>((bin.z * gridSize.y + bin.y) * gridSize.x + bin.x)].push_back(i);
  }

  // The next number, passing over the one that stands for nothing.
  grid.lastBuild = (grid.lastBuild + 1 == 0) ? 1 : grid.lastBuild + 1;
  published.store(grid.lastBuild, std::memory_order_release);
  return maxRadius;
}
}  // namespace

double SampledStructure::density() const
{
  return 1e-3 * this->mass /
         (this->unitCell.volume * Units::Angstrom * Units::Angstrom * Units::Angstrom *
          Units::AvogadroConstant);
}

double SampledStructure::gravimetricFactor() const
{
  return Units::Angstrom * Units::Angstrom * Units::AvogadroConstant / this->mass;
}

double SampledStructure::gravimetricVolumeFactor() const
{
  return (Units::Angstrom * Units::Angstrom * Units::Angstrom * 1.0e6) * Units::AvogadroConstant / this->mass;
}

bool SampledStructure::overlaps(const double3 &position, std::size_t skip) const
{
  if (this->positions.empty()) return false;

  const double maximumRadius = ensureNeighborGrid(*this);
  const SampledNeighborGrid &grid = this->neighborGrid;
  const int3 gridSize = grid.gridSize;
  const double3x3 &cell = this->unitCell.cell;

  // Stay in the caller's frame (MC SA places samples on the stored centre). Wrapping the query
  // first would rename the lattice image of `skip` that coincides with that centre.
  const double3 fractional = this->unitCell.inverseCell * position;
  const int3 pointBin = binOfUnwrappedFractional(fractional, gridSize);

  // Only atoms within the largest contact radius can contain the point. Walk bin shells outward
  // (wrapping periodically); all atoms in shell k are at least (k-1)·(minimum bin width) away.
  for (int k = 0;; ++k)
  {
    double lowerBound = static_cast<double>(k - 1) * grid.minimumBinWidth;
    if (lowerBound > maximumRadius) break;

    for (int ox = -k; ox <= k; ++ox)
    {
      for (int oy = -k; oy <= k; ++oy)
      {
        for (int oz = -k; oz <= k; ++oz)
        {
          if (std::max({std::abs(ox), std::abs(oy), std::abs(oz)}) != k) continue;

          auto [bx, lx] = binAndImage(pointBin.x + ox, gridSize.x);
          auto [by, ly] = binAndImage(pointBin.y + oy, gridSize.y);
          auto [bz, lz] = binAndImage(pointBin.z + oz, gridSize.z);

          double3 translation = cell * double3(static_cast<double>(lx), static_cast<double>(ly), static_cast<double>(lz));

          for (std::size_t j : grid.bins[static_cast<std::size_t>((bz * gridSize.y + by) * gridSize.x + bx)])
          {
            double3 atomImage = grid.wrappedPositions[j] + translation;
            // Ignore the image of `skip` that is the sphere centre the caller sampled on (the stored
            // position); other images of that atom must still be tested when 2r exceeds a cell edge.
            if (j == skip)
            {
              double3 fromCentre = atomImage - this->positions[skip];
              if (double3::dot(fromCentre, fromCentre) < 1.0e-12) continue;
            }
            double3 delta = position - atomImage;
            if (double3::dot(delta, delta) < this->radii[j] * this->radii[j]) return true;
          }
        }
      }
    }
  }
  return false;
}

std::optional<double> SampledStructure::freeRadius(const double3 &position) const
{
  double value = this->clearance(position);
  if (value < 0.0) return std::nullopt;
  return value;
}

// The clearance of 'position', walked shell by shell outward from its bin until no further shell can hold a
// nearer atom. With a 'floor' the walk ends as soon as an atom shows the clearance to be below it: the value
// then returned is below the floor (what it is exactly is not worked out), and a value at or above the floor
// is the exact clearance.
static double clearanceWalk(const SampledStructure &structure, const double3 &position, double floor)
{
  if (structure.positions.empty()) return std::numeric_limits<double>::max();

  const double maximumRadius = ensureNeighborGrid(structure);
  const SampledNeighborGrid &grid = structure.neighborGrid;
  const int3 gridSize = grid.gridSize;
  const double3x3 &cell = structure.unitCell.cell;

  const double3 fractional = structure.unitCell.inverseCell * position;
  const int3 pointBin = binOfUnwrappedFractional(fractional, gridSize);

  double smallest = std::numeric_limits<double>::max();
  for (int k = 0;; ++k)
  {
    double lowerBound = static_cast<double>(k - 1) * grid.minimumBinWidth;
    // Clearance to any atom in this shell is at least lowerBound - maxRadius.
    if (k > 0 && lowerBound - maximumRadius > smallest) break;
    if (smallest < floor) break;

    for (int ox = -k; ox <= k; ++ox)
    {
      for (int oy = -k; oy <= k; ++oy)
      {
        for (int oz = -k; oz <= k; ++oz)
        {
          if (std::max({std::abs(ox), std::abs(oy), std::abs(oz)}) != k) continue;

          auto [bx, lx] = binAndImage(pointBin.x + ox, gridSize.x);
          auto [by, ly] = binAndImage(pointBin.y + oy, gridSize.y);
          auto [bz, lz] = binAndImage(pointBin.z + oz, gridSize.z);

          double3 translation = cell * double3(static_cast<double>(lx), static_cast<double>(ly), static_cast<double>(lz));

          for (std::size_t j : grid.bins[static_cast<std::size_t>((bz * gridSize.y + by) * gridSize.x + bx)])
          {
            double3 delta = position - (grid.wrappedPositions[j] + translation);
            smallest = std::min(smallest, std::sqrt(double3::dot(delta, delta)) - structure.radii[j]);
          }
        }
      }
    }
  }
  return smallest;
}

double SampledStructure::clearance(const double3 &position) const
{
  return clearanceWalk(*this, position, std::numeric_limits<double>::lowest());
}

std::optional<double> SampledStructure::clearanceAtLeast(const double3 &position, double floor) const
{
  const double value = clearanceWalk(*this, position, floor);
  if (value < floor) return std::nullopt;
  return value;
}

SegmentBottleneck SampledStructure::segmentBottleneck(const double3 &position, const double3 &displacement) const
{
  SegmentBottleneck bottleneck{.radius = std::numeric_limits<double>::max(), .position = position};
  if (this->positions.empty()) return bottleneck;

  const double maximumRadius = ensureNeighborGrid(*this);
  const SampledNeighborGrid &grid = this->neighborGrid;
  const int3 gridSize = grid.gridSize;
  const double3x3 &cell = this->unitCell.cell;

  const double length = std::sqrt(double3::dot(displacement, displacement));
  const double length_squared = length * length;
  // Search from the midpoint: an atom farther than this cannot beat the current best clearance.
  const double3 midpoint = position + 0.5 * displacement;
  const double3 fractional = this->unitCell.inverseCell * midpoint;
  const int3 pointBin = binOfUnwrappedFractional(fractional, gridSize);

  for (int k = 0;; ++k)
  {
    double lowerBound = static_cast<double>(k - 1) * grid.minimumBinWidth;
    // dist(atom, segment) ≥ dist(atom, midpoint) - length/2, so clearance ≥ that − maxRadius.
    if (k > 0 && lowerBound - 0.5 * length - maximumRadius > bottleneck.radius) break;

    for (int ox = -k; ox <= k; ++ox)
    {
      for (int oy = -k; oy <= k; ++oy)
      {
        for (int oz = -k; oz <= k; ++oz)
        {
          if (std::max({std::abs(ox), std::abs(oy), std::abs(oz)}) != k) continue;

          auto [bx, lx] = binAndImage(pointBin.x + ox, gridSize.x);
          auto [by, ly] = binAndImage(pointBin.y + oy, gridSize.y);
          auto [bz, lz] = binAndImage(pointBin.z + oz, gridSize.z);

          double3 translation = cell * double3(static_cast<double>(lx), static_cast<double>(ly), static_cast<double>(lz));

          for (std::size_t j : grid.bins[static_cast<std::size_t>((bz * gridSize.y + by) * gridSize.x + bx)])
          {
            double3 atomImage = grid.wrappedPositions[j] + translation;
            double3 toAtom = atomImage - position;
            double t =
                length_squared > 0.0 ? std::clamp(double3::dot(toAtom, displacement) / length_squared, 0.0, 1.0)
                                     : 0.0;
            double3 closest = toAtom - t * displacement;
            double segment_clearance = std::sqrt(double3::dot(closest, closest)) - this->radii[j];
            if (segment_clearance < bottleneck.radius)
            {
              bottleneck.radius = segment_clearance;
              bottleneck.position = position + t * displacement;
            }
          }
        }
      }
    }
  }

  return bottleneck;
}

double SampledStructure::segmentClearance(const double3 &position, const double3 &displacement) const
{
  return this->segmentBottleneck(position, displacement).radius;
}

void SampledStructure::writeHeader(std::ostream &stream) const
{
  const SKSpaceGroupSetting& spaceGroup = SKSpaceGroupDataBase::setting(this->spaceGroupHallNumber);

  std::print(stream, "# Crystal: {}\n", this->name);
  std::print(stream, "# Space-group Hall-number: {}\n", this->spaceGroupHallNumber);
  std::print(stream, "# Space-group Hall-symbol: {}\n", spaceGroup.HallString());
  std::print(stream, "# Space-group HM-symbol: {}\n", spaceGroup.HMString());
  std::print(stream, "# Space-group IT number: {}\n", spaceGroup.number());
  std::print(stream, "# Number of framework atoms: {}\n", this->positions.size());
  std::print(stream, "# Crystal volume: {} [Å³]\n", this->unitCell.volume);
  std::print(stream, "# Crystal mass: {} [g/mol]\n", this->mass);
  std::print(stream, "# Crystal density: {} [kg/m³]\n", this->density());
}

void SampledProbe::writeHeader(std::ostream &stream) const
{
  std::print(stream, "# Probe atom: {} well-depth-factor: {} sigma: {}\n", this->name, this->wellDepthFactor,
             this->sizeParameter);
}
