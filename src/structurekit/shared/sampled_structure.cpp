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

void ensureNeighborGrid(const SampledStructure &structure)
{
  SampledNeighborGrid &grid = structure.neighborGrid;
  const double maxRadius = maximumRadiusOf(structure.radii);
  if (grid.atomCount == structure.positions.size() && grid.maximumRadius == maxRadius && !grid.bins.empty())
  {
    return;
  }

  grid.atomCount = structure.positions.size();
  grid.maximumRadius = maxRadius;

  const UnitCell &unitCell = structure.unitCell;
  const double3x3 &cell = unitCell.cell;
  const double volume = unitCell.volume;
  double3 a(cell[0][0], cell[0][1], cell[0][2]);
  double3 b(cell[1][0], cell[1][1], cell[1][2]);
  double3 c(cell[2][0], cell[2][1], cell[2][2]);
  double3 perpendicularWidths(volume / double3::cross(b, c).length(), volume / double3::cross(c, a).length(),
                              volume / double3::cross(a, b).length());

  // Roughly four atoms per bin; bin counts proportional to the perpendicular widths so the bins are
  // approximately metrically cubic (same construction as PoreAccessibility / SKVoronoi).
  const double targetBinSize =
      std::cbrt(volume / std::max(1.0, static_cast<double>(grid.atomCount) / 4.0));
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

  ensureNeighborGrid(*this);
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
    if (lowerBound > grid.maximumRadius) break;

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

double SampledStructure::clearance(const double3 &position) const
{
  if (this->positions.empty()) return std::numeric_limits<double>::max();

  ensureNeighborGrid(*this);
  const SampledNeighborGrid &grid = this->neighborGrid;
  const int3 gridSize = grid.gridSize;
  const double3x3 &cell = this->unitCell.cell;

  const double3 fractional = this->unitCell.inverseCell * position;
  const int3 pointBin = binOfUnwrappedFractional(fractional, gridSize);

  double smallest = std::numeric_limits<double>::max();
  for (int k = 0;; ++k)
  {
    double lowerBound = static_cast<double>(k - 1) * grid.minimumBinWidth;
    // Clearance to any atom in this shell is at least lowerBound - maxRadius.
    if (k > 0 && lowerBound - grid.maximumRadius > smallest) break;

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
            smallest = std::min(smallest, std::sqrt(double3::dot(delta, delta)) - this->radii[j]);
          }
        }
      }
    }
  }
  return smallest;
}

SegmentBottleneck SampledStructure::segmentBottleneck(const double3 &position, const double3 &displacement) const
{
  SegmentBottleneck bottleneck{.radius = std::numeric_limits<double>::max(), .position = position};
  if (this->positions.empty()) return bottleneck;

  ensureNeighborGrid(*this);
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
    if (k > 0 && lowerBound - 0.5 * length - grid.maximumRadius > bottleneck.radius) break;

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
  const SKSpaceGroupSetting &spaceGroup = SKSpaceGroupDataBase::spaceGroupData[this->spaceGroupHallNumber];

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
