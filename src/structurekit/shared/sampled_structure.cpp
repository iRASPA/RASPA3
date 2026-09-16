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
// How many cells out along each axis an overlapping image of a sphere of radius `maximumRadius` can
// still sit. Matches the brute-force structure's shell so the two agree on small cells.
int3 imageShell(const UnitCell &unitCell, double maximumRadius)
{
  const double3x3 &cell = unitCell.cell;
  double3 a(cell[0][0], cell[0][1], cell[0][2]);
  double3 b(cell[1][0], cell[1][1], cell[1][2]);
  double3 c(cell[2][0], cell[2][1], cell[2][2]);
  double spread = 0.5 * (a.length() + b.length() + c.length());
  double reach = 2.0 * spread + maximumRadius;

  double3 widths = unitCell.perpendicularWidths();
  auto along = [&](double width)
  { return static_cast<std::int32_t>(std::clamp(std::ceil(reach / std::max(width, 1.0e-9)), 1.0, 8.0)); };

  return int3(along(widths.x), along(widths.y), along(widths.z));
}

double maximumRadiusOf(const std::vector<double> &radii)
{
  return radii.empty() ? 0.0 : *std::ranges::max_element(radii);
}

// True when every sphere fits inside half the shortest lattice vector, so the ordinary minimum-image
// wrap cannot miss an overlapping copy.
bool minimumImageSufficient(const UnitCell &unitCell, double maximumRadius)
{
  const double3x3 &cell = unitCell.cell;
  double shortest = std::min({double3(cell[0][0], cell[0][1], cell[0][2]).length(),
                              double3(cell[1][0], cell[1][1], cell[1][2]).length(),
                              double3(cell[2][0], cell[2][1], cell[2][2]).length()});
  return 2.0 * maximumRadius <= shortest;
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
  const double maxRadius = maximumRadiusOf(this->radii);

  // Large cells: one minimum-image test per atom, plus the non-primary images of `skip` when that
  // atom's own sphere could reach across the cell.
  if (minimumImageSufficient(this->unitCell, maxRadius))
  {
    for (std::size_t index = 0; index < this->positions.size(); ++index)
    {
      if (index == skip) continue;
      double3 dr = this->unitCell.applyPeriodicBoundaryConditions(position - this->positions[index]);
      if (double3::dot(dr, dr) < this->radii[index] * this->radii[index]) return true;
    }
    return false;
  }

  const int3 shell = imageShell(this->unitCell, maxRadius);
  const double3x3 &cell = this->unitCell.cell;

  for (std::size_t index = 0; index < this->positions.size(); ++index)
  {
    const double radiusSquared = this->radii[index] * this->radii[index];
    for (std::int32_t nc = -shell.z; nc <= shell.z; ++nc)
    {
      for (std::int32_t nb = -shell.y; nb <= shell.y; ++nb)
      {
        for (std::int32_t na = -shell.x; na <= shell.x; ++na)
        {
          // A sample on atom `skip`'s sphere must ignore that atom's primary image, but not its
          // periodic copies: when 2r exceeds a cell edge the sphere overlaps itself across the boundary.
          if (index == skip && na == 0 && nb == 0 && nc == 0) continue;

          double3 translation =
              cell * double3(static_cast<double>(na), static_cast<double>(nb), static_cast<double>(nc));
          double3 dr = position - (this->positions[index] + translation);
          if (double3::dot(dr, dr) < radiusSquared) return true;
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
  const double maxRadius = maximumRadiusOf(this->radii);

  if (minimumImageSufficient(this->unitCell, maxRadius))
  {
    double smallest = std::numeric_limits<double>::max();
    for (std::size_t index = 0; index < this->positions.size(); ++index)
    {
      double3 dr = this->unitCell.applyPeriodicBoundaryConditions(position - this->positions[index]);
      smallest = std::min(smallest, std::sqrt(double3::dot(dr, dr)) - this->radii[index]);
    }
    return smallest;
  }

  const int3 shell = imageShell(this->unitCell, maxRadius);
  const double3x3 &cell = this->unitCell.cell;
  double smallest = std::numeric_limits<double>::max();

  for (std::size_t index = 0; index < this->positions.size(); ++index)
  {
    for (std::int32_t nc = -shell.z; nc <= shell.z; ++nc)
    {
      for (std::int32_t nb = -shell.y; nb <= shell.y; ++nb)
      {
        for (std::int32_t na = -shell.x; na <= shell.x; ++na)
        {
          double3 translation =
              cell * double3(static_cast<double>(na), static_cast<double>(nb), static_cast<double>(nc));
          double3 dr = position - (this->positions[index] + translation);
          smallest = std::min(smallest, std::sqrt(double3::dot(dr, dr)) - this->radii[index]);
        }
      }
    }
  }

  return smallest;
}

SegmentBottleneck SampledStructure::segmentBottleneck(const double3 &position, const double3 &displacement) const
{
  double length_squared = double3::dot(displacement, displacement);
  const double maxRadius = maximumRadiusOf(this->radii) + std::sqrt(length_squared);

  if (minimumImageSufficient(this->unitCell, maxRadius))
  {
    double3 midpoint = position + 0.5 * displacement;
    SegmentBottleneck bottleneck{.radius = std::numeric_limits<double>::max(), .position = position};

    for (std::size_t index = 0; index < this->positions.size(); ++index)
    {
      double3 dr = this->unitCell.applyPeriodicBoundaryConditions(this->positions[index] - midpoint) +
                   0.5 * displacement;
      double t =
          length_squared > 0.0 ? std::clamp(double3::dot(dr, displacement) / length_squared, 0.0, 1.0) : 0.0;
      double3 closest = dr - t * displacement;
      double segment_clearance = std::sqrt(double3::dot(closest, closest)) - this->radii[index];
      if (segment_clearance < bottleneck.radius)
      {
        bottleneck.radius = segment_clearance;
        bottleneck.position = position + t * displacement;
      }
    }
    return bottleneck;
  }

  const int3 shell = imageShell(this->unitCell, maxRadius);
  const double3x3 &cell = this->unitCell.cell;
  SegmentBottleneck bottleneck{.radius = std::numeric_limits<double>::max(), .position = position};

  for (std::size_t index = 0; index < this->positions.size(); ++index)
  {
    for (std::int32_t nc = -shell.z; nc <= shell.z; ++nc)
    {
      for (std::int32_t nb = -shell.y; nb <= shell.y; ++nb)
      {
        for (std::int32_t na = -shell.x; na <= shell.x; ++na)
        {
          double3 translation =
              cell * double3(static_cast<double>(na), static_cast<double>(nb), static_cast<double>(nc));
          double3 toAtom = this->positions[index] + translation - position;
          double t = length_squared > 0.0
                         ? std::clamp(double3::dot(toAtom, displacement) / length_squared, 0.0, 1.0)
                         : 0.0;
          double3 closest = toAtom - t * displacement;

          double segment_clearance = std::sqrt(double3::dot(closest, closest)) - this->radii[index];
          if (segment_clearance < bottleneck.radius)
          {
            bottleneck.radius = segment_clearance;
            bottleneck.position = position + t * displacement;
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
