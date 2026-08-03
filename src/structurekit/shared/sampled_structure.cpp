module;

module sampled_structure;

import std;

import double3;
import unit_cell;
import skspacegroupsetting;
import skspacegroupdatabase;
import units;

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
  for (std::size_t index = 0; index < this->positions.size(); ++index)
  {
    if (index == skip) continue;

    double3 dr = this->unitCell.applyPeriodicBoundaryConditions(position - this->positions[index]);
    if (double3::dot(dr, dr) < this->radii[index] * this->radii[index]) return true;
  }

  return false;
}

std::optional<double> SampledStructure::freeRadius(const double3 &position) const
{
  double smallest_radius = std::numeric_limits<double>::max();

  for (std::size_t index = 0; index < this->positions.size(); ++index)
  {
    double3 dr = this->unitCell.applyPeriodicBoundaryConditions(position - this->positions[index]);
    double rr = double3::dot(dr, dr);

    if (rr < this->radii[index] * this->radii[index]) return std::nullopt;

    smallest_radius = std::min(smallest_radius, std::sqrt(rr) - this->radii[index]);
  }

  return smallest_radius;
}

double SampledStructure::clearance(const double3 &position) const
{
  double smallest_radius = std::numeric_limits<double>::max();

  for (std::size_t index = 0; index < this->positions.size(); ++index)
  {
    double3 dr = this->unitCell.applyPeriodicBoundaryConditions(position - this->positions[index]);
    smallest_radius = std::min(smallest_radius, std::sqrt(double3::dot(dr, dr)) - this->radii[index]);
  }

  return smallest_radius;
}

SegmentBottleneck SampledStructure::segmentBottleneck(const double3 &position, const double3 &displacement) const
{
  double3 midpoint = position + 0.5 * displacement;
  double length_squared = double3::dot(displacement, displacement);

  SegmentBottleneck bottleneck{.radius = std::numeric_limits<double>::max(), .position = position};

  for (std::size_t index = 0; index < this->positions.size(); ++index)
  {
    // From the start of the segment to the nearest image of the atom, by way of the midpoint so that the
    // image chosen is the one nearest the segment as a whole rather than the one nearest an end of it.
    double3 dr = this->unitCell.applyPeriodicBoundaryConditions(this->positions[index] - midpoint) +
                 0.5 * displacement;

    // Where along the segment that atom comes closest, kept inside it.
    double t = length_squared > 0.0 ? std::clamp(double3::dot(dr, displacement) / length_squared, 0.0, 1.0) : 0.0;
    double3 closest = dr - t * displacement;

    double clearance = std::sqrt(double3::dot(closest, closest)) - this->radii[index];
    if (clearance < bottleneck.radius)
    {
      bottleneck.radius = clearance;
      bottleneck.position = position + t * displacement;
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
