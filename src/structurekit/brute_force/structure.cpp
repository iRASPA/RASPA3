module;

module brute_force_structure;

import std;

import int3;
import double3;
import double3x3;
import unit_cell;

BruteForceStructure BruteForceStructure::create(std::string name, const UnitCell &unitCell,
                                                std::vector<double3> positions, std::vector<double> radii)
{
  BruteForceStructure structure;
  structure.name = std::move(name);
  structure.unitCell = unitCell;
  structure.positions = std::move(positions);
  structure.radii = std::move(radii);
  structure.prepare();
  return structure;
}

void BruteForceStructure::prepare()
{
  const double3x3 &cell = this->unitCell.cell;

  this->fractionalPositions.clear();
  this->fractionalPositions.reserve(this->positions.size());
  for (const double3 &position : this->positions)
  {
    this->fractionalPositions.push_back(this->wrappedFractional(position));
  }

  this->maximumRadius = this->radii.empty() ? 0.0 : *std::ranges::max_element(this->radii);

  double3 a(cell[0][0], cell[0][1], cell[0][2]);
  double3 b(cell[1][0], cell[1][1], cell[1][2]);
  double3 c(cell[2][0], cell[2][1], cell[2][2]);

  // A difference of fractional coordinates brought into [-1/2, 1/2) along each axis is a Cartesian vector
  // no longer than this. It is the largest the ordinary nearest-image displacement can be, and so also a
  // bound on the distance to the nearest atom.
  this->spread = 0.5 * (a.length() + b.length() + c.length());

  // Which images can hold the nearest atom at all. The nearest one is at most `spread` away, so an image
  // translated by L puts every atom at least |L| - spread away and is worth searching only while that,
  // less the largest radius, is under `spread` as well.
  double reach = 2.0 * this->spread + this->maximumRadius;

  double3 widths = this->unitCell.perpendicularWidths();
  auto along = [&](double width)
  { return static_cast<std::int32_t>(std::clamp(std::ceil(reach / std::max(width, 1.0e-9)), 1.0, 8.0)); };

  this->shell = int3(along(widths.x), along(widths.y), along(widths.z));

  this->imageList.clear();
  for (std::int32_t nc = -this->shell.z; nc <= this->shell.z; ++nc)
  {
    for (std::int32_t nb = -this->shell.y; nb <= this->shell.y; ++nb)
    {
      for (std::int32_t na = -this->shell.x; na <= this->shell.x; ++na)
      {
        double3 translation =
            cell * double3(static_cast<double>(na), static_cast<double>(nb), static_cast<double>(nc));

        this->imageList.push_back(Image{.translation = translation, .distance = translation.length()});
      }
    }
  }

  std::ranges::sort(this->imageList, {}, &Image::distance);
}

double3 BruteForceStructure::wrappedFractional(const double3 &position) const
{
  double3 s = this->unitCell.inverseCell * position;
  return double3(s.x - std::floor(s.x), s.y - std::floor(s.y), s.z - std::floor(s.z));
}

// The displacement from `position` to atom `j`, with the fractional difference brought into [-1/2, 1/2)
// along each axis. Every image of that atom is this plus a lattice translation.
double3 BruteForceStructure::baseDisplacement(const double3 &fractional, std::size_t j) const
{
  double3 ds = fractional - this->fractionalPositions[j];
  ds = double3(ds.x - std::rint(ds.x), ds.y - std::rint(ds.y), ds.z - std::rint(ds.z));
  return this->unitCell.cell * ds;
}

double BruteForceStructure::clearance(const double3 &position) const
{
  double3 s = this->unitCell.inverseCell * position;

  double smallest = std::numeric_limits<double>::max();

  for (std::size_t j = 0; j < this->positions.size(); ++j)
  {
    double3 base = this->baseDisplacement(s, j);
    double reach = base.length() + this->radii[j];

    for (const Image &image : this->imageList)
    {
      // Every image from here on is at least this far off, and they are in order.
      if (image.distance - reach >= smallest) break;

      double3 dr = base + image.translation;
      smallest = std::min(smallest, dr.length() - this->radii[j]);
    }
  }

  return smallest;
}

bool BruteForceStructure::hasRoomFor(const double3 &position, double bound) const
{
  double3 s = this->unitCell.inverseCell * position;

  for (std::size_t j = 0; j < this->positions.size(); ++j)
  {
    double3 base = this->baseDisplacement(s, j);
    double reach = base.length() + this->radii[j];
    double limit = this->radii[j] + bound;

    for (const Image &image : this->imageList)
    {
      if (image.distance - reach >= bound) break;

      double3 dr = base + image.translation;
      if (double3::dot(dr, dr) < limit * limit) return false;
    }
  }

  return true;
}

double3 BruteForceStructure::nearestImage(const double3 &from, const double3 &to) const
{
  double3 ds = this->unitCell.inverseCell * (to - from);
  ds = double3(ds.x - std::rint(ds.x), ds.y - std::rint(ds.y), ds.z - std::rint(ds.z));

  double3 base = this->unitCell.cell * ds;
  double reach = base.length();

  double3 best = base;
  double bestSquared = double3::dot(best, best);

  for (const Image &image : this->imageList)
  {
    if (image.distance - reach >= std::sqrt(bestSquared)) break;

    double3 dr = base + image.translation;
    double squared = double3::dot(dr, dr);
    if (squared < bestSquared)
    {
      bestSquared = squared;
      best = dr;
    }
  }

  return best;
}

double BruteForceStructure::segmentClearance(const double3 &origin, const double3 &displacement) const
{
  double3 s = this->unitCell.inverseCell * origin;

  double lengthSquared = double3::dot(displacement, displacement);
  double length = std::sqrt(lengthSquared);

  double smallest = std::numeric_limits<double>::max();

  for (std::size_t j = 0; j < this->positions.size(); ++j)
  {
    double3 base = this->baseDisplacement(s, j);

    // The segment reaches `length` beyond its start, so an image is out of reach only by that much more.
    double reach = base.length() + this->radii[j] + length;

    for (const Image &image : this->imageList)
    {
      if (image.distance - reach >= smallest) break;

      // The atom is at `-(base + translation)` relative to the start of the segment.
      double3 toAtom = -1.0 * (base + image.translation);

      double t = lengthSquared > 0.0
                     ? std::clamp(double3::dot(toAtom, displacement) / lengthSquared, 0.0, 1.0)
                     : 0.0;
      double3 closest = toAtom - t * displacement;

      smallest = std::min(smallest, closest.length() - this->radii[j]);
    }
  }

  return smallest;
}

std::optional<double> BruteForceStructure::firstExit(const double3 &origin, const double3 &direction,
                                                     double limit) const
{
  double3 s = this->unitCell.inverseCell * origin;

  // Every atom the ray passes through gives the interval of the ray that is inside it. Sorted and merged,
  // the first gap is the first place the ray is outside everything.
  std::vector<std::pair<double, double>> intervals;

  for (std::size_t j = 0; j < this->positions.size(); ++j)
  {
    double3 base = this->baseDisplacement(s, j);
    double reach = base.length() + this->radii[j] + limit;

    for (const Image &image : this->imageList)
    {
      if (image.distance >= reach) break;

      double3 toOrigin = -1.0 * (base + image.translation);  // from the atom to the ray's origin

      // |toOrigin + t d|² = r², solved for t.
      double b = double3::dot(toOrigin, direction);
      double c = double3::dot(toOrigin, toOrigin) - this->radii[j] * this->radii[j];

      double discriminant = b * b - c;
      if (discriminant <= 0.0) continue;

      double root = std::sqrt(discriminant);
      double enter = -b - root;
      double leave = -b + root;

      if (leave <= 0.0 || enter >= limit) continue;

      intervals.emplace_back(std::max(enter, 0.0), std::min(leave, limit));
    }
  }

  if (intervals.empty()) return 0.0;

  std::ranges::sort(intervals);

  double reached = 0.0;
  for (const auto &[enter, leave] : intervals)
  {
    if (enter > reached) return reached;
    reached = std::max(reached, leave);
    if (reached >= limit) return std::nullopt;
  }

  return reached < limit ? std::optional<double>{reached} : std::nullopt;
}
