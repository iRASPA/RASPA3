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

  // Which images can be in contact with an atom, which is a question the same list has to answer for whatever
  // walks the pairs: an atom's nearest image of another is within half a cell along each axis or a little more on an
  // oblique cell, and from there the contact reaches on by the two radii. That is not limited the way the search for
  // the nearest atom is. A cell thin along one axis has an image in contact every step along it for as far as the
  // radii reach, and it can have plenty of void around all of them.
  const double contact = 2.0 * this->maximumRadius;

  double3 widths = this->unitCell.perpendicularWidths();
  auto along = [&](double width)
  {
    const double safeWidth = std::max(width, 1.0e-9);
    return static_cast<std::int32_t>(std::max(std::clamp(std::ceil(reach / safeWidth), 1.0, 8.0),
                                              std::ceil(contact / safeWidth) + 1.0));
  };

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

// At each turn the twenty-six lattice directions about the current point are tried at the current step, and
// so are the directions the nearest atoms single out. The lattice directions alone are not enough. The
// clearance is the least of the distances to the atoms, so where two or three of them are equally near it
// rises only along the crease between them: for two, the plane they are equidistant from, and for three, the
// line. A crease that lies between the lattice directions is one no lattice step climbs, since every such
// step brings one of the atoms nearer faster than it takes the walk along the crease, and the walk stops on
// the crease short of the vertex at its top, reporting a sphere smaller than the one that fits. So the
// bisector of each pair of the nearest atoms and the line along which each triple of them stays equidistant
// are tried as well, each with the sign that moves away from the atoms; along those the clearance rises until
// the next atom joins in, and at a vertex, where four do, nothing rises and the walk is done. Deterministic,
// so that two runs of the check agree with each other as well as with what they are checking.
std::pair<double3, double> BruteForceStructure::walkUphill(const double3 &from, double startingStep) const
{
  double3 best = from;
  double bestClearance = this->clearance(best);
  double step = startingStep;

  std::vector<double3> directions;
  std::vector<std::pair<double, double3>> nearest;  // the room past an atom, and the unit vector away from it

  while (step > 1.0e-10)
  {
    directions.clear();
    for (std::int32_t k = -1; k <= 1; ++k)
    {
      for (std::int32_t j = -1; j <= 1; ++j)
      {
        for (std::int32_t i = -1; i <= 1; ++i)
        {
          if (i == 0 && j == 0 && k == 0) continue;
          double3 direction(static_cast<double>(i), static_cast<double>(j), static_cast<double>(k));
          directions.push_back((1.0 / direction.length()) * direction);
        }
      }
    }

    // The atoms within a step of being the nearest, the four nearest of them.
    nearest.clear();
    for (std::size_t j = 0; j < this->positions.size(); ++j)
    {
      double3 away = -1.0 * this->nearestImage(best, this->positions[j]);
      double length = away.length();
      if (length <= 0.0) continue;
      double room = length - this->radii[j];
      if (room > bestClearance + step) continue;
      nearest.emplace_back(room, (1.0 / length) * away);
    }
    std::ranges::sort(nearest, [](const auto &a, const auto &b) { return a.first < b.first; });
    if (nearest.size() > 4) nearest.resize(4);

    for (std::size_t a = 0; a < nearest.size(); ++a)
    {
      for (std::size_t b = a + 1; b < nearest.size(); ++b)
      {
        double3 bisector = nearest[a].second + nearest[b].second;
        double bisectorLength = bisector.length();
        if (bisectorLength > 1.0e-12) directions.push_back((1.0 / bisectorLength) * bisector);

        for (std::size_t c = b + 1; c < nearest.size(); ++c)
        {
          double3 crease = double3::cross(nearest[a].second - nearest[b].second, nearest[a].second - nearest[c].second);
          double creaseLength = crease.length();
          if (creaseLength <= 1.0e-12) continue;
          crease = (1.0 / creaseLength) * crease;
          if (double3::dot(crease, nearest[a].second) < 0.0) crease = -1.0 * crease;
          directions.push_back(crease);
        }
      }
    }

    bool moved = false;
    for (const double3 &direction : directions)
    {
      double3 trial = best + step * direction;
      double clearance = this->clearance(trial);
      if (clearance > bestClearance)
      {
        bestClearance = clearance;
        best = trial;
        moved = true;
      }
    }
    if (!moved) step *= 0.5;
  }

  return {best, bestClearance};
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
