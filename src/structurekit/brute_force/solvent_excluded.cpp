module;

module brute_force_solvent_excluded;

import std;

import int3;
import double3;
import double3x3;
import unit_cell;
import randomnumbers;
import brute_force_structure;
import brute_force_voxels;
import brute_force_surface_area;

namespace
{
// The sampled accessible surface, in bins, so that "is there a probe centre nearer than the probe's radius"
// is a question about a handful of neighbouring bins rather than about every point of the cloud.
class ProbeCentreCloud
{
 public:
  ProbeCentreCloud(const BruteForceStructure &structure, const std::vector<double3> &points, double reach)
      : unitCell(structure.unitCell)
  {
    double3 widths = this->unitCell.perpendicularWidths();

    auto along = [&](double width)
    { return static_cast<std::int32_t>(std::clamp(std::floor(width / std::max(reach, 1.0e-6)), 1.0, 256.0)); };
    this->counts = int3(along(widths.x), along(widths.y), along(widths.z));

    std::size_t numberOfBins = static_cast<std::size_t>(this->counts.x) * static_cast<std::size_t>(this->counts.y) *
                               static_cast<std::size_t>(this->counts.z);

    this->wrapped.reserve(points.size());
    std::vector<std::size_t> bin;
    bin.reserve(points.size());

    std::vector<std::size_t> tally(numberOfBins, 0);
    for (const double3 &point : points)
    {
      double3 s = structure.wrappedFractional(point);
      this->wrapped.push_back(this->unitCell.cell * s);

      std::size_t index = this->binOf(s);
      bin.push_back(index);
      ++tally[index];
    }

    this->start.assign(numberOfBins + 1, 0);
    std::partial_sum(tally.begin(), tally.end(), this->start.begin() + 1);

    std::vector<std::size_t> cursor(this->start.begin(), this->start.end() - 1);
    this->contents.resize(points.size());
    for (std::size_t i = 0; i < points.size(); ++i) this->contents[cursor[bin[i]]++] = i;
  }

  // Whether some sampled probe centre is nearer to `position` than `distance`. A point of a generated patch
  // with one of those in reach is inside the space the probe can occupy, and so is not on its surface.
  bool anyWithin(const BruteForceStructure &structure, const double3 &position, double distance) const
  {
    double3 s = structure.wrappedFractional(position);
    double3 home = this->unitCell.cell * s;

    auto axis = [](double coordinate, std::int32_t count)
    {
      return std::clamp(static_cast<std::int32_t>(coordinate * static_cast<double>(count)), std::int32_t{0},
                        count - 1);
    };
    int3 centre(axis(s.x, this->counts.x), axis(s.y, this->counts.y), axis(s.z, this->counts.z));

    for (std::int32_t dk = -1; dk <= 1; ++dk)
    {
      for (std::int32_t dj = -1; dj <= 1; ++dj)
      {
        for (std::int32_t di = -1; di <= 1; ++di)
        {
          std::int32_t i = ((centre.x + di) % this->counts.x + this->counts.x) % this->counts.x;
          std::int32_t j = ((centre.y + dj) % this->counts.y + this->counts.y) % this->counts.y;
          std::int32_t k = ((centre.z + dk) % this->counts.z + this->counts.z) % this->counts.z;

          std::size_t index =
              (static_cast<std::size_t>(k) * static_cast<std::size_t>(this->counts.y) + static_cast<std::size_t>(j)) *
                  static_cast<std::size_t>(this->counts.x) +
              static_cast<std::size_t>(i);

          for (std::size_t slot = this->start[index]; slot < this->start[index + 1]; ++slot)
          {
            double3 dr = structure.nearestImage(home, this->wrapped[this->contents[slot]]);
            if (double3::dot(dr, dr) < distance * distance) return true;
          }
        }
      }
    }

    return false;
  }

 private:
  std::size_t binOf(const double3 &fractional) const
  {
    auto axis = [](double coordinate, std::int32_t count)
    {
      return std::clamp(static_cast<std::int32_t>(coordinate * static_cast<double>(count)), std::int32_t{0},
                        count - 1);
    };

    std::int32_t i = axis(fractional.x, this->counts.x);
    std::int32_t j = axis(fractional.y, this->counts.y);
    std::int32_t k = axis(fractional.z, this->counts.z);

    return (static_cast<std::size_t>(k) * static_cast<std::size_t>(this->counts.y) + static_cast<std::size_t>(j)) *
               static_cast<std::size_t>(this->counts.x) +
           static_cast<std::size_t>(i);
  }

  UnitCell unitCell;
  int3 counts{1, 1, 1};
  std::vector<double3> wrapped;
  std::vector<std::size_t> start;
  std::vector<std::size_t> contents;
};

// Where the corners found so far are, in bins over the cell, so that "have we been here already" is a
// question about a handful of neighbouring bins rather than about every corner found.
//
// It has to be a question about distance and not about a rounded coordinate. A corner is reached once from
// each atom it rests on and once for each triple of them, and the arithmetic that places it differs in the
// last bit or two between those routes; rounding the position and comparing the rounded value keeps a
// corner twice whenever the two routes fall either side of a rounding boundary. A crystal is exactly where
// that happens, because it puts its corners on the symmetry planes, and rounding boundaries lie on planes
// too. A corner kept twice is not a small error either: two probe centres a hair apart do not bury each
// other, so the patch there is laid down twice over.
class CornerBins
{
 public:
  CornerBins(const BruteForceStructure &structure, double binSize)
  {
    double3 widths = structure.unitCell.perpendicularWidths();

    auto along = [&](double width)
    { return static_cast<std::int32_t>(std::clamp(std::floor(width / binSize), 1.0, 128.0)); };
    this->counts = int3(along(widths.x), along(widths.y), along(widths.z));

    this->contents.resize(static_cast<std::size_t>(this->counts.x) * static_cast<std::size_t>(this->counts.y) *
                          static_cast<std::size_t>(this->counts.z));
  }

  // Whether `centre` is a corner not seen before, remembering it when it is.
  bool addIfNew(const BruteForceStructure &structure, const double3 &centre, double tolerance)
  {
    double3 s = structure.wrappedFractional(centre);
    int3 home = this->cellOf(s);

    for (std::int32_t dk = -1; dk <= 1; ++dk)
    {
      for (std::int32_t dj = -1; dj <= 1; ++dj)
      {
        for (std::int32_t di = -1; di <= 1; ++di)
        {
          for (const double3 &kept : this->contents[this->indexOf(home.x + di, home.y + dj, home.z + dk)])
          {
            double3 dr = structure.nearestImage(centre, kept);
            if (double3::dot(dr, dr) < tolerance * tolerance) return false;
          }
        }
      }
    }

    this->contents[this->indexOf(home.x, home.y, home.z)].push_back(centre);
    return true;
  }

 private:
  int3 cellOf(const double3 &fractional) const
  {
    auto axis = [](double coordinate, std::int32_t count)
    {
      return std::clamp(static_cast<std::int32_t>(coordinate * static_cast<double>(count)), std::int32_t{0},
                        count - 1);
    };
    return int3(axis(fractional.x, this->counts.x), axis(fractional.y, this->counts.y),
                axis(fractional.z, this->counts.z));
  }

  std::size_t indexOf(std::int32_t i, std::int32_t j, std::int32_t k) const
  {
    i = (i % this->counts.x + this->counts.x) % this->counts.x;
    j = (j % this->counts.y + this->counts.y) % this->counts.y;
    k = (k % this->counts.z + this->counts.z) % this->counts.z;

    return (static_cast<std::size_t>(k) * static_cast<std::size_t>(this->counts.y) +
            static_cast<std::size_t>(j)) *
               static_cast<std::size_t>(this->counts.x) +
           static_cast<std::size_t>(i);
  }

  int3 counts{1, 1, 1};
  std::vector<std::vector<double3>> contents;
};

// An atom of some image, as the pair and triple sweeps see it.
struct Neighbour
{
  std::size_t atom;
  double3 centre;   // Cartesian, the image's own position
  double radius;    // Å, inflated
  int3 translation;
};

// Whichever way round, one of the two vectors perpendicular to `axis`.
double3 somePerpendicular(const double3 &axis)
{
  double3 trial = std::abs(axis.x) < 0.9 ? double3(1.0, 0.0, 0.0) : double3(0.0, 1.0, 0.0);
  double3 out = trial - double3::dot(trial, axis) * axis;
  return (1.0 / out.length()) * out;
}

// How much surface one turn of the probe about a crease lays down, per unit angle turned: the integral of
// the area element across the arc from one contact to the other,
//
//   int (rho - r cos theta) r dtheta  =  r (rho theta - r sin theta) ,
//
// which is worth doing in closed form rather than in steps. What the steps were costing was not the
// smooth part of the sweep but the ends of it, where the probe stops being free partway through a step,
// and that is handled by finding the ends themselves instead.
double sweptAcrossArc(double circleRadius, double probeRadius, double from, double to)
{
  auto integral = [&](double theta)
  { return probeRadius * (circleRadius * theta - probeRadius * std::sin(theta)); };

  if (to < from) std::swap(from, to);

  // A circle of probe centres narrower than the probe folds the swept surface through the axis. The
  // surface stops where it reaches the axis, which is a gap in the middle of the arc.
  if (circleRadius >= probeRadius) return integral(to) - integral(from);

  double cusp = std::acos(std::clamp(circleRadius / probeRadius, -1.0, 1.0));

  double swept = 0.0;
  if (from < -cusp) swept += integral(std::min(to, -cusp)) - integral(from);
  if (to > cusp) swept += integral(to) - integral(std::max(from, cusp));
  return swept;
}
}  // namespace

BruteForceSolventExcluded BruteForceSolventExcluded::compute(const BruteForceStructure &structure,
                                                             const BruteForceStructure &inflated,
                                                             const BruteForceSurfaceArea &surface,
                                                             double probeRadius, std::size_t samplesPerAtom,
                                                             std::size_t creaseSteps, std::size_t cornerSamples)
{
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  BruteForceSolventExcluded self;
  self.probeRadius = probeRadius;

  std::size_t numberOfAtoms = structure.size();
  const double3x3 &cell = structure.unitCell.cell;

  // A probe centre is free when it is clear of every inflated atom. Generated centres sit exactly on one
  // inflated sphere, so the test is against rounding rather than against zero.
  auto freeProbeCentre = [&](const double3 &centre) { return inflated.clearance(centre) >= -1.0e-9; };

  // Who is near enough to whom for a probe to touch both, and for three of them to meet. Everything that
  // follows is local, so this is what keeps the sweeps below from being over every pair and triple there
  // is.
  std::vector<std::vector<Neighbour>> neighboursOf(numberOfAtoms);
  if (probeRadius > 0.0)
  {
    for (std::size_t atom = 0; atom < numberOfAtoms; ++atom)
    {
      double3 home = cell * structure.fractional()[atom];

      for (std::size_t other = 0; other < numberOfAtoms; ++other)
      {
        double3 elsewhere = cell * structure.fractional()[other];

        // Every image of the other atom is its nearest one plus a lattice translation, and the nearest is
        // no further off than that translation is long, so the walk stops where it can no longer reach.
        double3 base = structure.nearestImage(home, elsewhere);
        double reach = inflated.radii[atom] + inflated.radii[other];

        for (const BruteForceStructure::Image &image : inflated.images())
        {
          if (image.distance - base.length() >= reach) break;

          double3 dr = base + image.translation;

          double distanceSquared = double3::dot(dr, dr);
          if (distanceSquared >= reach * reach || distanceSquared < 1.0e-12) continue;

          double3 centre = home + dr;

          double3 s = structure.unitCell.inverseCell * (centre - elsewhere);
          int3 translation(static_cast<std::int32_t>(std::lround(s.x)),
                           static_cast<std::int32_t>(std::lround(s.y)),
                           static_cast<std::int32_t>(std::lround(s.z)));

          neighboursOf[atom].push_back(Neighbour{
              .atom = other, .centre = centre, .radius = inflated.radii[other], .translation = translation});
        }
      }
    }
  }

  // A corner and every direction in which the probe there is touching something. Three of them come from
  // the triple that found it; a fourth atom arriving at exactly the same distance is a corner the probe
  // rests on more heavily, and the patch there is bounded by all of them and not by three.
  struct Corner
  {
    double3 centre;
    std::vector<double3> contact;
  };

  std::vector<Corner> corners;
  if (probeRadius > 0.0)
  {
    CornerBins seen(structure, 0.5);

    for (std::size_t atom = 0; atom < numberOfAtoms; ++atom)
    {
      double3 first = cell * structure.fractional()[atom];
      double firstRadius = inflated.radii[atom];

      const std::vector<Neighbour> &around = neighboursOf[atom];

      for (std::size_t a = 0; a < around.size(); ++a)
      {
        for (std::size_t b = a + 1; b < around.size(); ++b)
        {
          const Neighbour &second = around[a];
          const Neighbour &third = around[b];

          // Trilateration in the frame of the three centres: the probe centre lies where three spheres meet.
          double3 ex = second.centre - first;
          double distance = ex.length();
          if (distance <= 0.0) continue;
          ex = (1.0 / distance) * ex;

          double3 toThird = third.centre - first;
          double i = double3::dot(ex, toThird);

          double3 ey = toThird - i * ex;
          double eyLength = ey.length();
          if (eyLength <= 1.0e-9) continue;  // three centres in a line have no corner
          ey = (1.0 / eyLength) * ey;

          double3 ez = double3::cross(ex, ey);
          double j = double3::dot(ey, toThird);

          double x = (firstRadius * firstRadius - second.radius * second.radius + distance * distance) /
                     (2.0 * distance);
          double y = (firstRadius * firstRadius - third.radius * third.radius + i * i + j * j - 2.0 * i * x) /
                     (2.0 * j);

          double zSquared = firstRadius * firstRadius - x * x - y * y;
          if (zSquared <= 0.0) continue;

          double z = std::sqrt(zSquared);

          for (double sign : {1.0, -1.0})
          {
            double3 centre = first + x * ex + y * ey + (sign * z) * ez;

            if (!freeProbeCentre(centre)) continue;

            // The same corner is reached from each of its atoms and from each triple of them, so it is kept
            // once, by where it is rather than by any rounding of where it is.
            if (!seen.addIfNew(structure, centre, 1.0e-5)) continue;

            // Everything the probe is touching there, not only the three that placed it. Missing a fourth
            // would leave the patch bounded by too few great circles and so too large.
            Corner corner;
            corner.centre = centre;

            for (const Neighbour &touching : around)
            {
              double3 to = touching.centre - centre;
              if (std::abs(to.length() - touching.radius) > 1.0e-7) continue;
              corner.contact.push_back((1.0 / touching.radius) * to);
            }
            {
              double3 to = first - centre;
              if (std::abs(to.length() - firstRadius) <= 1.0e-7)
              {
                corner.contact.push_back((1.0 / firstRadius) * to);
              }
            }

            if (corner.contact.size() < 3) continue;
            corners.push_back(std::move(corner));
          }
        }
      }
    }
  }

  self.numberOfCorners = corners.size();
  for (const Corner &corner : corners)
  {
    std::size_t touching = std::min(corner.contact.size(), self.cornersByContacts.size() - 1);
    ++self.cornersByContacts[touching];
  }

  // How near two corners come to each other. A wedge the arithmetic reached by two nearly equal routes
  // would appear twice and lay down its patch twice, so this is measured rather than assumed away.
  {
    double closest = std::numeric_limits<double>::max();
    std::size_t crowded = 0;

#pragma omp parallel for schedule(static) reduction(min : closest) reduction(+ : crowded)
    for (std::int64_t index = 0; index < static_cast<std::int64_t>(corners.size()); ++index)
    {
      double nearest = std::numeric_limits<double>::max();

      for (std::size_t other = 0; other < corners.size(); ++other)
      {
        if (other == static_cast<std::size_t>(index)) continue;

        double3 dr = structure.nearestImage(corners[static_cast<std::size_t>(index)].centre,
                                            corners[other].centre);
        nearest = std::min(nearest, dr.length());
      }

      closest = std::min(closest, nearest);
      if (nearest < 1.0e-3) ++crowded;
    }

    self.crowdedCorners = crowded;
    self.closestCorners = corners.size() > 1 ? closest : 0.0;
  }

  // Where the probe can be, as a cloud of points to hold every generated patch against.
  //
  // Two kinds of place go into it. The first is the accessible surface, sampled already by the surface-area
  // check, which is where the probe sits when it is touching the framework and rolling. The second is the
  // corners, which is where it sits when it is wedged and cannot roll at all -- and those have to be put in
  // by hand, because a corner is a single point of that surface and no sampling of an area will ever land
  // on one. Leaving them out is not a small omission: a concave patch is cut back mostly by the probes at
  // the corners next to it, so a cloud without them finds almost none of the burial there is.
  std::vector<double3> reachable = surface.exposedPoints;
  reachable.reserve(reachable.size() + corners.size());
  for (const Corner &corner : corners) reachable.push_back(corner.centre);

  ProbeCentreCloud cloud(structure, reachable, probeRadius);

  // A generated point is buried when a probe can be put over it. The comparison is made a whisker inside
  // the probe's radius: the centre that generates the point is itself at exactly that distance, and is not
  // evidence that the point is covered.
  const double burialDistance = probeRadius - 1.0e-7;
  auto buried = [&](const double3 &point)
  { return probeRadius > 0.0 && cloud.anyWithin(structure, point, burialDistance); };

  // The convex patches: each atom's own surface, wherever a probe resting against it there is free.
  {
    std::vector<double> convexOfAtom(numberOfAtoms, 0.0);
    std::vector<double> buriedOfAtom(numberOfAtoms, 0.0);
    std::vector<double> varianceOfAtom(numberOfAtoms, 0.0);

#pragma omp parallel for schedule(dynamic, 1)
    for (std::int64_t index = 0; index < static_cast<std::int64_t>(numberOfAtoms); ++index)
    {
      std::size_t atom = static_cast<std::size_t>(index);

      double radius = structure.radii[atom];
      double sphereArea = 4.0 * std::numbers::pi * radius * radius;
      double perSample = sphereArea / static_cast<double>(samplesPerAtom);

      RandomNumber random{numberOfAtoms + atom};

      std::size_t kept = 0;
      std::size_t lost = 0;

      for (std::size_t sample = 0; sample < samplesPerAtom; ++sample)
      {
        double3 direction = random.randomVectorOnUnitSphere();

        double3 centre = structure.positions[atom] + (radius + probeRadius) * direction;
        if (!freeProbeCentre(centre)) continue;

        double3 point = structure.positions[atom] + radius * direction;
        if (buried(point))
        {
          ++lost;
          continue;
        }

        ++kept;
      }

      convexOfAtom[atom] = static_cast<double>(kept) * perSample;
      buriedOfAtom[atom] = static_cast<double>(lost) * perSample;

      double fraction = static_cast<double>(kept) / static_cast<double>(samplesPerAtom);
      varianceOfAtom[atom] =
          sphereArea * sphereArea * fraction * (1.0 - fraction) / static_cast<double>(samplesPerAtom);
    }

    double variance = 0.0;
    for (std::size_t atom = 0; atom < numberOfAtoms; ++atom)
    {
      self.convexArea += convexOfAtom[atom];
      self.buriedConvexArea += buriedOfAtom[atom];
      variance += varianceOfAtom[atom];
    }
    self.convexAreaError = std::sqrt(variance);
  }

  if (probeRadius <= 0.0)
  {
    // With no probe there is nothing to roll: the surface is the atoms' own and has no other kind of patch.
    std::chrono::duration<double> timing = std::chrono::steady_clock::now() - begin;
    self.seconds = timing.count();
    return self;
  }

  // The saddle patches: a probe rolling along a pair of atoms sweeps part of a torus about the line between
  // them. Each unordered pair is taken once, by fixing the first atom in the home cell and ordering what
  // may stand as the second.
  {
    const std::size_t aroundTheCircle = std::max<std::size_t>(16, creaseSteps);
    const std::size_t acrossTheArc = std::max<std::size_t>(8, creaseSteps / 4);

    std::vector<double> saddleOfAtom(numberOfAtoms, 0.0);
    std::vector<double> buriedOfAtom(numberOfAtoms, 0.0);
    std::vector<std::size_t> pairsOfAtom(numberOfAtoms, 0);

#pragma omp parallel for schedule(dynamic, 1)
    for (std::int64_t index = 0; index < static_cast<std::int64_t>(numberOfAtoms); ++index)
    {
      std::size_t atom = static_cast<std::size_t>(index);

      double3 first = cell * structure.fractional()[atom];
      double firstInflated = inflated.radii[atom];

      for (const Neighbour &neighbour : neighboursOf[atom])
      {
        // Once per pair: the higher-numbered atom, or the same atom in a later image.
        bool ordered = neighbour.atom > atom;
        if (neighbour.atom == atom)
        {
          const int3 &t = neighbour.translation;
          ordered = t.x > 0 || (t.x == 0 && (t.y > 0 || (t.y == 0 && t.z > 0)));
        }
        if (!ordered) continue;

        double3 along = neighbour.centre - first;
        double distance = along.length();
        if (distance <= 0.0) continue;

        double3 axis = (1.0 / distance) * along;

        // Where along the line the circle of probe centres sits, and how wide it is.
        double offset =
            (distance * distance + firstInflated * firstInflated - neighbour.radius * neighbour.radius) /
            (2.0 * distance);
        double squared = firstInflated * firstInflated - offset * offset;
        if (squared <= 0.0) continue;

        double circleRadius = std::sqrt(squared);
        double3 circleCentre = first + offset * axis;

        double3 u = somePerpendicular(axis);
        double3 v = double3::cross(axis, u);

        ++pairsOfAtom[atom];

        // The arc of the probe's own surface between the two contacts, in the plane through the axis: the
        // angle is measured from the direction straight in towards the axis, so that the point at angle
        // theta stands off the axis by circleRadius - probeRadius cos theta.
        double toFirst = std::atan2(-offset, circleRadius);
        double toSecond = std::atan2(distance - offset, circleRadius);

        // Where the swept surface would fold through the axis, it stops instead. A circle of probe centres
        // narrower than the probe leaves a gap about the middle of the arc.
        double cusp = circleRadius < probeRadius ? std::acos(std::clamp(circleRadius / probeRadius, -1.0, 1.0))
                                                 : 0.0;

        double perTurn = sweptAcrossArc(circleRadius, probeRadius, toFirst, toSecond);
        if (perTurn <= 0.0) continue;

        double stepPhi = 2.0 * std::numbers::pi / static_cast<double>(aroundTheCircle);

        auto centreAt = [&](double phi)
        { return circleCentre + circleRadius * (std::cos(phi) * u + std::sin(phi) * v); };

        // Where round the circle the probe can go is a set of arcs, and the ends of them are where the
        // probe first meets a third atom. Stepping round and taking whole steps as free or not would put
        // each end wrong by up to half a step, which on a crystal is a systematic error rather than a
        // wobble, the ends of the crease sitting at the same place on every symmetry-related pair. So the
        // ends are found instead, by halving the interval the change happened in.
        auto boundaryBetween = [&](double open, double shut)
        {
          for (int refine = 0; refine < 40; ++refine)
          {
            double middle = 0.5 * (open + shut);
            (freeProbeCentre(centreAt(middle)) ? open : shut) = middle;
          }
          return 0.5 * (open + shut);
        };

        // How much of the turn the probe is free for, taken interval by interval so that the walk cannot
        // wrap past itself or count a stretch twice. Each step contributes all of itself, none of itself,
        // or the part of itself on the free side of the crossing found within it.
        double turned = 0.0;

        bool freeHere = freeProbeCentre(centreAt(0.0));
        bool freeAtStart = freeHere;

        for (std::size_t step = 0; step < aroundTheCircle; ++step)
        {
          double from = static_cast<double>(step) * stepPhi;
          double to = from + stepPhi;

          bool freeThere = (step + 1 == aroundTheCircle) ? freeAtStart : freeProbeCentre(centreAt(to));

          if (freeHere && freeThere)
          {
            turned += stepPhi;
          }
          else if (freeHere != freeThere)
          {
            double crossing = boundaryBetween(freeHere ? from : to, freeHere ? to : from);
            turned += freeHere ? crossing - from : to - crossing;
          }

          freeHere = freeThere;
        }

        turned = std::clamp(turned, 0.0, 2.0 * std::numbers::pi);

        double swept = turned * perTurn;
        if (swept <= 0.0) continue;

        // What is left to sample for is not how much surface there is but how much of it is buried under a
        // probe that can reach it from elsewhere, which is a fraction and needs far fewer points than the
        // area itself would. The points are weighted by how much surface each stands for, the area element
        // varying across the arc, so that the fraction is of area and not of points.
        double weighed = 0.0;
        double lost = 0.0;

        for (std::size_t step = 0; step < aroundTheCircle; ++step)
        {
          double phi = (static_cast<double>(step) + 0.5) * stepPhi;
          double3 radial = std::cos(phi) * u + std::sin(phi) * v;
          double3 centre = circleCentre + circleRadius * radial;

          if (!freeProbeCentre(centre)) continue;

          for (std::size_t across = 0; across < acrossTheArc; ++across)
          {
            double theta = toFirst + (static_cast<double>(across) + 0.5) * (toSecond - toFirst) /
                                         static_cast<double>(acrossTheArc);
            if (cusp > 0.0 && std::abs(theta) < cusp) continue;

            double standOff = circleRadius - probeRadius * std::cos(theta);
            if (standOff <= 0.0) continue;

            weighed += standOff;
            if (buried(centre + probeRadius * (-std::cos(theta) * radial + std::sin(theta) * axis)))
            {
              lost += standOff;
            }
          }
        }

        double buriedShare = weighed > 0.0 ? lost / weighed : 0.0;

        saddleOfAtom[atom] += swept * (1.0 - buriedShare);
        buriedOfAtom[atom] += swept * buriedShare;
      }
    }

    for (std::size_t atom = 0; atom < numberOfAtoms; ++atom)
    {
      self.saddleArea += saddleOfAtom[atom];
      self.buriedSaddleArea += buriedOfAtom[atom];
      self.numberOfPairs += pairsOfAtom[atom];
    }

    // The sweep is a quadrature rather than a sample, so it has no statistical error; what it has is a
    // discretisation error, and halving the steps is how that is seen.
    self.saddleAreaError = 0.0;
  }

  // The concave patches: where a free probe settles against three atoms at once it stops, and the piece of
  // its own surface facing into the corner is part of the boundary. Those positions are found by putting a
  // sphere at the right distance from three atoms at a time, which has at most two solutions.
  {
    // A corner and every direction in which the probe there is touching something. Three of them come from
    // the triple that found it; a fourth atom arriving at exactly the same distance is a corner the probe
    // rests on more heavily, and the patch there is bounded by all of them and not by three.
    const std::size_t samplesPerCorner = std::max<std::size_t>(256, cornerSamples);

    std::vector<double> concaveOfCorner(corners.size(), 0.0);
    std::vector<double> buriedOfCorner(corners.size(), 0.0);
    std::vector<double> varianceOfCorner(corners.size(), 0.0);

#pragma omp parallel for schedule(dynamic, 8)
    for (std::int64_t index = 0; index < static_cast<std::int64_t>(corners.size()); ++index)
    {
      const Corner &corner = corners[static_cast<std::size_t>(index)];

      // The patch is bounded by the great circles through pairs of contact directions -- but only by those
      // pairs that leave every other contact on one side, which for three contacts is all of them and for
      // four or more is the edges of their spherical convex hull. Bounding it by every pair would cut the
      // patch down to nothing wherever the probe rests on more than three atoms.
      std::vector<double3> edges;
      for (std::size_t a = 0; a + 1 < corner.contact.size(); ++a)
      {
        for (std::size_t b = a + 1; b < corner.contact.size(); ++b)
        {
          double3 normal = double3::cross(corner.contact[a], corner.contact[b]);
          double length = normal.length();
          if (length < 1.0e-9) continue;
          normal = (1.0 / length) * normal;

          bool anyAbove = false;
          bool anyBelow = false;
          for (std::size_t c = 0; c < corner.contact.size(); ++c)
          {
            if (c == a || c == b) continue;
            double side = double3::dot(normal, corner.contact[c]);
            if (side > 1.0e-9) anyAbove = true;
            if (side < -1.0e-9) anyBelow = true;
          }
          if (anyAbove && anyBelow) continue;  // the pair spans the hull rather than bounding it

          edges.push_back(anyBelow ? -1.0 * normal : normal);
        }
      }
      if (edges.size() < 3) continue;

      // Drawing over the whole sphere to find a small patch is mostly waste, so the draw is confined to a
      // cap that holds it and the cap's own area is what the fraction is taken of.
      double3 middle(0.0, 0.0, 0.0);
      for (const double3 &contact : corner.contact) middle += contact;

      double middleLength = middle.length();
      if (middleLength <= 1.0e-12) continue;
      middle = (1.0 / middleLength) * middle;

      double cosine = 1.0;
      for (const double3 &contact : corner.contact)
      {
        cosine = std::min(cosine, double3::dot(middle, contact));
      }
      cosine = std::clamp(cosine - 1.0e-6, -1.0, 1.0);

      double capArea = 2.0 * std::numbers::pi * probeRadius * probeRadius * (1.0 - cosine);

      double3 u = somePerpendicular(middle);
      double3 v = double3::cross(middle, u);

      RandomNumber random{static_cast<std::size_t>(index)};

      std::size_t kept = 0;
      std::size_t lost = 0;

      for (std::size_t sample = 0; sample < samplesPerCorner; ++sample)
      {
        // Uniform on the cap: the height is uniform and the angle about the axis is uniform.
        double height = cosine + (1.0 - cosine) * random.uniform();
        double ring = std::sqrt(std::max(1.0 - height * height, 0.0));
        double angle = 2.0 * std::numbers::pi * random.uniform();

        double3 direction = height * middle + (ring * std::cos(angle)) * u + (ring * std::sin(angle)) * v;

        bool inside = std::ranges::all_of(edges, [&](const double3 &normal)
                                          { return double3::dot(normal, direction) >= 0.0; });
        if (!inside) continue;

        if (buried(corner.centre + probeRadius * direction))
          ++lost;
        else
          ++kept;
      }

      double perSample = capArea / static_cast<double>(samplesPerCorner);
      concaveOfCorner[static_cast<std::size_t>(index)] = static_cast<double>(kept) * perSample;
      buriedOfCorner[static_cast<std::size_t>(index)] = static_cast<double>(lost) * perSample;

      double fraction = static_cast<double>(kept) / static_cast<double>(samplesPerCorner);
      varianceOfCorner[static_cast<std::size_t>(index)] =
          capArea * capArea * fraction * (1.0 - fraction) / static_cast<double>(samplesPerCorner);
    }

    double variance = 0.0;
    for (std::size_t corner = 0; corner < corners.size(); ++corner)
    {
      self.concaveArea += concaveOfCorner[corner];
      self.buriedConcaveArea += buriedOfCorner[corner];
      variance += varianceOfCorner[corner];

      std::size_t touching = std::min(corners[corner].contact.size(), self.concaveByContacts.size() - 1);
      self.concaveByContacts[touching] += concaveOfCorner[corner];
    }
    self.concaveAreaError = std::sqrt(variance);
  }

  std::chrono::duration<double> timing = std::chrono::steady_clock::now() - begin;
  self.seconds = timing.count();

  return self;
}
