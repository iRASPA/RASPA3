module;

module exact_solvent_excluded;

import std;

import double2;
import double3;
import int3;
import simulationbox;
import voronoi_accessibility;
import exact_boundary_components;
import exact_surface_patches;
import exact_union_volume;

constexpr std::size_t gaussOrder = exactQuadratureOrder;

// Gaps in a circle of latitude shorter than this are dropped: they are the seams where two caps meet almost
// tangentially and carry no solid angle worth the trigonometry.
constexpr double gapTolerance = 1.0e-12;

// How near a sphere has to pass to a vertex to count as touching it. The three spheres that put the vertex
// there pass through it to round-off, and a fourth either passes through it exactly, which a framework's own
// symmetry arranges often enough, or misses it by a length the structure sets. There is no scale in between,
// so the tolerance is not a parameter of the answer.
constexpr double touchTolerance = 1.0e-7;

// Two vertices nearer than this are the same vertex, seen from two of the atoms it lies on or in two lifts of
// the cell. Same argument: the copies agree to round-off and distinct vertices do not come this close.
constexpr double vertexTolerance = 1.0e-7;


// Gauss-Legendre nodes and weights on the unit interval, by Newton's method on the Legendre polynomial.
struct GaussRule
{
  std::array<double, gaussOrder> nodes{};
  std::array<double, gaussOrder> weights{};
};

const GaussRule& unitIntervalGaussRule()
{
  static const GaussRule rule = []()
  {
    GaussRule constructed;
    for (std::size_t i = 0; i < gaussOrder; ++i)
    {
      double abscissa =
          std::cos(std::numbers::pi * (static_cast<double>(i) + 0.75) / (static_cast<double>(gaussOrder) + 0.5));
      double derivative = 0.0;
      for (std::size_t iteration = 0; iteration < 100; ++iteration)
      {
        double previous = 1.0;
        double current = abscissa;
        for (std::size_t k = 2; k <= gaussOrder; ++k)
        {
          double next =
              ((2.0 * static_cast<double>(k) - 1.0) * abscissa * current - (static_cast<double>(k) - 1.0) * previous) /
              static_cast<double>(k);
          previous = current;
          current = next;
        }
        derivative = static_cast<double>(gaussOrder) * (abscissa * current - previous) / (abscissa * abscissa - 1.0);
        double step = current / derivative;
        abscissa -= step;
        if (std::abs(step) < 1.0e-15) break;
      }
      constructed.nodes[i] = 0.5 * (1.0 - abscissa);
      constructed.weights[i] = 1.0 / ((1.0 - abscissa * abscissa) * derivative * derivative);
    }
    return constructed;
  }();
  return rule;
}


double foldedPolarAngle(double angle)
{
  double wrapped = std::fmod(std::abs(angle), 2.0 * std::numbers::pi);
  return (wrapped > std::numbers::pi) ? 2.0 * std::numbers::pi - wrapped : wrapped;
}


// The frame the sphere is swept in, chosen so that no cap's axis sits on the polar axis, where latitude
// slicing degenerates. Out of a fixed set of candidates, so that the answer does not depend on the run.
std::array<double3, 3> sweepFrame(const std::vector<double3>& axes)
{
  static const std::array<double3, 6> candidates = {double3(1.0, 2.0, 3.0),   double3(-3.0, 1.0, 2.0),
                                                    double3(2.0, -3.0, 1.0),  double3(3.0, 2.0, -1.0),
                                                    double3(1.0, -3.0, -2.0), double3(-2.0, 3.0, -1.0)};

  double3 polarAxis = candidates.front() * (1.0 / candidates.front().length());
  double bestSeparation = -1.0;
  for (const double3& candidate : candidates)
  {
    double3 direction = candidate * (1.0 / candidate.length());
    double separation = 1.0;
    for (const double3& axis : axes) separation = std::min(separation, 1.0 - std::abs(double3::dot(direction, axis)));
    if (separation > bestSeparation)
    {
      bestSeparation = separation;
      polarAxis = direction;
    }
    if (bestSeparation > 0.01) break;
  }

  double3 helper(1.0, 0.0, 0.0);
  if (std::abs(polarAxis.y) < std::abs(polarAxis.x)) helper = double3(0.0, 1.0, 0.0);
  if (std::abs(polarAxis.z) < std::min(std::abs(polarAxis.x), std::abs(polarAxis.y))) helper = double3(0.0, 0.0, 1.0);
  double3 first = double3::cross(helper, polarAxis);
  first = first * (1.0 / first.length());
  double3 second = double3::cross(polarAxis, first);
  second = second * (1.0 / second.length());

  return {first, second, polarAxis};
}


// One cap, in the frame the sweep is done in.
struct SweepCircle
{
  double3 axis;
  double cosineHalfAngle{0.0};
  double halfAngle{0.0};
  double polarAngle{0.0};
  double cosinePolar{0.0};
  double sinePolar{0.0};
  double azimuth{0.0};
  double lowestLatitude{0.0};
  double highestLatitude{0.0};
  bool reachesOverPole{false};
  bool reachesOverAntipole{false};
};


SphericalRegion regionOutsideCaps(const std::vector<SphericalCap>& caps, std::size_t subdivisions)
{
  SphericalRegion region;

  std::vector<SweepCircle> circles;
  std::vector<double3> axes;
  for (const SphericalCap& cap : caps)
  {
    if (cap.cosineHalfAngle >= 1.0) continue;        // covers nothing
    if (cap.cosineHalfAngle <= -1.0) return region;  // covers the whole sphere, and the region is empty

    SweepCircle circle;
    circle.axis = cap.axis;
    circle.cosineHalfAngle = cap.cosineHalfAngle;
    circle.halfAngle = std::acos(cap.cosineHalfAngle);
    circles.push_back(circle);
    axes.push_back(circle.axis);
  }

  if (circles.empty())
  {
    region.solidAngle = 4.0 * std::numbers::pi;
    return region;
  }

  // A cap inside another covers nothing of its own and only adds latitudes at which nothing happens.
  if (circles.size() > 1)
  {
    std::vector<bool> redundant(circles.size(), false);
    for (std::size_t i = 0; i < circles.size(); ++i)
    {
      for (std::size_t j = 0; j < circles.size(); ++j)
      {
        if (i == j || redundant[j]) continue;
        double separation = std::acos(std::clamp(double3::dot(circles[i].axis, circles[j].axis), -1.0, 1.0));
        if (separation + circles[i].halfAngle <= circles[j].halfAngle)
        {
          redundant[i] = true;
          break;
        }
      }
    }
    std::size_t kept = 0;
    for (std::size_t i = 0; i < circles.size(); ++i)
    {
      if (!redundant[i]) circles[kept++] = circles[i];
    }
    circles.resize(kept);
    axes.clear();
    for (const SweepCircle& circle : circles) axes.push_back(circle.axis);
  }

  const std::array<double3, 3> frame = sweepFrame(axes);
  const double3 firstAxis = frame[0];
  const double3 secondAxis = frame[1];
  const double3 polarAxis = frame[2];

  std::vector<double> breakpoints{0.0, std::numbers::pi};
  for (SweepCircle& circle : circles)
  {
    circle.polarAngle = std::acos(std::clamp(double3::dot(circle.axis, polarAxis), -1.0, 1.0));
    circle.cosinePolar = std::cos(circle.polarAngle);
    circle.sinePolar = std::sin(circle.polarAngle);
    circle.azimuth = std::atan2(double3::dot(circle.axis, secondAxis), double3::dot(circle.axis, firstAxis));
    if (circle.azimuth < 0.0) circle.azimuth += 2.0 * std::numbers::pi;

    circle.lowestLatitude = foldedPolarAngle(circle.polarAngle - circle.halfAngle);
    circle.highestLatitude = foldedPolarAngle(circle.polarAngle + circle.halfAngle);
    circle.reachesOverPole = circle.polarAngle < circle.halfAngle;
    circle.reachesOverAntipole = std::numbers::pi - circle.polarAngle < circle.halfAngle;

    breakpoints.push_back(circle.lowestLatitude);
    breakpoints.push_back(circle.highestLatitude);
  }

  // Where two caps cross uncovered the boundary of the region turns a corner, and the uncovered length of a
  // latitude is analytic only between such latitudes. A crossing a third cap covers is not a corner.
  for (std::size_t j = 0; j + 1 < circles.size(); ++j)
  {
    for (std::size_t k = j + 1; k < circles.size(); ++k)
    {
      const SweepCircle& first = circles[j];
      const SweepCircle& second = circles[k];
      double alignment = double3::dot(first.axis, second.axis);
      double denominator = 1.0 - alignment * alignment;
      if (denominator < 1.0e-14) continue;

      double alongFirst = (first.cosineHalfAngle - alignment * second.cosineHalfAngle) / denominator;
      double alongSecond = (second.cosineHalfAngle - alignment * first.cosineHalfAngle) / denominator;
      double outOfPlaneSquared =
          (1.0 - alongFirst * first.cosineHalfAngle - alongSecond * second.cosineHalfAngle) / denominator;
      if (outOfPlaneSquared <= 0.0) continue;

      double3 inPlane = first.axis * alongFirst + second.axis * alongSecond;
      double3 outOfPlane = double3::cross(first.axis, second.axis) * std::sqrt(outOfPlaneSquared);
      for (std::size_t side = 0; side < 2; ++side)
      {
        double3 corner = (side == 0) ? inPlane + outOfPlane : inPlane - outOfPlane;
        bool covered = false;
        for (std::size_t l = 0; l < circles.size() && !covered; ++l)
        {
          if (l == j || l == k) continue;
          covered = double3::dot(corner, circles[l].axis) > circles[l].cosineHalfAngle + 1.0e-9;
        }
        if (!covered) breakpoints.push_back(std::acos(std::clamp(double3::dot(corner, polarAxis), -1.0, 1.0)));
      }
    }
  }

  std::sort(breakpoints.begin(), breakpoints.end());

  const GaussRule& rule = unitIntervalGaussRule();
  const double twoPi = 2.0 * std::numbers::pi;
  const std::size_t parts = std::max<std::size_t>(1, subdivisions);

  std::vector<std::size_t> cutting;
  std::vector<std::pair<double, double>> covered;

  for (std::size_t piece = 0; piece + 1 < breakpoints.size(); ++piece)
  {
    double pieceBegin = breakpoints[piece];
    double pieceEnd = breakpoints[piece + 1];
    if (pieceEnd - pieceBegin < 1.0e-14) continue;

    // Which caps cut the latitudes of this piece, and whether one of them covers them outright. Both can
    // change only at the ends of the piece, so both are settled once here.
    double interior = 0.5 * (pieceBegin + pieceEnd);
    cutting.clear();
    bool buried = false;
    for (std::size_t i = 0; i < circles.size(); ++i)
    {
      const SweepCircle& circle = circles[i];
      if (interior > circle.lowestLatitude && interior < circle.highestLatitude)
        cutting.push_back(i);
      else if ((interior <= circle.lowestLatitude && circle.reachesOverPole) ||
               (interior >= circle.highestLatitude && circle.reachesOverAntipole))
      {
        buried = true;
        break;
      }
    }
    if (buried) continue;

    for (std::size_t part = 0; part < parts; ++part)
    {
      double begin = pieceBegin + (pieceEnd - pieceBegin) * static_cast<double>(part) / static_cast<double>(parts);
      double end = pieceBegin + (pieceEnd - pieceBegin) * static_cast<double>(part + 1) / static_cast<double>(parts);
      double middle = 0.5 * (begin + end);

      // The uncovered length leaves the end of a piece like the square root of the distance to it, a cap
      // appearing there with that shape, so each half is anchored at its end and a square substituted.
      for (std::size_t half = 0; half < 2; ++half)
      {
        double anchor = (half == 0) ? begin : end;
        double span = (half == 0) ? middle - begin : end - middle;
        double direction = (half == 0) ? 1.0 : -1.0;

        for (std::size_t node = 0; node < gaussOrder; ++node)
        {
          double parameter = rule.nodes[node];
          double latitude = anchor + direction * span * parameter * parameter;
          double weight = 2.0 * span * parameter * rule.weights[node];
          double sineLatitude = std::sin(latitude);
          double cosineLatitude = std::cos(latitude);
          if (sineLatitude <= 0.0) continue;

          covered.clear();
          for (std::size_t i : cutting)
          {
            const SweepCircle& circle = circles[i];
            double cosineHalfWidth =
                (circle.cosineHalfAngle - cosineLatitude * circle.cosinePolar) / (sineLatitude * circle.sinePolar);
            if (cosineHalfWidth >= 1.0) continue;
            double halfWidth = (cosineHalfWidth <= -1.0) ? std::numbers::pi : std::acos(cosineHalfWidth);

            double arcBegin = circle.azimuth - halfWidth;
            double arcEnd = circle.azimuth + halfWidth;
            if (arcBegin < 0.0)
            {
              covered.emplace_back(arcBegin + twoPi, twoPi);
              covered.emplace_back(0.0, arcEnd);
            }
            else if (arcEnd > twoPi)
            {
              covered.emplace_back(arcBegin, twoPi);
              covered.emplace_back(0.0, arcEnd - twoPi);
            }
            else
            {
              covered.emplace_back(arcBegin, arcEnd);
            }
          }
          std::sort(covered.begin(), covered.end());

          double cursor = 0.0;
          for (std::size_t arc = 0; arc <= covered.size(); ++arc)
          {
            double gapEnd = (arc < covered.size()) ? covered[arc].first : twoPi;
            double gap = gapEnd - cursor;
            if (gap > gapTolerance)
            {
              region.solidAngle += sineLatitude * gap * weight;
              double3 direction3 = firstAxis * (sineLatitude * (std::sin(gapEnd) - std::sin(cursor))) +
                                   secondAxis * (sineLatitude * (std::cos(cursor) - std::cos(gapEnd))) +
                                   polarAxis * (cosineLatitude * gap);
              region.moment += direction3 * (sineLatitude * weight);
            }
            if (arc < covered.size()) cursor = std::max(cursor, covered[arc].second);
          }
        }
      }
    }
  }

  return region;
}


// One connected stretch of the circle in which two inflated spheres cut one another, that no third sphere
// covers. It is the curve the probe's centre follows as it rolls in the crease between the two atoms, so it
// carries one toroidal patch of the excluded surface.
struct ExposedArc
{
  double radius{0.0};          // s, the radius of the rolling circle
  double firstOffset{0.0};     // from the first atom's centre to the circle's plane, along the axis
  double secondOffset{0.0};    // from that plane on to the second atom's centre
  double separation{0.0};      // between the two centres
  double firstInflated{0.0};   // R_i
  double firstBare{0.0};       // r_i
  double secondBare{0.0};      // r_j
  double extent{0.0};          // the angle the stretch subtends about the axis
  std::int32_t component{-1};  // the connected surface it belongs to
};


// Whether the sphere of `atom` is the side of the circle toward (`neighbour`, `image`) that owns it, the two
// sides finding the same circle: the one that comes first as an atom of the home cell and a lift of it. Read
// from either side the answer is the same, which is what makes the circle counted once.
bool ownsCircle(std::size_t atom, std::size_t neighbour, const int3& image)
{
  if (atom != neighbour) return atom < neighbour;
  if (image.x != 0) return image.x > 0;
  if (image.y != 0) return image.y > 0;
  return image.z > 0;
}


std::vector<ExposedArc> exposedArcs(const VoronoiAccessibility& accessibility, double probeRadius,
                                    const BoundaryComponents& components)
{
  std::vector<ExposedArc> arcs;
  const double twoPi = 2.0 * std::numbers::pi;

  for (std::size_t i = 0; i < accessibility.atomPositions.size(); ++i)
  {
    const SphereBoundary& boundary = components.atoms[i];
    if (boundary.buried) continue;

    const double firstInflated = accessibility.atomRadii[i];
    const double3 centre = accessibility.atomPositions[i];

    for (const SphereCircle& circle : boundary.circles)
    {
      // The circle lies on both spheres it cuts, and each sphere's sweep finds it, so it is taken from the
      // side that owns it. That the two sides agree on which stretches are exposed is not an assumption --
      // covered by a third sphere is a statement about a point of space, not about which sphere it is looked
      // at from.
      if (!ownsCircle(i, circle.neighbourIndex, circle.neighbourImage)) continue;

      const double3 neighbourCentre =
          accessibility.atomPositions[circle.neighbourIndex] +
          accessibility.simulationBox.cell * double3(static_cast<double>(circle.neighbourImage.x),
                                                     static_cast<double>(circle.neighbourImage.y),
                                                     static_cast<double>(circle.neighbourImage.z));

      ExposedArc common;
      common.radius = firstInflated * circle.sineHalfAngle;
      common.firstOffset = firstInflated * circle.cosineHalfAngle;
      common.separation = (neighbourCentre - centre).length();
      common.secondOffset = common.separation - common.firstOffset;
      common.firstInflated = firstInflated;
      common.firstBare = firstInflated - probeRadius;
      common.secondBare = accessibility.atomRadii[circle.neighbourIndex] - probeRadius;
      if (common.radius <= 1.0e-12 || common.separation <= 1.0e-12) continue;

      const std::size_t count = circle.numberOfArcs();
      for (std::size_t k = 0; k < count && k < circle.arcPatch.size(); ++k)
      {
        std::int32_t patch = circle.arcPatch[k];
        if (patch < 0) continue;  // the stretch is covered and carries no surface

        double extent = twoPi;
        if (!circle.cornerAngles.empty())
        {
          double begin = circle.cornerAngles[k];
          double end = (k + 1 < circle.cornerAngles.size()) ? circle.cornerAngles[k + 1] : circle.cornerAngles[0] + twoPi;
          extent = end - begin;
        }
        if (extent <= 1.0e-14) continue;

        ExposedArc arc = common;
        arc.extent = extent;
        arc.component = components.componentOfPatch[i][static_cast<std::size_t>(patch)];
        arcs.push_back(arc);
      }
    }
  }

  return arcs;
}


// The shell between the accessible and the excluded surface, over the piece belonging to one arc. It is the
// solid swept by the segment from a point of the arc to the excluded surface as the point runs along the arc
// and the segment turns through the crease, and the integral over the turn is elementary. Where the rolling
// circle is smaller than the probe the toroidal patch folds back on itself and meets itself at a cusp: the
// turn then stops at the cusp on either side, and what the two stopped turns leave out is the wedge between
// the cusp circle and the axis, added back by the second term.
double shellArcVolume(const ExposedArc& arc, double probeRadius)
{
  const double s = arc.radius;
  const double r = probeRadius;
  const double first = std::atan2(arc.firstOffset, s);
  const double second = std::atan2(arc.secondOffset, s);

  if (s > r)
  {
    return arc.extent * r * r * (0.5 * s * (first + second) - r / 3.0 * (std::sin(first) + std::sin(second)));
  }

  const double cusp = std::acos(std::clamp(s / r, -1.0, 1.0));
  return arc.extent * r * r *
             (0.5 * s * (first + second - 2.0 * cusp) -
              r / 3.0 * (std::sin(first) + std::sin(second) - 2.0 * std::sin(cusp))) +
         arc.extent * s * s * std::sqrt(std::max(r * r - s * s, 0.0)) / 3.0;
}


// One stretch of the turn through the crease that carries excluded surface.
struct TurnStretch
{
  double begin{0.0};
  double end{0.0};
};

// The turn the toroidal patch of one arc is traced through, in the angle psi measured from the point of the
// rolling probe nearest the axis, and clipped to where the patch is a surface at all.
//
// The probe's surface is at a distance s - r cos psi from the axis, and the toroidal patch is what it sweeps
// between the two tangency directions, psi from -first to second. When the rolling circle is wider than the
// probe that whole turn is a surface. When it is not, the probe reaches across the axis and the sweep folds
// back through it: the distance above goes negative for |psi| < acos(s/r), and there the sweep traces the
// inside of the spindle, which is solid and not boundary. So the turn is cut at that cusp on either side, and
// what is left is one stretch, or two, or none.
std::size_t toroidalTurn(const ExposedArc& arc, double probeRadius, std::array<TurnStretch, 2>& turns)
{
  const double begin = -std::atan2(arc.firstOffset, arc.radius);
  const double end = std::atan2(arc.secondOffset, arc.radius);
  if (end <= begin) return 0uz;

  if (arc.radius > probeRadius)
  {
    turns[0] = {begin, end};
    return 1uz;
  }

  const double cusp = std::acos(std::clamp(arc.radius / probeRadius, -1.0, 1.0));
  std::size_t pieces = 0uz;
  if (begin < -cusp) turns[pieces++] = {begin, std::min(-cusp, end)};
  if (end > cusp) turns[pieces++] = {std::max(cusp, begin), end};
  return pieces;
}


// The area of the toroidal patch of one arc, which is the same turn with the surface element (s - r cos psi) r
// dphi dpsi integrated over it rather than a speed or a volume.
double toroidalPatchArea(const ExposedArc& arc, double probeRadius)
{
  std::array<TurnStretch, 2> turns;
  const std::size_t pieces = toroidalTurn(arc, probeRadius, turns);

  double total = 0.0;
  for (std::size_t piece = 0; piece < pieces; ++piece)
  {
    total += arc.radius * (turns[piece].end - turns[piece].begin) -
             probeRadius * (std::sin(turns[piece].end) - std::sin(turns[piece].begin));
  }
  return arc.extent * probeRadius * total;
}


// The arc's share of the derivative of the excluded volume, which is the integral of the normal speed of its
// toroidal patch. With the rolling circle at radius s(r) and offset t(r) from the first centre,
//
//   x = c + (t + r sin psi) e + (s - r cos psi) rho ,   n = cos psi rho - sin psi e ,
//   v_n = s' cos psi - t' sin psi - 1 ,                 dA = (s - r cos psi) r dphi dpsi ,
//
// and t' = (r_i - r_j)/d is constant while s' = (R_i - t t')/s, so the integral over the turn is elementary.
// The same cusp cuts the turn into the same pieces as the volume above.
double toroidalNormalSpeed(const ExposedArc& arc, double probeRadius)
{
  const double s = arc.radius;
  const double r = probeRadius;
  if (r <= 0.0) return 0.0;

  const double offsetSpeed = (arc.firstBare - arc.secondBare) / arc.separation;
  const double radiusSpeed = (arc.firstInflated - arc.firstOffset * offsetSpeed) / s;

  auto antiderivative = [&](double turn)
  {
    return s * radiusSpeed * std::sin(turn) - 0.5 * r * radiusSpeed * (turn + std::sin(turn) * std::cos(turn)) +
           s * offsetSpeed * std::cos(turn) + 0.5 * r * offsetSpeed * std::sin(turn) * std::sin(turn) - s * turn +
           r * std::sin(turn);
  };

  std::array<TurnStretch, 2> turns;
  const std::size_t pieces = toroidalTurn(arc, r, turns);

  double total = 0.0;
  for (std::size_t piece = 0; piece < pieces; ++piece)
  {
    total += antiderivative(turns[piece].end) - antiderivative(turns[piece].begin);
  }
  return arc.extent * r * total;
}


// One vertex of the accessible surface: a position at which the probe's centre touches three or more atoms at
// once, and so a corner where three or more toroidal patches meet. It carries the concave patch of the
// excluded surface.
struct SasVertex
{
  double3 position;
  std::size_t owner{0};  // the atom of the home cell whose sweep found it, and whose lift it is in
  std::vector<double3> tangents;     // unit, from the vertex toward each touching atom's centre
  std::vector<double3> hullNormals;  // the region on the probe's sphere is where all of these are met
  double3 velocity{0.0, 0.0, 0.0};   // how the vertex moves as the probe grows
  std::int32_t component{-1};
};


// Solves a three by three system by Cramer's rule, returning false where the matrix is too near singular for
// the answer to mean anything.
bool solveThreeByThree(const std::array<double3, 3>& rows, const double3& right, double3& answer)
{
  double determinant = double3::dot(rows[0], double3::cross(rows[1], rows[2]));
  double scale = rows[0].length() * rows[1].length() * rows[2].length();
  if (scale <= 0.0 || std::abs(determinant) < 1.0e-10 * scale) return false;

  std::array<double3, 3> columns = {double3(rows[0].x, rows[1].x, rows[2].x), double3(rows[0].y, rows[1].y, rows[2].y),
                                    double3(rows[0].z, rows[1].z, rows[2].z)};
  answer = double3(double3::dot(right, double3::cross(columns[1], columns[2])),
                   double3::dot(columns[0], double3::cross(right, columns[2])),
                   double3::dot(columns[0], double3::cross(columns[1], right))) *
           (1.0 / determinant);
  return true;
}


// The vertices of the accessible surface, one per position rather than one per description of it.
//
// The corners the decomposition already found are the vertices: a corner is a crossing of two of the circles
// bounding one sphere that no third sphere covers, which is exactly a point on three spheres and inside none.
// Each is found once from every sphere it lies on and every circle through it, in as many lifts of the cell,
// so the descriptions are gathered by position, and the atoms the surviving vertex touches are then asked for
// afresh. That last step is what makes a vertex on four spheres or more come out right without being a case:
// it is one vertex with four tangent directions, and the concave patch there is bounded by the hull of them.
std::vector<SasVertex> surfaceVertices(const VoronoiAccessibility& accessibility, const BoundaryComponents& components,
                                       std::size_t& degenerate, std::size_t& rejected)
{
  const SimulationBox& simulationBox = accessibility.simulationBox;
  const int3 gridSize = accessibility.gridSize;

  auto binIndex = [&](const double3& fractional)
  {
    int3 bin(std::min(gridSize.x - 1, static_cast<int>(fractional.x * static_cast<double>(gridSize.x))),
             std::min(gridSize.y - 1, static_cast<int>(fractional.y * static_cast<double>(gridSize.y))),
             std::min(gridSize.z - 1, static_cast<int>(fractional.z * static_cast<double>(gridSize.z))));
    return std::make_pair(bin, static_cast<std::size_t>((bin.z * gridSize.y + bin.y) * gridSize.x + bin.x));
  };

  const std::size_t numberOfBins =
      static_cast<std::size_t>(gridSize.x) * static_cast<std::size_t>(gridSize.y) * static_cast<std::size_t>(gridSize.z);
  std::vector<std::vector<std::size_t>> bins(numberOfBins);

  std::vector<SasVertex> vertices;
  std::vector<double3> wrapped;  // the same positions folded into the cell, for the comparison

  for (std::size_t i = 0; i < accessibility.atomPositions.size(); ++i)
  {
    const SphereBoundary& boundary = components.atoms[i];
    if (boundary.buried) continue;
    const double radius = accessibility.atomRadii[i];
    const double3 centre = accessibility.atomPositions[i];

    for (const SphereCircle& circle : boundary.circles)
    {
      for (double angle : circle.cornerAngles)
      {
        double3 position = centre + circle.direction(angle) * radius;
        double3 fractional = double3::fract(simulationBox.inverseCell * position);
        double3 folded = simulationBox.cell * fractional;
        auto [bin, index] = binIndex(fractional);

        // Already known? The copies of one vertex agree to round-off, so only the bins around it are
        // searched, and the search wraps because a vertex on the face of the cell folds to either side.
        bool seen = false;
        for (int ox = -1; ox <= 1 && !seen; ++ox)
        {
          for (int oy = -1; oy <= 1 && !seen; ++oy)
          {
            for (int oz = -1; oz <= 1 && !seen; ++oz)
            {
              int bx = (bin.x + ox + gridSize.x) % gridSize.x;
              int by = (bin.y + oy + gridSize.y) % gridSize.y;
              int bz = (bin.z + oz + gridSize.z) % gridSize.z;
              for (std::size_t other : bins[static_cast<std::size_t>((bz * gridSize.y + by) * gridSize.x + bx)])
              {
                double3 delta = simulationBox.applyPeriodicBoundaryConditions(folded - wrapped[other]);
                if (delta.length() < vertexTolerance)
                {
                  seen = true;
                  break;
                }
              }
            }
          }
        }
        if (seen) continue;

        SasVertex vertex;
        vertex.position = position;
        vertex.owner = i;
        bins[index].push_back(vertices.size());
        vertices.push_back(vertex);
        wrapped.push_back(folded);
      }
    }
  }

  // Which atoms each surviving vertex touches, asked of the structure rather than remembered from the
  // corner it came from, and the hull of the directions to them.
  std::vector<SasVertex> kept;
  kept.reserve(vertices.size());
  for (SasVertex& vertex : vertices)
  {
    bool inside = false;
    for (const NeighbourImage& image :
         accessibility.neighbourAtomImages(vertex.position, accessibility.maximumAtomRadius + touchTolerance))
    {
      double distance = image.delta.length();
      if (distance < image.radius - touchTolerance)
      {
        inside = true;  // the point is buried after all, and is on no surface
        break;
      }
      if (distance < image.radius + touchTolerance && distance > 1.0e-12)
      {
        vertex.tangents.push_back(image.delta * (1.0 / distance));
      }
    }
    if (inside || vertex.tangents.size() < 3)
    {
      ++rejected;
      continue;
    }
    if (vertex.tangents.size() > 3) ++degenerate;

    // The concave patch is bounded by the great circles through pairs of tangent directions, and for more
    // than three of them by those pairs that are edges of their spherical convex hull: the pair whose plane
    // leaves every other direction on one side.
    const std::size_t count = vertex.tangents.size();
    for (std::size_t a = 0; a + 1 < count; ++a)
    {
      for (std::size_t b = a + 1; b < count; ++b)
      {
        double3 normal = double3::cross(vertex.tangents[a], vertex.tangents[b]);
        double length = normal.length();
        if (length < 1.0e-9) continue;
        normal = normal * (1.0 / length);

        bool anyPositive = false;
        bool anyNegative = false;
        for (std::size_t c = 0; c < count; ++c)
        {
          if (c == a || c == b) continue;
          double side = double3::dot(normal, vertex.tangents[c]);
          if (side > 1.0e-9) anyPositive = true;
          if (side < -1.0e-9) anyNegative = true;
        }
        if (anyPositive && anyNegative) continue;  // not an edge of the hull
        vertex.hullNormals.push_back(anyNegative ? -normal : normal);
      }
    }
    if (vertex.hullNormals.size() < 3)
    {
      ++rejected;
      continue;
    }

    // How the vertex moves as the probe grows, from |x - c_m|^2 = (r_m + r)^2 differentiated: the tangency
    // conditions (x - c_m).x' = R_m, three of which fix it. Where more than three atoms meet they are
    // consistent, so the best conditioned triple is taken.
    bool moves = false;
    double bestScale = 0.0;
    for (std::size_t a = 0; a + 2 < count; ++a)
    {
      for (std::size_t b = a + 1; b + 1 < count; ++b)
      {
        for (std::size_t c = b + 1; c < count; ++c)
        {
          std::array<double3, 3> rows = {vertex.tangents[a], vertex.tangents[b], vertex.tangents[c]};
          double scale = std::abs(double3::dot(rows[0], double3::cross(rows[1], rows[2])));
          if (scale <= bestScale) continue;
          // Dividing (x - c_m).x' = R_m through by R_m turns each row into the unit direction from the vertex
          // toward that centre, which points the other way, so the right-hand side is minus one.
          double3 velocity;
          if (!solveThreeByThree(rows, double3(-1.0, -1.0, -1.0), velocity)) continue;
          bestScale = scale;
          vertex.velocity = velocity;
          moves = true;
        }
      }
    }
    if (!moves)
    {
      ++rejected;
      continue;
    }

    // Which connected surface the vertex is on. It lies on one of the bounding circles of its own atom, and
    // the arcs of those circles carry the patch they bound, so this is a search along a circle.
    const double radius = accessibility.atomRadii[vertex.owner];
    double3 direction = (vertex.position - accessibility.atomPositions[vertex.owner]) * (1.0 / radius);
    std::int32_t patch = components.patchOfCirclePoint(vertex.owner, direction);
    if (patch >= 0) vertex.component = components.componentOfPatch[vertex.owner][static_cast<std::size_t>(patch)];

    kept.push_back(std::move(vertex));
  }

  return kept;
}


// A cell list over the vertices, so that the probes near one probe are found without a scan. It is laid over
// the same bins as the atoms, which are already sized to the structure.
struct VertexGrid
{
  const SimulationBox* simulationBox{nullptr};
  int3 gridSize{1, 1, 1};
  double minimumBinWidth{0.0};
  std::vector<std::vector<std::size_t>> bins;
  std::vector<double3> folded;

  // Every image of every vertex, its own images included but not itself, whose position lies within `reach` of
  // `point`, as positions in the lift the point is in.
  std::vector<double3> within(const double3& point, double reach) const
  {
    double3 fractional = double3::fract(simulationBox->inverseCell * point);
    double3 wrappedPoint = simulationBox->cell * fractional;
    int3 pointBin(std::min(gridSize.x - 1, static_cast<int>(fractional.x * static_cast<double>(gridSize.x))),
                  std::min(gridSize.y - 1, static_cast<int>(fractional.y * static_cast<double>(gridSize.y))),
                  std::min(gridSize.z - 1, static_cast<int>(fractional.z * static_cast<double>(gridSize.z))));

    std::vector<double3> found;
    for (int k = 0;; ++k)
    {
      double lowerBound = static_cast<double>(k - 1) * minimumBinWidth;
      if (k > 0 && lowerBound > reach) break;

      for (int ox = -k; ox <= k; ++ox)
      {
        for (int oy = -k; oy <= k; ++oy)
        {
          for (int oz = -k; oz <= k; ++oz)
          {
            if (std::max({std::abs(ox), std::abs(oy), std::abs(oz)}) != k) continue;

            auto wrap = [](int coordinate, int extent)
            {
              int image = (coordinate >= 0) ? coordinate / extent : -((-coordinate + extent - 1) / extent);
              return std::make_pair(coordinate - image * extent, image);
            };
            auto [bx, lx] = wrap(pointBin.x + ox, gridSize.x);
            auto [by, ly] = wrap(pointBin.y + oy, gridSize.y);
            auto [bz, lz] = wrap(pointBin.z + oz, gridSize.z);

            double3 imageShift = simulationBox->cell * double3(static_cast<double>(lx), static_cast<double>(ly),
                                                               static_cast<double>(lz)) -
                                 wrappedPoint;

            for (std::size_t j : bins[static_cast<std::size_t>((bz * gridSize.y + by) * gridSize.x + bx)])
            {
              double3 delta = folded[j] + imageShift;
              double distance = delta.length();
              if (distance > vertexTolerance && distance <= reach) found.push_back(point + delta);
            }
          }
        }
      }
    }
    return found;
  }
};


VertexGrid buildVertexGrid(const VoronoiAccessibility& accessibility, const std::vector<SasVertex>& vertices)
{
  VertexGrid grid;
  grid.simulationBox = &accessibility.simulationBox;
  grid.gridSize = accessibility.gridSize;
  grid.minimumBinWidth = accessibility.minimumBinWidth;
  grid.bins.assign(static_cast<std::size_t>(grid.gridSize.x) * static_cast<std::size_t>(grid.gridSize.y) *
                       static_cast<std::size_t>(grid.gridSize.z),
                   {});
  grid.folded.reserve(vertices.size());

  for (std::size_t i = 0; i < vertices.size(); ++i)
  {
    double3 fractional = double3::fract(accessibility.simulationBox.inverseCell * vertices[i].position);
    grid.folded.push_back(accessibility.simulationBox.cell * fractional);
    int3 bin(std::min(grid.gridSize.x - 1, static_cast<int>(fractional.x * static_cast<double>(grid.gridSize.x))),
             std::min(grid.gridSize.y - 1, static_cast<int>(fractional.y * static_cast<double>(grid.gridSize.y))),
             std::min(grid.gridSize.z - 1, static_cast<int>(fractional.z * static_cast<double>(grid.gridSize.z))));
    grid.bins[static_cast<std::size_t>((bin.z * grid.gridSize.y + bin.y) * grid.gridSize.x + bin.x)].push_back(i);
  }
  return grid;
}


// What one vertex contributes to the shell and to the derivative of the excluded volume.
//
// Both come from the concave patch there, the part of the probe's own sphere at the vertex that is not reached
// into by a probe at any other vertex. That is a region of a sphere cut out by caps twice over: by the three
// or more great circles through pairs of tangent directions, which are where the concave patch runs into the
// toroidal ones, and by one cap per neighbouring vertex within two probe radii. The second is Theorem 5.1 of
// Quan and Stamm, and it is a cap because both points are at the probe's radius from their own vertex: the
// probe at x' reaches x_m + r w exactly where w leans further toward x' than the plane halfway between them,
// so clipping against the ball and taking the nearer vertex are the same clip.
//
// The shell piece is the cone from the vertex over its own boundary. Its lateral faces all pass through the
// vertex and contribute nothing, which leaves the concave patch at the probe's radius and the flat faces the
// neighbouring probes leave on those halfway planes, each a disc cut by half planes.
struct VertexTerms
{
  double shell{0.0};
  double normalSpeed{0.0};

  // The concave patch's own area, which is the solid angle left to it on the probe's sphere at the radius the
  // probe is resting at. The clipping is already done for the shell, so this costs a multiplication.
  double area{0.0};

  bool clipped{false};
  bool vanished{false};
};

VertexTerms vertexTerms(const SasVertex& vertex, double probeRadius, const std::vector<double3>& neighbours,
                        std::size_t subdivisions)
{
  VertexTerms terms;
  const double radius = probeRadius;
  if (radius <= 0.0) return terms;

  std::vector<SphericalCap> caps;
  caps.reserve(vertex.hullNormals.size() + neighbours.size());
  for (const double3& normal : vertex.hullNormals) caps.push_back({-normal, 0.0});
  for (const double3& other : neighbours)
  {
    double3 delta = other - vertex.position;
    double separation = delta.length();
    if (separation >= 2.0 * radius) continue;
    caps.push_back({delta * (1.0 / separation), separation / (2.0 * radius)});
    terms.clipped = true;
  }

  SphericalRegion region = regionOutsideCaps(caps, subdivisions);
  if (region.solidAngle <= 0.0)
  {
    terms.vanished = true;
    return terms;
  }

  terms.shell = region.solidAngle * radius * radius * radius / 3.0;
  terms.normalSpeed = -radius * radius * (double3::dot(vertex.velocity, region.moment) + region.solidAngle);
  terms.area = region.solidAngle * radius * radius;

  // The flat faces. Each lies halfway to a neighbouring probe, is a disc of what the probe's sphere leaves
  // there, and is cut back by the same hull planes as the concave patch and by the halfway planes to the other
  // neighbours.
  std::vector<double2> normals;
  std::vector<double> offsets;
  for (const double3& other : neighbours)
  {
    double3 delta = other - vertex.position;
    double separation = delta.length();
    if (separation >= 2.0 * radius) continue;

    double half = 0.5 * separation;
    double discRadius = std::sqrt(std::max(radius * radius - half * half, 0.0));
    if (discRadius <= 0.0) continue;

    double3 axis = delta * (1.0 / separation);
    double3 helper(1.0, 0.0, 0.0);
    if (std::abs(axis.y) < std::abs(axis.x)) helper = double3(0.0, 1.0, 0.0);
    if (std::abs(axis.z) < std::min(std::abs(axis.x), std::abs(axis.y))) helper = double3(0.0, 0.0, 1.0);
    double3 first = double3::cross(helper, axis);
    first = first * (1.0 / first.length());
    double3 second = double3::cross(axis, first);

    const double3 discCentre = vertex.position + axis * half;

    normals.clear();
    offsets.clear();
    for (const double3& normal : vertex.hullNormals)
    {
      // The disc keeps the side of the hull plane the patch is on, which is `(p - x_m).normal >= 0`.
      normals.emplace_back(-double3::dot(first, normal), -double3::dot(second, normal));
      offsets.push_back(half * double3::dot(axis, normal));
    }
    for (const double3& third : neighbours)
    {
      double3 towards = third - vertex.position;
      if ((towards - delta).length() < vertexTolerance) continue;  // this face's own neighbour
      if (towards.length() >= 2.0 * radius) continue;

      // Nearer to this vertex than to that one: `(p - midpoint).towards <= 0`.
      double3 midpoint = vertex.position + towards * 0.5;
      normals.emplace_back(double3::dot(first, towards), double3::dot(second, towards));
      offsets.push_back(double3::dot(midpoint - discCentre, towards));
    }

    terms.shell += clippedDiscArea(discRadius, normals, offsets) * half / 3.0;
  }

  return terms;
}


// Which side of the structure each connected surface faces, by the same argument the area is divided by: a
// surface that closes on a translate of itself has the void running away along it, one that closes on itself
// around void seals it, and only a surface standing round a cluster of atoms leaves the network anything to
// answer. Plus one is reachable, minus one is not, zero is undecided.
std::vector<int> componentSides(const BoundaryComponents& components, const ExactSurfaceAreaSample& patches,
                                const std::vector<ComponentLabel>& labels)
{
  std::vector<int> sides(components.numberOfComponents, 0);
  for (std::size_t label = 0; label < components.numberOfComponents; ++label)
  {
    const PoreBoundaryMoments& moments = patches.components[label];
    if (moments.area <= 0.0) continue;

    if (components.componentPercolates[label] != 0)
      sides[label] = 1;
    else if (-(moments.radiusWeightedArea + moments.originWeighted) > 0.0)
      sides[label] = -1;
    else if (labels[label].decided)
      sides[label] = labels[label].accessible ? 1 : -1;
  }
  return sides;
}


SolventExcludedGeometry solventExcludedGeometry(const VoronoiAccessibility& accessibility, double probeRadius,
                                                const BoundaryComponents& components,
                                                const std::vector<ComponentLabel>& labels, std::size_t subdivisions)
{
  // One sweep of the spheres serves four purposes: the area of each atom's own patch, which the shell wants; the
  // radius weighted area, which the volume of the union wants; the same patches on the bare spheres, which are
  // the convex part of the excluded surface; and the moments of each surface, which say which side it faces.
  return solventExcludedGeometry(
      accessibility, probeRadius, components, labels,
      exactAccessibleSurfaceAreaByComponent(accessibility, components, labels, subdivisions), subdivisions);
}


SolventExcludedGeometry solventExcludedGeometry(const VoronoiAccessibility& accessibility, double probeRadius,
                                                const BoundaryComponents& components,
                                                const std::vector<ComponentLabel>& labels,
                                                const ExactSurfaceAreaSample& patches, std::size_t subdivisions)
{
  SolventExcludedGeometry result;
  result.probeRadius = probeRadius;
  result.accessibleVolume = unionOfBallsVolume(accessibility, patches);

  const std::vector<int> sides = componentSides(components, patches, labels);
  auto sideOf = [&](std::int32_t component)
  {
    return (component >= 0 && static_cast<std::size_t>(component) < sides.size())
               ? sides[static_cast<std::size_t>(component)]
               : 0;
  };
  auto addToSide = [&](std::int32_t component, double amount)
  {
    result.distribution += amount;
    int side = sideOf(component);
    if (side > 0)
      result.accessibleDistribution += amount;
    else if (side < 0)
      result.inaccessibleDistribution += amount;
    else
      result.undecidedDistribution += amount;
  };

  // The same routing for area, given which of the three kinds of patch it is on.
  auto addAreaToSide = [&](std::int32_t component, double ExcludedSurfaceAreas::*kind, double amount)
  {
    result.area.*kind += amount;
    int side = sideOf(component);
    if (side > 0)
      result.accessibleArea.*kind += amount;
    else if (side < 0)
      result.inaccessibleArea.*kind += amount;
    else
      result.undecidedArea.*kind += amount;
  };

  // The patches. Each contributes the shell between the bare sphere and the inflated one over its own solid
  // angle, and nothing to the derivative: the bare spheres do not move as the probe grows. Its area on the
  // excluded surface is its own accessible area brought down onto the bare sphere, which the sweep has already
  // done atom by atom and surface by surface.
  double shell = 0.0;
  for (std::size_t i = 0; i < accessibility.atomPositions.size(); ++i)
  {
    const double inflated = accessibility.atomRadii[i];
    if (inflated <= 0.0 || i >= patches.atomArea.size()) continue;
    const double bare = std::max(0.0, inflated - probeRadius);
    const double ratio = bare / inflated;
    shell += patches.atomArea[i] * inflated / 3.0 * (1.0 - ratio * ratio * ratio);
  }
  for (std::size_t label = 0; label < patches.components.size(); ++label)
  {
    addAreaToSide(static_cast<std::int32_t>(label), &ExcludedSurfaceAreas::convex,
                  patches.components[label].convexArea);
  }

  // The arcs, which carry the toroidal patches.
  const std::vector<ExposedArc> arcs = exposedArcs(accessibility, probeRadius, components);
  result.numberOfArcs = arcs.size();
  for (const ExposedArc& arc : arcs)
  {
    if (arc.radius <= probeRadius) ++result.cuspedArcs;
    shell += shellArcVolume(arc, probeRadius);
    double speed = toroidalNormalSpeed(arc, probeRadius);
    result.toroidalDistribution += speed;
    addToSide(arc.component, speed);
    addAreaToSide(arc.component, &ExcludedSurfaceAreas::saddle, toroidalPatchArea(arc, probeRadius));
  }

  // The vertices, which carry the concave patches.
  const std::vector<SasVertex> vertices =
      surfaceVertices(accessibility, components, result.degenerateVertices, result.vanishedVertices);
  result.numberOfVertices = vertices.size();
  const VertexGrid grid = buildVertexGrid(accessibility, vertices);

  for (const SasVertex& vertex : vertices)
  {
    std::vector<double3> neighbours = grid.within(vertex.position, 2.0 * probeRadius);
    VertexTerms terms = vertexTerms(vertex, probeRadius, neighbours, subdivisions);
    if (terms.clipped) ++result.clippedVertices;
    if (terms.vanished) ++result.vanishedVertices;

    shell += terms.shell;
    result.concaveDistribution += terms.normalSpeed;
    addToSide(vertex.component, terms.normalSpeed);
    addAreaToSide(vertex.component, &ExcludedSurfaceAreas::concave, terms.area);
  }

  result.shellVolume = shell;
  result.excludedVolume = result.accessibleVolume - shell;
  result.poreVolume = accessibility.simulationBox.volume - result.excludedVolume;
  return result;
}


SolventExcludedGeometry solventExcludedGeometry(const VoronoiAccessibility& accessibility, double probeRadius,
                                                std::size_t subdivisions)
{
  BoundaryComponents components = boundaryComponents(accessibility);
  std::vector<ComponentLabel> labels = labelBoundaryComponents(accessibility, components);
  return solventExcludedGeometry(accessibility, probeRadius, components, labels, subdivisions);
}
