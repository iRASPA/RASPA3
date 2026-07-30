module;

module exact_solvent_excluded;

import std;

import double2;
import double3;
import int3;
import simulationbox;
import pore_accessibility;
import exact_boundary_components;
import exact_sphere_sweep;
import exact_surface_patches;
import exact_union_volume;

// How near a sphere has to pass to a vertex to count as touching it. The three spheres that put the vertex
// there pass through it to round-off, and a fourth either passes through it exactly, which a framework's own
// symmetry arranges often enough, or misses it by a length the structure sets. There is no scale in between,
// so the tolerance is not a parameter of the answer.
constexpr double touchTolerance = 1.0e-7;

// Two vertices nearer than this are the same vertex, seen from two of the atoms it lies on or in two lifts of
// the cell. Same argument: the copies agree to round-off and distinct vertices do not come this close.
constexpr double vertexTolerance = 1.0e-7;


SphericalRegion regionOutsideCaps(const std::vector<SphericalCap>& caps, std::size_t subdivisions)
{
  SphericalRegion region;

  std::vector<SweepCircle> circles;
  std::vector<double3> axes;
  for (const SphericalCap& cap : caps)
  {
    std::optional<SweepCircle> circle = makeSweepCircle(cap.axis, cap.cosineHalfAngle);
    if (!circle.has_value()) continue;                        // the cap covers nothing
    if (circle->cosineHalfAngle <= -1.0) return region;       // it covers the whole sphere, leaving nothing

    circles.push_back(circle.value());
    axes.push_back(circle->axis);
  }

  if (circles.empty())
  {
    region.solidAngle = 4.0 * std::numbers::pi;
    return region;
  }

  // A cap inside another covers nothing of its own and only adds latitudes at which nothing happens.
  if (circles.size() > 1)
  {
    pruneContainedDiscs(circles);
    axes.clear();
    for (const SweepCircle& circle : circles) axes.push_back(circle.axis);
  }

  const std::array<double3, 3> frame = sweepFrame(axes);
  const double3 firstAxis = frame[0];
  const double3 secondAxis = frame[1];
  const double3 polarAxis = frame[2];

  // Nothing has been over these caps before: they are the hull planes of one vertex and the probes around it,
  // assembled here and thrown away again, so the crossings among them have to be looked for.
  SweepWorkspace work;
  sweepExposedLatitudes(circles, frame, nullptr, subdivisions, work,
                        [&](const LatitudeGap& gap)
                        {
                          region.solidAngle += gap.sineLatitude * gap.span * gap.weight;
                          double3 normalIntegral =
                              firstAxis * (gap.sineLatitude * (gap.sineEnd - gap.sineBegin)) +
                              secondAxis * (gap.sineLatitude * (gap.cosineBegin - gap.cosineEnd)) +
                              polarAxis * (gap.cosineLatitude * gap.span);
                          region.moment += normalIntegral * (gap.sineLatitude * gap.weight);
                        });

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


std::vector<ExposedArc> exposedArcs(const PoreAccessibility& accessibility, double probeRadius,
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


// Which connected surface a vertex is on, as far as one of the atoms it touches can say, added to `votes`.
//
// `direction` is the vertex as a unit direction from that atom's centre, and the vertex is a corner of the
// atom's exposed region: a crossing of two of its bounding circles that no third sphere covers. Of the two
// stretches of a circle meeting at such a crossing exactly one is exposed, since crossing the second circle
// takes the first into that sphere or out of it, and the exposed one carries the patch the corner belongs to.
// That is index arithmetic on the stretches of one circle, with no angle compared against another.
//
// The reason to ask every atom and every circle through the point rather than one of them is symmetry. Where
// four or more spheres meet in a point, as a framework's own symmetry arranges over and over, several
// crossings coincide there, and a stretch on either side of one of them can be covered by one of the others:
// that corner then names no patch at all. The probe resting in a six-ring of LTA touches six atoms and is
// such a corner on thirty circles, of which twelve answer.
//
// The answers agree, and not by luck: every patch meeting at a vertex is joined to the others across the arcs
// running into it, so they are all one surface and any of them names it. Placing the vertex from its position
// instead, as the fallback at the caller does, is not the same thing and not as good. That asks which stretch
// of a circle holds an angle which is the endpoint of two of them, so round-off decides, and an atom in the
// wall of a cage has one patch facing the cage and another facing the pore outside it: the decision can put a
// vertex of the one on the other, carrying the whole of that vertex's shell through the wall.
void gatherComponentVotes(const BoundaryComponents& components, std::size_t atomIndex, const double3& direction,
                          std::vector<std::pair<std::int32_t, std::size_t>>& votes)
{
  const SphereBoundary& boundary = components.atoms[atomIndex];
  if (boundary.buried) return;

  for (const SphereCircle& circle : boundary.circles)
  {
    if (circle.arcPatch.empty()) continue;

    // On this circle, and then at a corner of it. Both comparisons rest on the argument the tolerances are
    // named for: a point of the surface is on a circle of an atom it touches, and at a corner of it, either
    // exactly to round-off or not at all.
    if (std::abs(double3::dot(direction, circle.axis) - circle.cosineHalfAngle) > touchTolerance) continue;

    const double angle = circle.angleOf(direction);
    for (std::size_t corner = 0; corner < circle.cornerAngles.size(); ++corner)
    {
      double difference = std::abs(circle.cornerAngles[corner] - angle);
      difference = std::min(difference, 2.0 * std::numbers::pi - difference);
      if (difference > vertexTolerance) continue;

      // Stretch `corner` begins at this angle and the one before it ends there.
      const std::size_t stretches = circle.arcPatch.size();
      const std::int32_t after = circle.arcPatch[corner];
      const std::int32_t patch = (after >= 0) ? after : circle.arcPatch[(corner + stretches - 1) % stretches];
      if (patch < 0) continue;

      const std::int32_t component = components.componentOfPatch[atomIndex][static_cast<std::size_t>(patch)];
      if (component < 0) continue;
      auto found = std::find_if(votes.begin(), votes.end(), [&](const auto& vote) { return vote.first == component; });
      if (found != votes.end())
        ++found->second;
      else
        votes.emplace_back(component, 1uz);
    }
  }
}


// Every corner of the decomposition, once each rather than once per description of it, with nothing filled in
// but where it is and which atom's sweep found it.
//
// A corner is a crossing of two of the circles bounding one sphere that no third sphere covers, which is
// exactly a point on three spheres and inside none: a vertex of the accessible surface. Each is found once
// from every sphere it lies on and every circle through it, in as many lifts of the cell, so the copies are
// gathered here by position. They agree to round-off, so only the bins around a candidate are searched, and
// the search wraps because a corner on the face of the cell folds to either side of it.
std::vector<SasVertex> distinctCorners(const PoreAccessibility& accessibility,
                                       const BoundaryComponents& components)
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

  std::vector<SasVertex> corners;
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
        bins[index].push_back(corners.size());
        corners.push_back(vertex);
        wrapped.push_back(folded);
      }
    }
  }

  return corners;
}


// The directions from a vertex to the centres of the atoms it touches, asked of the structure rather than
// remembered from the corner the vertex came from. That is what makes a vertex on four spheres or more come
// out right without being a case of its own: it is one vertex with four tangent directions.
//
// False where the point is not on the surface after all --- buried under an atom, or touching fewer than
// three. `votes` collects what each touching atom can say about which connected surface the vertex is on.
bool gatherTangents(const PoreAccessibility& accessibility, const BoundaryComponents& components,
                    SasVertex& vertex, std::vector<std::pair<std::int32_t, std::size_t>>& votes)
{
  votes.clear();
  for (const NeighbourImage& image :
       accessibility.neighbourAtomImages(vertex.position, accessibility.maximumAtomRadius + touchTolerance))
  {
    double distance = image.delta.length();
    if (distance < image.radius - touchTolerance) return false;  // buried after all, and on no surface

    if (distance < image.radius + touchTolerance && distance > 1.0e-12)
    {
      vertex.tangents.push_back(image.delta * (1.0 / distance));
      gatherComponentVotes(components, image.index, image.delta * (-1.0 / distance), votes);
    }
  }
  return vertex.tangents.size() >= 3;
}


// The planes bounding the concave patch on the probe's sphere. It is bounded by the great circles through
// pairs of tangent directions, and for more than three of them by those pairs that are edges of their
// spherical convex hull: the pair whose plane leaves every other direction on one side. False where fewer
// than three such planes are found, which leaves no patch to bound.
bool gatherHullNormals(SasVertex& vertex)
{
  constexpr double coplanarTolerance = 1.0e-9;

  const std::size_t touching = vertex.tangents.size();
  for (std::size_t a = 0; a + 1 < touching; ++a)
  {
    for (std::size_t b = a + 1; b < touching; ++b)
    {
      double3 normal = double3::cross(vertex.tangents[a], vertex.tangents[b]);
      double length = normal.length();
      if (length < coplanarTolerance) continue;
      normal = normal * (1.0 / length);

      bool anyPositive = false;
      bool anyNegative = false;
      for (std::size_t c = 0; c < touching; ++c)
      {
        if (c == a || c == b) continue;
        double side = double3::dot(normal, vertex.tangents[c]);
        if (side > coplanarTolerance) anyPositive = true;
        if (side < -coplanarTolerance) anyNegative = true;
      }
      if (anyPositive && anyNegative) continue;  // not an edge of the hull
      vertex.hullNormals.push_back(anyNegative ? -normal : normal);
    }
  }
  return vertex.hullNormals.size() >= 3;
}


// How the vertex moves as the probe grows, from |x - c_m|^2 = (r_m + r)^2 differentiated: the tangency
// conditions (x - c_m).x' = R_m, three of which fix it. Where more than three atoms meet they are consistent,
// so the best conditioned triple is taken. False where no triple is far enough from coplanar to solve.
bool solveVertexVelocity(SasVertex& vertex)
{
  const std::size_t touching = vertex.tangents.size();
  bool moves = false;
  double bestScale = 0.0;
  for (std::size_t a = 0; a + 2 < touching; ++a)
  {
    for (std::size_t b = a + 1; b + 1 < touching; ++b)
    {
      for (std::size_t c = b + 1; c < touching; ++c)
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
  return moves;
}


// Which connected surface the vertex is on: the surface the atoms it touches agree that it is on, and failing
// that the one a search from its own position lands on. See `gatherComponentVotes` for why the first is both
// the usual answer and the better one.
void placeOnSurface(const PoreAccessibility& accessibility, const BoundaryComponents& components,
                    const std::vector<std::pair<std::int32_t, std::size_t>>& votes, SasVertex& vertex)
{
  std::size_t mostVotes = 0;
  for (const auto& [component, votesFor] : votes)
  {
    if (votesFor > mostVotes)
    {
      mostVotes = votesFor;
      vertex.component = component;
    }
  }
  if (vertex.component >= 0) return;

  const double radius = accessibility.atomRadii[vertex.owner];
  double3 direction = (vertex.position - accessibility.atomPositions[vertex.owner]) * (1.0 / radius);
  std::int32_t patch = components.patchOfCirclePoint(vertex.owner, direction);
  if (patch >= 0) vertex.component = components.componentOfPatch[vertex.owner][static_cast<std::size_t>(patch)];
}


// The vertices of the accessible surface, one per position rather than one per description of it, each with
// the concave patch of the excluded surface it carries.
struct SurfaceVertices
{
  std::vector<SasVertex> vertices;

  // Vertices where more than three inflated spheres meet, which a framework's own symmetry produces rather
  // than a coincidence. Handled like any other and counted because a structure full of them is one to look
  // at twice.
  std::size_t onManySpheres{0};

  // Corners that turned out not to be vertices of the surface at all: buried under an atom after all, or
  // touching too few, or with the directions to them too nearly coplanar for the patch or for the way the
  // vertex moves to mean anything.
  std::size_t discarded{0};
};


SurfaceVertices surfaceVertices(const PoreAccessibility& accessibility, const BoundaryComponents& components)
{
  SurfaceVertices result;
  std::vector<SasVertex> corners = distinctCorners(accessibility, components);
  result.vertices.reserve(corners.size());

  std::vector<std::pair<std::int32_t, std::size_t>> votes;
  for (SasVertex& vertex : corners)
  {
    if (!gatherTangents(accessibility, components, vertex, votes))
    {
      ++result.discarded;
      continue;
    }

    // Counted here rather than below, so that the count is of the vertices a structure has and not of the
    // ones that also came through the two steps after this.
    if (vertex.tangents.size() > 3) ++result.onManySpheres;

    if (!gatherHullNormals(vertex) || !solveVertexVelocity(vertex))
    {
      ++result.discarded;
      continue;
    }

    placeOnSurface(accessibility, components, votes, vertex);
    result.vertices.push_back(std::move(vertex));
  }
  return result;
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


VertexGrid buildVertexGrid(const PoreAccessibility& accessibility, const std::vector<SasVertex>& vertices)
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
    double3 first = perpendicularTo(axis);
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


namespace
{

// Where a piece of the excluded surface puts what it measures.
//
// The surface is made of three kinds of piece --- the convex patches on the atoms, the toroidal ones along the
// arcs, the concave ones at the vertices --- and each of them contributes an area, a volume of shell, and a
// derivative in the probe radius. Every one of those has to land in three places at once: the total, the
// surface the piece belongs to, and the side of the structure that surface faces. This is that routing, in one
// place rather than repeated at each of the three sites.
//
// A piece whose surface is not known, which is what a negative component means, still counts towards the
// totals. It is the per-surface accounts it cannot be entered in.
struct SurfaceTotals
{
  SolventExcludedGeometry& result;
  const std::vector<SurfaceSide>& sides;

  int sideOf(std::int32_t component) const
  {
    return (component >= 0 && static_cast<std::size_t>(component) < sides.size())
               ? sides[static_cast<std::size_t>(component)].side
               : 0;
  }

  // The derivative of the excluded volume in the probe radius.
  void addDistribution(std::int32_t component, double amount)
  {
    result.distribution += amount;
    if (component >= 0 && static_cast<std::size_t>(component) < result.componentDistribution.size())
    {
      result.componentDistribution[static_cast<std::size_t>(component)] += amount;
    }
    const int side = sideOf(component);
    if (side > 0)
      result.accessibleDistribution += amount;
    else if (side < 0)
      result.inaccessibleDistribution += amount;
    else
      result.undecidedDistribution += amount;
  }

  // The volume between the accessible surface and the excluded one over this piece.
  void addShell(std::int32_t component, double amount)
  {
    result.shellVolume += amount;
    if (component >= 0 && static_cast<std::size_t>(component) < result.componentShellVolume.size())
    {
      result.componentShellVolume[static_cast<std::size_t>(component)] += amount;
    }
  }

  // Area, of whichever of the three kinds this piece is.
  void addArea(std::int32_t component, double ExcludedSurfaceAreas::* kind, double amount)
  {
    result.area.*kind += amount;
    const int side = sideOf(component);
    if (side > 0)
      result.accessibleArea.*kind += amount;
    else if (side < 0)
      result.inaccessibleArea.*kind += amount;
    else
      result.undecidedArea.*kind += amount;
  }
};

}  // namespace


SolventExcludedGeometry solventExcludedGeometry(const PoreAccessibility& accessibility, double probeRadius,
                                                const BoundaryComponents& components,
                                                const std::vector<ComponentVerdict>& verdicts, std::size_t subdivisions)
{
  // One sweep of the spheres serves four purposes: the area of each atom's own patch, which the shell wants; the
  // radius weighted area, which the volume of the union wants; the same patches on the bare spheres, which are
  // the convex part of the excluded surface; and the moments of each surface, which say which side it faces.
  return solventExcludedGeometry(
      accessibility, probeRadius, components, verdicts,
      exactAccessibleSurfaceAreaByComponent(accessibility, components, verdicts, subdivisions), subdivisions);
}


SolventExcludedGeometry solventExcludedGeometry(const PoreAccessibility& accessibility, double probeRadius,
                                                const BoundaryComponents& components,
                                                const std::vector<ComponentVerdict>& verdicts,
                                                const MeasuredPatches& patches, std::size_t subdivisions)
{
  SolventExcludedGeometry result;
  result.probeRadius = probeRadius;
  result.accessibleVolume = unionOfBallsVolume(accessibility, patches, components);

  // The same division of the surfaces the area was split by, so that a volume and an area of this structure
  // never disagree about which side a surface is on.
  const std::vector<SurfaceSide> sides = surfaceSides(components, patches, verdicts);
  result.componentSide.assign(components.numberOfComponents, 0);
  for (std::size_t component = 0; component < components.numberOfComponents; ++component)
  {
    result.componentSide[component] = sides[component].side;
  }
  result.componentDistribution.assign(components.numberOfComponents, 0.0);
  result.componentEnclosedVolume.assign(components.numberOfComponents, 0.0);

  // The shell over each surface's own patches, which the sweep accumulated patch by patch. The arcs and the
  // vertices add theirs below, where they are visited anyway.
  result.componentShellVolume.assign(components.numberOfComponents, 0.0);
  const std::size_t measured = std::min(components.numberOfComponents, patches.components.size());
  for (std::size_t component = 0; component < measured; ++component)
  {
    result.componentShellVolume[component] = patches.components[component].shellVolume;
  }

  SurfaceTotals totals{result, sides};

  // The patches. Each contributes the shell between the bare sphere and the inflated one over its own solid
  // angle, and nothing to the derivative: the bare spheres do not move as the probe grows. Its area on the
  // excluded surface is its own accessible area brought down onto the bare sphere, which the sweep has already
  // done atom by atom and surface by surface.
  //
  // Summed over the atoms rather than over the arcs. The bare radius is the one the sweep itself carried the
  // patches down by, taken from `accessibility` rather than from `probeRadius`: the two say the same thing
  // for any caller that inflated these atoms by this probe, and reading it from the same place as the sweep
  // is what leaves this total and the per-surface ones the sweep accumulated in step by construction.
  double patchShell = 0.0;
  for (std::size_t i = 0; i < accessibility.atomPositions.size(); ++i)
  {
    const double inflated = accessibility.atomRadii[i];
    if (inflated <= 0.0 || i >= patches.atomArea.size()) continue;
    const double ratio = accessibility.bareRadius(i) / inflated;
    patchShell += patches.atomArea[i] * inflated / 3.0 * (1.0 - ratio * ratio * ratio);
  }
  result.shellVolume = patchShell;
  for (std::size_t component = 0; component < patches.components.size(); ++component)
  {
    totals.addArea(static_cast<std::int32_t>(component), &ExcludedSurfaceAreas::convex,
                   patches.components[component].convexArea);
  }

  // The arcs, which carry the toroidal patches.
  const std::vector<ExposedArc> arcs = exposedArcs(accessibility, probeRadius, components);
  result.diagnostics.numberOfArcs = arcs.size();
  for (const ExposedArc& arc : arcs)
  {
    if (arc.radius <= probeRadius) ++result.diagnostics.cuspedArcs;
    totals.addShell(arc.component, shellArcVolume(arc, probeRadius));
    double speed = toroidalNormalSpeed(arc, probeRadius);
    result.toroidalDistribution += speed;
    totals.addDistribution(arc.component, speed);
    totals.addArea(arc.component, &ExcludedSurfaceAreas::saddle, toroidalPatchArea(arc, probeRadius));
  }

  // The vertices, which carry the concave patches.
  const SurfaceVertices found = surfaceVertices(accessibility, components);
  const std::vector<SasVertex>& vertices = found.vertices;
  result.diagnostics.degenerateVertices = found.onManySpheres;
  result.diagnostics.discardedCorners = found.discarded;
  result.diagnostics.numberOfVertices = vertices.size();
  const VertexGrid grid = buildVertexGrid(accessibility, vertices);

  for (const SasVertex& vertex : vertices)
  {
    std::vector<double3> neighbours = grid.within(vertex.position, 2.0 * probeRadius);
    VertexTerms terms = vertexTerms(vertex, probeRadius, neighbours, subdivisions);
    if (terms.clipped) ++result.diagnostics.clippedVertices;
    if (terms.vanished) ++result.diagnostics.vanishedVertices;

    totals.addShell(vertex.component, terms.shell);
    result.concaveDistribution += terms.normalSpeed;
    totals.addDistribution(vertex.component, terms.normalSpeed);
    totals.addArea(vertex.component, &ExcludedSurfaceAreas::concave, terms.area);
  }

  result.excludedVolume = result.accessibleVolume - result.shellVolume;
  result.poreVolume = accessibility.simulationBox.volume - result.excludedVolume;

  // What each bounded surface holds, and with it the division of the pore volume. A surface that closes
  // encloses room for the probe's centre, which the divergence theorem gives from the same arcs, and the probe
  // reaches out from that room as far as the shell over the surface: the two together are the volume this
  // probe opens up behind that surface. A surface that runs away encloses nothing and is left at zero, its
  // pore's volume being what the sealed ones leave of the total.
  for (std::size_t component = 0; component < components.numberOfComponents; ++component)
  {
    if (component >= patches.components.size()) continue;
    const BoundaryMoments& moments = patches.components[component];
    if (moments.area <= 0.0 || components.componentPercolates[component] != 0) continue;

    const double enclosed = -(moments.radiusWeightedArea + moments.originWeighted) / 3.0;
    result.componentEnclosedVolume[component] = enclosed;

    if (sides[component].side < 0)
      result.inaccessiblePoreVolume += enclosed + result.componentShellVolume[component];
    else if (sides[component].side == 0)
      result.undecidedPoreVolume += enclosed + result.componentShellVolume[component];
  }
  result.accessiblePoreVolume = result.poreVolume - result.inaccessiblePoreVolume - result.undecidedPoreVolume;
  return result;
}


SolventExcludedGeometry solventExcludedGeometry(const PoreAccessibility& accessibility, double probeRadius,
                                                std::size_t subdivisions)
{
  BoundaryComponents components = boundaryComponents(accessibility);
  std::vector<ComponentVerdict> verdicts = boundaryComponentVerdicts(accessibility, components);
  return solventExcludedGeometry(accessibility, probeRadius, components, verdicts, subdivisions);
}
