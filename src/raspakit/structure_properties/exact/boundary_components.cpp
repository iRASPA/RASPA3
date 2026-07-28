module;

module exact_boundary_components;

import std;

import int3;
import double3;
import simulationbox;
import voronoi_accessibility;
import voronoi_channels;

namespace
{

// A crossing is only hidden by a third sphere if it is hidden with room to spare. A framework is
// symmetric, so three spheres meeting in one point is the ordinary case rather than a coincidence, and
// deciding on the sign of a rounding error whether the third one hides the crossing loses an edge of the
// patch. The same allowance is made in the sweep, and for the same reason.
constexpr double coverTolerance = 1.0e-9;

// Crossings closer together than this along a circle are the same crossing. Several circles through one
// point is again the ordinary case, and left alone it would leave arcs of no length between them.
constexpr double angleTolerance = 1.0e-9;

// How far in from an edge a patch is represented. Far enough that the point is strictly inside and the
// arithmetic below has a margin, near enough that it cannot have left the patch.
constexpr double insetAngle = 1.0e-6;


bool coveredDirection(const std::vector<SphereCircle>& circles, const double3& direction, std::size_t skipFirst,
                      std::size_t skipSecond)
{
  for (std::size_t c = 0; c < circles.size(); ++c)
  {
    if (c == skipFirst || c == skipSecond) continue;
    if (double3::dot(direction, circles[c].axis) > circles[c].cosineHalfAngle + coverTolerance) return true;
  }
  return false;
}


// A unit vector perpendicular to `axis`, chosen so that the cross product behind it is well conditioned.
double3 perpendicularTo(const double3& axis)
{
  double3 helper(1.0, 0.0, 0.0);
  if (std::abs(axis.y) < std::abs(axis.x)) helper = double3(0.0, 1.0, 0.0);
  if (std::abs(axis.z) < std::min(std::abs(axis.x), std::abs(axis.y))) helper = double3(0.0, 0.0, 1.0);
  double3 perpendicular = double3::cross(helper, axis);
  return perpendicular * (1.0 / perpendicular.length());
}


// A little way into the exposed region from the point of `circle` at angle `t`: the same point rotated
// away from the circle's axis, which is the direction the region lies in.
double3 insetFromCircle(const SphereCircle& circle, double t)
{
  double3 outward = circle.first * std::cos(t) + circle.second * std::sin(t);
  double angle = std::min(circle.halfAngle + insetAngle, std::numbers::pi);
  return circle.axis * std::cos(angle) + outward * std::sin(angle);
}


struct UnionFind
{
  std::vector<std::size_t> parent;

  explicit UnionFind(std::size_t count) : parent(count)
  {
    for (std::size_t i = 0; i < count; ++i) parent[i] = i;
  }

  std::size_t find(std::size_t i)
  {
    while (parent[i] != i)
    {
      parent[i] = parent[parent[i]];
      i = parent[i];
    }
    return i;
  }

  void join(std::size_t i, std::size_t j)
  {
    std::size_t a = find(i);
    std::size_t b = find(j);
    if (a != b) parent[b] = a;
  }
};


// The circles cutting one atom's sphere, with the neighbour each belongs to. Returns false when some
// neighbour swallows the sphere whole, in which case the atom carries no surface at all.
bool boundingCircles(const VoronoiAccessibility& accessibility, std::size_t atomIndex,
                     std::vector<SphereCircle>& circles)
{
  const double radius = accessibility.atomRadii[atomIndex];
  const double3 centre = accessibility.atomPositions[atomIndex];

  circles.clear();
  for (const NeighbourImage& neighbour :
       accessibility.neighbourAtomImages(centre, radius + accessibility.maximumAtomRadius))
  {
    double distance = neighbour.delta.length();
    if (distance < 1.0e-12)
    {
      // This sphere itself, or another sitting exactly on it: only a strictly larger one covers anything,
      // and then it covers all of it.
      if (neighbour.radius > radius) return false;
      continue;
    }

    double cosineHalfAngle =
        (radius * radius + distance * distance - neighbour.radius * neighbour.radius) / (2.0 * radius * distance);
    if (cosineHalfAngle >= 1.0) continue;   // reaches nothing of this sphere
    if (cosineHalfAngle <= -1.0) return false;  // swallows it

    SphereCircle circle;
    circle.axis = neighbour.delta * (1.0 / distance);
    circle.cosineHalfAngle = cosineHalfAngle;
    circle.halfAngle = std::acos(cosineHalfAngle);
    circle.sineHalfAngle = std::sin(circle.halfAngle);
    circle.neighbourIndex = neighbour.index;
    circle.neighbourImage = neighbour.image;
    circle.first = perpendicularTo(circle.axis);
    circle.second = double3::cross(circle.axis, circle.first);
    circles.push_back(circle);
  }

  // A circle whose disc lies inside another's bounds nothing and covers nothing of its own. Dropping it is
  // no loss of an edge: the surface the two atoms share is buried under the third, on both of their spheres
  // alike, so there is no piece of boundary there to join anything across.
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
  }
  return true;
}


// One uncovered crossing of two of the circles: a corner of the patches that meet there.
struct Crossing
{
  std::size_t firstCircle{0};
  std::size_t secondCircle{0};
  double3 direction{0.0, 0.0, 0.0};
};


std::vector<Crossing> uncoveredCrossings(const std::vector<SphereCircle>& circles)
{
  std::vector<Crossing> crossings;
  for (std::size_t j = 0; j + 1 < circles.size(); ++j)
  {
    for (std::size_t k = j + 1; k < circles.size(); ++k)
    {
      const SphereCircle& first = circles[j];
      const SphereCircle& second = circles[k];
      double alignment = double3::dot(first.axis, second.axis);
      double denominator = 1.0 - alignment * alignment;
      if (denominator < 1.0e-14) continue;  // parallel axes never cross

      double alongFirst = (first.cosineHalfAngle - alignment * second.cosineHalfAngle) / denominator;
      double alongSecond = (second.cosineHalfAngle - alignment * first.cosineHalfAngle) / denominator;
      double outOfPlaneSquared =
          (1.0 - alongFirst * first.cosineHalfAngle - alongSecond * second.cosineHalfAngle) / denominator;
      if (outOfPlaneSquared <= 0.0) continue;  // the circles miss one another

      double3 inPlane = first.axis * alongFirst + second.axis * alongSecond;
      double3 outOfPlane = double3::cross(first.axis, second.axis) * std::sqrt(outOfPlaneSquared);
      for (std::size_t side = 0; side < 2; ++side)
      {
        double3 direction = (side == 0) ? inPlane + outOfPlane : inPlane - outOfPlane;
        direction = direction * (1.0 / direction.length());
        if (coveredDirection(circles, direction, j, k)) continue;
        crossings.push_back(Crossing{j, k, direction});
      }
    }
  }
  return crossings;
}


// Where a corner sits in a circle's own list of corner angles, the angles having been merged and sorted.
std::optional<std::size_t> cornerIndex(const SphereCircle& circle, const double3& direction)
{
  double angle = circle.angleOf(direction);
  for (std::size_t i = 0; i < circle.cornerAngles.size(); ++i)
  {
    double difference = std::abs(circle.cornerAngles[i] - angle);
    difference = std::min(difference, 2.0 * std::numbers::pi - difference);
    if (difference < 1.0e-7) return i;
  }
  return std::nullopt;
}


// The arc of `circle` holding the angle `t`.
std::size_t arcHolding(const SphereCircle& circle, double t)
{
  if (circle.cornerAngles.empty()) return 0;
  for (std::size_t k = 0; k < circle.cornerAngles.size(); ++k)
  {
    double begin = circle.cornerAngles[k];
    double end = circle.cornerAngles[(k + 1) % circle.cornerAngles.size()];
    double width = end - begin;
    if (width <= 0.0) width += 2.0 * std::numbers::pi;
    double offset = t - begin;
    if (offset < 0.0) offset += 2.0 * std::numbers::pi;
    if (offset <= width) return k;
  }
  return 0;
}


// How far a direction is from being covered: the least room it has to spare against any of the circles.
// Positive is exposed, and the larger it is the further into the open the direction points.
double exposureMargin(const std::vector<SphereCircle>& circles, const double3& direction)
{
  double margin = std::numeric_limits<double>::max();
  for (const SphereCircle& circle : circles)
  {
    margin = std::min(margin, circle.cosineHalfAngle - double3::dot(direction, circle.axis));
  }
  return margin;
}


// Directions in the exposed region to try as a way point between two of its points. Away from each
// neighbour is where the room is, so the reversed axes are the natural candidates, together with the
// directions between pairs of them and the coordinate axes for the cases where neither serves. Only the
// most open of them are kept: they are tried in turn, and the cost of that has to stay bounded.
std::vector<double3> exposedWayPoints(const std::vector<SphereCircle>& circles)
{
  constexpr std::size_t limit = 16;

  std::vector<double3> candidates;
  candidates.reserve(circles.size() * 2 + 6);
  for (const SphereCircle& circle : circles) candidates.push_back(circle.axis * -1.0);
  for (std::size_t i = 0; i < circles.size(); ++i)
  {
    for (std::size_t j = i + 1; j < circles.size(); ++j)
    {
      double3 between = circles[i].axis * -1.0 + circles[j].axis * -1.0;
      double length = between.length();
      if (length > 1.0e-6) candidates.push_back(between * (1.0 / length));
    }
  }
  candidates.push_back(double3(1.0, 0.0, 0.0));
  candidates.push_back(double3(-1.0, 0.0, 0.0));
  candidates.push_back(double3(0.0, 1.0, 0.0));
  candidates.push_back(double3(0.0, -1.0, 0.0));
  candidates.push_back(double3(0.0, 0.0, 1.0));
  candidates.push_back(double3(0.0, 0.0, -1.0));

  std::vector<std::pair<double, double3>> ranked;
  for (const double3& candidate : candidates)
  {
    double margin = exposureMargin(circles, candidate);
    if (margin > coverTolerance) ranked.emplace_back(margin, candidate);
  }
  std::sort(ranked.begin(), ranked.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
  if (ranked.size() > limit) ranked.resize(limit);

  std::vector<double3> wayPoints;
  wayPoints.reserve(ranked.size());
  for (const auto& [margin, direction] : ranked) wayPoints.push_back(direction);
  return wayPoints;
}


// Going from `direction` along the great circle toward the axis of `circles[reference]`, the point where
// the first of the circles is crossed.
//
// That point is on the edge of the patch holding `direction`, and it is the whole of what a lookup needs:
// the path up to it crosses nothing, so it stays in the exposed region, so the patch it lands on the edge of
// is the one it started in. There is always such a point, the axis of a circle being inside its own cap, and
// finding it is a root of a sinusoid per circle rather than a search over paths.
std::optional<double3> firstCrossingToward(const std::vector<SphereCircle>& circles, const double3& direction,
                                          std::size_t reference)
{
  const double3 axis = circles[reference].axis;
  double alignment = double3::dot(direction, axis);
  double3 tangent = axis - direction * alignment;
  double length = tangent.length();
  if (length < 1.0e-9) return std::nullopt;  // the direction is that axis: this circle says nothing
  tangent = tangent * (1.0 / length);

  // Along the path the height above any axis is `a cos t + b sin t`, so a crossing is where that sinusoid
  // reaches the circle's own height. The path is followed only as far as the reference axis, which is
  // covered, so the first crossing is at or before it.
  const double limit = std::acos(std::clamp(alignment, -1.0, 1.0));
  double first = limit;
  for (const SphereCircle& circle : circles)
  {
    double a = double3::dot(direction, circle.axis);
    double b = double3::dot(tangent, circle.axis);
    double amplitude = std::hypot(a, b);
    if (amplitude < circle.cosineHalfAngle) continue;  // this sinusoid never reaches that height

    double phase = std::atan2(b, a);
    double half = std::acos(std::clamp(circle.cosineHalfAngle / amplitude, -1.0, 1.0));
    for (double root : {phase - half, phase + half, phase - half + 2.0 * std::numbers::pi,
                        phase + half + 2.0 * std::numbers::pi})
    {
      if (root > 1.0e-12 && root < first) first = root;
    }
  }
  return direction * std::cos(first) + tangent * std::sin(first);
}


// The way points of one sphere, grouped by which of them a single arc joins to which. A group is then a
// piece of the exposed region reachable within itself, and two points that both reach one group are joined
// through it however many legs that takes.
void groupWayPoints(SphereBoundary& boundary)
{
  const std::size_t count = boundary.wayPoints.size();
  boundary.wayPointGroup.assign(count, -1);

  UnionFind together(count);
  for (std::size_t i = 0; i < count; ++i)
  {
    for (std::size_t j = i + 1; j < count; ++j)
    {
      if (exposedGreatCircleArc(boundary.circles, boundary.wayPoints[i], boundary.wayPoints[j]))
      {
        together.join(i, j);
      }
    }
  }

  std::unordered_map<std::size_t, std::int32_t> numbering;
  for (std::size_t i = 0; i < count; ++i)
  {
    auto [entry, inserted] = numbering.try_emplace(together.find(i), static_cast<std::int32_t>(numbering.size()));
    boundary.wayPointGroup[i] = entry->second;
  }
}


// Whether the point `x` of the exposed region lies on the interior side of a loop, in the sense of Quan and
// Stamm: the side the exposed surface the loop bounds is on.
//
// Among the circles carrying the loop's arcs, take the one whose circle is nearest to `x` along the sphere,
// and the point of it the shortest path lands on. That path crosses none of those circles --- crossing one
// would have reached it sooner --- so `x` and the landing point are on the same side of the loop, the loop
// being made of arcs of those very circles. So it is enough to see where the path lands: on an arc of the
// loop is the interior side, on a stretch of the same circle that some other sphere covers is the exterior.
// Nothing is searched for and nothing can fail to be found.
bool insideLoop(const SphereBoundary& boundary, const std::vector<std::size_t>& arcBase,
                const std::vector<std::int32_t>& loopOfArc, const std::vector<std::size_t>& loopCircles,
                std::int32_t loop, const double3& x)
{
  double nearest = std::numeric_limits<double>::max();
  std::size_t landedOn = 0;
  double3 landing(0.0, 0.0, 0.0);

  for (std::size_t c : loopCircles)
  {
    const SphereCircle& circle = boundary.circles[c];
    double alignment = std::clamp(double3::dot(x, circle.axis), -1.0, 1.0);
    double distance = std::abs(std::acos(alignment) - circle.halfAngle);
    if (distance >= nearest) continue;

    double3 perpendicular = x - circle.axis * alignment;
    double length = perpendicular.length();
    if (length < 1.0e-12) continue;  // `x` is on that axis, where every point of the circle is as near

    perpendicular = perpendicular * (1.0 / length);
    nearest = distance;
    landedOn = c;
    landing = circle.axis * circle.cosineHalfAngle + perpendicular * circle.sineHalfAngle;
  }
  if (nearest == std::numeric_limits<double>::max()) return false;

  const SphereCircle& circle = boundary.circles[landedOn];
  std::size_t arc = arcHolding(circle, circle.angleOf(landing));
  return loopOfArc[arcBase[landedOn] + arc] == loop;
}


// The middle of arc `k` of `circle`, as an angle.
double arcMidAngle(const SphereCircle& circle, std::size_t k)
{
  if (circle.cornerAngles.empty()) return 0.0;
  double begin = circle.cornerAngles[k];
  double end = circle.cornerAngles[(k + 1) % circle.cornerAngles.size()];
  double width = end - begin;
  if (width <= 0.0) width += 2.0 * std::numbers::pi;
  return begin + 0.5 * width;
}

}  // namespace


bool exposedGreatCircleArc(const std::vector<SphereCircle>& circles, const double3& from, const double3& to)
{
  double alignment = std::clamp(double3::dot(from, to), -1.0, 1.0);
  double span = std::acos(alignment);
  if (span < 1.0e-12) return true;

  double3 perpendicular = to - from * alignment;
  double length = perpendicular.length();
  if (length < 1.0e-14) return false;  // antipodal: no arc is determined
  perpendicular = perpendicular * (1.0 / length);

  for (const SphereCircle& circle : circles)
  {
    // Along the arc the height above the circle's axis is `atStart cos t + atPerpendicular sin t`, one
    // sinusoid, so its largest value over the arc is at an end or at the single interior turning point.
    double atStart = double3::dot(from, circle.axis);
    double atPerpendicular = double3::dot(perpendicular, circle.axis);
    double highest = std::max(atStart, atStart * std::cos(span) + atPerpendicular * std::sin(span));
    double phase = std::atan2(atPerpendicular, atStart);
    if (phase >= 0.0 && phase <= span)
    {
      highest = std::max(highest, std::sqrt(atStart * atStart + atPerpendicular * atPerpendicular));
    }
    if (highest > circle.cosineHalfAngle) return false;
  }
  return true;
}


bool connectedOnSphere(const SphereBoundary& boundary, const double3& from, const double3& to)
{
  if (exposedGreatCircleArc(boundary.circles, from, to)) return true;

  // Otherwise through the way points, which have already been grouped by the arcs among themselves: a group
  // both ends reach is a path between them of as many legs as that group takes.
  std::vector<std::uint8_t> fromEnd, toEnd;
  for (std::size_t w = 0; w < boundary.wayPoints.size(); ++w)
  {
    std::size_t group = static_cast<std::size_t>(boundary.wayPointGroup[w]);
    if (group >= fromEnd.size())
    {
      fromEnd.resize(group + 1, 0);
      toEnd.resize(group + 1, 0);
    }
    if (exposedGreatCircleArc(boundary.circles, from, boundary.wayPoints[w])) fromEnd[group] = 1;
    if (exposedGreatCircleArc(boundary.circles, boundary.wayPoints[w], to)) toEnd[group] = 1;
  }
  for (std::size_t group = 0; group < fromEnd.size(); ++group)
  {
    if (fromEnd[group] != 0 && toEnd[group] != 0) return true;
  }
  return false;
}


std::int32_t BoundaryComponents::patchOfDirection(std::size_t atomIndex, const double3& unitDirection) const
{
  const SphereBoundary& boundary = atoms[atomIndex];
  if (boundary.buried) return -1;
  if (boundary.numberOfPatches == 1) return 0;
  if (boundary.circles.empty()) return -1;

  // Walk to the edge of the patch and read the patch off there. The circles are tried as the direction to
  // walk in from the nearest outward, which is where the edge is likely to be closest; each gives an edge
  // point of the same patch, so the first that yields one is as good as any.
  std::vector<std::pair<double, std::size_t>> byNearest;
  byNearest.reserve(boundary.circles.size());
  for (std::size_t c = 0; c < boundary.circles.size(); ++c)
  {
    byNearest.emplace_back(-double3::dot(unitDirection, boundary.circles[c].axis), c);
  }
  std::sort(byNearest.begin(), byNearest.end());

  for (const auto& [nearness, c] : byNearest)
  {
    std::optional<double3> edge = firstCrossingToward(boundary.circles, unitDirection, c);
    if (!edge.has_value()) continue;
    std::int32_t patch = patchOfCirclePoint(atomIndex, edge.value());
    if (patch >= 0) return patch;
  }

  // Nothing was found on an edge. What is left is to join the direction to a patch by an arc of a great
  // circle, which proves the same thing where it works and is all that is left where the walk did not.
  for (std::size_t patch = 0; patch < boundary.patchRepresentative.size(); ++patch)
  {
    if (connectedOnSphere(boundary, unitDirection, boundary.patchRepresentative[patch]))
    {
      return static_cast<std::int32_t>(patch);
    }
  }
  return -1;
}


std::int32_t BoundaryComponents::patchOfCirclePoint(std::size_t atomIndex, const double3& unitDirection) const
{
  const SphereBoundary& boundary = atoms[atomIndex];
  if (boundary.buried) return -1;
  if (boundary.circles.empty() || boundary.numberOfPatches == 1) return boundary.numberOfPatches == 1 ? 0 : -1;

  // Which circle the point is on. Several may pass through it, so they are tried nearest first: the point
  // comes from the sweep's own arithmetic rather than from the circle's, and agrees with it only to
  // rounding.
  std::vector<std::pair<double, std::size_t>> onCircle;
  for (std::size_t c = 0; c < boundary.circles.size(); ++c)
  {
    double residual = std::abs(double3::dot(unitDirection, boundary.circles[c].axis) -
                               boundary.circles[c].cosineHalfAngle);
    if (residual < 1.0e-6) onCircle.emplace_back(residual, c);
  }
  std::sort(onCircle.begin(), onCircle.end());

  for (const auto& [residual, c] : onCircle)
  {
    const SphereCircle& circle = boundary.circles[c];
    if (circle.arcPatch.empty()) continue;
    std::size_t arc = arcHolding(circle, circle.angleOf(unitDirection));

    // The point is a corner of the exposed region as often as not, and there either of the two arcs
    // meeting there may be the one holding it. Whichever is exposed answers, both being edges of the same
    // patch where both are.
    const std::size_t count = circle.arcPatch.size();
    for (std::size_t attempt = 0; attempt < 3; ++attempt)
    {
      std::size_t k = (attempt == 0) ? arc : ((attempt == 1) ? (arc + count - 1) % count : (arc + 1) % count);
      if (circle.arcPatch[k] >= 0) return circle.arcPatch[k];
    }
  }
  return -1;
}


std::vector<ComponentLabel> labelBoundaryComponents(const VoronoiAccessibility& accessibility,
                                                    const BoundaryComponents& components)
{
  // The walk out from the surface. Each step is kept inside the free ball about the point it starts from,
  // so the whole of it is void: a node ball reached this way holds a point connected to the surface through
  // the void, and that settles which pore the surface faces without appealing to any nearness heuristic.
  constexpr double firstStep = 1.0e-3;   // to get off the surface, where the clearance is zero
  constexpr double stepFraction = 0.9;   // of the room at each point, leaving the step strictly inside
  constexpr std::size_t stepLimit = 64;
  constexpr double reachLimit = 32.0;    // no cavity in a framework is wider than this

  std::vector<ComponentLabel> labels(components.numberOfComponents);
  for (std::size_t label = 0; label < components.numberOfComponents; ++label)
  {
    ComponentLabel& answer = labels[label];
    double3 bestSurface(0.0, 0.0, 0.0);
    bool haveSurface = false;

    for (const auto& [atomIndex, patch] : components.componentCandidates[label])
    {
      const SphereBoundary& boundary = components.atoms[atomIndex];
      if (patch >= boundary.patchRepresentative.size()) continue;

      const double3 outward = boundary.patchRepresentative[patch];
      const double3 surface = accessibility.atomPositions[atomIndex] + outward * accessibility.atomRadii[atomIndex];
      if (!haveSurface)
      {
        bestSurface = surface;
        haveSurface = true;
      }

      double walked = firstStep;
      for (std::size_t step = 0; step < stepLimit && walked < reachLimit; ++step)
      {
        double3 point = surface + outward * walked;
        double room = accessibility.clearance(point);
        if (room <= 0.0) break;  // the ray has run into another atom; the walk cannot be continued

        answer.steps = step + 1;
        answer.walked = walked;
        if (std::optional<std::size_t> holder = accessibility.containingNode(point); holder.has_value())
        {
          std::size_t node = holder.value();
          answer.decided = true;
          answer.proved = true;
          answer.poreId = accessibility.channels.nodePoreId[node];
          answer.accessible = accessibility.nodeAccessible[node] != 0;
          break;
        }
        walked += stepFraction * room;
      }
      if (answer.proved) break;
    }

    if (answer.proved || !haveSurface) continue;

    // No ball was reached from any of them: fall back on what the classifier makes of a point on the
    // surface. This is the per-arc route, but taken once for the whole surface, so it can still be wrong
    // about which pore and no longer wrong about only part of a pocket's boundary.
    PointClassification classification = accessibility.classify(bestSurface);
    if (classification.inside || classification.resample) continue;
    answer.decided = true;
    answer.accessible = classification.accessible;
    answer.poreId = classification.poreId;
  }
  return labels;
}


BoundaryComponents boundaryComponents(const VoronoiAccessibility& accessibility, LoopMerge rule)
{
  const std::size_t numberOfAtoms = accessibility.atomPositions.size();
  BoundaryComponents result;
  result.atoms.assign(numberOfAtoms, SphereBoundary{});

  // ---- the patches of each sphere ------------------------------------------------------------------
  for (std::size_t atomIndex = 0; atomIndex < numberOfAtoms; ++atomIndex)
  {
    SphereBoundary& boundary = result.atoms[atomIndex];
    if (!boundingCircles(accessibility, atomIndex, boundary.circles))
    {
      boundary.buried = true;
      continue;
    }

    // A sphere no neighbour reaches is one patch with no edges at all.
    if (boundary.circles.empty())
    {
      boundary.numberOfPatches = 1;
      boundary.patchRepresentative.push_back(double3(0.0, 0.0, 1.0));
      continue;
    }

    std::vector<Crossing> crossings = uncoveredCrossings(boundary.circles);

    // The crossings on each circle, merged where they coincide and sorted, cut it into arcs.
    for (SphereCircle& circle : boundary.circles) circle.cornerAngles.clear();
    for (const Crossing& crossing : crossings)
    {
      boundary.circles[crossing.firstCircle].cornerAngles.push_back(
          boundary.circles[crossing.firstCircle].angleOf(crossing.direction));
      boundary.circles[crossing.secondCircle].cornerAngles.push_back(
          boundary.circles[crossing.secondCircle].angleOf(crossing.direction));
    }
    for (SphereCircle& circle : boundary.circles)
    {
      std::sort(circle.cornerAngles.begin(), circle.cornerAngles.end());
      circle.cornerAngles.erase(std::unique(circle.cornerAngles.begin(), circle.cornerAngles.end(),
                                            [](double a, double b) { return std::abs(a - b) < angleTolerance; }),
                                circle.cornerAngles.end());
      // The first and last can also be the same crossing, the circle closing on itself.
      if (circle.cornerAngles.size() > 1 &&
          2.0 * std::numbers::pi - (circle.cornerAngles.back() - circle.cornerAngles.front()) < angleTolerance)
      {
        circle.cornerAngles.pop_back();
      }
    }

    // Which arcs are exposed, judged in the middle: between two consecutive crossings the answer cannot
    // change, a crossing being the only place where the boundary of the covered part can cross the circle.
    std::vector<std::size_t> arcBase(boundary.circles.size(), 0);
    std::size_t totalArcs = 0;
    for (std::size_t c = 0; c < boundary.circles.size(); ++c)
    {
      arcBase[c] = totalArcs;
      totalArcs += boundary.circles[c].numberOfArcs();
    }

    std::vector<std::uint8_t> exposed(totalArcs, 0);
    for (std::size_t c = 0; c < boundary.circles.size(); ++c)
    {
      SphereCircle& circle = boundary.circles[c];
      circle.arcPatch.assign(circle.numberOfArcs(), -1);
      for (std::size_t k = 0; k < circle.numberOfArcs(); ++k)
      {
        double3 middle = circle.direction(arcMidAngle(circle, k));
        if (!coveredDirection(boundary.circles, middle, c, c)) exposed[arcBase[c] + k] = 1;
      }
    }

    // At an uncovered crossing the exposed region takes exactly one of the four quadrants, so of the two
    // arcs of each circle that end there exactly one is exposed, and those two bound the same patch.
    UnionFind patches(totalArcs);
    for (const Crossing& crossing : crossings)
    {
      std::array<std::optional<std::size_t>, 2> incident;
      const std::array<std::size_t, 2> circleIndices = {crossing.firstCircle, crossing.secondCircle};
      for (std::size_t side = 0; side < 2; ++side)
      {
        const SphereCircle& circle = boundary.circles[circleIndices[side]];
        std::optional<std::size_t> corner = cornerIndex(circle, crossing.direction);
        if (!corner.has_value()) continue;

        std::size_t after = corner.value();
        std::size_t before = (after + circle.cornerAngles.size() - 1) % circle.cornerAngles.size();
        bool afterExposed = exposed[arcBase[circleIndices[side]] + after] != 0;
        bool beforeExposed = exposed[arcBase[circleIndices[side]] + before] != 0;
        if (afterExposed && !beforeExposed) incident[side] = arcBase[circleIndices[side]] + after;
        else if (beforeExposed && !afterExposed) incident[side] = arcBase[circleIndices[side]] + before;
      }
      if (incident[0].has_value() && incident[1].has_value())
      {
        patches.join(incident[0].value(), incident[1].value());
      }
    }

    // What the walk round the edges leaves separate is either a patch of its own or another loop of the
    // same patch, a patch with a hole in it having more than one. Which of the two it is is what `rule`
    // decides; where it cannot, the patch is left cut in two, which costs a second classification of the
    // same surface but cannot mislabel it.
    std::vector<std::size_t> roots;
    std::unordered_map<std::size_t, std::size_t> loopRepresentative;
    for (std::size_t arc = 0; arc < totalArcs; ++arc)
    {
      if (exposed[arc] == 0) continue;
      std::size_t root = patches.find(arc);
      if (loopRepresentative.emplace(root, arc).second) roots.push_back(root);
    }

    auto insetOfArc = [&](std::size_t arc)
    {
      std::size_t c =
          static_cast<std::size_t>(std::upper_bound(arcBase.begin(), arcBase.end(), arc) - arcBase.begin() - 1);
      return insetFromCircle(boundary.circles[c], arcMidAngle(boundary.circles[c], arc - arcBase[c]));
    };

    std::vector<double3> loopPoint;
    loopPoint.reserve(roots.size());
    for (std::size_t root : roots) loopPoint.push_back(insetOfArc(loopRepresentative[root]));

    auto circleOfArc = [&](std::size_t arc)
    {
      return static_cast<std::size_t>(std::upper_bound(arcBase.begin(), arcBase.end(), arc) - arcBase.begin() - 1);
    };

    const std::chrono::steady_clock::time_point mergeBegan = std::chrono::steady_clock::now();

    if (roots.size() > 1 && rule == LoopMerge::paths)
    {
      // A little way in from every exposed edge, along with the most open directions of the region. The edge
      // points are what a region that winds through a crevice needs: it may hold none of the open directions
      // at all, but it is bounded by its own edges everywhere, so consecutive points along them are a leg
      // apart. These are found only where there is something to merge, the walk round the edges having
      // already settled a region bounded by a single loop.
      boundary.wayPoints = exposedWayPoints(boundary.circles);
      for (std::size_t arc = 0; arc < totalArcs; ++arc)
      {
        if (exposed[arc] != 0) boundary.wayPoints.push_back(insetOfArc(arc));
      }
      groupWayPoints(boundary);

      for (std::size_t i = 0; i < roots.size(); ++i)
      {
        for (std::size_t j = i + 1; j < roots.size(); ++j)
        {
          if (patches.find(loopRepresentative[roots[i]]) == patches.find(loopRepresentative[roots[j]])) continue;
          if (connectedOnSphere(boundary, loopPoint[i], loopPoint[j]))
          {
            patches.join(loopRepresentative[roots[i]], loopRepresentative[roots[j]]);
          }
        }
      }
    }
    else if (roots.size() > 1)
    {
      // The loops as the walk round the edges left them, which is what the nesting test is about and has to
      // be taken down before any merging moves the roots.
      constexpr std::size_t sampleLimit = 8;

      std::unordered_map<std::size_t, std::int32_t> loopOfRoot;
      for (std::size_t i = 0; i < roots.size(); ++i) loopOfRoot.emplace(roots[i], static_cast<std::int32_t>(i));

      std::vector<std::int32_t> loopOfArc(totalArcs, -1);
      std::vector<std::vector<std::size_t>> loopCircles(roots.size());
      std::vector<std::vector<double3>> loopSamples(roots.size());
      for (std::size_t arc = 0; arc < totalArcs; ++arc)
      {
        if (exposed[arc] == 0) continue;
        std::int32_t loop = loopOfRoot.at(patches.find(arc));
        loopOfArc[arc] = loop;
        loopCircles[static_cast<std::size_t>(loop)].push_back(circleOfArc(arc));
        loopSamples[static_cast<std::size_t>(loop)].push_back(insetOfArc(arc));
      }
      for (std::size_t i = 0; i < roots.size(); ++i)
      {
        std::sort(loopCircles[i].begin(), loopCircles[i].end());
        loopCircles[i].erase(std::unique(loopCircles[i].begin(), loopCircles[i].end()), loopCircles[i].end());
      }

      // One point of a loop decides it, the two loops not crossing, except where they run along a circle they
      // share: there the nearest point of it is the asking loop's own arc rather than the asked loop's, and
      // the test says no to a pair it has nothing to say about. Asking at a few points round the loop rather
      // than at one costs a few more of a cheap test and leaves that case to the points elsewhere on it.
      auto anyPointInside = [&](std::int32_t loop, std::int32_t of)
      {
        const std::vector<double3>& samples = loopSamples[static_cast<std::size_t>(of)];
        const std::size_t stride = std::max<std::size_t>(1, samples.size() / sampleLimit);
        for (std::size_t s = 0; s < samples.size(); s += stride)
        {
          if (insideLoop(boundary, arcBase, loopOfArc, loopCircles[static_cast<std::size_t>(loop)], loop,
                         samples[s]))
          {
            return true;
          }
        }
        return false;
      };

      for (std::size_t i = 0; i < roots.size(); ++i)
      {
        for (std::size_t j = i + 1; j < roots.size(); ++j)
        {
          if (patches.find(loopRepresentative[roots[i]]) == patches.find(loopRepresentative[roots[j]])) continue;
          if (anyPointInside(static_cast<std::int32_t>(i), static_cast<std::int32_t>(j)) &&
              anyPointInside(static_cast<std::int32_t>(j), static_cast<std::int32_t>(i)))
          {
            patches.join(loopRepresentative[roots[i]], loopRepresentative[roots[j]]);
          }
        }
      }
    }

    result.mergeSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - mergeBegan).count();
    result.loopsToMerge += (roots.size() > 1) ? roots.size() : 0;

    // Number the patches and record where each is to be asked about.
    std::unordered_map<std::size_t, std::int32_t> patchOfRoot;
    for (std::size_t arc = 0; arc < totalArcs; ++arc)
    {
      if (exposed[arc] == 0) continue;
      std::size_t root = patches.find(arc);
      auto [entry, inserted] = patchOfRoot.try_emplace(root, static_cast<std::int32_t>(boundary.numberOfPatches));
      if (inserted)
      {
        ++boundary.numberOfPatches;
        std::size_t c = static_cast<std::size_t>(std::upper_bound(arcBase.begin(), arcBase.end(), arc) -
                                                 arcBase.begin() - 1);
        boundary.patchRepresentative.push_back(
            insetFromCircle(boundary.circles[c], arcMidAngle(boundary.circles[c], arc - arcBase[c])));
      }
      std::size_t c =
          static_cast<std::size_t>(std::upper_bound(arcBase.begin(), arcBase.end(), arc) - arcBase.begin() - 1);
      boundary.circles[c].arcPatch[arc - arcBase[c]] = entry->second;
    }
  }

  // ---- the graph of patches, joined across the circles they share ----------------------------------
  std::vector<std::size_t> patchBase(numberOfAtoms + 1, 0);
  for (std::size_t atomIndex = 0; atomIndex < numberOfAtoms; ++atomIndex)
  {
    patchBase[atomIndex + 1] = patchBase[atomIndex] + result.atoms[atomIndex].numberOfPatches;
  }
  const std::size_t totalPatches = patchBase[numberOfAtoms];
  result.numberOfPatches = totalPatches;

  struct Edge
  {
    std::size_t to{0};
    int3 offset{0, 0, 0};
  };
  std::vector<std::vector<Edge>> adjacency(totalPatches);

  for (std::size_t atomIndex = 0; atomIndex < numberOfAtoms; ++atomIndex)
  {
    const SphereBoundary& boundary = result.atoms[atomIndex];
    if (boundary.buried) continue;
    const double radius = accessibility.atomRadii[atomIndex];
    const double3 centre = accessibility.atomPositions[atomIndex];

    for (const SphereCircle& circle : boundary.circles)
    {
      const std::size_t neighbourIndex = circle.neighbourIndex;
      const int3 image = circle.neighbourImage;
      const SphereBoundary& other = result.atoms[neighbourIndex];
      if (other.buried) continue;

      // The same circle on the neighbour's own sphere: it sees this atom in the opposite image.
      const SphereCircle* mirror = nullptr;
      for (const SphereCircle& candidate : other.circles)
      {
        if (candidate.neighbourIndex == atomIndex && candidate.neighbourImage.x == -image.x &&
            candidate.neighbourImage.y == -image.y && candidate.neighbourImage.z == -image.z)
        {
          mirror = &candidate;
          break;
        }
      }

      for (std::size_t k = 0; k < circle.numberOfArcs(); ++k)
      {
        if (circle.arcPatch[k] < 0) continue;
        std::size_t here = patchBase[atomIndex] + static_cast<std::size_t>(circle.arcPatch[k]);

        if (mirror == nullptr)
        {
          ++result.looseEdges;
          continue;
        }

        // Where the middle of this arc is, seen from the neighbour's home sphere. The point is on both
        // spheres, so its direction there lands on the mirror circle, and the arc of that circle holding it
        // is the patch on the other side of this edge.
        double3 point = centre + circle.direction(arcMidAngle(circle, k)) * radius;
        double3 shifted = point - accessibility.simulationBox.cell * double3(static_cast<double>(image.x),
                                                                            static_cast<double>(image.y),
                                                                            static_cast<double>(image.z));
        double3 direction = shifted - accessibility.atomPositions[neighbourIndex];
        double length = direction.length();
        if (length < 1.0e-12)
        {
          ++result.looseEdges;
          continue;
        }
        direction = direction * (1.0 / length);

        std::size_t farArc = arcHolding(*mirror, mirror->angleOf(direction));
        if (mirror->arcPatch.empty() || mirror->arcPatch[farArc] < 0)
        {
          ++result.looseEdges;
          continue;
        }
        std::size_t there = patchBase[neighbourIndex] + static_cast<std::size_t>(mirror->arcPatch[farArc]);

        // The neighbour's home sphere carries this piece of surface a translation `image` away from where
        // this atom carries it, so bringing the two into one frame is that translation.
        adjacency[here].push_back(Edge{there, image});
        adjacency[there].push_back(Edge{here, int3(-image.x, -image.y, -image.z)});
      }
    }
  }

  // ---- the connected surfaces, and the translations around them ------------------------------------
  result.componentOfPatch.assign(numberOfAtoms, {});
  result.offsetOfPatch.assign(numberOfAtoms, {});
  for (std::size_t atomIndex = 0; atomIndex < numberOfAtoms; ++atomIndex)
  {
    result.componentOfPatch[atomIndex].assign(result.atoms[atomIndex].numberOfPatches, -1);
    result.offsetOfPatch[atomIndex].assign(result.atoms[atomIndex].numberOfPatches, int3(0, 0, 0));
  }

  auto atomOfPatch = [&](std::size_t patch)
  {
    return static_cast<std::size_t>(std::upper_bound(patchBase.begin(), patchBase.end(), patch) - patchBase.begin() -
                                    1);
  };

  std::vector<std::int32_t> component(totalPatches, -1);
  std::vector<int3> offset(totalPatches, int3(0, 0, 0));
  for (std::size_t seed = 0; seed < totalPatches; ++seed)
  {
    if (component[seed] >= 0) continue;
    std::int32_t label = static_cast<std::int32_t>(result.numberOfComponents++);
    result.componentPercolates.push_back(0);
    result.componentTranslations.emplace_back();

    // The translations this surface closes on itself by, kept down to a basis of them. A framework of any size
    // returns to an already lifted patch thousands of times, and all but the first few say the same thing, so
    // one is kept only when it says something the ones before it did not.
    std::vector<int3>& translations = result.componentTranslations.back();
    auto remember = [&translations](const int3& translation)
    {
      if (latticeVectorRank(translations) == 3) return;
      std::vector<int3> extended = translations;
      extended.push_back(translation);
      if (latticeVectorRank(extended) > latticeVectorRank(translations)) translations.push_back(translation);
    };

    component[seed] = label;
    offset[seed] = int3(0, 0, 0);
    std::vector<std::size_t> queue{seed};
    while (!queue.empty())
    {
      std::size_t patch = queue.back();
      queue.pop_back();
      for (const Edge& edge : adjacency[patch])
      {
        int3 carried(offset[patch].x + edge.offset.x, offset[patch].y + edge.offset.y,
                     offset[patch].z + edge.offset.z);
        if (component[edge.to] < 0)
        {
          component[edge.to] = label;
          offset[edge.to] = carried;
          queue.push_back(edge.to);
        }
        else if (offset[edge.to].x != carried.x || offset[edge.to].y != carried.y ||
                 offset[edge.to].z != carried.z)
        {
          // The surface has come back to a patch it had already reached, in a different copy of the cell: it
          // closes on a translate of itself and so runs away in that direction. The disagreement is that
          // translation, and over all such edges these generate every translation the surface has.
          result.componentPercolates[static_cast<std::size_t>(label)] = 1;
          remember(int3(carried.x - offset[edge.to].x, carried.y - offset[edge.to].y,
                        carried.z - offset[edge.to].z));
        }
      }
    }
  }

  result.componentDimensionality.reserve(result.numberOfComponents);
  for (const std::vector<int3>& translations : result.componentTranslations)
  {
    result.componentDimensionality.push_back(latticeVectorRank(translations));
  }

  for (std::size_t patch = 0; patch < totalPatches; ++patch)
  {
    std::size_t atomIndex = atomOfPatch(patch);
    std::size_t local = patch - patchBase[atomIndex];
    result.componentOfPatch[atomIndex][local] = component[patch];
    result.offsetOfPatch[atomIndex][local] = offset[patch];
    if (adjacency[patch].empty()) ++result.unjoinedPatches;
  }

  // ---- where to ask each surface which pore it faces ------------------------------------------------
  //
  // Any patch of the surface would answer for all of it, but not equally well: a patch in a crevice has
  // little room in front of it and the room is what a classification can be proved from. So the patch with
  // the most room is taken, measured by stepping out along the normal and asking how large a free ball fits.
  // A handful are kept rather than one: the room in front of a patch says how likely a free ball is to be
  // reachable from it, but not that it is, and trying the next costs one more walk.
  constexpr std::size_t candidateLimit = 6;

  std::vector<std::vector<std::pair<double, std::pair<std::size_t, std::size_t>>>> ranked(
      result.numberOfComponents);
  for (std::size_t atomIndex = 0; atomIndex < numberOfAtoms; ++atomIndex)
  {
    const SphereBoundary& boundary = result.atoms[atomIndex];
    for (std::size_t patch = 0; patch < boundary.numberOfPatches; ++patch)
    {
      std::int32_t label = result.componentOfPatch[atomIndex][patch];
      if (label < 0) continue;

      double3 outward = boundary.patchRepresentative[patch];
      double3 point = accessibility.atomPositions[atomIndex] + outward * accessibility.atomRadii[atomIndex];
      double room = 0.0;
      for (double step : {0.25, 0.5, 1.0, 2.0})
      {
        room = std::max(room, accessibility.clearance(point + outward * step));
      }

      auto& candidates = ranked[static_cast<std::size_t>(label)];
      candidates.emplace_back(room, std::pair<std::size_t, std::size_t>{atomIndex, patch});

      // Kept sorted and short, so that this stays linear in the patches rather than in their number times
      // their number.
      if (candidates.size() > 2 * candidateLimit)
      {
        std::partial_sort(candidates.begin(), candidates.begin() + candidateLimit, candidates.end(),
                          [](const auto& a, const auto& b) { return a.first > b.first; });
        candidates.resize(candidateLimit);
      }
    }
  }

  result.componentRepresentative.assign(result.numberOfComponents, {0, 0});
  result.componentCandidates.assign(result.numberOfComponents, {});
  for (std::size_t label = 0; label < result.numberOfComponents; ++label)
  {
    auto& candidates = ranked[label];
    std::sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
    if (candidates.size() > candidateLimit) candidates.resize(candidateLimit);
    for (const auto& [room, where] : candidates) result.componentCandidates[label].push_back(where);
    if (!candidates.empty()) result.componentRepresentative[label] = candidates.front().second;
  }

  return result;
}
