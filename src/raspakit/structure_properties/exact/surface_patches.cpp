module;

module exact_surface_patches;

import std;

import double3;
import int3;
import voronoi_channels;
import voronoi_accessibility;
import exact_boundary_components;

constexpr std::size_t gaussOrder = exactQuadratureOrder;

// Arcs shorter than this are dropped. They are the seams between two neighbouring circles that meet
// almost tangentially, carry no area worth the classifier call, and are the one place where the closed
// form for the arc loses digits.
constexpr double arcTolerance = 1.0e-12;


// One of the circles in which a neighbour cuts the sphere being measured, described by the direction of
// that neighbour and the polar half angle of the disc it covers.
//
// The two latitudes are where the circle turns back on itself: between them it crosses each circle of
// latitude in a genuine arc, outside them the latitude is either wholly inside the disc or wholly
// outside it. Which of the two it is depends only on whether the disc reaches over the pole on that
// side, so the pair of flags settles it. This is what lets a whole piece of the integral be skipped
// when some neighbour buries it, and lets the rest be integrated against a handful of circles rather
// than against all of them.
struct BoundingCircle
{
  double3 axis;
  double cosineHalfAngle{0.0};
  double halfAngle{0.0};
  double polarAngle{0.0};   // of the axis, in the frame the sphere is swept in
  double cosinePolar{0.0};
  double sinePolar{0.0};
  double azimuth{0.0};      // of the axis, in the same frame
  double lowestLatitude{0.0};
  double highestLatitude{0.0};
  bool reachesOverPole{false};
  bool reachesOverAntipole{false};
};


// Gauss-Legendre nodes and weights on the unit interval, found once by Newton's method on the Legendre
// polynomial with the usual Chebyshev starting guess.
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
      double abscissa = std::cos(std::numbers::pi * (static_cast<double>(i) + 0.75) /
                                 (static_cast<double>(gaussOrder) + 0.5));
      double derivative = 0.0;
      for (std::size_t iteration = 0; iteration < 100; ++iteration)
      {
        double previous = 1.0;
        double current = abscissa;
        for (std::size_t k = 2; k <= gaussOrder; ++k)
        {
          double next = ((2.0 * static_cast<double>(k) - 1.0) * abscissa * current -
                         (static_cast<double>(k) - 1.0) * previous) /
                        static_cast<double>(k);
          previous = current;
          current = next;
        }
        derivative = static_cast<double>(gaussOrder) * (abscissa * current - previous) /
                     (abscissa * abscissa - 1.0);
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


// A polar angle folded back into [0, pi], which is where the extreme latitudes of a circle live: a
// circle whose axis lies at polar angle beta turns back at |beta - theta| and at beta + theta, and the
// second of those reaches past the far pole and comes back when beta + theta exceeds pi.
double foldedPolarAngle(double angle)
{
  double wrapped = std::fmod(std::abs(angle), 2.0 * std::numbers::pi);
  return (wrapped > std::numbers::pi) ? 2.0 * std::numbers::pi - wrapped : wrapped;
}


// The frame the sphere is swept in. Latitude slicing degenerates for a circle whose axis sits on the
// polar axis, so the polar axis is chosen to be as far as possible from every one of them, out of a
// fixed set of candidates: the answer must not depend on which run computed it.
std::array<double3, 3> sweepFrame(const std::vector<double3>& axes)
{
  // Normalised where they are used rather than written out as unit vectors: the frame has to be
  // orthonormal to a rounding error, since a point placed with a frame that is not lands off the sphere
  // it was meant to be on, and then a neighbour's surface is the wrong side of it.
  static const std::array<double3, 6> candidates = {double3(1.0, 2.0, 3.0),  double3(-3.0, 1.0, 2.0),
                                                    double3(2.0, -3.0, 1.0), double3(3.0, 2.0, -1.0),
                                                    double3(1.0, -3.0, -2.0), double3(-2.0, 3.0, -1.0)};

  double3 polarAxis = candidates.front() * (1.0 / candidates.front().length());
  double bestSeparation = -1.0;
  for (const double3& candidate : candidates)
  {
    double3 direction = candidate * (1.0 / candidate.length());
    double separation = 1.0;
    for (const double3& axis : axes)
    {
      separation = std::min(separation, 1.0 - std::abs(double3::dot(direction, axis)));
    }
    if (separation > bestSeparation)
    {
      bestSeparation = separation;
      polarAxis = direction;
    }
    if (bestSeparation > 0.01) break;
  }

  // Any two directions perpendicular to the polar axis will do for the azimuth; taking the smallest
  // component of the axis to build the first of them keeps the cross product well conditioned.
  double3 helper(1.0, 0.0, 0.0);
  if (std::abs(polarAxis.y) < std::abs(polarAxis.x)) helper = double3(0.0, 1.0, 0.0);
  if (std::abs(polarAxis.z) < std::min(std::abs(polarAxis.x), std::abs(polarAxis.y)))
    helper = double3(0.0, 0.0, 1.0);
  double3 first = double3::cross(helper, polarAxis);
  first = first * (1.0 / first.length());
  double3 second = double3::cross(polarAxis, first);
  second = second * (1.0 / second.length());

  return {first, second, polarAxis};
}


// Everything needed to integrate over one sphere, reused across the spheres so that the buffers are
// allocated once for a structure rather than once per atom.
struct SphereWorkspace
{
  std::vector<BoundingCircle> circles;
  std::vector<double3> axes;
  std::vector<double> breakpoints;
  std::vector<std::size_t> cuttingCircles;  // the circles that cut the piece being integrated
  std::vector<std::pair<double, double>> coveredArcs;
};


// Where an arc's side is read off, when it is read off the connected surface it belongs to rather than
// asked about arc by arc. `origins` is the point each surface's moments are taken about, in that surface's
// own frame.
struct ComponentRoute
{
  const BoundaryComponents* components{nullptr};
  const std::vector<ComponentLabel>* labels{nullptr};
  const std::vector<double3>* origins{nullptr};
};


// The exposed area of one probe-inflated sphere, added to `sample` split by the classifier.
void measureSphere(const VoronoiAccessibility& accessibility, std::size_t atomIndex, std::size_t subdivisions,
                   bool classifyArcs, const ComponentRoute* route, SphereWorkspace& work,
                   ExactSurfaceAreaSample& sample)
{
  const double radius = accessibility.atomRadii[atomIndex];
  const double3 centre = accessibility.atomPositions[atomIndex];

  // Every sphere that can reach this one, its own periodic images included.
  work.circles.clear();
  work.axes.clear();
  for (const auto& [delta, neighbourRadius] :
       accessibility.neighbourAtoms(centre, radius + accessibility.maximumAtomRadius))
  {
    double distance = delta.length();
    if (distance < 1.0e-12)
    {
      // The sphere itself, which the query returns along with the rest, or another sitting on top of it.
      // Only a strictly larger one covers anything, and then it covers everything.
      if (neighbourRadius > radius) return;
      continue;
    }

    // The circle of intersection sits at this polar angle from the neighbour's direction. At least one
    // means the neighbour reaches nothing of this sphere; at most minus one means it swallows it whole.
    double cosineHalfAngle =
        (radius * radius + distance * distance - neighbourRadius * neighbourRadius) / (2.0 * radius * distance);
    if (cosineHalfAngle >= 1.0) continue;
    if (cosineHalfAngle <= -1.0) return;

    BoundingCircle circle;
    circle.axis = delta * (1.0 / distance);
    circle.cosineHalfAngle = cosineHalfAngle;
    circle.halfAngle = std::acos(cosineHalfAngle);
    work.circles.push_back(circle);
    work.axes.push_back(circle.axis);
  }

  // A circle whose disc lies inside another's carries no boundary and no coverage of its own, and
  // leaving it in would only add latitudes at which nothing happens.
  if (work.circles.size() > 1)
  {
    std::vector<bool> redundant(work.circles.size(), false);
    for (std::size_t i = 0; i < work.circles.size(); ++i)
    {
      for (std::size_t j = 0; j < work.circles.size(); ++j)
      {
        if (i == j || redundant[j]) continue;
        double separation = std::acos(std::clamp(double3::dot(work.circles[i].axis, work.circles[j].axis), -1.0, 1.0));
        if (separation + work.circles[i].halfAngle <= work.circles[j].halfAngle)
        {
          redundant[i] = true;
          break;
        }
      }
    }
    std::size_t kept = 0;
    for (std::size_t i = 0; i < work.circles.size(); ++i)
    {
      if (!redundant[i]) work.circles[kept++] = work.circles[i];
    }
    work.circles.resize(kept);
    work.axes.clear();
    for (const BoundingCircle& circle : work.circles) work.axes.push_back(circle.axis);
  }

  const std::array<double3, 3> frame = sweepFrame(work.axes);
  const double3 firstAxis = frame[0];
  const double3 secondAxis = frame[1];
  const double3 polarAxis = frame[2];

  work.breakpoints.clear();
  work.breakpoints.push_back(0.0);
  work.breakpoints.push_back(std::numbers::pi);
  for (BoundingCircle& circle : work.circles)
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

    work.breakpoints.push_back(circle.lowestLatitude);
    work.breakpoints.push_back(circle.highestLatitude);
  }

  // The corners of the patch are latitudes too. Where two of the circles cross, the arcs they cover run
  // into one another and the uncovered length turns a corner; between such latitudes it is an analytic
  // function of latitude, which is what the quadrature needs it to be. A crossing that a third circle
  // covers is not a corner of anything -- the boundary there belongs to that third circle and passes
  // through smoothly -- so those are dropped, and what is left is the handful of corners the patch
  // really has rather than one for every pair.
  for (std::size_t j = 0; j + 1 < work.circles.size(); ++j)
  {
    for (std::size_t k = j + 1; k < work.circles.size(); ++k)
    {
      const BoundingCircle& first = work.circles[j];
      const BoundingCircle& second = work.circles[k];
      double alignment = double3::dot(first.axis, second.axis);
      double denominator = 1.0 - alignment * alignment;
      if (denominator < 1.0e-14) continue;  // parallel axes: concentric circles never cross

      double alongFirst = (first.cosineHalfAngle - alignment * second.cosineHalfAngle) / denominator;
      double alongSecond = (second.cosineHalfAngle - alignment * first.cosineHalfAngle) / denominator;
      double outOfPlaneSquared =
          (1.0 - alongFirst * first.cosineHalfAngle - alongSecond * second.cosineHalfAngle) / denominator;
      if (outOfPlaneSquared <= 0.0) continue;  // the circles miss each other

      double3 inPlane = first.axis * alongFirst + second.axis * alongSecond;
      double3 outOfPlane = double3::cross(first.axis, second.axis) * std::sqrt(outOfPlaneSquared);
      for (std::size_t side = 0; side < 2; ++side)
      {
        double3 corner = (side == 0) ? inPlane + outOfPlane : inPlane - outOfPlane;

        // Covered has to mean covered with room to spare. A framework is symmetric, so three spheres
        // meeting in one point is the ordinary case rather than a coincidence, and a corner that a third
        // circle happens to pass exactly through is still a corner. Deciding that on the sign of a
        // rounding error loses it, and with it the latitude at which the integrand turns, so the test is
        // given room: a corner wrongly kept costs one more panel, a corner wrongly dropped costs digits.
        bool covered = false;
        for (std::size_t l = 0; l < work.circles.size() && !covered; ++l)
        {
          if (l == j || l == k) continue;
          covered = double3::dot(corner, work.circles[l].axis) > work.circles[l].cosineHalfAngle + 1.0e-9;
        }
        if (!covered)
        {
          work.breakpoints.push_back(std::acos(std::clamp(double3::dot(corner, polarAxis), -1.0, 1.0)));
        }
      }
    }
  }

  std::sort(work.breakpoints.begin(), work.breakpoints.end());

  const GaussRule& rule = unitIntervalGaussRule();
  const double twoPi = 2.0 * std::numbers::pi;
  const std::size_t parts = std::max<std::size_t>(1, subdivisions);

  for (std::size_t piece = 0; piece + 1 < work.breakpoints.size(); ++piece)
  {
    double pieceBegin = work.breakpoints[piece];
    double pieceEnd = work.breakpoints[piece + 1];
    if (pieceEnd - pieceBegin < 1.0e-14) continue;

    // Which circles cut the latitudes of this piece, and whether any of them covers them entirely. Both
    // are settled anywhere inside the piece, its ends being exactly the latitudes at which either can
    // change, so they are settled once here rather than at every node of the quadrature.
    double interior = 0.5 * (pieceBegin + pieceEnd);
    work.cuttingCircles.clear();
    bool buried = false;
    for (std::size_t i = 0; i < work.circles.size(); ++i)
    {
      const BoundingCircle& circle = work.circles[i];
      if (interior > circle.lowestLatitude && interior < circle.highestLatitude)
      {
        work.cuttingCircles.push_back(i);
      }
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

      // The uncovered length of a latitude leaves a piece end like the square root of the distance to
      // it, a circle appearing or vanishing there with a width of that shape. Anchoring the two halves
      // at the ends and substituting a square makes the integrand smooth in the new variable, which is
      // what leaves the quadrature error at round-off rather than at the square root of it.
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

          // Each crossing circle covers one arc of this latitude, centred on its own azimuth. An arc
          // that runs over zero is entered as its two pieces, so that the sweep below stays on [0, 2pi).
          work.coveredArcs.clear();
          for (std::size_t i : work.cuttingCircles)
          {
            const BoundingCircle& circle = work.circles[i];
            double cosineHalfWidth = (circle.cosineHalfAngle - cosineLatitude * circle.cosinePolar) /
                                     (sineLatitude * circle.sinePolar);
            if (cosineHalfWidth >= 1.0) continue;
            double halfWidth = (cosineHalfWidth <= -1.0) ? std::numbers::pi : std::acos(cosineHalfWidth);

            double arcBegin = circle.azimuth - halfWidth;
            double arcEnd = circle.azimuth + halfWidth;
            if (arcBegin < 0.0)
            {
              work.coveredArcs.emplace_back(arcBegin + twoPi, twoPi);
              work.coveredArcs.emplace_back(0.0, arcEnd);
            }
            else if (arcEnd > twoPi)
            {
              work.coveredArcs.emplace_back(arcBegin, twoPi);
              work.coveredArcs.emplace_back(0.0, arcEnd - twoPi);
            }
            else
            {
              work.coveredArcs.emplace_back(arcBegin, arcEnd);
            }
          }
          std::sort(work.coveredArcs.begin(), work.coveredArcs.end());

          // What is left of the latitude, arc by arc: each gap is one connected piece of exposed
          // surface, so its length is area and its midpoint is where that area is asked about.
          double cursor = 0.0;
          for (std::size_t arc = 0; arc <= work.coveredArcs.size(); ++arc)
          {
            double gapEnd = (arc < work.coveredArcs.size()) ? work.coveredArcs[arc].first : twoPi;

            double gap = gapEnd - cursor;
            if (gap > arcTolerance)
            {
              double area = radius * radius * sineLatitude * gap * weight;
              ++sample.numberOfArcs;
              sample.area += area;
              sample.radiusWeightedArea += radius * area;

              if (classifyArcs)
              {
                // The integral of the outward normal over this arc. Its azimuthal part is elementary, being
                // the difference of a sine and of a cosine between the ends of the same gap whose length
                // gave the area, so the moments cost three sums and no new geometry.
                double3 normalIntegral = firstAxis * (sineLatitude * (std::sin(gapEnd) - std::sin(cursor))) +
                                         secondAxis * (sineLatitude * (std::cos(cursor) - std::cos(gapEnd))) +
                                         polarAxis * (cosineLatitude * gap);
                double3 arcVectorArea = normalIntegral * (radius * radius * sineLatitude * weight);

                auto directionAt = [&](double azimuth)
                {
                  return firstAxis * (sineLatitude * std::cos(azimuth)) +
                         secondAxis * (sineLatitude * std::sin(azimuth)) + polarAxis * cosineLatitude;
                };

                // The azimuthal integrals of the products of two of the normal's components, which are what
                // the first moment of the enclosed region needs beyond the normal itself. They are the same
                // two endpoints once more, at twice the angle.
                double spread = 0.25 * (std::sin(2.0 * gapEnd) - std::sin(2.0 * cursor));
                double integralCosine = std::sin(gapEnd) - std::sin(cursor);
                double integralSine = std::cos(cursor) - std::cos(gapEnd);
                double integralCosineCosine = 0.5 * gap + spread;
                double integralSineSine = 0.5 * gap - spread;
                double integralSineCosine = 0.25 * (std::cos(2.0 * cursor) - std::cos(2.0 * gapEnd));

                // Everything one arc adds to the surface it belongs to, given where that surface's moments
                // are taken from. Both routes below add the same thing and differ only in which surface it
                // goes to and which copy of the cell the arc is carried into.
                auto addArcTo = [&](PoreBoundaryMoments& moments, const double3& shifted, const double3& origin)
                {
                  const double3 delta = shifted - origin;

                  // The tensor of the arc applied to `delta`, in the frame the sweep is done in.
                  double alongFirst = double3::dot(delta, firstAxis);
                  double alongSecond = double3::dot(delta, secondAxis);
                  double alongPolar = double3::dot(delta, polarAxis);
                  double3 tensorTimesDelta =
                      firstAxis * (sineLatitude * sineLatitude * alongFirst * integralCosineCosine +
                                   sineLatitude * sineLatitude * alongSecond * integralSineCosine +
                                   sineLatitude * cosineLatitude * alongPolar * integralCosine) +
                      secondAxis * (sineLatitude * sineLatitude * alongFirst * integralSineCosine +
                                    sineLatitude * sineLatitude * alongSecond * integralSineSine +
                                    sineLatitude * cosineLatitude * alongPolar * integralSine) +
                      polarAxis * (cosineLatitude * sineLatitude * alongFirst * integralCosine +
                                   cosineLatitude * sineLatitude * alongSecond * integralSine +
                                   cosineLatitude * cosineLatitude * alongPolar * gap);
                  tensorTimesDelta = tensorTimesDelta * (radius * radius * sineLatitude * weight);

                  moments.area += area;
                  moments.radiusWeightedArea += radius * area;
                  moments.originWeighted += double3::dot(delta, arcVectorArea);
                  moments.vectorArea += arcVectorArea;
                  moments.enclosedFirstMoment +=
                      arcVectorArea * (-0.5 * (double3::dot(delta, delta) + radius * radius)) -
                      tensorTimesDelta * radius;
                };

                if (route != nullptr)
                {
                  // Which patch the arc lies on. Either end of the gap is a point of that patch's own edge,
                  // where a bounding circle crosses this latitude, and the arcs of the circles carry the
                  // patch they bound: so the lookup is exact wherever there is such an end. A latitude no
                  // circle cuts has none, and only there is a path on the sphere looked for instead.
                  std::int32_t patch = -1;
                  if (cursor > 0.0)
                  {
                    patch = route->components->patchOfCirclePoint(atomIndex, directionAt(cursor));
                  }
                  if (patch < 0 && gapEnd < twoPi)
                  {
                    patch = route->components->patchOfCirclePoint(atomIndex, directionAt(gapEnd));
                  }
                  if (patch < 0)
                  {
                    patch = route->components->patchOfDirection(atomIndex, directionAt(0.5 * (cursor + gapEnd)));
                  }

                  std::int32_t label =
                      (patch < 0) ? -1 : route->components->componentOfPatch[atomIndex][static_cast<std::size_t>(patch)];
                  if (label < 0)
                  {
                    sample.undecided += area;
                    ++sample.unplacedArcs;
                    sample.unplacedArea += area;
                  }
                  else
                  {
                    // Which side the arc is on is not settled here. It is settled for the whole surface once
                    // the surface has been measured, since what decides it is the volume the surface encloses
                    // and that is not known until then.

                    // The arc is on an atom of the home cell, but the surface it belongs to runs through
                    // several copies of the cell, and it closes only once the pieces are brought into one
                    // frame. Which translation does that for this patch was settled when the surface was
                    // assembled, so nothing here has to be inferred from where the arc happens to be.
                    int3 offset = route->components->offsetOfPatch[atomIndex][static_cast<std::size_t>(patch)];
                    double3 shifted =
                        centre + accessibility.simulationBox.cell * double3(static_cast<double>(offset.x),
                                                                           static_cast<double>(offset.y),
                                                                           static_cast<double>(offset.z));

                    addArcTo(sample.components[static_cast<std::size_t>(label)], shifted,
                             (*route->origins)[static_cast<std::size_t>(label)]);
                  }
                }
                else
                {
                  double3 outward = directionAt(0.5 * (cursor + gapEnd));
                  PointClassification classification = accessibility.classify(centre + outward * radius);

                  if (classification.resample || classification.inside)
                    sample.undecided += area;
                  else if (classification.accessible)
                    sample.accessible += area;
                  else
                    sample.inaccessible += area;

                  if (classification.poreId >= 0)
                  {
                    // The arc faces one lift of its pore, and the pore's boundary closes only in the frame
                    // its own nodes are in, so the arc is carried there before it is added. Taking the
                    // origin at the pore's own first node rather than at some fixed point keeps |x - c|
                    // small, which matters: what survives of the choice of origin is the closure defect
                    // times the distance to it.
                    const VoronoiPore& pore =
                        accessibility.channels.pores[static_cast<std::size_t>(classification.poreId)];
                    double3 origin = accessibility.network.nodes[pore.nodeIndices.front()].position;
                    double3 shifted = centre + accessibility.simulationBox.cell *
                                                   double3(static_cast<double>(classification.latticeOffset.x),
                                                           static_cast<double>(classification.latticeOffset.y),
                                                           static_cast<double>(classification.latticeOffset.z));

                    addArcTo(sample.pores[static_cast<std::size_t>(classification.poreId)], shifted, origin);
                  }
                }
              }
            }
            if (arc < work.coveredArcs.size()) cursor = std::max(cursor, work.coveredArcs[arc].second);
          }
        }
      }
    }
  }
}


ExactSurfaceAreaSample exactAccessibleSurfaceArea(const VoronoiAccessibility& accessibility,
                                                  std::size_t subdivisions, bool classifyArcs)
{
  ExactSurfaceAreaSample sample;
  if (classifyArcs) sample.pores.assign(accessibility.channels.pores.size(), PoreBoundaryMoments{});

  SphereWorkspace work;
  for (std::size_t i = 0; i < accessibility.atomPositions.size(); ++i)
  {
    measureSphere(accessibility, i, subdivisions, classifyArcs, nullptr, work, sample);
  }
  return sample;
}


std::vector<double3> surfaceMomentOrigins(const VoronoiAccessibility& accessibility,
                                          const BoundaryComponents& components)
{
  // Each surface's moments are taken about a point of the surface itself, carried into the surface's own
  // frame. The choice cannot change a closed surface's volume; what it changes is how much of the closure
  // defect shows up in it, and a near point is what keeps that small.
  std::vector<double3> origins(components.numberOfComponents, double3(0.0, 0.0, 0.0));
  for (std::size_t label = 0; label < components.numberOfComponents; ++label)
  {
    const auto [atomIndex, patch] = components.componentRepresentative[label];
    if (patch >= components.offsetOfPatch[atomIndex].size()) continue;
    int3 offset = components.offsetOfPatch[atomIndex][patch];
    origins[label] = accessibility.atomPositions[atomIndex] +
                     accessibility.simulationBox.cell * double3(static_cast<double>(offset.x),
                                                                static_cast<double>(offset.y),
                                                                static_cast<double>(offset.z));
  }
  return origins;
}


ExactSurfaceAreaSample exactAccessibleSurfaceAreaByComponent(const VoronoiAccessibility& accessibility,
                                                             const BoundaryComponents& components,
                                                             const std::vector<ComponentLabel>& labels,
                                                             std::size_t subdivisions)
{
  ExactSurfaceAreaSample sample;
  sample.components.assign(components.numberOfComponents, PoreBoundaryMoments{});

  std::vector<double3> origins = surfaceMomentOrigins(accessibility, components);

  ComponentRoute route;
  route.components = &components;
  route.labels = &labels;
  route.origins = &origins;

  SphereWorkspace work;
  for (std::size_t i = 0; i < accessibility.atomPositions.size(); ++i)
  {
    measureSphere(accessibility, i, subdivisions, true, &route, work, sample);
  }

  // Which side each surface is on, now that each has been measured whole. Two of the three cases are settled
  // by the surface itself. A surface that closes on a translate of itself runs away through the crystal, and
  // so does the void along it, so its area is reachable. A surface that closes on itself has the void on one
  // side and the solid on the other, and the sign of the volume it encloses says which: positive is void, and
  // void enclosed by a closed surface can reach nothing outside it, so that area is not reachable. Only a
  // surface enclosing solid --- a cluster of atoms with the void outside it --- leaves a question for the
  // network, which is whether the void around the cluster is itself sealed.
  for (std::size_t label = 0; label < components.numberOfComponents; ++label)
  {
    const PoreBoundaryMoments& moments = sample.components[label];
    if (moments.area <= 0.0) continue;

    ++sample.numberOfSurfaces;

    bool reachable;
    if (components.componentPercolates[label] != 0)
    {
      ++sample.runawaySurfaces;
      reachable = true;
    }
    else if (-(moments.radiusWeightedArea + moments.originWeighted) > 0.0)
    {
      ++sample.sealedSurfaces;
      reachable = false;
    }
    else
    {
      ++sample.clusterSurfaces;
      const ComponentLabel& answer = labels[label];
      if (!answer.decided)
      {
        sample.undecided += moments.area;
        continue;
      }
      reachable = answer.accessible;
    }

    if (reachable)
      sample.accessible += moments.area;
    else
      sample.inaccessible += moments.area;
  }
  return sample;
}


ExactSurfaceAreaSample exactAccessibleSurfaceAreaByComponent(const VoronoiAccessibility& accessibility,
                                                             std::size_t subdivisions)
{
  BoundaryComponents components = boundaryComponents(accessibility);
  std::vector<ComponentLabel> labels = labelBoundaryComponents(accessibility, components);
  return exactAccessibleSurfaceAreaByComponent(accessibility, components, labels, subdivisions);
}
