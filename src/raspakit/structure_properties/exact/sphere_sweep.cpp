module;

module exact_sphere_sweep;

import std;

import double3;

double3 perpendicularTo(const double3& axis)
{
  double3 helper(1.0, 0.0, 0.0);
  if (std::abs(axis.y) < std::abs(axis.x)) helper = double3(0.0, 1.0, 0.0);
  if (std::abs(axis.z) < std::min(std::abs(axis.x), std::abs(axis.y))) helper = double3(0.0, 0.0, 1.0);
  double3 perpendicular = double3::cross(helper, axis);
  return perpendicular * (1.0 / perpendicular.length());
}


double foldedPolarAngle(double angle)
{
  double wrapped = std::fmod(std::abs(angle), 2.0 * std::numbers::pi);
  return (wrapped > std::numbers::pi) ? 2.0 * std::numbers::pi - wrapped : wrapped;
}


std::array<double3, 3> sweepFrame(const std::vector<double3>& axes)
{
  // Normalised where they are used rather than written out as unit vectors: the frame has to be orthonormal
  // to a rounding error, since a point placed with a frame that is not lands off the sphere it was meant to
  // be on, and then a neighbour's surface is the wrong side of it.
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

  double3 first = perpendicularTo(polarAxis);
  double3 second = double3::cross(polarAxis, first);
  second = second * (1.0 / second.length());

  return {first, second, polarAxis};
}


std::optional<SweepCircle> makeSweepCircle(const double3& axis, double cosineHalfAngle)
{
  if (cosineHalfAngle >= 1.0) return std::nullopt;  // the cap covers nothing of the sphere

  SweepCircle circle;
  circle.axis = axis;
  circle.cosineHalfAngle = cosineHalfAngle;
  if (cosineHalfAngle <= -1.0)
  {
    circle.halfAngle = std::numbers::pi;
    circle.sineHalfAngle = 0.0;
    return circle;
  }
  circle.halfAngle = std::acos(cosineHalfAngle);
  circle.sineHalfAngle = std::sin(circle.halfAngle);
  return circle;
}


void prepareSweep(std::vector<SweepCircle>& circles, const std::array<double3, 3>& frame,
                  const std::vector<double3>* knownCrossings, SweepWorkspace& work)
{
  const double3& firstAxis = frame[0];
  const double3& secondAxis = frame[1];
  const double3& polarAxis = frame[2];

  std::vector<double>& breakpoints = work.breakpoints;
  breakpoints.clear();
  breakpoints.push_back(0.0);
  breakpoints.push_back(std::numbers::pi);

  for (SweepCircle& circle : circles)
  {
    circle.polarAngle = std::acos(std::clamp(double3::dot(circle.axis, polarAxis), -1.0, 1.0));
    circle.cosinePolar = std::cos(circle.polarAngle);
    circle.sinePolar = std::sin(circle.polarAngle);
    circle.azimuth = std::atan2(double3::dot(circle.axis, secondAxis), double3::dot(circle.axis, firstAxis));
    if (circle.azimuth < 0.0) circle.azimuth += 2.0 * std::numbers::pi;
    circle.cosineAzimuth = std::cos(circle.azimuth);
    circle.sineAzimuth = std::sin(circle.azimuth);

    circle.lowestLatitude = foldedPolarAngle(circle.polarAngle - circle.halfAngle);
    circle.highestLatitude = foldedPolarAngle(circle.polarAngle + circle.halfAngle);
    circle.reachesOverPole = circle.polarAngle < circle.halfAngle;
    circle.reachesOverAntipole = std::numbers::pi - circle.polarAngle < circle.halfAngle;

    breakpoints.push_back(circle.lowestLatitude);
    breakpoints.push_back(circle.highestLatitude);
  }

  // In order of their own azimuth, which is what the sweep wants of them: a cap covers an arc centred on its
  // azimuth, so taking them in this order leaves the covered arcs of every circle of latitude sorted but for
  // the one or two that run over zero. Ordered once here rather than at every piece of the integral.
  std::sort(circles.begin(), circles.end(),
            [](const SweepCircle& a, const SweepCircle& b) { return a.azimuth < b.azimuth; });

  // The corners of the region are latitudes too. Where two of the caps cross, the arcs they cover run into
  // one another and the uncovered length turns a corner; between such latitudes it is an analytic function of
  // latitude, which is what the quadrature needs it to be. A crossing that a third cap covers is not a corner
  // of anything --- the boundary there belongs to that third cap and passes through smoothly --- so what is
  // left is the handful of corners the region really has rather than one for every pair.
  if (knownCrossings == nullptr)
  {
    uncoveredCrossings(circles, work.crossings);
    for (const CapCrossing& crossing : work.crossings)
    {
      breakpoints.push_back(std::acos(std::clamp(double3::dot(crossing.direction, polarAxis), -1.0, 1.0)));
    }
  }
  else
  {
    for (const double3& crossing : *knownCrossings)
    {
      breakpoints.push_back(std::acos(std::clamp(double3::dot(crossing, polarAxis), -1.0, 1.0)));
    }
  }

  std::sort(breakpoints.begin(), breakpoints.end());
}


void panelBoundaries(double begin, double end, std::size_t subdivisions, bool cut, std::vector<double>& panels)
{
  const std::size_t parts = std::max<std::size_t>(1, subdivisions);

  panels.clear();
  for (std::size_t part = 0; part <= parts; ++part)
  {
    panels.push_back(begin + (end - begin) * static_cast<double>(part) / static_cast<double>(parts));
  }
  panels.front() = begin;
  panels.back() = end;

  if (!cut) return;

  // However little room is left, halving reaches it in a few dozen steps; the bound is there so that a
  // degenerate piece ending on a pole exactly cannot spin here.
  constexpr std::size_t halvingLimit = 60uz;

  const double roomBelow = begin;
  for (std::size_t halving = 0; halving < halvingLimit; ++halving)
  {
    if (panels[1] - panels[0] <= poleClearance * roomBelow) break;
    panels.insert(panels.begin() + 1, 0.5 * (panels[0] + panels[1]));
  }

  const double roomAbove = std::numbers::pi - end;
  for (std::size_t halving = 0; halving < halvingLimit; ++halving)
  {
    const std::size_t last = panels.size() - 1;
    if (panels[last] - panels[last - 1] <= poleClearance * roomAbove) break;
    panels.insert(panels.end() - 1, 0.5 * (panels[last - 1] + panels[last]));
  }
}


const GaussRule& unitIntervalGaussRule()
{
  static const GaussRule rule = []()
  {
    GaussRule constructed;
    for (std::size_t i = 0; i < exactQuadratureOrder; ++i)
    {
      double abscissa = std::cos(std::numbers::pi * (static_cast<double>(i) + 0.75) /
                                 (static_cast<double>(exactQuadratureOrder) + 0.5));
      double derivative = 0.0;
      for (std::size_t iteration = 0; iteration < 100; ++iteration)
      {
        double previous = 1.0;
        double current = abscissa;
        for (std::size_t k = 2; k <= exactQuadratureOrder; ++k)
        {
          double next =
              ((2.0 * static_cast<double>(k) - 1.0) * abscissa * current - (static_cast<double>(k) - 1.0) * previous) /
              static_cast<double>(k);
          previous = current;
          current = next;
        }
        derivative =
            static_cast<double>(exactQuadratureOrder) * (abscissa * current - previous) / (abscissa * abscissa - 1.0);
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
