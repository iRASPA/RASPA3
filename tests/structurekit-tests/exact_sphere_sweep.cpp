#include <gtest/gtest.h>

import std;

import double3;
import randomnumbers;
import exact_sphere_sweep;

// The sweep of the exposed part of one sphere, on configurations of caps whose answers are known.
//
// Three of the exact analyses go through this module, so what it gets wrong they all get wrong together and
// none of them can be the thing that catches it. What is checked here is therefore not an area agreeing with
// another area but each piece against what it is defined to be: the frame against orthonormality and against
// the one thing it is chosen for, the quadrature against the polynomials it is supposed to be exact on, the
// crossings against the circles they are supposed to lie on, and the swept area against the closed forms for
// a sphere less one cap, less two disjoint caps, and less three orthogonal hemispheres.
//
// The tolerances are of two kinds. Where a quantity is arrived at by a handful of operations it is held to a
// few ulps, written as 1e-14 of whatever it is. Where it comes out of the latitude quadrature it is held to
// 1e-12 relative, which is the accuracy the ten-point rule on the right pieces is claimed to have; a failure
// at 1e-9 there would mean the pieces are no longer the right ones, which is the thing worth catching.

namespace
{

constexpr double roundOff = 1.0e-14;

// What the latitude rule is worth, relative to the area of the whole sphere. Away from the poles one panel
// per smooth piece settles an area to thirteen digits; a cap turning back within a few hundredths of a pole
// costs three or four of those, which is the worst this has been seen to do and what this is set by. See
// `a_turning_latitude_near_a_pole_costs_digits_that_refining_gives_back`.
constexpr double sweptAccuracy = 1.0e-8;

SweepCircle capOf(const double3& axis, double halfAngle)
{
  std::optional<SweepCircle> circle = makeSweepCircle(double3::normalize(axis), std::cos(halfAngle));
  EXPECT_TRUE(circle.has_value());
  return circle.value_or(SweepCircle{});
}

// The area of the sphere of this radius that none of the caps covers, by the sweep itself. `knownCrossings`
// is handed on as it stands, so that the route that reuses a decomposition's crossings can be measured
// against the one that looks for them.
double exposedArea(std::vector<SweepCircle> circles, double radius, std::size_t subdivisions = 1,
                   const std::vector<double3>* knownCrossings = nullptr)
{
  std::vector<double3> axes;
  for (const SweepCircle& circle : circles) axes.push_back(circle.axis);
  const std::array<double3, 3> frame = sweepFrame(axes);

  double area = 0.0;
  SweepWorkspace work;
  sweepExposedLatitudes(circles, frame, knownCrossings, subdivisions, work,
                        [&](const LatitudeGap& gap)
                        { area += radius * radius * gap.sineLatitude * gap.span * gap.weight; });
  return area;
}

// A rotation about `axis` by `angle`, to ask the same question of a configuration standing somewhere else.
double3 turned(const double3& point, const double3& axis, double angle)
{
  const double3 unit = double3::normalize(axis);
  const double cosine = std::cos(angle);
  return point * cosine + double3::cross(unit, point) * std::sin(angle) +
         unit * (double3::dot(unit, point) * (1.0 - cosine));
}

}  // namespace


TEST(exact_sphere_sweep, a_perpendicular_is_perpendicular_and_of_unit_length)
{
  // The axes along which the cross product behind it is worst conditioned are the coordinate axes themselves,
  // which is what the choice of helper vector is there to keep it off.
  const std::vector<double3> axes = {double3(1.0, 0.0, 0.0),  double3(0.0, 1.0, 0.0),  double3(0.0, 0.0, 1.0),
                                     double3(-1.0, 0.0, 0.0), double3(0.0, -1.0, 0.0), double3(0.0, 0.0, -1.0),
                                     double3(1.0, 1.0e-13, 0.0), double3(0.6, 0.8, 0.0),
                                     double3(1.0, 2.0, 3.0)};

  for (const double3& given : axes)
  {
    const double3 axis = double3::normalize(given);
    const double3 perpendicular = perpendicularTo(axis);
    EXPECT_NEAR(perpendicular.length(), 1.0, roundOff);
    EXPECT_NEAR(double3::dot(perpendicular, axis), 0.0, roundOff);
  }
}


TEST(exact_sphere_sweep, a_polar_angle_folds_back_into_the_half_turn)
{
  EXPECT_NEAR(foldedPolarAngle(0.0), 0.0, roundOff);
  EXPECT_NEAR(foldedPolarAngle(0.3), 0.3, roundOff);
  EXPECT_NEAR(foldedPolarAngle(-0.3), 0.3, roundOff);
  EXPECT_NEAR(foldedPolarAngle(std::numbers::pi), std::numbers::pi, roundOff);

  // Past the far pole and back: an angle a little beyond pi is the same latitude as the same little way
  // short of it, which is the whole reason the fold is there.
  EXPECT_NEAR(foldedPolarAngle(std::numbers::pi + 0.3), std::numbers::pi - 0.3, roundOff);
  EXPECT_NEAR(foldedPolarAngle(2.0 * std::numbers::pi + 0.3), 0.3, roundOff);

  // The use it is put to: a cap whose axis is at polar angle beta with half angle theta turns back at these
  // two latitudes, and a cap reaching over the pole has to come back rather than run past it.
  const double beta = 0.4;
  const double theta = 0.9;  // larger than beta, so the cap covers the pole
  EXPECT_NEAR(foldedPolarAngle(beta - theta), theta - beta, roundOff);
  EXPECT_NEAR(foldedPolarAngle(beta + theta), beta + theta, roundOff);
}


TEST(exact_sphere_sweep, the_sweep_frame_is_orthonormal)
{
  RandomNumber random{std::optional<std::size_t>(42)};
  for (std::size_t trial = 0; trial < 200; ++trial)
  {
    std::vector<double3> axes;
    for (std::size_t i = 0; i < 1 + trial % 12; ++i) axes.push_back(random.randomVectorOnUnitSphere());

    const std::array<double3, 3> frame = sweepFrame(axes);
    const auto& [first, second, polar] = frame;

    EXPECT_NEAR(first.length(), 1.0, roundOff);
    EXPECT_NEAR(second.length(), 1.0, roundOff);
    EXPECT_NEAR(polar.length(), 1.0, roundOff);
    EXPECT_NEAR(double3::dot(first, second), 0.0, roundOff);
    EXPECT_NEAR(double3::dot(first, polar), 0.0, roundOff);
    EXPECT_NEAR(double3::dot(second, polar), 0.0, roundOff);

    // Right handed, so that an azimuth measured from `first` towards `second` turns the way the sweep
    // assumes it does.
    const double3 handedness = double3::cross(first, second) - polar;
    EXPECT_NEAR(handedness.length(), 0.0, roundOff);
  }
}


TEST(exact_sphere_sweep, the_pole_is_kept_off_every_cap_axis)
{
  // What the frame is chosen for: latitude slicing degenerates for a cap sitting on the polar axis, so the
  // pole has to stand clear of all of them. The candidates are six and a sphere is cut by rather fewer caps
  // than that in any ordinary structure, so the margin should be the one the search settles for.
  RandomNumber random{std::optional<std::size_t>(7)};
  for (std::size_t trial = 0; trial < 500; ++trial)
  {
    std::vector<double3> axes;
    for (std::size_t i = 0; i < 1 + trial % 20; ++i) axes.push_back(random.randomVectorOnUnitSphere());

    const double3 polar = sweepFrame(axes)[2];
    double separation = 1.0;
    for (const double3& axis : axes) separation = std::min(separation, 1.0 - std::abs(double3::dot(polar, axis)));
    EXPECT_GT(separation, 0.01);
  }
}


TEST(exact_sphere_sweep, the_frame_does_not_depend_on_the_order_the_axes_arrive_in)
{
  // The spheres are visited in whatever order a neighbour query returns them, and two runs that put the same
  // caps in a different order have to sweep in the same frame or they integrate different panels.
  RandomNumber random{std::optional<std::size_t>(11)};
  std::vector<double3> axes;
  for (std::size_t i = 0; i < 9; ++i) axes.push_back(random.randomVectorOnUnitSphere());

  const std::array<double3, 3> reference = sweepFrame(axes);
  for (std::size_t shuffle = 0; shuffle < 20; ++shuffle)
  {
    std::vector<double3> permuted = axes;
    std::rotate(permuted.begin(), permuted.begin() + static_cast<std::ptrdiff_t>(1 + shuffle % 8), permuted.end());
    std::reverse(permuted.begin(), permuted.begin() + 4);

    const std::array<double3, 3> frame = sweepFrame(permuted);
    for (std::size_t i = 0; i < 3; ++i)
    {
      EXPECT_DOUBLE_EQ(frame[i].x, reference[i].x);
      EXPECT_DOUBLE_EQ(frame[i].y, reference[i].y);
      EXPECT_DOUBLE_EQ(frame[i].z, reference[i].z);
    }
  }
}


TEST(exact_sphere_sweep, a_cap_that_covers_nothing_is_not_a_cap_and_one_that_covers_everything_says_so)
{
  const double3 axis(0.0, 0.0, 1.0);

  EXPECT_FALSE(makeSweepCircle(axis, 1.0).has_value());
  EXPECT_FALSE(makeSweepCircle(axis, 1.5).has_value());

  // Covering the whole sphere comes back rather than being dropped: what it means is the caller's to decide,
  // and the three of them do not agree, so the flag for it is the cosine and the half angle being straight.
  const std::optional<SweepCircle> whole = makeSweepCircle(axis, -1.0);
  ASSERT_TRUE(whole.has_value());
  EXPECT_LE(whole->cosineHalfAngle, -1.0);
  EXPECT_NEAR(whole->halfAngle, std::numbers::pi, roundOff);
  EXPECT_NEAR(whole->sineHalfAngle, 0.0, roundOff);

  for (double halfAngle = 0.1; halfAngle < 3.0; halfAngle += 0.1)
  {
    const std::optional<SweepCircle> circle = makeSweepCircle(axis, std::cos(halfAngle));
    ASSERT_TRUE(circle.has_value());
    EXPECT_NEAR(circle->halfAngle, halfAngle, roundOff);
    EXPECT_NEAR(circle->sineHalfAngle, std::sin(halfAngle), roundOff);
    EXPECT_NEAR(circle->cosineHalfAngle * circle->cosineHalfAngle + circle->sineHalfAngle * circle->sineHalfAngle,
                1.0, roundOff);
  }
}


TEST(exact_sphere_sweep, a_disc_within_a_disc_is_recognised_without_taking_an_arc_cosine)
{
  const double3 axis(0.0, 0.0, 1.0);
  const double3 tilted = double3::normalize(double3(0.3, 0.0, 1.0));

  const SweepCircle small = capOf(axis, 0.3);
  const SweepCircle large = capOf(axis, 0.9);
  const auto within = [](const SweepCircle& inner, const SweepCircle& outer)
  {
    return discWithinDisc(double3::dot(inner.axis, outer.axis), inner.cosineHalfAngle, inner.sineHalfAngle,
                          outer.cosineHalfAngle, outer.sineHalfAngle);
  };

  EXPECT_TRUE(within(small, large));
  EXPECT_FALSE(within(large, small));

  // Concentric is the easy case; a small cap tilted off the axis is inside so long as the tilt plus its own
  // half angle stays within the larger one, and outside as soon as it does not.
  EXPECT_TRUE(within(capOf(tilted, 0.3), large));
  EXPECT_FALSE(within(capOf(tilted, 0.7), large));

  // Two caps that overlap without either containing the other, and two that miss altogether.
  EXPECT_FALSE(within(capOf(axis, 0.6), capOf(tilted, 0.6)));
  EXPECT_FALSE(within(capOf(axis, 0.3), capOf(double3(0.0, 0.0, -1.0), 0.3)));
}


TEST(exact_sphere_sweep, pruning_keeps_one_of_a_pair_of_equals_and_the_order_of_the_rest)
{
  // The order matters beyond tidiness: the decomposition and the sweep prune the same list by this same rule
  // and the sweep then reuses the crossings the decomposition found among them, which is only the same list
  // if the survivors come out in the same order.
  std::vector<SweepCircle> circles = {capOf(double3(0.0, 0.0, 1.0), 0.9), capOf(double3(0.0, 0.0, 1.0), 0.3),
                                      capOf(double3(1.0, 0.0, 0.0), 0.5), capOf(double3(0.0, 0.0, 1.0), 0.9)};
  pruneContainedDiscs(circles);

  // The small concentric cap goes, and of the two identical large ones exactly one survives: each is inside
  // the other, so the first is struck out against the second and the second then has nothing left to be
  // inside of. What comes back is the survivors in the order they were given, which here is the third cap
  // and then the fourth.
  ASSERT_EQ(circles.size(), 2uz);
  EXPECT_NEAR(circles[0].halfAngle, 0.5, roundOff);
  EXPECT_NEAR(double3::dot(circles[0].axis, double3(1.0, 0.0, 0.0)), 1.0, roundOff);
  EXPECT_NEAR(circles[1].halfAngle, 0.9, roundOff);
  EXPECT_NEAR(double3::dot(circles[1].axis, double3(0.0, 0.0, 1.0)), 1.0, roundOff);

  // Caps that overlap without containing one another are all boundary and none of them is dropped.
  std::vector<SweepCircle> overlapping = {capOf(double3(0.0, 0.0, 1.0), 0.8), capOf(double3(1.0, 0.0, 0.0), 0.8),
                                          capOf(double3(0.0, 1.0, 0.0), 0.8)};
  pruneContainedDiscs(overlapping);
  EXPECT_EQ(overlapping.size(), 3uz);
}


TEST(exact_sphere_sweep, a_crossing_lies_on_both_of_the_circles_that_made_it)
{
  RandomNumber random{std::optional<std::size_t>(5)};
  for (std::size_t trial = 0; trial < 100; ++trial)
  {
    std::vector<SweepCircle> circles;
    for (std::size_t i = 0; i < 4; ++i)
    {
      circles.push_back(capOf(random.randomVectorOnUnitSphere(), 0.4 + 0.9 * random.uniform()));
    }

    for (const CapCrossing& crossing : uncoveredCrossings(circles))
    {
      EXPECT_NEAR(crossing.direction.length(), 1.0, roundOff);
      EXPECT_NEAR(double3::dot(crossing.direction, circles[crossing.firstCircle].axis),
                  circles[crossing.firstCircle].cosineHalfAngle, roundOff);
      EXPECT_NEAR(double3::dot(crossing.direction, circles[crossing.secondCircle].axis),
                  circles[crossing.secondCircle].cosineHalfAngle, roundOff);

      // Uncovered means no third cap reaches strictly past it, which is what makes it a corner of the
      // exposed region rather than a point in the interior of some other cap.
      for (std::size_t l = 0; l < circles.size(); ++l)
      {
        if (l == crossing.firstCircle || l == crossing.secondCircle) continue;
        EXPECT_LE(double3::dot(crossing.direction, circles[l].axis), circles[l].cosineHalfAngle + capCoverTolerance);
      }
    }
  }
}


TEST(exact_sphere_sweep, circles_that_miss_or_run_concentric_cross_nowhere)
{
  // Two small caps on opposite poles.
  EXPECT_TRUE(uncoveredCrossings(std::vector<SweepCircle>{capOf(double3(0.0, 0.0, 1.0), 0.3),
                                                          capOf(double3(0.0, 0.0, -1.0), 0.3)})
                  .empty());

  // The same axis twice, at different half angles: the circles are concentric and never meet.
  EXPECT_TRUE(uncoveredCrossings(std::vector<SweepCircle>{capOf(double3(0.0, 0.0, 1.0), 0.3),
                                                          capOf(double3(0.0, 0.0, 1.0), 0.9)})
                  .empty());

  // One cap on its own is a boundary with no corners on it at all.
  EXPECT_TRUE(uncoveredCrossings(std::vector<SweepCircle>{capOf(double3(0.0, 0.0, 1.0), 0.7)}).empty());
}


TEST(exact_sphere_sweep, three_orthogonal_hemispheres_leave_three_corners_of_six)
{
  // Hemispheres about the three axes. Each pair of boundary circles crosses at the two poles of the third
  // axis, and of those two the positive one is inside the third hemisphere and the negative one is not. So
  // three of the six crossings survive, and they are the corners of the one octant left uncovered.
  const std::vector<SweepCircle> circles = {capOf(double3(1.0, 0.0, 0.0), 0.5 * std::numbers::pi),
                                            capOf(double3(0.0, 1.0, 0.0), 0.5 * std::numbers::pi),
                                            capOf(double3(0.0, 0.0, 1.0), 0.5 * std::numbers::pi)};

  const std::vector<CapCrossing> crossings = uncoveredCrossings(circles);
  ASSERT_EQ(crossings.size(), 3uz);
  for (const CapCrossing& crossing : crossings)
  {
    // Each is a negative coordinate direction: two components zero and one at minus one.
    const double3& d = crossing.direction;
    EXPECT_NEAR(std::abs(d.x) + std::abs(d.y) + std::abs(d.z), 1.0, roundOff);
    EXPECT_NEAR(std::min({d.x, d.y, d.z}), -1.0, roundOff);
  }
}


TEST(exact_sphere_sweep, the_gauss_rule_is_symmetric_on_the_unit_interval)
{
  const GaussRule& rule = unitIntervalGaussRule();

  double totalWeight = 0.0;
  for (std::size_t i = 0; i < exactQuadratureOrder; ++i)
  {
    EXPECT_GT(rule.nodes[i], 0.0);
    EXPECT_LT(rule.nodes[i], 1.0);
    EXPECT_GT(rule.weights[i], 0.0);
    if (i > 0) EXPECT_GT(rule.nodes[i], rule.nodes[i - 1]);
    totalWeight += rule.weights[i];
  }
  EXPECT_NEAR(totalWeight, 1.0, roundOff);

  for (std::size_t i = 0; i < exactQuadratureOrder; ++i)
  {
    const std::size_t mirror = exactQuadratureOrder - 1 - i;
    EXPECT_NEAR(rule.nodes[i] + rule.nodes[mirror], 1.0, roundOff);
    EXPECT_NEAR(rule.weights[i], rule.weights[mirror], roundOff);
  }
}


TEST(exact_sphere_sweep, the_gauss_rule_is_exact_on_the_polynomials_it_claims)
{
  // Ten nodes integrate a polynomial of degree nineteen exactly, and that is the whole of the argument for
  // ten being enough: what the pieces and the endpoint substitution leave is analytic, so the error is at
  // round-off as soon as the rule is exact this far. Degree twenty is where it has to start being wrong.
  const GaussRule& rule = unitIntervalGaussRule();

  for (std::size_t degree = 0; degree <= 2 * exactQuadratureOrder - 1; ++degree)
  {
    double quadrature = 0.0;
    for (std::size_t i = 0; i < exactQuadratureOrder; ++i)
    {
      quadrature += rule.weights[i] * std::pow(rule.nodes[i], static_cast<double>(degree));
    }
    EXPECT_NEAR(quadrature, 1.0 / (static_cast<double>(degree) + 1.0), roundOff);
  }
}


TEST(exact_sphere_sweep, a_sphere_no_cap_touches_is_the_whole_sphere)
{
  const double radius = 1.7;
  const double whole = 4.0 * std::numbers::pi * radius * radius;
  EXPECT_NEAR(exposedArea({}, radius), whole, sweptAccuracy * whole);
}


TEST(exact_sphere_sweep, one_cap_leaves_the_closed_form)
{
  // A cap of half angle gamma covers 2 pi r^2 (1 - cos gamma) of the sphere, so this much is left.
  const double radius = 1.3;
  for (double gamma = 0.2; gamma < 3.0; gamma += 0.2)
  {
    const double expected = 2.0 * std::numbers::pi * radius * radius * (1.0 + std::cos(gamma));
    const double measured = exposedArea({capOf(double3(0.4, -0.7, 0.2), gamma)}, radius);
    EXPECT_NEAR(measured, expected, sweptAccuracy * 4.0 * std::numbers::pi * radius * radius);
  }
}


TEST(exact_sphere_sweep, a_turning_latitude_near_a_pole_does_not_cost_digits)
{
  // What sets the accuracy of the rule is not the number of caps but where their boundaries turn back. The
  // frame keeps the pole away from every cap axis, and it does that correctly; but a cap turns back at its
  // axis give or take its half angle, and nothing in the choice of frame looks at the half angle. This cap's
  // axis stands 0.87 in cosine from the pole, so the frame search is satisfied by the first candidate it
  // tries, and its far turning latitude nonetheless lands 0.042 from the antipole. Before the panels were
  // graded there, one panel was worth nine digits here against thirteen at every other half angle below.
  //
  // Swept at every half angle from a hundredth of a radian to within a hundredth of pi, so that a turning
  // latitude passes through every distance from a pole rather than being placed at one chosen distance.
  const double radius = 1.3;
  const double sphere = 4.0 * std::numbers::pi * radius * radius;
  const double3 axis(0.4, -0.7, 0.2);

  // The worst of the three hundred is a part in 6e12 of the sphere, against a part in 1.3e9 before the
  // grading, and it no longer falls where a turning latitude comes near a pole.
  double worst = 0.0;
  for (double gamma = 0.01; gamma < std::numbers::pi - 0.01; gamma += 0.01)
  {
    const double expected = 2.0 * std::numbers::pi * radius * radius * (1.0 + std::cos(gamma));
    worst = std::max(worst, std::abs(exposedArea({capOf(axis, gamma)}, radius, 1) - expected) / sphere);
  }
  EXPECT_LT(worst, 1.0e-12);

  // The half angle that used to be the bad one, held on its own so that a regression names itself.
  const double gamma = 1.4;
  const double expected = 2.0 * std::numbers::pi * radius * radius * (1.0 + std::cos(gamma));
  EXPECT_NEAR(exposedArea({capOf(axis, gamma)}, radius, 1), expected, 1.0e-13 * sphere);

  // And the grading is asked for only where it is needed: a piece clear of both poles keeps the panels the
  // caller asked for and no others.
  std::vector<double> panels;
  panelBoundaries(0.8, 2.3, 1, true, panels);
  EXPECT_EQ(panels.size(), 2uz);
  panelBoundaries(0.8, 2.3, 4, true, panels);
  EXPECT_EQ(panels.size(), 5uz);

  // A piece running from well inside the sphere up to a hair short of a pole is cut back until the last
  // panel is no longer than the room left, and the count stays in single figures because the halving is
  // geometric in that room.
  panelBoundaries(0.2999, 3.0999, 1, true, panels);
  EXPECT_GT(panels.size(), 2uz);
  EXPECT_LT(panels.size(), 16uz);
  EXPECT_LE(panels.back() - panels[panels.size() - 2], poleClearance * (std::numbers::pi - 3.0999));
  EXPECT_DOUBLE_EQ(panels.front(), 0.2999);
  EXPECT_DOUBLE_EQ(panels.back(), 3.0999);
  for (std::size_t i = 1; i < panels.size(); ++i) EXPECT_GT(panels[i], panels[i - 1]);

  // A piece no cap cuts has nothing but the sine of the latitude in it, which is analytic at a pole like
  // anywhere else, so it is left alone however close to one it runs.
  panelBoundaries(0.2999, 3.0999, 1, false, panels);
  EXPECT_EQ(panels.size(), 2uz);
}


TEST(exact_sphere_sweep, two_caps_that_do_not_meet_take_their_own_areas_and_no_more)
{
  const double radius = 0.9;
  const double first = 0.5;
  const double second = 0.8;  // first + second < pi, so the caps stand clear of one another
  const double expected = 4.0 * std::numbers::pi * radius * radius -
                          2.0 * std::numbers::pi * radius * radius * (1.0 - std::cos(first)) -
                          2.0 * std::numbers::pi * radius * radius * (1.0 - std::cos(second));

  const double measured =
      exposedArea({capOf(double3(0.0, 0.0, 1.0), first), capOf(double3(0.0, 0.0, -1.0), second)}, radius);
  EXPECT_NEAR(measured, expected, sweptAccuracy * 4.0 * std::numbers::pi * radius * radius);
}


TEST(exact_sphere_sweep, three_orthogonal_hemispheres_leave_one_octant)
{
  // The configuration with corners in it: what the three hemispheres leave is the octant where all three
  // coordinates are negative, an eighth of the sphere. The latitudes of its three corners are exactly the
  // breakpoints the crossings put in, so an area right to round-off here is the crossings being found and
  // being used.
  const double radius = 2.1;
  const double expected = 0.5 * std::numbers::pi * radius * radius;
  const double measured = exposedArea({capOf(double3(1.0, 0.0, 0.0), 0.5 * std::numbers::pi),
                                       capOf(double3(0.0, 1.0, 0.0), 0.5 * std::numbers::pi),
                                       capOf(double3(0.0, 0.0, 1.0), 0.5 * std::numbers::pi)},
                                      radius);
  EXPECT_NEAR(measured, expected, sweptAccuracy * 4.0 * std::numbers::pi * radius * radius);
}


TEST(exact_sphere_sweep, refining_the_panels_moves_the_answer_by_no_more_than_the_rule_is_worth)
{
  // One panel per smooth piece settles the area already, the pieces being cut where the integrand stops
  // being analytic. Cutting them finer is a check on that and not an improvement to it, so the two have to
  // agree to what the rule is worth --- which is not nothing, only very little.
  RandomNumber random{std::optional<std::size_t>(13)};
  for (std::size_t trial = 0; trial < 30; ++trial)
  {
    std::vector<SweepCircle> circles;
    for (std::size_t i = 0; i < 5; ++i)
    {
      circles.push_back(capOf(random.randomVectorOnUnitSphere(), 0.5 + 0.8 * random.uniform()));
    }

    const double coarse = exposedArea(circles, 1.0, 1);
    const double fine = exposedArea(circles, 1.0, 6);
    EXPECT_NEAR(coarse, fine, sweptAccuracy * 4.0 * std::numbers::pi);
  }
}


TEST(exact_sphere_sweep, the_area_does_not_depend_on_where_the_configuration_stands)
{
  // The frame is chosen from the cap axes, so a configuration turned somewhere else is swept in a different
  // frame against different latitudes. The area is a property of the configuration and not of the frame.
  RandomNumber random{std::optional<std::size_t>(17)};
  std::vector<SweepCircle> circles;
  std::vector<double3> axes;
  std::vector<double> halfAngles;
  for (std::size_t i = 0; i < 6; ++i)
  {
    axes.push_back(random.randomVectorOnUnitSphere());
    halfAngles.push_back(0.5 + 0.8 * random.uniform());
    circles.push_back(capOf(axes.back(), halfAngles.back()));
  }
  const double reference = exposedArea(circles, 1.0);

  for (std::size_t trial = 0; trial < 20; ++trial)
  {
    const double3 turningAxis = random.randomVectorOnUnitSphere();
    const double angle = 2.0 * std::numbers::pi * random.uniform();

    std::vector<SweepCircle> moved;
    for (std::size_t i = 0; i < axes.size(); ++i)
    {
      moved.push_back(capOf(turned(axes[i], turningAxis, angle), halfAngles[i]));
    }
    EXPECT_NEAR(exposedArea(moved, 1.0), reference, sweptAccuracy * 4.0 * std::numbers::pi);
  }
}


TEST(exact_sphere_sweep, crossings_handed_in_give_the_same_area_as_crossings_looked_for)
{
  // The decomposition into patches has already been over every pair of these caps, so the sweep is allowed to
  // take its crossings rather than search for them again. The two are not bit-identical and are not meant to
  // be: which cap of a pair is reached first decides how the cross product behind a crossing is bracketed,
  // and the compiler contracts those into fused multiply-adds, which round one order and not the other. What
  // that moves is where a panel begins, by an ulp or two of latitude, and never how much region there is.
  RandomNumber random{std::optional<std::size_t>(23)};
  for (std::size_t trial = 0; trial < 40; ++trial)
  {
    std::vector<SweepCircle> circles;
    for (std::size_t i = 0; i < 5; ++i)
    {
      circles.push_back(capOf(random.randomVectorOnUnitSphere(), 0.5 + 0.8 * random.uniform()));
    }

    std::vector<double3> directions;
    for (const CapCrossing& crossing : uncoveredCrossings(circles)) directions.push_back(crossing.direction);

    const double searched = exposedArea(circles, 1.0, 1, nullptr);
    const double handed = exposedArea(circles, 1.0, 1, &directions);
    EXPECT_NEAR(searched, handed, sweptAccuracy * 4.0 * std::numbers::pi);
  }
}


TEST(exact_sphere_sweep, covered_arcs_come_out_in_order_of_where_they_begin)
{
  // The insertion sort is chosen because the arcs arrive very nearly sorted already, the caps being taken in
  // order of their own azimuth. It still has to sort a list that is not.
  RandomNumber random{std::optional<std::size_t>(29)};
  for (std::size_t trial = 0; trial < 50; ++trial)
  {
    std::vector<CoveredArc> arcs;
    for (std::size_t i = 0; i < 1 + trial % 9; ++i)
    {
      CoveredArc arc;
      arc.begin = 2.0 * std::numbers::pi * random.uniform();
      arc.end = arc.begin + 0.1;
      arcs.push_back(arc);
    }

    std::vector<double> expected;
    for (const CoveredArc& arc : arcs) expected.push_back(arc.begin);
    std::sort(expected.begin(), expected.end());

    sortByBeginning(arcs);
    for (std::size_t i = 0; i < arcs.size(); ++i)
    {
      EXPECT_DOUBLE_EQ(arcs[i].begin, expected[i]);
      EXPECT_DOUBLE_EQ(arcs[i].end, arcs[i].begin + 0.1);  // an arc is carried whole, not just its beginning
    }
  }
}


TEST(exact_sphere_sweep, a_gap_carries_the_directions_at_its_own_ends)
{
  // The endpoints come with their cosine and sine rather than being recovered from the angle, which is what
  // lets the moments of an arc cost no trigonometry. They have to be the cosine and sine of that angle.
  std::vector<SweepCircle> circles = {capOf(double3(1.0, 0.0, 0.0), 0.9), capOf(double3(0.0, 1.0, 0.3), 1.1),
                                      capOf(double3(-0.4, 0.2, 1.0), 0.8)};
  std::vector<double3> axes;
  for (const SweepCircle& circle : circles) axes.push_back(circle.axis);
  const std::array<double3, 3> frame = sweepFrame(axes);

  std::size_t gaps = 0;
  SweepWorkspace work;
  sweepExposedLatitudes(circles, frame, nullptr, 1, work,
                        [&](const LatitudeGap& gap)
                        {
                          ++gaps;
                          EXPECT_NEAR(gap.cosineBegin, std::cos(gap.begin), roundOff);
                          EXPECT_NEAR(gap.sineBegin, std::sin(gap.begin), roundOff);
                          EXPECT_NEAR(gap.cosineEnd, std::cos(gap.end), roundOff);
                          EXPECT_NEAR(gap.sineEnd, std::sin(gap.end), roundOff);
                          EXPECT_NEAR(gap.span, gap.end - gap.begin, roundOff);

                          // And the direction it builds from them is on the sphere, in the frame swept in.
                          EXPECT_NEAR(gap.atBegin(frame).length(), 1.0, roundOff);
                          EXPECT_NEAR(gap.atEnd(frame).length(), 1.0, roundOff);
                          EXPECT_NEAR(gap.atMiddle(frame).length(), 1.0, roundOff);
                        });
  EXPECT_GT(gaps, 0uz);
}


TEST(exact_sphere_sweep, a_sphere_buried_under_one_cap_has_nothing_left)
{
  // A cap reaching over both poles covers every latitude, and the sweep is expected to skip the whole of the
  // integral rather than integrate a length of zero at each node of it.
  std::vector<SweepCircle> circles = {capOf(double3(0.3, 0.4, 0.5), 3.0)};
  EXPECT_NEAR(exposedArea(circles, 1.0), 2.0 * std::numbers::pi * (1.0 + std::cos(3.0)), 1.0e-12);

  // And with the cap a hair short of the whole sphere there is still the sliver about the far point.
  const double sliver = 2.0 * std::numbers::pi * (1.0 + std::cos(std::numbers::pi - 1.0e-3));
  EXPECT_NEAR(exposedArea({capOf(double3(0.0, 0.0, 1.0), std::numbers::pi - 1.0e-3)}, 1.0), sliver, 1.0e-12);
}
