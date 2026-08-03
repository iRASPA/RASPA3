#include <gtest/gtest.h>

import std;

import uint3;
import double3;
import double3x3;
import unit_cell;
import grid_area_curve;

// The coarea estimator, on fields whose level sets have areas that can be written down.
//
// The identity being relied on is that the integral of |grad A| over the points whose value falls in a band
// equals the integral of the area of the level set over that band. It is exact for a smooth field, so anything
// the estimator gets wrong here is either the arithmetic or the treatment of the cell, and the three cases
// below separate those: a field that varies along one axis tests the normalisation, a radial field tests the
// geometry, and the same radial field in an oblique cell tests that the gradient is carried from fractional
// coordinates to positions properly rather than by dividing each component by its own spacing.

namespace
{

double3 positionOf(const UnitCell &box, uint3 gridSize, std::size_t i, std::size_t j, std::size_t k)
{
  double3 fractional(static_cast<double>(i) / static_cast<double>(gridSize.x),
                     static_cast<double>(j) / static_cast<double>(gridSize.y),
                     static_cast<double>(k) / static_cast<double>(gridSize.z));
  return box.cell * fractional;
}

template <typename Function>
std::vector<float> sample(const UnitCell &box, uint3 gridSize, Function value)
{
  std::vector<float> field(static_cast<std::size_t>(gridSize.x) * gridSize.y * gridSize.z);
  for (std::size_t k = 0; k < gridSize.z; ++k)
  {
    for (std::size_t j = 0; j < gridSize.y; ++j)
    {
      for (std::size_t i = 0; i < gridSize.x; ++i)
      {
        field[(k * gridSize.y + j) * gridSize.x + i] =
            static_cast<float>(value(positionOf(box, gridSize, i, j, k)));
      }
    }
  }
  return field;
}

// The nearest image of a displacement, for any cell.
double3 nearestImage(const UnitCell &box, double3 displacement)
{
  double3 s = box.inverseCell * displacement;
  s.x -= std::round(s.x);
  s.y -= std::round(s.y);
  s.z -= std::round(s.z);
  return box.cell * s;
}

}  // namespace


// A field that varies along one axis only, as a sine so that it is periodic and the wrap costs nothing. Each
// level between the extremes is met twice per period, so every level set is a pair of planes and the area is
// twice the cross-section, the same at every level. A curve that is not flat here has its normalisation wrong.
TEST(grid_area_curve, a_sine_along_one_axis_gives_two_planes_at_every_level)
{
  const double side = 20.0;
  const UnitCell box(side, side, side);
  const uint3 gridSize(64, 64, 64);

  std::vector<float> field =
      sample(box, gridSize, [&](double3 position) { return std::sin(2.0 * std::numbers::pi * position.x / side); });

  AreaCurve curve = areaAgainstLevel(gridSize, box, field, -0.9, 0.9, 36);

  const double expected = 2.0 * side * side;
  ASSERT_EQ(curve.points.size(), 36u);
  for (const AreaCurvePoint &point : curve.points)
  {
    EXPECT_NEAR(point.area, expected, 0.02 * expected)
        << "at level " << point.level << " the two planes should carry " << expected << " square Angstrom";
  }
}


// A field that is the distance to a point. Its level sets are spheres, so the area at level a must be the area
// of a sphere of radius a, over more than an order of magnitude in a. This is the estimator's geometry rather
// than its bookkeeping: the gradient of a distance has magnitude one everywhere, so getting this wrong means
// the volume element or the bin width is wrong.
TEST(grid_area_curve, a_distance_to_a_point_gives_spheres_of_the_right_area)
{
  const double side = 24.0;
  const UnitCell box(side, side, side);
  const uint3 gridSize(120, 120, 120);
  const double3 centre(0.5 * side, 0.5 * side, 0.5 * side);

  std::vector<float> field =
      sample(box, gridSize, [&](double3 position) { return nearestImage(box, position - centre).length(); });

  // Kept clear of half the cell, past which the spheres begin to meet the faces and are no longer spheres.
  AreaCurve curve = areaAgainstLevel(gridSize, box, field, 1.0, 10.0, 90);

  for (const AreaCurvePoint &point : curve.points)
  {
    const double expected = 4.0 * std::numbers::pi * point.level * point.level;
    EXPECT_NEAR(point.area, expected, 0.04 * expected)
        << "the level set at " << point.level << " is a sphere of that radius";
  }

  // The same statement read the other way: the share of the cell below a level is the ball of that radius.
  const double atFive = 4.0 / 3.0 * std::numbers::pi * 125.0 / (side * side * side);
  auto below = [&](double level)
  {
    for (const AreaCurvePoint &point : curve.points)
    {
      if (point.level >= level) return point.fractionBelow;
    }
    return 1.0;
  };
  EXPECT_NEAR(below(5.0), atFive, 0.01);
}


// The same spheres in a cell with no right angles in it. The answer cannot change, a sphere being a sphere, so
// this isolates the one step that knows about the shape of the cell. Dividing each component of the gradient
// by its own grid spacing, which is right for a cube, is wrong here and is wrong by enough to see.
TEST(grid_area_curve, an_oblique_cell_does_not_change_the_area_of_a_sphere)
{
  const double side = 24.0;
  const double alpha = 75.0 * std::numbers::pi / 180.0;
  const double beta = 80.0 * std::numbers::pi / 180.0;
  const double gamma = 70.0 * std::numbers::pi / 180.0;

  // The standard upper-triangular cell for the three lengths and the three angles.
  const double ax = side;
  const double bx = side * std::cos(gamma);
  const double by = side * std::sin(gamma);
  const double cx = side * std::cos(beta);
  const double cy = side * (std::cos(alpha) - std::cos(beta) * std::cos(gamma)) / std::sin(gamma);
  const double cz = std::sqrt(side * side - cx * cx - cy * cy);

  const UnitCell box(double3x3(double3(ax, 0.0, 0.0), double3(bx, by, 0.0), double3(cx, cy, cz)));
  const uint3 gridSize(120, 120, 120);
  const double3 centre = box.cell * double3(0.5, 0.5, 0.5);

  std::vector<float> field =
      sample(box, gridSize, [&](double3 position) { return nearestImage(box, position - centre).length(); });

  // The perpendicular widths of an oblique cell are shorter than its edges, so the spheres run into the faces
  // sooner and the range has to be kept shorter than in the cube.
  AreaCurve curve = areaAgainstLevel(gridSize, box, field, 1.0, 8.0, 70);

  for (const AreaCurvePoint &point : curve.points)
  {
    const double expected = 4.0 * std::numbers::pi * point.level * point.level;
    EXPECT_NEAR(point.area, expected, 0.05 * expected)
        << "the cell being oblique cannot change the area of a sphere of radius " << point.level;
  }
}


// Reading the curve between its bins, which is what a comparison against a single marching-cubes area needs.
TEST(grid_area_curve, the_curve_can_be_read_off_between_its_bins)
{
  const double side = 24.0;
  const UnitCell box(side, side, side);
  const uint3 gridSize(96, 96, 96);
  const double3 centre(0.5 * side, 0.5 * side, 0.5 * side);

  std::vector<float> field =
      sample(box, gridSize, [&](double3 position) { return nearestImage(box, position - centre).length(); });

  AreaCurve curve = areaAgainstLevel(gridSize, box, field, 1.0, 9.0, 40);

  for (const double radius : {3.0, 5.5, 7.25})
  {
    const double expected = 4.0 * std::numbers::pi * radius * radius;
    EXPECT_NEAR(curve.areaAt(radius), expected, 0.05 * expected);
  }

  EXPECT_EQ(curve.areaAt(0.5), 0.0) << "below the range there is nothing to report";
  EXPECT_EQ(curve.areaAt(20.0), 0.0) << "nor above it";
}
