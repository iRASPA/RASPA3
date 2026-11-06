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


// An energy field at a wall: a slab of pore voxels at one value beside voxels held at a ceiling orders of
// magnitude above the range. Every level between the pore's value and the ceiling is crossed within the voxel
// beside the wall, so each such level is two planes, the slab having two faces, and no level below the pore's
// value is crossed anywhere. Measured by a band symmetric about the voxel's value, the voxel beside the wall was
// credited with every level under the range as well, half its weight was counted as lying below the range, and
// the resolution reported was the height of the wall.
TEST(grid_area_curve, a_wall_is_credited_only_with_the_levels_it_climbs_through)
{
  const double side = 20.0;
  const UnitCell box(side, side, side);
  const uint3 gridSize(32, 32, 32);

  std::vector<float> field(32 * 32 * 32, -10.0f);
  for (std::size_t k = 16; k < 32; ++k)
  {
    for (std::size_t j = 0; j < 32; ++j)
    {
      for (std::size_t i = 0; i < 32; ++i) field[(k * 32 + j) * 32 + i] = 1.0e7f;
    }
  }

  AreaCurve curve = areaAgainstLevel(gridSize, box, field, -500.0, 500.0, 100);

  const double twoPlanes = 2.0 * side * side;
  EXPECT_NEAR(curve.areaAt(200.0), twoPlanes, 0.05 * twoPlanes);
  EXPECT_NEAR(curve.areaAt(-300.0), 0.0, 1.0e-9) << "no level below the pore's value is crossed anywhere";
  EXPECT_NEAR(curve.fractionBelowRange, 0.0, 1.0e-12);
  EXPECT_LE(curve.meanSpanPerVoxel, 3.0 * 1000.0) << "the resolution is read within the range, not off the wall";
}


// The same wall held at whatever a float can hold. What a voxel beside the wall gives each level does not depend
// on how high the wall stands: its gradient and the band of values it covers both grow with the wall, and what a
// bin gets is the one over the other. So the curve has to be the same for a wall at ten million as for one at
// the largest float there is.
//
// Which bins the band reaches is worked out by dividing how far the band extends by the width of a bin, and a
// band that reaches 1e38 past a range a thousand wide is 1e37 bins long, which no integer holds. Turned into one
// regardless, the number is whatever the machine makes of it: the largest integer on some, which happens to be
// right, and the most negative on others, which puts the last bin before the first and leaves the voxel out of
// every bin. The wall is then there in the field and gone from the curve.
TEST(grid_area_curve, a_wall_as_high_as_a_float_goes_is_the_same_wall)
{
  const double side = 20.0;
  const UnitCell box(side, side, side);
  const uint3 gridSize(32, 32, 32);
  const double twoPlanes = 2.0 * side * side;

  auto curveWithWallAt = [&](float ceiling)
  {
    std::vector<float> field(32 * 32 * 32, -10.0f);
    for (std::size_t k = 16; k < 32; ++k)
    {
      for (std::size_t j = 0; j < 32; ++j)
      {
        for (std::size_t i = 0; i < 32; ++i) field[(k * 32 + j) * 32 + i] = ceiling;
      }
    }
    return areaAgainstLevel(gridSize, box, field, -500.0, 500.0, 100);
  };

  const AreaCurve reference = curveWithWallAt(1.0e7f);
  ASSERT_NEAR(reference.areaAt(200.0), twoPlanes, 0.05 * twoPlanes);

  for (const float ceiling : {1.0e12f, 1.0e20f, 1.0e30f, std::numeric_limits<float>::max()})
  {
    const AreaCurve curve = curveWithWallAt(ceiling);
    ASSERT_EQ(curve.points.size(), reference.points.size());
    for (const double level : {-5.0, 5.0, 200.0, 495.0})
    {
      EXPECT_NEAR(curve.areaAt(level), reference.areaAt(level), 1.0e-4 * twoPlanes)
          << "wall at " << ceiling << ", level " << level;
    }
    EXPECT_NEAR(curve.areaAt(-300.0), 0.0, 1.0e-9) << "wall at " << ceiling;
    EXPECT_NEAR(curve.fractionBelowRange, 0.0, 1.0e-12) << "wall at " << ceiling;
    EXPECT_NEAR(curve.fractionAboveRange, reference.fractionAboveRange, 1.0e-4) << "wall at " << ceiling;
    EXPECT_NEAR(curve.points.back().fractionBelow, reference.points.back().fractionBelow, 1.0e-4)
        << "wall at " << ceiling;
  }
}

// And the other way down: a pit as deep as a float goes beside voxels in the range. Every level between the
// pit and the voxels' own value is crossed in the voxel beside it, so each of those is two planes, and none
// above that value is crossed anywhere. Here it is the first bin that is worked out from a number no integer
// holds.
TEST(grid_area_curve, a_pit_as_deep_as_a_float_goes_is_the_same_pit)
{
  const double side = 20.0;
  const UnitCell box(side, side, side);
  const uint3 gridSize(32, 32, 32);
  const double twoPlanes = 2.0 * side * side;

  auto curveWithPitAt = [&](float floor)
  {
    std::vector<float> field(32 * 32 * 32, 10.0f);
    for (std::size_t k = 16; k < 32; ++k)
    {
      for (std::size_t j = 0; j < 32; ++j)
      {
        for (std::size_t i = 0; i < 32; ++i) field[(k * 32 + j) * 32 + i] = floor;
      }
    }
    return areaAgainstLevel(gridSize, box, field, -500.0, 500.0, 100);
  };

  const AreaCurve reference = curveWithPitAt(-1.0e7f);
  ASSERT_NEAR(reference.areaAt(-200.0), twoPlanes, 0.05 * twoPlanes);

  for (const float floor : {-1.0e12f, -1.0e20f, -1.0e30f, -std::numeric_limits<float>::max()})
  {
    const AreaCurve curve = curveWithPitAt(floor);
    for (const double level : {-495.0, -200.0, -5.0, 5.0})
    {
      EXPECT_NEAR(curve.areaAt(level), reference.areaAt(level), 1.0e-4 * twoPlanes)
          << "pit at " << floor << ", level " << level;
    }
    EXPECT_NEAR(curve.areaAt(300.0), 0.0, 1.0e-9) << "pit at " << floor;
    EXPECT_NEAR(curve.fractionAboveRange, 0.0, 1.0e-12) << "pit at " << floor;
    EXPECT_NEAR(curve.fractionBelowRange, reference.fractionBelowRange, 1.0e-4) << "pit at " << floor;
  }
}


// The bins span the whole of the range and their centres stop half a bin short of either end of it, but a level in
// that last half bin is in the range all the same. It is in the end bin, and has that bin's area, not none.
TEST(grid_area_curve, the_last_half_bin_at_either_end_of_the_range_is_read_from_the_end_bin)
{
  AreaCurve curve;
  curve.lowestLevel = 0.0;
  curve.highestLevel = 4.0;
  curve.binWidth = 1.0;
  curve.points = {AreaCurvePoint{0.5, 10.0, 0.0}, AreaCurvePoint{1.5, 20.0, 0.0}, AreaCurvePoint{2.5, 30.0, 0.0},
                  AreaCurvePoint{3.5, 40.0, 0.0}};

  EXPECT_DOUBLE_EQ(curve.areaAt(0.0), 10.0);
  EXPECT_DOUBLE_EQ(curve.areaAt(0.25), 10.0);
  EXPECT_DOUBLE_EQ(curve.areaAt(0.5), 10.0);
  EXPECT_DOUBLE_EQ(curve.areaAt(1.0), 15.0);
  EXPECT_DOUBLE_EQ(curve.areaAt(3.5), 40.0);
  EXPECT_DOUBLE_EQ(curve.areaAt(3.75), 40.0);
  EXPECT_DOUBLE_EQ(curve.areaAt(4.0), 40.0);

  // Outside the range there is still nothing to read.
  EXPECT_DOUBLE_EQ(curve.areaAt(-0.001), 0.0);
  EXPECT_DOUBLE_EQ(curve.areaAt(4.001), 0.0);
  EXPECT_DOUBLE_EQ(curve.areaAt(std::numeric_limits<double>::quiet_NaN()), 0.0);
}

// One bin is one number for the whole of the range, not a number at its middle and nothing either side.
TEST(grid_area_curve, a_single_bin_answers_for_the_whole_range)
{
  AreaCurve curve;
  curve.lowestLevel = 2.0;
  curve.highestLevel = 6.0;
  curve.binWidth = 4.0;
  curve.points = {AreaCurvePoint{4.0, 7.0, 0.0}};

  for (double level : {2.0, 3.0, 4.0, 5.0, 6.0}) EXPECT_DOUBLE_EQ(curve.areaAt(level), 7.0) << "at " << level;
  EXPECT_DOUBLE_EQ(curve.areaAt(1.999), 0.0);
  EXPECT_DOUBLE_EQ(curve.areaAt(6.001), 0.0);
}

// The same off a curve the pass built: a sine reads about the same at the very end of a range as just inside it.
TEST(grid_area_curve, a_real_curve_reads_at_the_ends_of_its_range)
{
  const uint3 gridSize(64, 64, 64);
  const double edge = 20.0;
  UnitCell box(edge, edge, edge);
  std::vector<float> field(64uz * 64uz * 64uz);
  for (std::size_t k = 0; k < 64; ++k)
    for (std::size_t j = 0; j < 64; ++j)
      for (std::size_t i = 0; i < 64; ++i)
        field[(k * 64 + j) * 64 + i] = static_cast<float>(std::sin(2.0 * std::numbers::pi * static_cast<double>(i) / 64.0));

  AreaCurve curve = areaAgainstLevel(gridSize, box, field, -0.5, 0.5, 20);
  ASSERT_EQ(curve.points.size(), 20uz);
  EXPECT_GT(curve.areaAt(-0.5), 0.0);
  EXPECT_GT(curve.areaAt(0.5), 0.0);
  EXPECT_DOUBLE_EQ(curve.areaAt(-0.5), curve.points.front().area);
  EXPECT_DOUBLE_EQ(curve.areaAt(0.5), curve.points.back().area);
}
