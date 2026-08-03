#include <gtest/gtest.h>

import std;

import uint3;
import double3;
import unit_cell;
import grid_pore_size;

// The two halves of a pore-size distribution, on fields whose answers are known in closed form.
//
// The energy route cannot be checked against the geometric one directly: the two bound their voids by
// different things, a level set of an energy against the surfaces of hard spheres, and there is no level at
// which a Lennard-Jones probe and a hard sphere of any radius occupy the same region. What can be checked is
// the arithmetic that sits between a field and a curve, and that is what is done here, by handing it fields
// that are distances already.
//
// A signed distance is the fixed point of `distanceToIsosurface`: hand it one and it must hand the same one
// back, since the distance from a point to the zero level of a distance is the distance itself. That is the
// sharpest test available, because it compares against an exact answer rather than against another
// approximation, and it is the sense in which the energy route reduces to the geometric one, whose clearance
// field is a signed distance to the surfaces it was built from.

namespace
{

// The field is sampled endpoint-exclusive, as every grid in the structure-property routes is.
double3 positionOf(const UnitCell &box, uint3 gridSize, std::size_t i, std::size_t j, std::size_t k)
{
  double3 fractional(static_cast<double>(i) / static_cast<double>(gridSize.x),
                     static_cast<double>(j) / static_cast<double>(gridSize.y),
                     static_cast<double>(k) / static_cast<double>(gridSize.z));
  return box.cell * fractional;
}

// Fills a field from a function of position, in the grid's own order.
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

// The nearest image of a displacement in a cell that is a cube of this side.
double3 nearestImage(double3 displacement, double side)
{
  return double3(displacement.x - side * std::round(displacement.x / side),
                 displacement.y - side * std::round(displacement.y / side),
                 displacement.z - side * std::round(displacement.z / side));
}

}  // namespace


// A slab of open space between two planes. The boundary is flat and lies along the grid, so interpolating
// along an edge is exact there, and nothing but the search itself is being tested: the distance has to come
// back to the precision of the arithmetic and not to the precision of the grid.
TEST(grid_pore_size, a_slab_gives_back_its_own_distance_to_within_rounding)
{
  const double side = 20.0;

  // Half a spacing off a whole number of steps, so that the boundary falls between grid points and the
  // interpolation is asked to do something. On a field that runs straight along the edge it crosses, which
  // this one does, it should still land on the boundary exactly.
  const double halfWidth = 4.25;
  const UnitCell box(side, side, side);
  const uint3 gridSize(40, 40, 40);

  // Positive within `halfWidth` of the plane z = side/2, negative outside it.
  auto slab = [&](double3 position) { return halfWidth - std::abs(position.z - 0.5 * side); };
  std::vector<float> openness = sample(box, gridSize, slab);

  std::vector<float> distance = distanceToIsosurface(gridSize, box, openness, 0.0);

  double worst = 0.0;
  for (std::size_t k = 0; k < gridSize.z; ++k)
  {
    for (std::size_t j = 0; j < gridSize.y; ++j)
    {
      for (std::size_t i = 0; i < gridSize.x; ++i)
      {
        const std::size_t voxel = (k * gridSize.y + j) * gridSize.x + i;
        const double expected = slab(positionOf(box, gridSize, i, j, k));
        if (expected < 0.0)
        {
          EXPECT_LT(distance[voxel], 0.0f) << "a point outside the slab was called void";
          continue;
        }
        worst = std::max(worst, std::abs(static_cast<double>(distance[voxel]) - expected));
      }
    }
  }

  EXPECT_LT(worst, 1.0e-5) << "the distance to a flat boundary should be exact, and is off by " << worst;
}


// A sphere of open space. The boundary curves away from the grid, so a crossing placed on an edge is a little
// outside the sphere and the distance a little too large, by an amount that goes as the square of the
// spacing. What is checked is that the error is of that size and one-signed, since a distance that came back
// too small would mean the crossings are being placed by interpolating over more than one edge, which on an
// energy field is the difference between a sensible answer and a meaningless one.
TEST(grid_pore_size, a_sphere_gives_back_its_own_distance_to_within_the_spacing)
{
  const double side = 20.0;
  const double radius = 6.0;
  const UnitCell box(side, side, side);
  const uint3 gridSize(60, 60, 60);
  const double3 centre(0.5 * side, 0.5 * side, 0.5 * side);
  const double spacing = side / 60.0;

  auto sphere = [&](double3 position) { return radius - nearestImage(position - centre, side).length(); };
  std::vector<float> openness = sample(box, gridSize, sphere);

  std::vector<float> distance = distanceToIsosurface(gridSize, box, openness, 0.0);

  double worstAbove = 0.0;
  double worstBelow = 0.0;
  for (std::size_t k = 0; k < gridSize.z; ++k)
  {
    for (std::size_t j = 0; j < gridSize.y; ++j)
    {
      for (std::size_t i = 0; i < gridSize.x; ++i)
      {
        const std::size_t voxel = (k * gridSize.y + j) * gridSize.x + i;
        const double expected = sphere(positionOf(box, gridSize, i, j, k));
        if (expected < 0.0) continue;

        const double error = static_cast<double>(distance[voxel]) - expected;
        worstAbove = std::max(worstAbove, error);
        worstBelow = std::max(worstBelow, -error);
      }
    }
  }

  EXPECT_LT(worstAbove, spacing) << "the distance inside a sphere is too large by " << worstAbove;
  EXPECT_LT(worstBelow, 0.1 * spacing) << "the distance inside a sphere is too small by " << worstBelow
                                       << ", which a crossing placed on a single edge cannot do";
}


// The largest sphere that fits in a spherical void and covers any point of it is that void itself, so the
// Gelb-Gubbins field is flat at the radius throughout, and the distribution is a single spike at the
// diameter. It is the one void whose whole curve is known without computing anything.
TEST(grid_pore_size, every_point_of_a_spherical_pore_is_given_the_size_of_the_sphere)
{
  const double side = 20.0;
  const double radius = 6.0;
  const UnitCell box(side, side, side);
  const uint3 gridSize(60, 60, 60);
  const double3 centre(0.5 * side, 0.5 * side, 0.5 * side);
  const double spacing = side / 60.0;

  auto sphere = [&](double3 position) { return radius - nearestImage(position - centre, side).length(); };
  std::vector<float> openness = sample(box, gridSize, sphere);

  std::vector<float> distance = distanceToIsosurface(gridSize, box, openness, 0.0);
  const double slack = coveringSlack(box, gridSize);
  std::vector<float> poreRadius = poreRadiusField(gridSize, box, distance, slack);

  std::size_t voidVoxels = 0;
  double smallest = std::numeric_limits<double>::max();
  double largest = 0.0;
  for (std::size_t voxel = 0; voxel < distance.size(); ++voxel)
  {
    if (distance[voxel] < 0.0f) continue;
    ++voidVoxels;
    smallest = std::min(smallest, static_cast<double>(poreRadius[voxel]));
    largest = std::max(largest, static_cast<double>(poreRadius[voxel]));
  }

  EXPECT_GT(voidVoxels, 0u);
  EXPECT_NEAR(largest, radius, spacing) << "the widest point of the pore is not the radius";
  EXPECT_NEAR(smallest, radius, 2.0 * spacing)
      << "some point of the pore was given a size smaller than the pore, which is what the slack is for";

  // The volume of the void is the volume of the sphere, which is what says the level set was found in the
  // right place to begin with.
  const double measured = static_cast<double>(voidVoxels) / static_cast<double>(distance.size());
  const double expected = (4.0 / 3.0) * std::numbers::pi * radius * radius * radius / (side * side * side);
  EXPECT_NEAR(measured, expected, 0.01 * expected);
}


// Two slabs of different widths, one twice the other, stacked in the same cell. The curve then has to have
// exactly two steps, at the two widths, carrying weight in proportion to the volumes. It is the smallest
// case in which the binning is doing anything, and it checks the cumulative curve is a share of the void
// rather than of the cell.
TEST(grid_pore_size, two_slabs_give_two_steps_in_proportion_to_their_volumes)
{
  const double side = 24.0;
  const UnitCell box(side, side, side);
  const uint3 gridSize(48, 48, 48);

  // One slab of half-width 2 about z = 6, another of half-width 4 about z = 18.
  auto slabs = [&](double3 position)
  {
    const double narrow = 2.0 - std::abs(position.z - 6.0);
    const double wide = 4.0 - std::abs(position.z - 18.0);
    return std::max(narrow, wide);
  };
  std::vector<float> openness = sample(box, gridSize, slabs);

  std::vector<float> distance = distanceToIsosurface(gridSize, box, openness, 0.0);
  const double slack = coveringSlack(box, gridSize);
  std::vector<float> poreRadius = poreRadiusField(gridSize, box, distance, slack);

  PoreSizeCurve curve = binPoreSizes(poreRadius, distance, {}, {}, 16.0, 320);

  // A slab of half-width w holds spheres of diameter 2w, so the two sizes are 4 and 8, and the wide slab
  // holds twice the volume of the narrow one.
  auto shareAtLeast = [&](double diameter)
  {
    for (const PoreSizeCurvePoint &point : curve.points)
    {
      if (point.diameter >= diameter) return point.cumulativeVoidFraction;
    }
    return 0.0;
  };

  EXPECT_NEAR(curve.largestDiameter, 8.0, 0.2);
  EXPECT_NEAR(shareAtLeast(1.0), 1.0, 0.02) << "every point of the void is in a pore at least 1 across";
  EXPECT_NEAR(shareAtLeast(5.0), 2.0 / 3.0, 0.03) << "the wide slab holds two thirds of the void";
  EXPECT_NEAR(shareAtLeast(9.0), 0.0, 0.02) << "nothing is wider than the wide slab";
  EXPECT_NEAR(curve.meanDiameter, (4.0 + 2.0 * 8.0) / 3.0, 0.2) << "the mean size is the volume-weighted one";
}


// The same two slabs, with the narrow one given four times the weight per point of the wide one. Counting
// volume, the wide slab holds two thirds; counting weight, the narrow one now does. This is the whole of what
// the Boltzmann weighting does on a real landscape, with the factor of four standing in for the exp(-A/kT)
// of a place the molecule likes, so it is worth checking on a case where the answer can be worked out by
// hand.
TEST(grid_pore_size, weighting_the_narrow_slab_moves_the_curve_onto_it)
{
  const double side = 24.0;
  const UnitCell box(side, side, side);
  const uint3 gridSize(48, 48, 48);

  auto slabs = [&](double3 position)
  {
    const double narrow = 2.0 - std::abs(position.z - 6.0);
    const double wide = 4.0 - std::abs(position.z - 18.0);
    return std::max(narrow, wide);
  };
  std::vector<float> openness = sample(box, gridSize, slabs);

  std::vector<float> distance = distanceToIsosurface(gridSize, box, openness, 0.0);
  const double slack = coveringSlack(box, gridSize);
  std::vector<float> poreRadius = poreRadiusField(gridSize, box, distance, slack);

  // Four to a point in the narrow slab, one to a point in the wide one, by which half of the cell the point
  // is in.
  std::vector<double> weight(distance.size(), 0.0);
  for (std::size_t k = 0; k < gridSize.z; ++k)
  {
    const double z = side * static_cast<double>(k) / static_cast<double>(gridSize.z);
    for (std::size_t j = 0; j < gridSize.y; ++j)
    {
      for (std::size_t i = 0; i < gridSize.x; ++i)
      {
        weight[i + gridSize.x * (j + gridSize.y * k)] = (z < 12.0) ? 4.0 : 1.0;
      }
    }
  }

  PoreSizeCurve curve = binPoreSizes(poreRadius, distance, {}, weight, 16.0, 320);

  auto shareAtLeast = [&](double diameter)
  {
    for (const PoreSizeCurvePoint &point : curve.points)
    {
      if (point.diameter >= diameter) return point.cumulativeVoidFraction;
    }
    return 0.0;
  };

  // The narrow slab is a third of the void at weight four, the wide one two thirds at weight one, so the
  // wide slab keeps 2 / (4 + 2) of the weight.
  EXPECT_NEAR(shareAtLeast(1.0), 1.0, 0.02) << "the curve is still a share of the whole, whatever the weights";
  EXPECT_NEAR(shareAtLeast(5.0), 1.0 / 3.0, 0.03) << "weighting the narrow slab takes the curve off the wide one";
  EXPECT_NEAR(curve.meanDiameter, (4.0 * 4.0 + 2.0 * 8.0) / 6.0, 0.2) << "and the mean size comes down with it";
}
