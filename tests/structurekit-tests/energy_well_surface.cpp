#include <gtest/gtest.h>

import std;

import uint3;
import double3;
import double3x3;
import unit_cell;
import crystal;
import surface_curvature;
import energy_isosurface;
import energy_shared_well_surface;

// The well surface, on landscapes whose wells can be written down.
//
// A Lennard-Jones pair crosses zero at sigma and bottoms out at 2^(1/6) sigma, so everything about this
// construction on an isolated atom is known in closed form: how far the walk goes, how deep it ends up, and
// what the map does to the area. The last of those is the sharpest test available, because the offset of a
// sphere is a similarity: every normal is radial and every walk is the same length, so the mapped mesh is the
// original mesh scaled, and the area ratio is exactly 2^(1/3) whatever marching cubes made of the sphere. Any
// error in the extraction cancels between numerator and denominator and what is left is the map.
//
// The slit tests are for the part that is not in the specification of the thing: telling a well a wall owns
// from one it shares with the wall facing it.

namespace
{

constexpr double ceiling = 1.0e7;

double lennardJones(double distance, double sigma, double epsilon)
{
  if (!(distance > 0.2 * sigma)) return ceiling;

  const double ratio = sigma / distance;
  const double six = ratio * ratio * ratio * ratio * ratio * ratio;

  return std::min(4.0 * epsilon * (six * six - six), ceiling);
}

// A grid filled in from a function of the Cartesian position of each node.
template <typename Landscape>
std::vector<float> fieldOf(const UnitCell &unitCell, uint3 gridSize, Landscape &&landscape)
{
  std::vector<float> field(static_cast<std::size_t>(gridSize.x) * gridSize.y * gridSize.z, 0.0f);

  for (std::size_t k = 0; k < gridSize.z; ++k)
  {
    for (std::size_t j = 0; j < gridSize.y; ++j)
    {
      for (std::size_t i = 0; i < gridSize.x; ++i)
      {
        const double3 fractional(static_cast<double>(i) / static_cast<double>(gridSize.x),
                                 static_cast<double>(j) / static_cast<double>(gridSize.y),
                                 static_cast<double>(k) / static_cast<double>(gridSize.z));

        field[(k * gridSize.y + j) * gridSize.x + i] = static_cast<float>(landscape(unitCell.cell * fractional));
      }
    }
  }

  return field;
}

// Where a decreasing function crosses zero, by bisection on the function itself rather than on the grid, so
// that a test of the walk starts exactly on the surface the walk is supposed to start on.
template <typename Landscape>
double zeroCrossing(Landscape &&landscape, double low, double high)
{
  for (int iteration = 0; iteration < 200; ++iteration)
  {
    const double middle = 0.5 * (low + high);
    if (landscape(middle) > 0.0)
    {
      low = middle;
    }
    else
    {
      high = middle;
    }
  }

  return 0.5 * (low + high);
}

}  // namespace

// The sampler reproduces a linear field exactly, value and gradient both, which is the least that can be asked
// of a trilinear interpolation and of the transform out of index space. Sampled away from the cell boundary,
// where a field that is linear rather than periodic has its jump.
TEST(energy_well_surface, samples_a_linear_field_exactly)
{
  const UnitCell unitCell(10.0, 10.0, 10.0);
  const uint3 gridSize{32, 32, 32};

  const double3 slope(3.0, 5.0, -2.0);
  std::vector<float> field =
      fieldOf(unitCell, gridSize, [&](double3 position) { return double3::dot(slope, position); });

  const PeriodicFieldSampler sampler{gridSize, field};
  ASSERT_TRUE(sampler.usable());

  const double3 fractional(0.41234, 0.52345, 0.46789);
  const double3 position = unitCell.cell * fractional;

  EXPECT_NEAR(sampler.value(fractional), double3::dot(slope, position), 1.0e-3);

  const double3 gradient =
      cartesianGradientOfGridGradient(unitCell.inverseCell, gridSize, sampler.gradient(fractional));

  EXPECT_NEAR(gradient.x, slope.x, 1.0e-4);
  EXPECT_NEAR(gradient.y, slope.y, 1.0e-4);
  EXPECT_NEAR(gradient.z, slope.z, 1.0e-4);
}

// The walk out of an isolated atom lands on the pair well: 2^(1/6) sigma from the centre, so
// (2^(1/6) - 1) sigma out from the zero crossing it started at, and one epsilon deep.
TEST(energy_well_surface, walks_to_the_lennard_jones_well)
{
  const double sigma = 3.0;
  const double epsilon = 100.0;
  const double edge = 12.0;

  const UnitCell unitCell(edge, edge, edge);
  const uint3 gridSize{192, 192, 192};
  const double3 centre = unitCell.cell * double3(0.5, 0.5, 0.5);

  std::vector<float> field = fieldOf(unitCell, gridSize,
                                     [&](double3 position)
                                     { return lennardJones((position - centre).length(), sigma, epsilon); });

  const PeriodicFieldSampler sampler{gridSize, field};
  ASSERT_TRUE(sampler.usable());

  const double3 outward(1.0, 0.0, 0.0);
  const double3 start = centre + outward * sigma;

  const WellWalk walk = walkToWell(sampler, unitCell, start, outward, edge / 192.0, 3.0, 0.0);

  ASSERT_TRUE(walk.found);
  EXPECT_FALSE(walk.atWall);
  EXPECT_FALSE(walk.shared);

  const double expected = (std::pow(2.0, 1.0 / 6.0) - 1.0) * sigma;
  EXPECT_NEAR(walk.distance, expected, 0.01);
  EXPECT_NEAR(walk.depth, -epsilon, 0.01 * epsilon);
}

// The whole map on an isolated atom. The zero surface is the sphere of radius sigma and the well surface is the
// sphere of radius 2^(1/6) sigma, so the area ratio is 2^(1/3) exactly --- the offset being a similarity, the
// mapped mesh is the extracted mesh scaled, and whatever marching cubes got wrong about the sphere divides out.
// The surface is convex everywhere and there is no second wall anywhere, so nothing is shared and nothing folds.
TEST(energy_well_surface, offsetting_an_isolated_atom_expands_it_by_the_cube_root_of_two)
{
  const double sigma = 3.0;
  const double epsilon = 100.0;
  const double edge = 12.0;
  const std::size_t points = 192;

  Crystal framework;
  framework.name = "one-atom";
  framework.unitCell = UnitCell(edge, edge, edge);
  framework.mass = 1.0;

  const uint3 gridSize{static_cast<std::uint32_t>(points), static_cast<std::uint32_t>(points),
                       static_cast<std::uint32_t>(points)};
  const double3 centre = framework.unitCell.cell * double3(0.5, 0.5, 0.5);

  std::vector<float> field = fieldOf(framework.unitCell, gridSize,
                                     [&](double3 position)
                                     { return lennardJones((position - centre).length(), sigma, epsilon); });

  std::vector<double3> corners = EnergyIsosurface::trianglesOfIsosurface(field, gridSize, 0.0);
  ASSERT_FALSE(corners.empty());

  // kT in the same units the landscape is in, epsilon being the only scale this test has.
  const double thermalEnergy = 0.5 * epsilon;
  const WellSurface surface = wellSurfaceOfField(framework, field, gridSize, corners, 0.0, thermalEnergy, 3.0);

  ASSERT_GT(surface.numberOfTriangles, 1000u);
  EXPECT_EQ(surface.numberOfUnmappedTriangles, 0u);
  EXPECT_EQ(surface.numberOfFoldedTriangles, 0u);
  EXPECT_EQ(surface.numberOfVerticesAtWall, 0u);
  EXPECT_LT(surface.sharedArea, 1.0e-6 * surface.area);

  // The zero surface itself, against 4 pi sigma^2, which is a check on the extraction rather than on the map.
  EXPECT_NEAR(surface.zeroArea, 4.0 * std::numbers::pi * sigma * sigma, 0.02 * surface.zeroArea);

  EXPECT_NEAR(surface.area / surface.zeroArea, std::pow(2.0, 1.0 / 3.0), 0.01);
  EXPECT_NEAR(surface.meanWalk, (std::pow(2.0, 1.0 / 6.0) - 1.0) * sigma, 0.01);
  EXPECT_NEAR(surface.meanDepth, -epsilon, 0.02 * epsilon);

  // A sphere seen from outside is convex over the whole of it, and the offset of a sphere is a sphere.
  ASSERT_GT(surface.curvature.classified(), 0.0);
  EXPECT_GT(surface.curvature.convexFraction(), 0.99);

  // Every well is one epsilon deep, so the weight is the same everywhere and the weighted area is the area
  // times that one number. This is the degenerate case the level-set route is stuck in permanently, and it is
  // worth pinning: on this landscape the weighting really does say nothing.
  const double weight = std::exp(epsilon / thermalEnergy);
  EXPECT_NEAR(surface.enhancement(), weight, 0.05 * weight);
}

// A slit narrow enough that the two walls have one well between them. The walk from either wall ends in the
// middle, and the well is reported as shared, which is what stops the pair of walls being counted as twice the
// surface that is really there.
TEST(energy_well_surface, a_narrow_slit_shares_one_well_between_its_walls)
{
  const double sigma = 3.0;
  const double epsilon = 100.0;
  const double gap = 2.1 * sigma;

  auto wall = [&](double x)
  {
    double total = 0.0;
    for (int image = -2; image <= 2; ++image)
    {
      total += lennardJones(std::abs(x - static_cast<double>(image) * gap), sigma, epsilon);
    }
    return std::min(total, ceiling);
  };

  const UnitCell unitCell(gap, 4.0, 4.0);
  const uint3 gridSize{512, 4, 4};

  std::vector<float> field = fieldOf(unitCell, gridSize, [&](double3 position) { return wall(position.x); });
  const PeriodicFieldSampler sampler{gridSize, field};
  ASSERT_TRUE(sampler.usable());

  const double start = zeroCrossing(wall, 0.5 * sigma, 0.5 * gap);
  const double3 outward(1.0, 0.0, 0.0);

  const WellWalk walk =
      walkToWell(sampler, unitCell, double3(start, 2.0, 2.0), outward, gap / 512.0, gap, 0.0);

  ASSERT_TRUE(walk.found);
  EXPECT_TRUE(walk.shared);

  // The gap is symmetric, so the one well is exactly halfway across it.
  EXPECT_NEAR(walk.distance, 0.5 * gap - start, 0.02);
}

// The same construction with the walls far enough apart to have a well each. The walk stops at this wall's own
// well, and carrying on from there crests a barrier before it finds anything else, so the well is this wall's
// own and its area is counted once.
TEST(energy_well_surface, a_wide_slit_gives_each_wall_a_well_of_its_own)
{
  const double sigma = 3.0;
  const double epsilon = 100.0;
  const double gap = 6.0 * sigma;

  auto wall = [&](double x)
  {
    double total = 0.0;
    for (int image = -2; image <= 2; ++image)
    {
      total += lennardJones(std::abs(x - static_cast<double>(image) * gap), sigma, epsilon);
    }
    return std::min(total, ceiling);
  };

  const UnitCell unitCell(gap, 4.0, 4.0);
  const uint3 gridSize{1024, 4, 4};

  std::vector<float> field = fieldOf(unitCell, gridSize, [&](double3 position) { return wall(position.x); });
  const PeriodicFieldSampler sampler{gridSize, field};
  ASSERT_TRUE(sampler.usable());

  const double start = zeroCrossing(wall, 0.5 * sigma, 0.5 * gap);
  const double3 outward(1.0, 0.0, 0.0);

  const WellWalk walk =
      walkToWell(sampler, unitCell, double3(start, 2.0, 2.0), outward, gap / 1024.0, 0.75 * gap, 0.0);

  ASSERT_TRUE(walk.found);
  EXPECT_FALSE(walk.shared);

  // Far from the other wall, so the well sits where an isolated wall would put it.
  EXPECT_NEAR(walk.distance, (std::pow(2.0, 1.0 / 6.0) - 1.0) * sigma, 0.05);
}

// A ray fired into a void with no wall within reach finds no well and says so, rather than stopping somewhere
// arbitrary. The field is flat, so there is no minimum anywhere along the ray.
TEST(energy_well_surface, a_ray_with_no_wall_within_reach_finds_nothing)
{
  const UnitCell unitCell(10.0, 10.0, 10.0);
  const uint3 gridSize{16, 16, 16};

  // Falling gently along x for ever, in the periodic sense: no minimum, and no wall to stop at.
  std::vector<float> field = fieldOf(unitCell, gridSize, [&](double3 position) { return -position.x; });

  const PeriodicFieldSampler sampler{gridSize, field};
  ASSERT_TRUE(sampler.usable());

  const WellWalk walk = walkToWell(sampler, unitCell, double3(1.0, 5.0, 5.0), double3(1.0, 0.0, 0.0), 0.1, 3.0, 0.0);

  EXPECT_FALSE(walk.found);
}
