#include <gtest/gtest.h>

#include <cstring>

import std;

import double3;
import double3x3;
import unit_cell;
import randomnumbers;
import opencl;
import sampled_structure;
import mc_surface_area;
import mc_opencl_surface_area;
import integration_surface_area;
import integration_opencl_surface_area;

// The areas of spheres, by points thrown at them and by a grid of points over them, on the processor and on the
// device. Arrangements of overlapping spheres have areas known in closed form, whichever route is asked.

namespace
{

SampledStructure structureOf(const UnitCell &cell, const std::vector<double3> &positions, double radius)
{
  SampledStructure structure;
  structure.name = "sampled_surface_area_test";
  structure.unitCell = cell;
  for (const double3 &position : positions)
  {
    structure.positions.push_back(position);
    structure.radii.push_back(radius);
  }
  return structure;
}

struct ReportCleanup
{
  ~ReportCleanup()
  {
    for (const std::filesystem::directory_entry &entry : std::filesystem::directory_iterator("."))
    {
      if (entry.path().filename().string().starts_with("sampled_surface_area_test.")) std::filesystem::remove(entry.path());
    }
  }
};

double integrated(const SampledStructure &structure, std::size_t slices)
{
  Integration_SurfaceArea route;
  route.run(structure, SampledProbe{}, slices);
  return route.surfaceArea;
}

double integratedOnDevice(const SampledStructure &structure, std::size_t slices)
{
  Integration_OpenCL_SurfaceArea route;
  route.run(structure, SampledProbe{}, slices);
  return route.surfaceArea;
}

}  // namespace


// Two spheres a radius apart along the axis of the grid each lose a cap of a quarter of the sphere, and that plane is
// a whole number of stacks from a pole. So a grid of equal shares gets the area exactly. Stacks placed on the poles
// give a pole a stack's worth of points for an area of nothing, and the answer is off by far more than the grid's
// own resolution.
TEST(sampled_surface_area, the_grid_gives_every_point_an_equal_share_of_the_sphere)
{
  ReportCleanup cleanup;
  const double radius = 2.0;
  SampledStructure pair = structureOf(UnitCell(30.0, 30.0, 30.0), {double3(15.0, 15.0, 14.0), double3(15.0, 15.0, 16.0)}, radius);
  const double expected = 2.0 * 4.0 * std::numbers::pi * radius * radius * 0.75;

  EXPECT_NEAR(integrated(pair, 64), expected, 1.0e-9 * expected);

  OpenCL::initialize();
  if (OpenCL::clContext.has_value()) EXPECT_NEAR(integratedOnDevice(pair, 64), expected, 1.0e-3 * expected);
}

// A chain of spheres through the cell, each cut by its neighbours in two planes a repeat apart. What is left of each
// is a band of area 2 pi r times the repeat, the same whether its nearest neighbours bury it or further ones do too,
// and none of those neighbours is anything but an image of the sphere itself.
TEST(sampled_surface_area, a_sphere_is_buried_by_its_own_images)
{
  ReportCleanup cleanup;
  OpenCL::initialize();

  for (double repeat : {4.0, 2.0})
  {
    const double radius = 3.0;
    SampledStructure chain = structureOf(UnitCell(repeat, 30.0, 30.0), {double3(0.0, 15.0, 15.0)}, radius);
    const double expected = 2.0 * std::numbers::pi * radius * repeat;
    EXPECT_NEAR(integrated(chain, 512), expected, 5.0e-3 * expected) << "repeat " << repeat;

    MC_SurfaceArea sampled;
    sampled.run(chain, SampledProbe{}, 20, 2000);
    EXPECT_NEAR(sampled.surfaceArea, expected, 0.03 * expected) << "repeat " << repeat;

    if (!OpenCL::clContext.has_value()) continue;
    EXPECT_NEAR(integratedOnDevice(chain, 512), expected, 5.0e-3 * expected) << "repeat " << repeat;

    MC_OpenCL_SurfaceArea sampledOnDevice;
    sampledOnDevice.run(chain, SampledProbe{}, 20, 2000);
    EXPECT_NEAR(sampledOnDevice.surfaceArea, expected, 0.03 * expected) << "repeat " << repeat;
  }
}

// The device against the processor, on the grid, where there is no scatter between them: cells square and skewed,
// including one whose edges are all long and whose width across is not.
TEST(sampled_surface_area, the_device_agrees_with_the_processor_in_any_cell)
{
  OpenCL::initialize();
  if (!OpenCL::clContext.has_value()) GTEST_SKIP() << "no OpenCL device";
  ReportCleanup cleanup;

  struct Case
  {
    const char *name;
    UnitCell cell;
    std::size_t atoms;
    double radius;
  };
  const std::vector<Case> cases = {
      {"cube", UnitCell(10.0, 10.0, 10.0), 10, 1.6},
      {"skew", UnitCell(double3x3(double3(10.0, 0.0, 0.0), double3(4.0, 9.0, 0.0), double3(3.0, 3.5, 11.0))), 10, 1.6},
      // The edges are 9 Å and more, and yet the shortest repeat is 4.2 Å, so this sphere is buried by its own images
      // though it is less than half as wide as any edge.
      {"long edges, narrow", UnitCell(double3x3(double3(9.0, 0.0, 0.0), double3(8.0, 5.0, 0.0), double3(7.5, 4.0, 4.0))), 1, 2.4}};

  RandomNumber random{std::optional<std::size_t>(21)};
  for (const Case &entry : cases)
  {
    std::vector<double3> positions;
    for (std::size_t i = 0; i < entry.atoms; ++i)
      positions.push_back(entry.cell.cell * double3(random.uniform(), random.uniform(), random.uniform()));
    SampledStructure structure = structureOf(entry.cell, positions, entry.radius);

    const double expected = integrated(structure, 128);
    const double whole = static_cast<double>(entry.atoms) * 4.0 * std::numbers::pi * entry.radius * entry.radius;
    ASSERT_GT(expected, 0.01 * whole) << entry.name;
    ASSERT_LT(expected, 0.99 * whole) << entry.name;
    EXPECT_NEAR(integratedOnDevice(structure, 128), expected, 3.0e-3 * expected) << entry.name;

    MC_OpenCL_SurfaceArea sampledOnDevice;
    sampledOnDevice.run(structure, SampledProbe{}, 20, 2000);
    EXPECT_NEAR(sampledOnDevice.surfaceArea, expected, std::max(4.0 * sampledOnDevice.surfaceAreaError, 0.02 * expected))
        << entry.name;
  }
}

// One sphere with nothing near it reads the same on every pass, and the same numbers added up and squared can come out
// with the variance a rounding below zero. There is no scatter, so the error is nothing, and it is a number.
TEST(sampled_surface_area, passes_that_all_read_the_same_have_no_error)
{
  ReportCleanup cleanup;
  OpenCL::initialize();

  // The build is entitled to assume that nothing is ever not a number, so `isnan` may be answered without looking, and
  // so may a mask of the bits once they have been through a double. So the bytes are copied out as they are stored.
  auto notFinite = [](const double &stored)
  {
    std::uint64_t bits = 0;
    std::memcpy(&bits, &stored, sizeof(bits));
    return (bits & 0x7ff0000000000000ull) == 0x7ff0000000000000ull;
  };

  std::size_t notANumber = 0;
  for (std::size_t i = 0; i < 500; ++i)
  {
    const double radius = 1.0 + 0.00737 * static_cast<double>(i);
    SampledStructure alone = structureOf(UnitCell(30.0, 30.0, 30.0), {double3(15.0, 15.0, 15.0)}, radius);

    MC_SurfaceArea route;
    route.run(alone, SampledProbe{}, 3 + i % 37, 1);
    const double area = 4.0 * std::numbers::pi * radius * radius;
    EXPECT_NEAR(route.surfaceArea, area, 1.0e-9 * area);
    if (notFinite(route.surfaceAreaError)) ++notANumber;
    else EXPECT_LT(route.surfaceAreaError, 1.0e-6 * area);
  }
  EXPECT_EQ(notANumber, 0uz);
}
