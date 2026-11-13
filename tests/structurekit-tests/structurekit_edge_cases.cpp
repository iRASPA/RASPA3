#include <gtest/gtest.h>

import std;

import uint3;
import double3;
import double3x3;
import unit_cell;
import brute_force_structure;
import energy_shared_bet_surface_area;
import energy_shared_nldft;

// A 1x1x1 grid has only the k = 0 wave vector, so the largest |k| is zero: the attraction table was built and read
// as 0 / 0. The uniform bulk gas of an empty box needs only k = 0 and is still recovered.
TEST(structurekit_edge_cases, nldft_on_a_single_voxel_grid_holds_the_bulk_gas)
{
  NLDFTOptions options;
  options.pressurePoints = 2;
  options.xLow = 0.1;
  options.xHigh = 0.3;

  const uint3 gridSize{1, 1, 1};
  const UnitCell cell(20.0, 20.0, 20.0);
  const std::vector<float> energy(1, 0.0f);

  const NLDFTGridIsotherm result = nldftIsothermOnGrid(energy, gridSize, cell, options);
  ASSERT_EQ(result.isotherm.size(), 2u);
  for (const IsothermPoint& point : result.isotherm)
  {
    // NaN is detected through its printed form: the project flags fold std::isnan away
    ASSERT_EQ(std::to_string(point.moleculesPerCell).find("nan"), std::string::npos);
    const double expected = nldftGasDensity(options, result.bulk, point.relativePressure) * cell.volume;
    EXPECT_NEAR(point.moleculesPerCell, expected, 0.05 * expected);
  }
}

// A badly reduced cell needs more than eight images along its thin axis before the nearest atom is surely found: the
// shell was capped at eight, and the clearance was that of a farther image.
TEST(structurekit_edge_cases, brute_force_clearance_finds_the_nearest_image_of_a_skewed_cell)
{
  // the simple cubic lattice of spacing one, spanned by a = x, b = 20 x + y, c = z
  const UnitCell cell(double3x3(double3(1.0, 0.0, 0.0), double3(20.0, 1.0, 0.0), double3(0.0, 0.0, 1.0)));
  const BruteForceStructure structure = BruteForceStructure::create("skewed", cell, {double3(0.0, 0.0, 0.0)}, {0.0});

  // nearest lattice point (9, 0, 0), reached from the base image only nine cells along a
  const double3 point(9.3, 0.45, 0.2);
  EXPECT_NEAR(structure.clearance(point), std::sqrt(0.3 * 0.3 + 0.45 * 0.45 + 0.2 * 0.2), 1e-12);
}

// A site smaller than the overlap a negative bound allows has a negative limit, within which nothing is: there is room
// (a sphere of radius 0.2, a point 0.1 from its centre and an allowed overlap of 0.5). The limit is no longer squared
// into a positive one; the early exit of the image walk happened to hide that for the nearest image.
TEST(structurekit_edge_cases, brute_force_room_with_a_bound_deeper_than_the_radius)
{
  const BruteForceStructure structure =
      BruteForceStructure::create("point", UnitCell(10.0, 10.0, 10.0), {double3(5.0, 5.0, 5.0)}, {0.2});
  EXPECT_TRUE(structure.hasRoomFor(double3(5.1, 5.0, 5.0), -0.5));
  EXPECT_FALSE(structure.hasRoomFor(double3(5.1, 5.0, 5.0), -0.05));
  EXPECT_FALSE(structure.hasRoomFor(double3(5.5, 5.0, 5.0), 1.0));
}
