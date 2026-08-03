#include <gtest/gtest.h>

import std;

import int3;
import uint3;
import double3;
import unit_cell;
import grid_connected_components;
import grid_blocking_cover;

// Covering pockets with spheres, on regions whose covering can be worked out by hand.
//
// The covering is shared between the geometric route and the energy one, which differ only in what they hand
// it: which regions are sealed, where the molecule is entitled to be, and what field says how much room there
// is at a point. So the cases here are stated in those terms rather than in terms of either route.

namespace
{

// A field of zero everywhere, into which regions are painted, plus the components of it above a level.
struct Region
{
  uint3 gridSize;
  UnitCell box;
  std::vector<float> openness;

  Region(uint3 size, double side) : gridSize(size), box(side, side, side)
  {
    this->openness.assign(static_cast<std::size_t>(size.x) * size.y * size.z, 0.0f);
  }

  void paint(int3 low, int3 high, float value)
  {
    for (std::size_t k = static_cast<std::size_t>(low.z); k <= static_cast<std::size_t>(high.z); ++k)
    {
      for (std::size_t j = static_cast<std::size_t>(low.y); j <= static_cast<std::size_t>(high.y); ++j)
      {
        for (std::size_t i = static_cast<std::size_t>(low.x); i <= static_cast<std::size_t>(high.x); ++i)
        {
          this->openness[(k * this->gridSize.y + j) * this->gridSize.x + i] = value;
        }
      }
    }
  }

  GridComponents componentsAbove(double level) const
  {
    return GridComponents::compute(this->gridSize, this->openness, level);
  }
};

// The nearest image of a fractional separation, in the units of the cell side.
double wrappedDistance(double3 a, double3 b, double side)
{
  double3 d = a - b;
  d.x -= std::round(d.x);
  d.y -= std::round(d.y);
  d.z -= std::round(d.z);
  return (d * side).length();
}

}  // namespace

TEST(grid_blocking_cover, one_cavity_gets_one_sphere_over_its_middle)
{
  // A cube of open points off on its own, with nothing reachable anywhere. The sphere goes at the middle of
  // it and reaches the corners, plus half a voxel diagonal for the corner point standing for its voxel.
  const double side = 20.0;
  const uint3 gridSize(20, 20, 20);
  Region region(gridSize, side);
  region.paint(int3(4, 4, 4), int3(8, 8, 8), 1.0f);

  GridComponents components = region.componentsAbove(0.5);
  ASSERT_EQ(components.pores.size(), 1u);
  ASSERT_EQ(components.numberOfPockets, 1u);
  ASSERT_EQ(components.pores[0].numberOfVoxels, 125u);

  std::vector<std::uint8_t> needsCover(1, 1);
  std::vector<std::uint8_t> reachable(region.openness.size(), 0);

  GridBlockingCover cover = coverPockets(gridSize, region.box, components.voxelPore, components.pores, needsCover,
                                         region.openness, reachable);

  ASSERT_EQ(cover.spheres.size(), 1u);
  EXPECT_EQ(cover.numberOfPocketVoxels, 125u);
  EXPECT_EQ(cover.numberOfClippedSpheres, 0u) << "there is nothing reachable to stop it";
  EXPECT_EQ(cover.numberOfRefusedPoints, 0u);

  // The centre is the middle of the cube, at grid point (6, 6, 6) of twenty.
  const double3 middle(6.0 / 20.0, 6.0 / 20.0, 6.0 / 20.0);
  EXPECT_LT(wrappedDistance(cover.spheres[0].centreFractional, middle, side), 1.0e-9);

  // The furthest point is a corner, two grid steps away in each direction, and the spacing is one Angstrom.
  const double toCorner = std::sqrt(3.0) * 2.0;
  const double halfDiagonal = 0.5 * std::sqrt(3.0);
  EXPECT_NEAR(cover.spheres[0].radius, toCorner + halfDiagonal, 1.0e-9);
  EXPECT_EQ(cover.spheres[0].voxelsCovered, 125u);
}

TEST(grid_blocking_cover, a_sphere_stops_at_the_reachable_room_beside_it)
{
  // The same cavity, with a rod of reachable points running past it four grid steps from its middle. The
  // sphere would have grown to about 3.5 and is cut back to 4 minus nothing, the rod being exactly four away.
  const double side = 20.0;
  const uint3 gridSize(20, 20, 20);
  Region region(gridSize, side);
  region.paint(int3(4, 4, 4), int3(8, 8, 8), 1.0f);

  GridComponents components = region.componentsAbove(0.5);
  ASSERT_EQ(components.pores.size(), 1u);

  std::vector<std::uint8_t> needsCover(1, 1);
  std::vector<std::uint8_t> reachable(region.openness.size(), 0);
  for (std::size_t i = 0; i < gridSize.x; ++i)
  {
    // Along x at y = 10, z = 6: four steps in y from the middle of the cavity at (6, 6, 6).
    reachable[(6 * gridSize.y + 10) * gridSize.x + i] = 1;
  }

  GridBlockingCover cover = coverPockets(gridSize, region.box, components.voxelPore, components.pores, needsCover,
                                         region.openness, reachable);

  ASSERT_FALSE(cover.spheres.empty());
  EXPECT_TRUE(cover.spheres[0].clipped) << "the rod is nearer than the corner of the cavity";
  EXPECT_NEAR(cover.spheres[0].radius, 4.0, 1.0e-9) << "cut back to exactly the four steps to the rod";
  EXPECT_GE(cover.numberOfClippedSpheres, 1u);

  // However many spheres it takes, the two things that must hold of them are that the cavity is accounted
  // for and that none of them reaches the rod.
  std::size_t covered = 0;
  for (const GridBlockingSphere &sphere : cover.spheres) covered += sphere.voxelsCovered;
  EXPECT_EQ(covered + cover.numberOfRefusedPoints, 125u) << "every point of the cavity is covered or refused";

  for (const GridBlockingSphere &sphere : cover.spheres)
  {
    const double3 centre = sphere.centreFractional;
    for (std::size_t i = 0; i < gridSize.x; ++i)
    {
      const double3 rod(static_cast<double>(i) / 20.0, 10.0 / 20.0, 6.0 / 20.0);
      EXPECT_GE(wrappedDistance(centre, rod, side), sphere.radius - 1.0e-9)
          << "a sphere reaching the rod would block room the molecule is owed";
    }
  }
}

TEST(grid_blocking_cover, a_cut_that_leaves_part_of_a_cavity_bare_is_followed_up)
{
  // The rod moved in to three steps, which is nearer than the corners of the cavity. One sphere can no
  // longer hold the whole of it, and the covering has to come back for what is left.
  const double side = 20.0;
  const uint3 gridSize(20, 20, 20);
  Region region(gridSize, side);
  region.paint(int3(4, 4, 4), int3(8, 8, 8), 1.0f);

  GridComponents components = region.componentsAbove(0.5);
  std::vector<std::uint8_t> needsCover(1, 1);
  std::vector<std::uint8_t> reachable(region.openness.size(), 0);
  for (std::size_t i = 0; i < gridSize.x; ++i)
  {
    reachable[(6 * gridSize.y + 9) * gridSize.x + i] = 1;
  }

  GridBlockingCover cover = coverPockets(gridSize, region.box, components.voxelPore, components.pores, needsCover,
                                         region.openness, reachable);

  std::size_t covered = 0;
  for (const GridBlockingSphere &sphere : cover.spheres) covered += sphere.voxelsCovered;
  EXPECT_EQ(covered + cover.numberOfRefusedPoints, 125u) << "the cavity is finished, however many it took";

  for (const GridBlockingSphere &sphere : cover.spheres)
  {
    for (std::size_t i = 0; i < gridSize.x; ++i)
    {
      const double3 rod(static_cast<double>(i) / 20.0, 9.0 / 20.0, 6.0 / 20.0);
      EXPECT_GE(wrappedDistance(sphere.centreFractional, rod, side), sphere.radius - 1.0e-9)
          << "no sphere may reach the rod, whatever it takes to cover the cavity";
    }
  }
}

TEST(grid_blocking_cover, only_the_regions_asked_for_are_covered)
{
  // Two cavities, of which one is marked as needing cover and the other is reachable. Nothing is drawn for
  // the second, and the first is kept away from it.
  const double side = 20.0;
  const uint3 gridSize(20, 20, 20);
  Region region(gridSize, side);
  region.paint(int3(3, 3, 3), int3(5, 5, 5), 1.0f);
  region.paint(int3(13, 13, 13), int3(15, 15, 15), 1.0f);

  GridComponents components = region.componentsAbove(0.5);
  ASSERT_EQ(components.pores.size(), 2u);

  // Whichever of the two holds the point at (4, 4, 4) is the one to cover.
  const std::int32_t first = components.voxelPore[(4 * gridSize.y + 4) * gridSize.x + 4];
  ASSERT_GE(first, 0);

  std::vector<std::uint8_t> needsCover(components.pores.size(), 0);
  needsCover[static_cast<std::size_t>(first)] = 1;

  std::vector<std::uint8_t> reachable(region.openness.size(), 0);
  for (std::size_t v = 0; v < region.openness.size(); ++v)
  {
    const std::int32_t pore = components.voxelPore[v];
    if (pore >= 0 && pore != first) reachable[v] = 1;
  }

  GridBlockingCover cover = coverPockets(gridSize, region.box, components.voxelPore, components.pores, needsCover,
                                         region.openness, reachable);

  EXPECT_EQ(cover.numberOfPocketVoxels, 27u) << "only the region asked for was walked";
  for (const GridBlockingSphere &sphere : cover.spheres)
  {
    EXPECT_EQ(sphere.pocket, static_cast<std::size_t>(first));
  }
}

TEST(grid_blocking_cover, a_cavity_across_the_cell_boundary_is_one_body)
{
  // A cube straddling the face at x = 0. Followed out of the cell it is one cavity with a centre inside it,
  // and not two halves with a centre halfway across the cell from either.
  const double side = 20.0;
  const uint3 gridSize(20, 20, 20);
  Region region(gridSize, side);
  region.paint(int3(18, 8, 8), int3(19, 10, 10), 1.0f);
  region.paint(int3(0, 8, 8), int3(1, 10, 10), 1.0f);

  GridComponents components = region.componentsAbove(0.5);
  ASSERT_EQ(components.pores.size(), 1u) << "the two halves are joined across the boundary";
  ASSERT_EQ(components.pores[0].numberOfVoxels, 36u);

  std::vector<std::uint8_t> needsCover(1, 1);
  std::vector<std::uint8_t> reachable(region.openness.size(), 0);

  GridBlockingCover cover = coverPockets(gridSize, region.box, components.voxelPore, components.pores, needsCover,
                                         region.openness, reachable);

  ASSERT_EQ(cover.spheres.size(), 1u);
  EXPECT_EQ(cover.spheres[0].voxelsCovered, 36u);

  // The cavity runs from x = 18 to x = 1 the short way, so its middle is at x = 19.5 of twenty, and a sphere
  // there reaches the whole of it. A centre placed at the mean of the wrapped coordinates would have landed
  // near the far side of the cell and needed a radius of about ten.
  EXPECT_LT(cover.spheres[0].radius, 4.0) << "a sphere over the true middle is small";
  const double x = cover.spheres[0].centreFractional.x;
  EXPECT_TRUE(x > 0.9 || x < 0.1) << "the centre lies in the cavity, at x = " << x;
}
