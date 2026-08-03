#include <gtest/gtest.h>

import std;

import uint3;
import grid_percolation;
import grid_connected_components;

// The escape level: for each point, how far the field has to be allowed to fall before that point is joined
// to a region running through the crystal.
//
// The cases below are built on a coarse grid whose values are written out by hand, so that the answer is a
// matter of reading the picture rather than of trusting the sweep that produced it.

namespace
{

// A field of `background` everywhere, into which named points and boxes are painted.
struct Field
{
  uint3 gridSize;
  std::vector<float> value;

  Field(uint3 size, float background)
      : gridSize(size), value(static_cast<std::size_t>(size.x) * size.y * size.z, background)
  {
  }

  float &at(std::size_t i, std::size_t j, std::size_t k)
  {
    return this->value[(k * this->gridSize.y + j) * this->gridSize.x + i];
  }

  float read(std::size_t i, std::size_t j, std::size_t k) const
  {
    return this->value[(k * this->gridSize.y + j) * this->gridSize.x + i];
  }
};

float escapeAt(const GridPercolation &sweep, uint3 gridSize, std::size_t i, std::size_t j, std::size_t k)
{
  return sweep.escapeOpenness[(k * gridSize.y + j) * gridSize.x + i];
}

}  // namespace

TEST(grid_percolation, a_channel_of_one_width_costs_its_points_nothing)
{
  // One rod of open points running the length of x, all of a width, everything else shut. The rod is a
  // channel at its own level, so a point of it is on the percolating region the moment it is switched on and
  // gives up nothing to leave.
  const uint3 gridSize(8, 8, 8);
  Field field(gridSize, 0.0f);
  for (std::size_t i = 0; i < gridSize.x; ++i) field.at(i, 4, 4) = 5.0f;

  GridPercolation sweep = sweepPercolation(gridSize, field.value, 1.0f, true);

  ASSERT_TRUE(sweep.percolates);
  for (std::size_t i = 0; i < gridSize.x; ++i)
  {
    EXPECT_FLOAT_EQ(escapeAt(sweep, gridSize, i, 4, 4), 5.0f) << "at " << i;
  }
}

TEST(grid_percolation, a_channel_is_worth_no_more_than_its_tightest_point)
{
  // The same rod, pinched. Nothing in it can get anywhere until the pinch is switched on, so the widest part
  // of the channel is no better off than the narrowest, and every point of it gives up the same.
  const uint3 gridSize(8, 8, 8);
  Field field(gridSize, 0.0f);
  for (std::size_t i = 0; i < gridSize.x; ++i) field.at(i, 4, 4) = 5.0f;
  field.at(6, 4, 4) = 2.0f;

  GridPercolation sweep = sweepPercolation(gridSize, field.value, 1.0f, true);

  ASSERT_TRUE(sweep.percolates);
  EXPECT_FLOAT_EQ(sweep.percolationOpenness, 2.0) << "the rod does not run anywhere until the pinch is on";
  for (std::size_t i = 0; i < gridSize.x; ++i)
  {
    EXPECT_FLOAT_EQ(escapeAt(sweep, gridSize, i, 4, 4), 2.0f) << "at " << i;
  }
}

TEST(grid_percolation, a_pocket_gives_up_the_neck_between_it_and_the_channel)
{
  // A channel along x at y = 4, and a cavity beside it at y = 6 joined to it by a single point of 2. The
  // cavity is wide open at 9, so a point in it has to come down to 2 to get out, and the one point of the
  // neck has to come down to its own value, being no better off than the neck itself.
  const uint3 gridSize(8, 8, 8);
  Field field(gridSize, 0.0f);
  for (std::size_t i = 0; i < gridSize.x; ++i) field.at(i, 4, 4) = 5.0f;
  field.at(3, 5, 4) = 2.0f;
  field.at(3, 6, 4) = 9.0f;
  field.at(4, 6, 4) = 8.0f;
  field.at(2, 6, 4) = 7.0f;

  GridPercolation sweep = sweepPercolation(gridSize, field.value, 1.0f, true);

  ASSERT_TRUE(sweep.percolates);
  EXPECT_FLOAT_EQ(sweep.percolationOpenness, 5.0) << "the channel closes on itself at its own level";

  for (std::size_t i = 0; i < gridSize.x; ++i)
  {
    EXPECT_FLOAT_EQ(escapeAt(sweep, gridSize, i, 4, 4), 5.0f) << "the channel is its own way out";
  }
  EXPECT_FLOAT_EQ(escapeAt(sweep, gridSize, 3, 6, 4), 2.0f) << "the cavity leaves by the neck";
  EXPECT_FLOAT_EQ(escapeAt(sweep, gridSize, 4, 6, 4), 2.0f);
  EXPECT_FLOAT_EQ(escapeAt(sweep, gridSize, 2, 6, 4), 2.0f);
  EXPECT_FLOAT_EQ(escapeAt(sweep, gridSize, 3, 5, 4), 2.0f) << "the neck is no better off than itself";
}

TEST(grid_percolation, the_deeper_of_two_necks_is_the_one_that_is_used)
{
  // The same cavity, joined to the channel twice: once by a neck of 2 and once by a neck of 4. The wider of
  // the two is the way out, so the cavity gives up 4 rather than 2, while the tighter neck itself still has
  // to come down to its own value before it is joined to anything.
  const uint3 gridSize(8, 8, 8);
  Field field(gridSize, 0.0f);
  for (std::size_t i = 0; i < gridSize.x; ++i) field.at(i, 4, 4) = 6.0f;
  field.at(2, 5, 4) = 2.0f;
  field.at(5, 5, 4) = 4.0f;
  for (std::size_t i = 2; i <= 5; ++i) field.at(i, 6, 4) = 9.0f;

  GridPercolation sweep = sweepPercolation(gridSize, field.value, 1.0f, true);

  ASSERT_TRUE(sweep.percolates);
  for (std::size_t i = 2; i <= 5; ++i)
  {
    EXPECT_FLOAT_EQ(escapeAt(sweep, gridSize, i, 6, 4), 4.0f) << "at " << i << ": the wider neck is the way out";
  }
  EXPECT_FLOAT_EQ(escapeAt(sweep, gridSize, 5, 5, 4), 4.0f) << "the wider neck leaves at its own level";
  EXPECT_FLOAT_EQ(escapeAt(sweep, gridSize, 2, 5, 4), 2.0f)
      << "the tighter neck is joined to nothing until it is switched on itself";
}

TEST(grid_percolation, a_sealed_cavity_never_gets_out)
{
  // A channel and a cavity with nothing between them. Nothing the sweep does joins the two, so the cavity
  // keeps the value that stands for having no way out at any level.
  const uint3 gridSize(8, 8, 8);
  Field field(gridSize, 0.0f);
  for (std::size_t i = 0; i < gridSize.x; ++i) field.at(i, 4, 4) = 5.0f;
  field.at(3, 6, 4) = 9.0f;
  field.at(4, 6, 4) = 8.0f;

  GridPercolation sweep = sweepPercolation(gridSize, field.value, 1.0f, true);

  ASSERT_TRUE(sweep.percolates);
  const float noWayOut = std::numeric_limits<float>::lowest();
  EXPECT_FLOAT_EQ(escapeAt(sweep, gridSize, 3, 6, 4), noWayOut);
  EXPECT_FLOAT_EQ(escapeAt(sweep, gridSize, 4, 6, 4), noWayOut);
  EXPECT_FLOAT_EQ(escapeAt(sweep, gridSize, 0, 4, 4), 5.0f) << "the channel is unaffected by it";
}

TEST(grid_percolation, recording_the_escape_changes_nothing_else)
{
  // Whatever the bookkeeping costs, it must not move any of the numbers the sweep already reported.
  const uint3 gridSize(12, 12, 12);
  Field field(gridSize, 0.0f);

  std::mt19937 engine(20260801);
  std::uniform_real_distribution<float> uniform(0.0f, 10.0f);
  for (float &value : field.value) value = uniform(engine);

  GridPercolation plain = sweepPercolation(gridSize, field.value, 3.0f, false);
  GridPercolation recorded = sweepPercolation(gridSize, field.value, 3.0f, true);

  EXPECT_EQ(plain.percolates, recorded.percolates);
  EXPECT_EQ(plain.numberOfVoxels, recorded.numberOfVoxels);
  EXPECT_EQ(plain.dimensionalityAtThreshold, recorded.dimensionalityAtThreshold);
  EXPECT_DOUBLE_EQ(plain.percolationOpenness, recorded.percolationOpenness);
  EXPECT_DOUBLE_EQ(plain.highestOpenness, recorded.highestOpenness);
  EXPECT_DOUBLE_EQ(plain.highestOpennessOnPath, recorded.highestOpennessOnPath);
  for (std::size_t d = 0; d < 3; ++d)
  {
    EXPECT_DOUBLE_EQ(plain.opennessByDimension[d], recorded.opennessByDimension[d]);
  }
  EXPECT_TRUE(plain.escapeOpenness.empty());
}

TEST(grid_percolation, the_escape_agrees_with_asking_at_every_level_separately)
{
  // The sweep answers every point in one pass, which is the whole point of it and also the thing that could
  // go wrong. The definition it is meant to be computing can be applied directly instead: at each level in
  // turn, take the components of the field held at that level and see whether the one holding the point runs
  // through the crystal. The highest level at which it does is the answer.
  //
  // That is a different piece of code, run once per level rather than once in total, so the two agreeing is
  // worth something. The field is coarse and takes few values to keep the direct version affordable.
  const uint3 gridSize(10, 10, 10);
  Field field(gridSize, 0.0f);

  // Two points in five take part, which is above the level at which a lattice like this starts running from
  // side to side and well below the level at which everything on it is joined to everything else. So the
  // field both percolates and has pieces sealed off from what percolates, which is the case worth checking.
  std::mt19937 engine(20260801);
  std::uniform_int_distribution<int> uniform(0, 9);
  for (float &value : field.value) value = static_cast<float>(uniform(engine));

  const float floorLevel = 6.0f;
  GridPercolation sweep = sweepPercolation(gridSize, field.value, floorLevel, true);
  ASSERT_TRUE(sweep.percolates);

  const float noWayOut = std::numeric_limits<float>::lowest();
  std::vector<float> expected(field.value.size(), noWayOut);
  for (int level = 9; level >= static_cast<int>(floorLevel); --level)
  {
    GridComponents components = GridComponents::compute(gridSize, field.value, static_cast<double>(level));
    for (std::size_t v = 0; v < field.value.size(); ++v)
    {
      if (expected[v] != noWayOut) continue;
      const std::int32_t pore = components.voxelPore[v];
      if (pore < 0) continue;
      if (!components.pores[static_cast<std::size_t>(pore)].isChannel) continue;
      expected[v] = static_cast<float>(level);
    }
  }

  std::size_t trapped = 0;
  for (std::size_t v = 0; v < field.value.size(); ++v)
  {
    if (field.value[v] < floorLevel) continue;
    ASSERT_FLOAT_EQ(sweep.escapeOpenness[v], expected[v]) << "at voxel " << v;
    if (expected[v] == noWayOut) ++trapped;
  }

  // The check is only worth having if the field had somewhere to be trapped in and somewhere to get out by.
  EXPECT_GT(trapped, 0u) << "nothing was sealed off, so the case being checked never arose";
}
