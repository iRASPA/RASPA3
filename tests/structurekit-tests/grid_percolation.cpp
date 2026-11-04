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


// Three rods of one openness along x, y and z, crossing at the centre, with one point of the y rod more open
// than the rest. Whichever rod the sweep happens to close first, all three close at the same level, so the
// pore system runs in three directions at its threshold; and the open point is reachable from any of them
// without passing anything less open than that level, so it is the most open point on the path. Both were
// read off at the first closing, when the other rods and the rest of the closing point's own neighbours were
// still to come, and came out as the order of ties fell.
TEST(grid_percolation, channels_that_open_together_are_counted_together_at_the_threshold)
{
  const uint3 gridSize(8, 8, 8);
  Field field(gridSize, 0.0f);
  for (std::size_t i = 0; i < 8; ++i)
  {
    field.at(i, 4, 4) = 5.0f;
    field.at(4, i, 4) = 5.0f;
    field.at(4, 4, i) = 5.0f;
  }
  field.at(4, 0, 4) = 9.0f;

  GridPercolation sweep = sweepPercolation(gridSize, field.value, 1.0f, true);

  ASSERT_TRUE(sweep.percolates);
  EXPECT_NEAR(sweep.percolationOpenness, 5.0, 1.0e-12);
  EXPECT_EQ(sweep.dimensionalityAtThreshold, 3);
  EXPECT_NEAR(sweep.highestOpennessOnPath, 9.0, 1.0e-12);
  for (int dimension = 0; dimension < 3; ++dimension)
  {
    EXPECT_NEAR(sweep.opennessByDimension[static_cast<std::size_t>(dimension)], 5.0, 1.0e-12);
  }
}


// A rod that closes at a higher level than the rest keeps the threshold to itself: the other rods join below
// it, and the dimensionality at the threshold stays one, as it should.
TEST(grid_percolation, a_channel_that_opens_alone_keeps_its_own_threshold)
{
  const uint3 gridSize(8, 8, 8);
  Field field(gridSize, 0.0f);
  for (std::size_t i = 0; i < 8; ++i)
  {
    field.at(i, 4, 4) = 6.0f;
    field.at(4, i, 4) = 5.0f;
    field.at(4, 4, i) = 5.0f;
  }
  field.at(4, 4, 4) = 6.0f;  // the crossing belongs to the x rod, which the other two rods overwrote
  field.at(4, 0, 4) = 9.0f;

  GridPercolation sweep = sweepPercolation(gridSize, field.value, 1.0f, true);

  ASSERT_TRUE(sweep.percolates);
  EXPECT_NEAR(sweep.percolationOpenness, 6.0, 1.0e-12);
  EXPECT_EQ(sweep.dimensionalityAtThreshold, 1);
  // The y rod holds the open point, but joins only at 5: below the threshold, so off the path.
  EXPECT_NEAR(sweep.highestOpennessOnPath, 6.0, 1.0e-12);
  EXPECT_NEAR(sweep.opennessByDimension[0], 6.0, 1.0e-12);
  EXPECT_NEAR(sweep.opennessByDimension[2], 5.0, 1.0e-12);
}


// Two rods that never meet, of one width, and a wide place in one of them. Whichever rod holds it, it is reachable
// without passing anything narrower than the threshold, so it is on the path; the sweep reaching the other rod's
// closing point first does not change that.
TEST(grid_percolation, a_wide_place_in_either_rod_tied_at_the_threshold_is_on_the_path)
{
  const uint3 gridSize(8, 8, 8);
  for (std::size_t wider : {2uz, 6uz})
  {
    Field field(gridSize, 0.0f);
    for (std::size_t i = 0; i < gridSize.x; ++i)
    {
      field.at(i, 2, 2) = 5.0f;
      field.at(i, 6, 6) = 5.0f;
    }
    field.at(3, wider, wider) = 9.0f;

    GridPercolation sweep = sweepPercolation(gridSize, field.value, 1.0f, false);
    ASSERT_TRUE(sweep.percolates);
    EXPECT_FLOAT_EQ(sweep.percolationOpenness, 5.0f);
    EXPECT_FLOAT_EQ(sweep.highestOpennessOnPath, 9.0f) << "wide place in the rod at " << wider;
    EXPECT_EQ(sweep.dimensionalityAtThreshold, 1);
  }
}

// A rod and a sheet tied at one level: the pore system runs two ways at it, wherever the two lie in the grid. Two
// rods along different directions that never touch still run one way each.
TEST(grid_percolation, the_dimensionality_at_a_tied_threshold_is_that_of_the_region_that_runs_most_ways)
{
  const uint3 gridSize(8, 8, 8);
  for (std::size_t rodAt : {1uz, 6uz})
  {
    const std::size_t sheetAt = rodAt == 1 ? 5 : 2;
    Field field(gridSize, 0.0f);
    for (std::size_t i = 0; i < gridSize.x; ++i) field.at(i, rodAt, rodAt) = 5.0f;
    for (std::size_t i = 0; i < gridSize.x; ++i)
      for (std::size_t j = 0; j < gridSize.y; ++j) field.at(i, j, sheetAt) = 5.0f;

    GridPercolation sweep = sweepPercolation(gridSize, field.value, 1.0f, false);
    ASSERT_TRUE(sweep.percolates);
    EXPECT_EQ(sweep.dimensionalityAtThreshold, 2) << "rod at " << rodAt;
  }

  Field crossed(gridSize, 0.0f);
  for (std::size_t i = 0; i < gridSize.x; ++i) crossed.at(i, 2, 2) = 5.0f;
  for (std::size_t j = 0; j < gridSize.y; ++j) crossed.at(6, j, 6) = 5.0f;
  GridPercolation sweep = sweepPercolation(gridSize, crossed.value, 1.0f, false);
  ASSERT_TRUE(sweep.percolates);
  EXPECT_EQ(sweep.dimensionalityAtThreshold, 1);
}


// Openness is whatever the caller sweeps by, and the energy analyses sweep by the negated energy, so all of it can be
// below zero. The most open point on the path is then a negative number and has to come back as one.
TEST(grid_percolation, the_most_open_point_on_the_path_can_be_negative)
{
  const uint3 gridSize(8, 8, 8);
  Field field(gridSize, -1000.0f);
  for (std::size_t i = 0; i < gridSize.x; ++i) field.at(i, 4, 4) = -5.0f;
  field.at(3, 4, 4) = -2.0f;

  GridPercolation sweep = sweepPercolation(gridSize, field.value, -100.0f, false);
  ASSERT_TRUE(sweep.percolates);
  EXPECT_FLOAT_EQ(sweep.percolationOpenness, -5.0f);
  EXPECT_FLOAT_EQ(sweep.highestOpennessOnPath, -2.0f);
}
