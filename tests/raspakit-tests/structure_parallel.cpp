#include <gtest/gtest.h>

import std;

import structure_parallel;

// The two loops the structural analyses hand their independent work to: `forEachIndex`, which deals the indices
// out interleaved, and `forEachBlock`, which gives each worker a contiguous stretch.
//
// The threaded route cannot be run from here: the thread pool is a singleton that cannot be put back the way it
// was found, so starting one would leave every test after this one running on it. What can be run from here,
// and is, is the arithmetic those routes are built out of. `laneCount`, `forEachIndexOfLane` and `blockOfLane`
// are the whole of the dealing-out and are functions of their arguments alone, so the promises that matter ---
// every index taken, none of them twice, and the same lane taking the same ones every time --- are checked
// below at every shape of loop the analyses can present, with no threads involved. What is left untested is the
// dispatch itself, which is a handful of lines around a pool that the analyses exercise end to end.
//
// The serial route through either loop is the one every run takes unless threads were asked for, and the
// promise it makes is the stronger one: every index is visited, exactly once, in order, by the calling thread.

namespace
{

// The indices one lane takes, in the order it takes them.
std::vector<std::size_t> lane(std::size_t count, std::size_t lanes, std::size_t worker)
{
  std::vector<std::size_t> taken;
  forEachIndexOfLane(count, lanes, worker, [&](std::size_t index) { taken.push_back(index); });
  return taken;
}

}  // namespace

TEST(structure_parallel, without_a_pool_there_is_one_worker) { EXPECT_EQ(workersAvailable(), 1uz); }

TEST(structure_parallel, a_request_for_lanes_is_held_to_what_can_run_them)
{
  EXPECT_EQ(laneCount(8, 8), 8uz);
  EXPECT_EQ(laneCount(3, 8), 3uz);  // fewer asked for than there are
  EXPECT_EQ(laneCount(8, 3), 3uz);  // more asked for than there are

  // Neither number can bring it below one: a loop always has the calling thread to run on.
  EXPECT_EQ(laneCount(0, 8), 1uz);
  EXPECT_EQ(laneCount(8, 0), 1uz);
  EXPECT_EQ(laneCount(0, 0), 1uz);
  EXPECT_EQ(laneCount(1, 1), 1uz);
}

TEST(structure_parallel, every_index_is_taken_by_exactly_one_lane)
{
  for (std::size_t count = 0; count <= 24; ++count)
  {
    for (std::size_t lanes = 1; lanes <= 8; ++lanes)
    {
      std::vector<std::size_t> visits(count, 0);
      std::size_t total = 0;
      for (std::size_t worker = 0; worker < lanes; ++worker)
      {
        for (std::size_t index : lane(count, lanes, worker))
        {
          ASSERT_LT(index, count) << "count " << count << " lanes " << lanes << " worker " << worker;
          ++visits[index];
          ++total;
        }
      }

      EXPECT_EQ(total, count) << "count " << count << " lanes " << lanes;
      for (std::size_t index = 0; index < count; ++index)
      {
        EXPECT_EQ(visits[index], 1uz) << "index " << index << " of " << count << " on " << lanes << " lanes";
      }
    }
  }
}

TEST(structure_parallel, a_lane_takes_its_indices_in_order_and_the_same_ones_every_time)
{
  for (std::size_t count = 0; count <= 24; ++count)
  {
    for (std::size_t lanes = 1; lanes <= 8; ++lanes)
    {
      for (std::size_t worker = 0; worker < lanes; ++worker)
      {
        const std::vector<std::size_t> taken = lane(count, lanes, worker);
        EXPECT_EQ(taken, lane(count, lanes, worker)) << "count " << count << " lanes " << lanes;
        EXPECT_TRUE(std::ranges::is_sorted(taken)) << "count " << count << " lanes " << lanes;
        if (!taken.empty()) EXPECT_EQ(taken.front(), worker);
      }
    }
  }
}

TEST(structure_parallel, the_lanes_are_within_one_index_of_one_another)
{
  // What interleaving buys over dealing out blocks is that no lane can be left holding a run of the expensive
  // indices while the others have finished, and the first half of that is that they hold the same number.
  for (std::size_t count = 0; count <= 24; ++count)
  {
    for (std::size_t lanes = 1; lanes <= 8; ++lanes)
    {
      std::size_t fewest = count;
      std::size_t most = 0;
      for (std::size_t worker = 0; worker < lanes; ++worker)
      {
        const std::size_t taken = lane(count, lanes, worker).size();
        fewest = std::min(fewest, taken);
        most = std::max(most, taken);
      }
      EXPECT_LE(most - fewest, 1uz) << "count " << count << " lanes " << lanes;
    }
  }
}

TEST(structure_parallel, more_lanes_than_indices_leaves_the_last_of_them_empty)
{
  EXPECT_EQ(lane(3, 8, 2), (std::vector<std::size_t>{2}));
  EXPECT_TRUE(lane(3, 8, 3).empty());
  EXPECT_TRUE(lane(3, 8, 7).empty());
  EXPECT_TRUE(lane(0, 4, 0).empty());
}

TEST(structure_parallel, a_lane_of_none_is_no_loop_at_all)
{
  // Not a shape `forEachIndex` can produce, `laneCount` never returning zero, but the loop is exported and
  // stepping by nothing would not end.
  EXPECT_TRUE(lane(5, 0, 0).empty());
}

TEST(structure_parallel, one_lane_takes_everything_in_order)
{
  const std::vector<std::size_t> taken = lane(6, 1, 0);
  EXPECT_EQ(taken, (std::vector<std::size_t>{0, 1, 2, 3, 4, 5}));
}

TEST(structure_parallel, the_blocks_cover_every_index_once_and_follow_one_another)
{
  for (std::size_t count = 0; count <= 24; ++count)
  {
    for (std::size_t lanes = 1; lanes <= 8; ++lanes)
    {
      std::size_t expected = 0;
      for (std::size_t worker = 0; worker < lanes; ++worker)
      {
        const auto [begin, end] = blockOfLane(count, lanes, worker);
        ASSERT_LE(begin, end) << "count " << count << " lanes " << lanes << " worker " << worker;
        ASSERT_LE(end, count) << "count " << count << " lanes " << lanes << " worker " << worker;

        // Each stretch starts where the one before it stopped, which is what makes them a partition and what
        // keeps two workers off the same cache line.
        EXPECT_EQ(begin, expected) << "count " << count << " lanes " << lanes << " worker " << worker;
        expected = end;
      }
      EXPECT_EQ(expected, count) << "count " << count << " lanes " << lanes;
    }
  }
}

TEST(structure_parallel, the_blocks_are_within_one_index_of_one_another)
{
  for (std::size_t count = 0; count <= 24; ++count)
  {
    for (std::size_t lanes = 1; lanes <= 8; ++lanes)
    {
      std::size_t fewest = count;
      std::size_t most = 0;
      for (std::size_t worker = 0; worker < lanes; ++worker)
      {
        const auto [begin, end] = blockOfLane(count, lanes, worker);
        fewest = std::min(fewest, end - begin);
        most = std::max(most, end - begin);
      }
      EXPECT_LE(most - fewest, 1uz) << "count " << count << " lanes " << lanes;
    }
  }
}

TEST(structure_parallel, there_is_no_block_beyond_the_last_lane)
{
  EXPECT_EQ(blockOfLane(10, 4, 0), (std::pair<std::size_t, std::size_t>{0, 3}));
  EXPECT_EQ(blockOfLane(10, 4, 1), (std::pair<std::size_t, std::size_t>{3, 6}));
  EXPECT_EQ(blockOfLane(10, 4, 2), (std::pair<std::size_t, std::size_t>{6, 8}));
  EXPECT_EQ(blockOfLane(10, 4, 3), (std::pair<std::size_t, std::size_t>{8, 10}));

  // Asked for a lane there is not, or for a lane of no lanes at all, what comes back is an empty stretch at the
  // end rather than one that would be walked.
  EXPECT_EQ(blockOfLane(10, 4, 4), (std::pair<std::size_t, std::size_t>{10, 10}));
  EXPECT_EQ(blockOfLane(10, 0, 0), (std::pair<std::size_t, std::size_t>{10, 10}));
}

TEST(structure_parallel, a_block_loop_with_no_pool_is_one_stretch_on_the_calling_thread)
{
  constexpr std::size_t count = 11;
  std::vector<std::pair<std::size_t, std::size_t>> stretches;
  forEachBlock(count, 8,
               [&](std::size_t worker, std::size_t begin, std::size_t end)
               {
                 EXPECT_EQ(worker, 0uz);
                 stretches.emplace_back(begin, end);
               });

  EXPECT_EQ(stretches, (std::vector<std::pair<std::size_t, std::size_t>>{{0, count}}));

  // An empty loop is still one call, over nothing.
  stretches.clear();
  forEachBlock(0, 8, [&](std::size_t, std::size_t begin, std::size_t end) { stretches.emplace_back(begin, end); });
  EXPECT_EQ(stretches, (std::vector<std::pair<std::size_t, std::size_t>>{{0, 0}}));
}

TEST(structure_parallel, an_exception_from_a_block_comes_back_out)
{
  EXPECT_THROW(forEachBlock(5, workersAvailable(),
                            [](std::size_t, std::size_t, std::size_t) { throw std::runtime_error("the stretch"); }),
               std::runtime_error);
}

TEST(structure_parallel, every_index_is_visited_exactly_once_and_in_order)
{
  constexpr std::size_t count = 37;
  std::vector<std::size_t> visits(count, 0);
  std::vector<std::size_t> order;

  forEachIndex(count, workersAvailable(),
               [&](std::size_t worker, std::size_t index)
               {
                 EXPECT_EQ(worker, 0uz);
                 ++visits[index];
                 order.push_back(index);
               });

  ASSERT_EQ(order.size(), count);
  for (std::size_t index = 0; index < count; ++index)
  {
    EXPECT_EQ(visits[index], 1uz);
    EXPECT_EQ(order[index], index);
  }
}

TEST(structure_parallel, an_empty_loop_and_a_loop_of_one_are_both_allowed)
{
  std::size_t visited = 0;
  forEachIndex(0, workersAvailable(), [&](std::size_t, std::size_t) { ++visited; });
  EXPECT_EQ(visited, 0uz);

  forEachIndex(1, workersAvailable(),
               [&](std::size_t worker, std::size_t index)
               {
                 EXPECT_EQ(worker, 0uz);
                 EXPECT_EQ(index, 0uz);
                 ++visited;
               });
  EXPECT_EQ(visited, 1uz);
}

TEST(structure_parallel, asking_for_lanes_there_is_no_pool_for_runs_them_serially)
{
  // A caller may ask for any number of workers; what it gets is what there is to run them on. With no pool
  // that is the calling thread, and the loop is the serial one rather than a queue nobody will drain.
  constexpr std::size_t count = 11;
  std::vector<std::size_t> order;
  forEachIndex(count, 8,
               [&](std::size_t worker, std::size_t index)
               {
                 EXPECT_EQ(worker, 0uz);
                 order.push_back(index);
               });

  ASSERT_EQ(order.size(), count);
  for (std::size_t index = 0; index < count; ++index) EXPECT_EQ(order[index], index);
}

TEST(structure_parallel, an_exception_from_one_index_comes_back_out)
{
  EXPECT_THROW(forEachIndex(5, workersAvailable(),
                            [](std::size_t, std::size_t index)
                            {
                              if (index == 3) throw std::runtime_error("index three");
                            }),
               std::runtime_error);
}
