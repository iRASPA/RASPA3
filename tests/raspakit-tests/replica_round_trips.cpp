#include <gtest/gtest.h>

import std;

import archive;
import replica_round_trips;

// The round-trip tracker follows configuration identities through accepted replica exchanges and
// counts a round trip when a configuration returns to the coldest replica after having visited the
// hottest one. These tests drive it with hand-made swap sequences on a four-replica ladder.

namespace
{
// moves the walker sitting at 'from' to 'to' through successive neighbour swaps, one sweep per step
void walk(ReplicaRoundTrips& tracker, std::size_t from, std::size_t to)
{
  while (from != to)
  {
    const std::size_t next = (to > from) ? from + 1 : from - 1;
    tracker.recordAcceptedSwap(from, next);
    tracker.endOfSweep();
    from = next;
  }
}
}  // namespace

TEST(REPLICA_ROUND_TRIPS, no_swaps_labels_ends_but_completes_nothing)
{
  ReplicaRoundTrips tracker;
  tracker.initialize(4);
  for (std::size_t sweep = 0; sweep < 10; ++sweep)
  {
    tracker.endOfSweep();
  }

  EXPECT_EQ(tracker.sweeps, 10uz);
  EXPECT_EQ(tracker.roundTrips, 0uz);
  EXPECT_EQ(tracker.meanRoundTripSweeps(), 0.0);

  // walker 0 sits at the bottom (Up), walker 3 at the top (Down); the middle ones are never labeled
  EXPECT_EQ(tracker.upFraction(0).value(), 1.0);
  EXPECT_EQ(tracker.upFraction(3).value(), 0.0);
  EXPECT_FALSE(tracker.upFraction(1).has_value());
  EXPECT_FALSE(tracker.upFraction(2).has_value());
}

TEST(REPLICA_ROUND_TRIPS, walker_traversing_the_ladder_and_back_counts_one_round_trip)
{
  ReplicaRoundTrips tracker;
  tracker.initialize(4);
  tracker.endOfSweep();  // sweep 1: labels walker 0 Up at the bottom, walker 3 Down at the top

  // carry walker 0 from replica 0 to replica 3 (3 sweeps) and back down (3 sweeps)
  walk(tracker, 0, 3);
  EXPECT_EQ(tracker.roundTrips, 0uz);  // reached the top, not yet back
  EXPECT_EQ(tracker.walkerAtReplica[3], 0uz);

  walk(tracker, 3, 0);
  EXPECT_EQ(tracker.walkerAtReplica[0], 0uz);
  EXPECT_EQ(tracker.roundTrips, 1uz);
  EXPECT_EQ(tracker.roundTripsPerWalker[0], 1uz);
  EXPECT_EQ(tracker.roundTripsPerWalker[1], 0uz);
  EXPECT_EQ(tracker.roundTripsPerWalker[2], 0uz);
  EXPECT_EQ(tracker.roundTripsPerWalker[3], 0uz);

  // left the bottom after sweep 1, back at the bottom at sweep 7
  EXPECT_DOUBLE_EQ(tracker.meanRoundTripSweeps(), 6.0);

  // once back at the bottom the walker is Up again; another excursion is a second round trip
  walk(tracker, 0, 3);
  walk(tracker, 3, 0);
  EXPECT_EQ(tracker.roundTrips, 2uz);
  EXPECT_DOUBLE_EQ(tracker.meanRoundTripSweeps(), 6.0);
}

TEST(REPLICA_ROUND_TRIPS, excursion_without_reaching_the_top_is_not_a_round_trip)
{
  ReplicaRoundTrips tracker;
  tracker.initialize(4);
  tracker.endOfSweep();

  walk(tracker, 0, 2);  // up to replica 2 only
  walk(tracker, 2, 0);  // and back
  EXPECT_EQ(tracker.roundTrips, 0uz);

  // the walker that started at the top is Down; bringing it to the bottom closes a round trip for it,
  // since it has visited the hottest replica (where it started)
  EXPECT_EQ(tracker.walkerAtReplica[3], 3uz);
  walk(tracker, 3, 0);
  EXPECT_EQ(tracker.roundTrips, 1uz);
  EXPECT_EQ(tracker.roundTripsPerWalker[3], 1uz);
}

TEST(REPLICA_ROUND_TRIPS, up_fraction_counts_labeled_walkers_per_temperature)
{
  ReplicaRoundTrips tracker;
  tracker.initialize(3);
  tracker.endOfSweep();  // walker 0 Up at replica 0, walker 2 Down at replica 2; walker 1 unlabeled

  // swap replicas 0 and 1: the Up walker moves to replica 1, the unlabeled walker arrives at the
  // bottom and becomes Up
  tracker.recordAcceptedSwap(0, 1);
  tracker.endOfSweep();

  EXPECT_EQ(tracker.upVisitsPerReplica[0], 2uz);
  EXPECT_EQ(tracker.downVisitsPerReplica[0], 0uz);
  EXPECT_EQ(tracker.upVisitsPerReplica[1], 1uz);
  EXPECT_EQ(tracker.downVisitsPerReplica[1], 0uz);
  EXPECT_EQ(tracker.upVisitsPerReplica[2], 0uz);
  EXPECT_EQ(tracker.downVisitsPerReplica[2], 2uz);
  EXPECT_EQ(tracker.upFraction(1).value(), 1.0);

  // now bring the Down walker from the top to the middle: replica 1 has seen one Up and one Down
  tracker.recordAcceptedSwap(1, 2);
  tracker.endOfSweep();
  EXPECT_DOUBLE_EQ(tracker.upFraction(1).value(), 0.5);
}

TEST(REPLICA_ROUND_TRIPS, archive_round_trip_restores_state)
{
  ReplicaRoundTrips original;
  original.initialize(5);
  original.endOfSweep();
  walk(original, 0, 4);
  walk(original, 4, 0);
  walk(original, 0, 2);
  ASSERT_EQ(original.roundTrips, 1uz);

  const std::filesystem::path path = std::filesystem::temp_directory_path() / "raspa3_replica_round_trips.bin";
  {
    std::ofstream stream(path, std::ios::binary);
    Archive<std::ofstream> archive(stream);
    archive << original;
  }
  ReplicaRoundTrips restored;
  {
    std::ifstream stream(path, std::ios::binary);
    Archive<std::ifstream> archive(stream);
    archive >> restored;
  }
  std::filesystem::remove(path);

  EXPECT_EQ(restored.numberOfReplicas, original.numberOfReplicas);
  EXPECT_EQ(restored.sweeps, original.sweeps);
  EXPECT_EQ(restored.walkerAtReplica, original.walkerAtReplica);
  EXPECT_EQ(restored.direction, original.direction);
  EXPECT_EQ(restored.leftBottomSweep, original.leftBottomSweep);
  EXPECT_EQ(restored.roundTrips, original.roundTrips);
  EXPECT_EQ(restored.roundTripSweepsTotal, original.roundTripSweepsTotal);
  EXPECT_EQ(restored.roundTripsPerWalker, original.roundTripsPerWalker);
  EXPECT_EQ(restored.upVisitsPerReplica, original.upVisitsPerReplica);
  EXPECT_EQ(restored.downVisitsPerReplica, original.downVisitsPerReplica);

  // continuing on the restored tracker gives the same result as continuing on the original
  walk(original, 2, 4);
  walk(original, 4, 0);
  walk(restored, 2, 4);
  walk(restored, 4, 0);
  EXPECT_EQ(restored.roundTrips, 2uz);
  EXPECT_EQ(restored.roundTrips, original.roundTrips);
  EXPECT_DOUBLE_EQ(restored.meanRoundTripSweeps(), original.meanRoundTripSweeps());
}
