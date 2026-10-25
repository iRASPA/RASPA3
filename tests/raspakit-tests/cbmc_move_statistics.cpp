#include <gtest/gtest.h>

import std;

import archive;
import move_statistics;
import cbmc_statistics;

// The CBMC step-size statistics are part of the binary restart file. A round trip must restore every
// counter and step size.

static void expectEqualStatistics(const MoveStatistics<double>& a, const MoveStatistics<double>& b)
{
  EXPECT_EQ(a.counts, b.counts);
  EXPECT_EQ(a.constructed, b.constructed);
  EXPECT_EQ(a.accepted, b.accepted);
  EXPECT_EQ(a.totalCounts, b.totalCounts);
  EXPECT_EQ(a.totalConstructed, b.totalConstructed);
  EXPECT_EQ(a.totalAccepted, b.totalAccepted);
  EXPECT_EQ(a.maxChange, b.maxChange);
  EXPECT_EQ(a.lowerLimit, b.lowerLimit);
  EXPECT_EQ(a.upperLimit, b.upperLimit);
}

TEST(CBMC_MOVE_STATISTICS, archive_round_trip_restores_all_step_sizes)
{
  CBMC::InternalMoveStatistics original;
  original.ringDisplacementChange.maxChange = 0.37;
  original.ringDisplacementChange.counts = 12.0;
  original.ringDisplacementChange.accepted = 5.0;
  original.ringRotationChange.maxChange = 0.61;
  original.ringRotationChange.totalCounts = 400.0;
  original.rigidTiltRotationChange.maxChange = 0.42;
  original.rigidTiltRotationChange.counts = 300.0;
  original.rigidTiltRotationChange.constructed = 300.0;
  original.rigidTiltRotationChange.accepted = 150.0;
  original.rigidTiltRotationChange.totalAccepted = 1500.0;

  const std::filesystem::path path =
      std::filesystem::temp_directory_path() / "raspa3_cbmc_move_statistics_round_trip.bin";
  {
    std::ofstream stream(path, std::ios::binary);
    Archive<std::ofstream> archive(stream);
    archive << original;
  }
  CBMC::InternalMoveStatistics restored;
  {
    std::ifstream stream(path, std::ios::binary);
    Archive<std::ifstream> archive(stream);
    archive >> restored;
  }
  std::filesystem::remove(path);

  expectEqualStatistics(restored.ringDisplacementChange, original.ringDisplacementChange);
  expectEqualStatistics(restored.ringRotationChange, original.ringRotationChange);
  expectEqualStatistics(restored.ringCrankshaftMove, original.ringCrankshaftMove);
  expectEqualStatistics(restored.rigidTiltRotationChange, original.rigidTiltRotationChange);
}

TEST(CBMC_MOVE_STATISTICS, optimize_adapts_rigid_tilt_like_ring_step_sizes)
{
  CBMC::InternalMoveStatistics statistics;
  const double initial = statistics.rigidTiltRotationChange.maxChange;

  // Acceptance far above the target: the step size must grow (clamped scaling 1.5).
  statistics.rigidTiltRotationChange.counts = 100.0;
  statistics.rigidTiltRotationChange.accepted = 100.0;
  statistics.optimize();
  EXPECT_GT(statistics.rigidTiltRotationChange.maxChange, initial);
  EXPECT_EQ(statistics.rigidTiltRotationChange.counts, 0.0);  // counters reset for the next window

  // Acceptance far below the target: the step size must shrink, but never below the lower limit.
  for (int i = 0; i != 50; ++i)
  {
    statistics.rigidTiltRotationChange.counts = 100.0;
    statistics.rigidTiltRotationChange.accepted = 0.0;
    statistics.optimize();
  }
  EXPECT_EQ(statistics.rigidTiltRotationChange.maxChange, statistics.rigidTiltRotationChange.lowerLimit);
}
