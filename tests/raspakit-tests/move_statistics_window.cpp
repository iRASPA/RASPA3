#include <gtest/gtest.h>

import std;

import archive;
import double3;
import move_statistics;
import mc_moves_move_types;
import mc_moves_statistics;

// The periodic status reports show the acceptance 'since the previous report' next to the cumulative
// one. The window is the difference between the cumulative counters and a snapshot taken when the
// previous report was written ('markReported'). The snapshot is part of the restart file (version 2 of
// the MoveStatistics record); a version-1 record must still read, with an empty snapshot.

TEST(MOVE_STATISTICS_WINDOW, window_is_difference_with_snapshot)
{
  MCMoveStatistics statistics;
  for (int i = 0; i != 10; ++i)
  {
    statistics.addTrial(Move::Types::ReinsertionCBMC);
    statistics.addConstructed(Move::Types::ReinsertionCBMC);
  }
  for (int i = 0; i != 4; ++i) statistics.addAccepted(Move::Types::ReinsertionCBMC);

  const MoveStatistics<double>& s = std::get<MoveStatistics<double>>(statistics[Move::Types::ReinsertionCBMC]);
  EXPECT_EQ(s.windowCounts(), 10.0);
  EXPECT_EQ(s.windowAccepted(), 4.0);
  EXPECT_EQ(s.windowConstructed(), 10.0);

  s.markReported();
  EXPECT_EQ(s.windowCounts(), 0.0);
  EXPECT_EQ(s.windowAccepted(), 0.0);
  EXPECT_EQ(s.totalCounts, 10.0);  // cumulative counters untouched
  EXPECT_EQ(s.totalAccepted, 4.0);

  for (int i = 0; i != 5; ++i) statistics.addTrial(Move::Types::ReinsertionCBMC);
  statistics.addAccepted(Move::Types::ReinsertionCBMC);
  EXPECT_EQ(s.windowCounts(), 5.0);
  EXPECT_EQ(s.windowAccepted(), 1.0);
  EXPECT_EQ(s.totalCounts, 15.0);

  // the step-size optimization resets the optimization counters but must not disturb the report window
  statistics.optimizeMCMoves();
  EXPECT_EQ(s.windowCounts(), 5.0);
  EXPECT_EQ(s.windowAccepted(), 1.0);
}

TEST(MOVE_STATISTICS_WINDOW, per_direction_window_for_double3_moves)
{
  MCMoveStatistics statistics;
  statistics.addTrial(Move::Types::Translation, 0);
  statistics.addTrial(Move::Types::Translation, 2);
  statistics.addAccepted(Move::Types::Translation, 2);

  const MoveStatistics<double3>& s = std::get<MoveStatistics<double3>>(statistics[Move::Types::Translation]);
  EXPECT_EQ(s.windowCounts().x, 1.0);
  EXPECT_EQ(s.windowCounts().y, 0.0);
  EXPECT_EQ(s.windowCounts().z, 1.0);
  EXPECT_EQ(s.windowAccepted().z, 1.0);

  s.markReported();
  statistics.addTrial(Move::Types::Translation, 1);
  EXPECT_EQ(s.windowCounts().x, 0.0);
  EXPECT_EQ(s.windowCounts().y, 1.0);
  EXPECT_EQ(s.windowCounts().z, 0.0);
}

TEST(MOVE_STATISTICS_WINDOW, archive_round_trip_keeps_snapshot)
{
  MoveStatistics<double> original{.maxChange = 0.3, .lowerLimit = 0.01, .upperLimit = 0.5};
  original.counts = 3.0;
  original.totalCounts = 120.0;
  original.totalConstructed = 110.0;
  original.totalAccepted = 40.0;
  original.reportedCounts = 100.0;
  original.reportedConstructed = 95.0;
  original.reportedAccepted = 30.0;

  const std::filesystem::path path = std::filesystem::temp_directory_path() / "raspa3_move_statistics_window.bin";
  {
    std::ofstream stream(path, std::ios::binary);
    Archive<std::ofstream> archive(stream);
    archive << original;
  }
  MoveStatistics<double> restored;
  {
    std::ifstream stream(path, std::ios::binary);
    Archive<std::ifstream> archive(stream);
    archive >> restored;
  }
  std::filesystem::remove(path);

  EXPECT_EQ(restored, original);
  EXPECT_EQ(restored.windowCounts(), 20.0);
  EXPECT_EQ(restored.windowAccepted(), 10.0);
}

TEST(MOVE_STATISTICS_WINDOW, version_1_record_reads_with_empty_snapshot)
{
  // hand-written version-1 record: everything up to and including 'optimize', no snapshot
  const std::filesystem::path path = std::filesystem::temp_directory_path() / "raspa3_move_statistics_window_v1.bin";
  {
    std::ofstream stream(path, std::ios::binary);
    Archive<std::ofstream> archive(stream);
    archive << std::uint64_t{1};
    archive << double{3.0};     // counts
    archive << double{3.0};     // constructed
    archive << double{1.0};     // accepted
    archive << std::size_t{7};  // allCounts
    archive << double{120.0};   // totalCounts
    archive << double{110.0};   // totalConstructed
    archive << double{40.0};    // totalAccepted
    archive << double{0.3};     // maxChange
    archive << double{0.5};     // targetAcceptance
    archive << double{0.01};    // lowerLimit
    archive << double{0.5};     // upperLimit
    archive << bool{true};      // optimize
  }
  MoveStatistics<double> restored;
  restored.reportedCounts = 99.0;  // must be reset by a version-1 read
  {
    std::ifstream stream(path, std::ios::binary);
    Archive<std::ifstream> archive(stream);
    archive >> restored;
  }
  std::filesystem::remove(path);

  EXPECT_EQ(restored.totalCounts, 120.0);
  EXPECT_EQ(restored.totalAccepted, 40.0);
  EXPECT_EQ(restored.maxChange, 0.3);
  EXPECT_EQ(restored.reportedCounts, 0.0);
  EXPECT_EQ(restored.windowCounts(), 120.0);  // first report after the restart covers the whole history
}
