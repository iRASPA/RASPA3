#include <gtest/gtest.h>

import std;

import double3;
import move_statistics;
import mc_moves_move_types;
import mc_moves_statistics;

// A double3 statistics holds three channels (the Cartesian directions of a translation, or the insertion,
// deletion and lambda-change sub-moves of a CFCMC swap) that are optimized on their own acceptance. A channel
// without trials in the optimization window has no acceptance to optimize on and keeps its step.

TEST(MOVE_STATISTICS_OPTIMIZE, channels_without_trials_keep_their_step)
{
  MoveStatistics<double3> s{.maxChange = double3(1.0), .lowerLimit = double3(0.01), .upperLimit = double3(1.5)};
  s.counts.x = 10.0;
  s.accepted.x = 5.0;
  s.optimizeAcceptance();
  EXPECT_NEAR(s.maxChange.x, (5.0 / 11.0) / 0.5, 1e-12);
  EXPECT_EQ(s.maxChange.y, 1.0);
  EXPECT_EQ(s.maxChange.z, 1.0);
  EXPECT_EQ(s.counts, double3(0.0));

  s.counts.z = 10.0;
  s.accepted.z = 10.0;
  s.optimizeAcceptance();
  EXPECT_NEAR(s.maxChange.x, (5.0 / 11.0) / 0.5, 1e-12);
  EXPECT_EQ(s.maxChange.y, 1.0);
  EXPECT_EQ(s.maxChange.z, 1.5);
}

TEST(MOVE_STATISTICS_OPTIMIZE, cfcmc_lambda_step_follows_the_lambda_changes_only)
{
  MCMoveStatistics statistics;
  const Move::Types move = Move::Types::GibbsConventionalCFCMC;
  const double initialStep = statistics.getMaxChange(move, 2);

  // a window of insertions and deletions, none accepted, and lambda changes accepted at the target rate
  for (int i = 0; i != 10; ++i)
  {
    statistics.addTrial(move, 0);
    statistics.addTrial(move, 1);
    statistics.addTrial(move, 2);
  }
  for (int i = 0; i != 5; ++i) statistics.addAccepted(move, 2);
  statistics.optimizeMCMoves();
  EXPECT_NEAR(statistics.getMaxChange(move, 2), initialStep * (5.0 / 11.0) / 0.5, 1e-12);

  // a window with lambda changes alone, accepted at the target rate: the step is kept
  const double step = statistics.getMaxChange(move, 2);
  for (int i = 0; i != 10; ++i) statistics.addTrial(move, 2);
  for (int i = 0; i != 5; ++i) statistics.addAccepted(move, 2);
  statistics.optimizeMCMoves();
  EXPECT_NEAR(statistics.getMaxChange(move, 2), step * (5.0 / 11.0) / 0.5, 1e-12);

  // a window with only insertions and deletions: the lambda step is not touched
  const double kept = statistics.getMaxChange(move, 2);
  for (int i = 0; i != 10; ++i) statistics.addTrial(move, 0);
  statistics.optimizeMCMoves();
  EXPECT_EQ(statistics.getMaxChange(move, 2), kept);
}
