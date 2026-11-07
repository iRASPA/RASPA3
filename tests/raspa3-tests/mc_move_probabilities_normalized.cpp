#include <gtest/gtest.h>

import std;

import mc_moves_move_types;
import mc_moves_probabilities;

// The normalized move probabilities of a component: fractions summing to one, or all zero for a component
// without any move (the division by the zero total gave not-a-number in the move-probability output).
TEST(MC_MOVE_PROBABILITIES_NORMALIZED, no_moves_normalize_to_zero_and_moves_to_fractions)
{
  const std::vector<double> none = MCMoveProbabilities().normalizedMap();
  EXPECT_TRUE(std::ranges::all_of(none, [](double p) { return p == 0.0; }));

  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Rotation, 3.0);
  const std::vector<double> normalized = probabilities.normalizedMap();
  EXPECT_DOUBLE_EQ(normalized[std::to_underlying(Move::Types::Translation)], 0.25);
  EXPECT_DOUBLE_EQ(normalized[std::to_underlying(Move::Types::Rotation)], 0.75);
  EXPECT_DOUBLE_EQ(std::accumulate(normalized.begin(), normalized.end(), 0.0), 1.0);
}
