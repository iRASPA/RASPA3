#include <gtest/gtest.h>

import std;

import randomnumbers;
import dctst_free_energy;

// The reinsertion support of a mirror-symmetric umbrella window and the proposal drawn from it. In mirror mode
// the coordinate is folded to |x|, and the draw is mirrored with probability one half; the support must then
// start at 0 at the lowest, or the part of it below 0 is reached by both signs and gets twice the proposal
// density, which the acceptance rule (no proposal correction) turns into a bias towards the barrier.

namespace
{

constexpr double gridBegin = -1.0;
constexpr double gridEnd = 1.0;

// the fraction of 'samples' draws whose folded coordinate falls in [a, b)
double foldedFraction(RandomNumber& rng, double low, double high, double a, double b, std::size_t samples)
{
  std::size_t hits = 0;
  for (std::size_t i = 0; i < samples; ++i)
  {
    const double x = std::abs(dctstSampleReinsertionCoordinate(rng, low, high, true));
    if (x >= a && x < b) ++hits;
  }
  return static_cast<double>(hits) / static_cast<double>(samples);
}

}  // namespace

TEST(DCTST_MIRROR_REINSERTION, a_window_clear_of_the_barrier_keeps_its_full_support)
{
  const auto [low, high] = dctstFoldedReinsertionSupport(0.6, 0.3, gridBegin, gridEnd, true);
  EXPECT_DOUBLE_EQ(low, 0.3);
  EXPECT_DOUBLE_EQ(high, 0.9);
}

TEST(DCTST_MIRROR_REINSERTION, a_near_barrier_window_has_its_support_clipped_at_the_fold)
{
  // center 0.2, half-width 0.5: unclipped, the support would run from -0.3 to 0.7
  const auto [low, high] = dctstFoldedReinsertionSupport(0.2, 0.5, gridBegin, gridEnd, true);
  EXPECT_DOUBLE_EQ(low, 0.0);
  EXPECT_DOUBLE_EQ(high, 0.7);

  // the window on the barrier itself
  const auto [lowBarrier, highBarrier] = dctstFoldedReinsertionSupport(0.0, 0.5, gridBegin, gridEnd, true);
  EXPECT_DOUBLE_EQ(lowBarrier, 0.0);
  EXPECT_DOUBLE_EQ(highBarrier, 0.5);
}

TEST(DCTST_MIRROR_REINSERTION, without_mirror_symmetry_the_support_is_clipped_to_the_grid_only)
{
  const auto [low, high] = dctstFoldedReinsertionSupport(0.2, 0.5, gridBegin, gridEnd, false);
  EXPECT_DOUBLE_EQ(low, -0.3);
  EXPECT_DOUBLE_EQ(high, 0.7);

  const auto [lowEdge, highEdge] = dctstFoldedReinsertionSupport(-0.9, 0.5, gridBegin, gridEnd, false);
  EXPECT_DOUBLE_EQ(lowEdge, gridBegin);
  EXPECT_DOUBLE_EQ(highEdge, -0.4);
}

TEST(DCTST_MIRROR_REINSERTION, the_mirrored_draw_is_uniform_in_the_folded_coordinate_over_the_support)
{
  // the near-barrier window of the test above: before the clip, |x| < 0.3 had twice the density of the rest
  const auto [low, high] = dctstFoldedReinsertionSupport(0.2, 0.5, gridBegin, gridEnd, true);
  RandomNumber rng(7);
  constexpr std::size_t samples = 400000;

  const double near = foldedFraction(rng, low, high, 0.0, 0.3, samples);
  const double far = foldedFraction(rng, low, high, 0.3, 0.6, samples);
  EXPECT_NEAR(near, 3.0 / 7.0, 0.005);
  EXPECT_NEAR(far, 3.0 / 7.0, 0.005);
  EXPECT_NEAR(near / far, 1.0, 0.02);

  // every draw lies in the support, with either sign equally often
  std::size_t negatives = 0;
  for (std::size_t i = 0; i < samples; ++i)
  {
    const double x = dctstSampleReinsertionCoordinate(rng, low, high, true);
    ASSERT_LT(std::abs(x), high);
    if (x < 0.0) ++negatives;
  }
  EXPECT_NEAR(static_cast<double>(negatives) / static_cast<double>(samples), 0.5, 0.005);
}

TEST(DCTST_MIRROR_REINSERTION, without_mirror_symmetry_the_draw_stays_on_the_interval)
{
  RandomNumber rng(11);
  for (std::size_t i = 0; i < 10000; ++i)
  {
    const double x = dctstSampleReinsertionCoordinate(rng, -0.3, 0.7, false);
    ASSERT_GE(x, -0.3);
    ASSERT_LT(x, 0.7);
  }
}
