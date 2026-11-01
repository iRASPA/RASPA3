#include <gtest/gtest.h>

import std;

import double3;
import simulationbox;
import randomnumbers;
import dctst_barrier;

TEST(dcTST, periodic_window_support_wraps_through_the_period)
{
  const SimulationBox box(10.0, 10.0, 10.0);
  const DcTSTBarrier barrier(box, double3(0.5, 0.5, 0.5), double3(1.0, 0.0, 0.0), 10.0);
  const double lambda = barrier.lambda;
  ASSERT_DOUBLE_EQ(lambda, 10.0);

  // a window centered on the identified endpoints ±λ/2: the support is the union of a cap at
  // −λ/2 and a cap at +λ/2. The old reinsertion sampler clipped to [center − w, center + w] ∩
  // [−λ/2, λ/2) and so never proposed the positive cap, while membership used the periodic
  // distance and accepted the positive cap as 'inside' — q(n → o) = 0 for those jumps.
  const double center = -0.5 * lambda;
  const double halfWidth = 0.2 * lambda;
  EXPECT_TRUE(barrier.inPeriodicWindow(center + 0.5 * halfWidth, center, halfWidth));
  EXPECT_TRUE(barrier.inPeriodicWindow(0.5 * lambda - 0.5 * halfWidth, center, halfWidth));
  EXPECT_FALSE(barrier.inPeriodicWindow(0.0, center, halfWidth));

  const double clipLow = std::max(-0.5 * lambda, center - halfWidth);
  const double clipHigh = std::min(0.5 * lambda, center + halfWidth);

  RandomNumber random(1);
  std::size_t negativeCap = 0, positiveCap = 0, outsideCartesianClip = 0;
  constexpr std::size_t samples = 8000;
  for (std::size_t i = 0; i != samples; ++i)
  {
    const double x = barrier.samplePeriodicWindow(random, center, halfWidth);
    EXPECT_TRUE(barrier.inPeriodicWindow(x, center, halfWidth));
    EXPECT_GE(x, -0.5 * lambda);
    EXPECT_LT(x, 0.5 * lambda);
    if (x < 0.0)
      ++negativeCap;
    else
      ++positiveCap;
    if (x < clipLow || x >= clipHigh) ++outsideCartesianClip;
  }
  EXPECT_GT(negativeCap, samples / 5);
  EXPECT_GT(positiveCap, samples / 5);
  EXPECT_GT(outsideCartesianClip, samples / 5) << "the reverse of a wrap must be proposable";

  // a window that does not wrap stays a single interval
  const double interior = 0.0;
  for (std::size_t i = 0; i != 200; ++i)
  {
    const double x = barrier.samplePeriodicWindow(random, interior, halfWidth);
    EXPECT_LE(std::abs(x - interior), halfWidth + 1.0e-12);
  }
}
