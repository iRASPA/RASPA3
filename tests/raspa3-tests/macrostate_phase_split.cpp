#include <gtest/gtest.h>

import std;

import isotherm_bet;

// Where a macrostate distribution ln P(N) is cut in two phases, and when it is not cut at all. An N that was
// never sampled takes the lower of its sampled neighbours and is never a valley on its own account.

namespace
{

constexpr double logZero = -1e300;
constexpr double logZeroThreshold = -1e299;
constexpr double depthThreshold = 0.2;

std::optional<std::pair<std::size_t, double>> splitOf(const std::vector<double>& logProbability)
{
  return findMacrostatePhaseSplit(logProbability, logZeroThreshold, depthThreshold);
}

// A peak of the given height at the given centre, falling off quadratically, over N = 0 ... size - 1.
std::vector<double> peak(std::size_t size, double centre, double height, double width)
{
  std::vector<double> values(size);
  for (std::size_t n = 0; n < size; ++n)
  {
    const double x = (static_cast<double>(n) - centre) / width;
    values[n] = height - x * x;
  }
  return values;
}

std::vector<double> logSumOf(const std::vector<double>& a, const std::vector<double>& b)
{
  std::vector<double> values(a.size());
  for (std::size_t n = 0; n < a.size(); ++n)
  {
    const double largest = std::max(a[n], b[n]);
    values[n] = largest + std::log(std::exp(a[n] - largest) + std::exp(b[n] - largest));
  }
  return values;
}

}  // namespace


TEST(MACROSTATE_PHASE_SPLIT, two_peaks_are_cut_at_the_valley_between_them)
{
  const std::vector<double> bimodal = logSumOf(peak(60, 10.0, 0.0, 4.0), peak(60, 45.0, -1.0, 6.0));
  const std::optional<std::pair<std::size_t, double>> split = splitOf(bimodal);
  ASSERT_TRUE(split.has_value());
  // the valley lies between the peaks, nearer the narrow one, and is as deep as the lower peak stands above it
  EXPECT_GT(split->first, 15uz);
  EXPECT_LT(split->first, 40uz);
  const double valley = *std::min_element(bimodal.begin() + 15, bimodal.begin() + 40);
  EXPECT_NEAR(split->second, -1.0 - valley, 1.0e-12);
  EXPECT_GT(split->second, 5.0);
}


TEST(MACROSTATE_PHASE_SPLIT, one_peak_is_not_cut)
{
  EXPECT_FALSE(splitOf(peak(60, 30.0, 0.0, 8.0)).has_value());
}


// A hole in the sampling inside a peak is read as the peak around it. Filled twenty below everything
// sampled, as it was, it cut the one peak in two at the hole with a depth of twenty.
TEST(MACROSTATE_PHASE_SPLIT, an_unsampled_state_inside_a_peak_is_not_a_valley)
{
  for (std::size_t hole : {30uz, 36uz, 44uz})
  {
    std::vector<double> holed = peak(60, 30.0, 0.0, 8.0);
    holed[hole] = logZero;
    EXPECT_FALSE(splitOf(holed).has_value()) << "hole at " << hole;
  }

  // several holes, as one block of many has them high up where it visits seldom
  std::vector<double> sparse = peak(60, 30.0, 0.0, 8.0);
  for (std::size_t n = 41; n < 60; n += 2) sparse[n] = logZero;
  EXPECT_FALSE(splitOf(sparse).has_value());

  // nor do the ends that were never reached make one
  std::vector<double> clipped = peak(60, 30.0, 0.0, 8.0);
  for (std::size_t n = 0; n < 8; ++n) clipped[n] = logZero;
  for (std::size_t n = 52; n < 60; ++n) clipped[n] = logZero;
  EXPECT_FALSE(splitOf(clipped).has_value());
}


// Two peaks with the valley between them never visited: the cut is as deep as the flanks have come down,
// no deeper, and it still finds the two phases where the flanks are sampled into the valley.
TEST(MACROSTATE_PHASE_SPLIT, a_gap_between_two_peaks_is_as_deep_as_its_sampled_flanks)
{
  std::vector<double> gapped = logSumOf(peak(60, 10.0, 0.0, 4.0), peak(60, 45.0, -1.0, 6.0));
  for (std::size_t n = 22; n <= 32; ++n) gapped[n] = logZero;
  const std::optional<std::pair<std::size_t, double>> split = splitOf(gapped);
  ASSERT_TRUE(split.has_value());
  EXPECT_GE(split->first, 21uz);
  EXPECT_LE(split->first, 33uz);
  EXPECT_NEAR(split->second, -1.0 - std::min(gapped[21], gapped[33]), 1.0e-12);

  // the same two peaks sampled at their tops only, with nothing in between that has come down: there is
  // no valley in evidence, and no cut
  std::vector<double> tops(60, logZero);
  for (std::size_t n = 8; n <= 12; ++n) tops[n] = 0.0;
  for (std::size_t n = 43; n <= 47; ++n) tops[n] = -1.0;
  EXPECT_FALSE(splitOf(tops).has_value());
}


TEST(MACROSTATE_PHASE_SPLIT, nothing_sampled_is_not_cut)
{
  EXPECT_FALSE(splitOf(std::vector<double>(60, logZero)).has_value());
  EXPECT_FALSE(splitOf(std::vector<double>{0.0, 0.0}).has_value());
}
