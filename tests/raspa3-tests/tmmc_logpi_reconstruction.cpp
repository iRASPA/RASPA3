#include <gtest/gtest.h>

import std;

import double3;
import parallel_tmmc;

// ln Pi from a pooled collection matrix with the Widom anchor laid over its low-N increments. The anchor
// replaces the collection-matrix increments inside its block, so a link the collection matrix cannot
// support there (no accepted deletion at low N, a starved λ hop) is bridged rather than cutting off every
// macrostate above it; outside the block the collection-matrix increments stand, shifted as one piece.

namespace
{

constexpr double logZero = -1e300;
constexpr double logZeroThreshold = -1e299;

bool sampled(double value) { return value > logZeroThreshold; }

// A row with the given accepted-up and accepted-down probability sums out of 'attempts' attempts each way.
double3 row(double down, double up, double attempts = 10.0) { return double3(down, 2.0 * attempts - down - up, up); }

// The collection-matrix increment of the link i -> i+1.
double linkIncrement(const std::vector<double3>& matrix, std::size_t i)
{
  const double forwardTotal = matrix[i].x + matrix[i].y + matrix[i].z;
  const double reverseTotal = matrix[i + 1].x + matrix[i + 1].y + matrix[i + 1].z;
  return std::log(matrix[i].z) - std::log(forwardTotal) - std::log(matrix[i + 1].x) + std::log(reverseTotal);
}

std::vector<double> noIncrements(std::size_t n) { return std::vector<double>(n, 0.0); }
std::vector<bool> noFlags(std::size_t n) { return std::vector<bool>(n, false); }

}  // namespace

TEST(TMMC_LOGPI_RECONSTRUCTION, without_an_anchor_the_chain_is_the_detailed_balance_recursion)
{
  const std::vector<double3> matrix = {row(0.0, 4.0), row(1.0, 3.0), row(2.0, 2.0), row(3.0, 1.0), row(4.0, 0.0)};
  const std::vector<double> logPi =
      tmmcReconstructLogPi(matrix, 5, 1, noIncrements(5), noFlags(5), 0, logZero, logZeroThreshold);

  ASSERT_EQ(logPi.size(), 5uz);
  EXPECT_DOUBLE_EQ(logPi[0], 0.0);
  for (std::size_t i = 0; i + 1 < 5; ++i)
  {
    EXPECT_NEAR(logPi[i + 1] - logPi[i], linkIncrement(matrix, i), 1e-12) << "link " << i;
  }
}

TEST(TMMC_LOGPI_RECONSTRUCTION, without_an_anchor_an_unsupported_link_truncates_the_chain)
{
  // no accepted deletion out of N = 2: the link 1 -> 2 has no reverse rate
  const std::vector<double3> matrix = {row(0.0, 4.0), row(1.0, 3.0), row(0.0, 2.0), row(3.0, 1.0), row(4.0, 0.0)};
  const std::vector<double> logPi =
      tmmcReconstructLogPi(matrix, 5, 1, noIncrements(5), noFlags(5), 0, logZero, logZeroThreshold);

  EXPECT_TRUE(sampled(logPi[0]));
  EXPECT_TRUE(sampled(logPi[1]));
  EXPECT_FALSE(sampled(logPi[2]));
  EXPECT_FALSE(sampled(logPi[3]));
  EXPECT_FALSE(sampled(logPi[4]));
}

TEST(TMMC_LOGPI_RECONSTRUCTION, the_anchor_bridges_an_unsupported_low_N_link_and_keeps_the_chain_above_it)
{
  // the strongly binding adsorbent: not one deletion accepted out of N = 1 or N = 2, the two links below
  // the anchor's end are unsupported, the links above it are fine
  const std::vector<double3> matrix = {row(0.0, 4.0), row(0.0, 3.0), row(0.0, 2.5), row(2.0, 2.0),
                                       row(3.0, 1.0), row(3.5, 0.5), row(4.0, 0.0)};
  std::vector<double> widomLogIncrement = noIncrements(7);
  std::vector<bool> widomHasIncrement = noFlags(7);
  widomLogIncrement[0] = 5.0;
  widomLogIncrement[1] = 3.0;
  widomLogIncrement[2] = 1.0;
  widomHasIncrement[0] = widomHasIncrement[1] = widomHasIncrement[2] = true;
  const std::size_t widomAnchorEnd = 3;

  const std::vector<double> logPi = tmmcReconstructLogPi(matrix, 7, 1, widomLogIncrement, widomHasIncrement,
                                                         widomAnchorEnd, logZero, logZeroThreshold);

  ASSERT_EQ(logPi.size(), 7uz);
  for (std::size_t n = 0; n < 7; ++n)
  {
    EXPECT_TRUE(sampled(logPi[n])) << "N = " << n << " was cut off";
  }
  // the anchored increments inside the block
  EXPECT_NEAR(logPi[1] - logPi[0], 5.0, 1e-12);
  EXPECT_NEAR(logPi[2] - logPi[1], 3.0, 1e-12);
  EXPECT_NEAR(logPi[3] - logPi[2], 1.0, 1e-12);
  // the collection-matrix increments above it, untouched
  for (std::size_t i = 3; i + 1 < 7; ++i)
  {
    EXPECT_NEAR(logPi[i + 1] - logPi[i], linkIncrement(matrix, i), 1e-12) << "link " << i;
  }
}

TEST(TMMC_LOGPI_RECONSTRUCTION, a_counted_increment_at_the_first_macrostate_gauges_an_unvisited_empty_box)
{
  // the walker never sat at N = 0 (row empty), but the test insertions at N = 0 were counted
  const std::vector<double3> matrix = {double3(0.0, 0.0, 0.0), row(1.0, 3.0), row(2.0, 2.0), row(3.0, 1.0)};
  std::vector<double> widomLogIncrement = noIncrements(4);
  std::vector<bool> widomHasIncrement = noFlags(4);
  widomLogIncrement[0] = 2.0;
  widomHasIncrement[0] = true;

  const std::vector<double> logPi =
      tmmcReconstructLogPi(matrix, 4, 1, widomLogIncrement, widomHasIncrement, 1, logZero, logZeroThreshold);

  EXPECT_DOUBLE_EQ(logPi[0], 0.0);
  EXPECT_NEAR(logPi[1], 2.0, 1e-12);
  EXPECT_NEAR(logPi[2] - logPi[1], linkIncrement(matrix, 1), 1e-12);
  EXPECT_NEAR(logPi[3] - logPi[2], linkIncrement(matrix, 2), 1e-12);
}

TEST(TMMC_LOGPI_RECONSTRUCTION, a_hole_in_the_anchor_block_falls_back_on_the_collection_matrix)
{
  // the increment out of N = 1 has no counted test-insertion mean but a supported link; the one out
  // of N = 2 is anchored again. The block spans 0..2 and the hole is filled from the collection matrix.
  const std::vector<double3> matrix = {row(0.0, 4.0), row(1.0, 3.0), row(2.0, 2.0), row(3.0, 1.0), row(4.0, 0.0)};
  std::vector<double> widomLogIncrement = noIncrements(5);
  std::vector<bool> widomHasIncrement = noFlags(5);
  widomLogIncrement[0] = 4.0;
  widomLogIncrement[2] = 1.5;
  widomHasIncrement[0] = widomHasIncrement[2] = true;

  const std::vector<double> logPi =
      tmmcReconstructLogPi(matrix, 5, 1, widomLogIncrement, widomHasIncrement, 3, logZero, logZeroThreshold);

  EXPECT_NEAR(logPi[1] - logPi[0], 4.0, 1e-12);
  EXPECT_NEAR(logPi[2] - logPi[1], linkIncrement(matrix, 1), 1e-12);
  EXPECT_NEAR(logPi[3] - logPi[2], 1.5, 1e-12);
  EXPECT_NEAR(logPi[4] - logPi[3], linkIncrement(matrix, 3), 1e-12);
}

TEST(TMMC_LOGPI_RECONSTRUCTION, an_N_lambda_chain_sums_its_lambda_links_and_the_anchor_bridges_a_starved_hop)
{
  // two λ bins per macrostate, three macrostates; the hop (1, 1) -> (2, 0) was never accepted downward
  const std::vector<double3> matrix = {row(0.0, 4.0), row(1.0, 3.0), row(1.5, 2.5),
                                       row(2.0, 2.0), row(0.0, 1.5), row(3.0, 0.0)};
  const std::size_t nLambda = 2;

  // unanchored: N = 0 -> 1 is the sum of two links, N = 2 is cut off
  {
    const std::vector<double> logPi =
        tmmcReconstructLogPi(matrix, 3, nLambda, noIncrements(3), noFlags(3), 0, logZero, logZeroThreshold);
    ASSERT_EQ(logPi.size(), 3uz);
    EXPECT_NEAR(logPi[1] - logPi[0], linkIncrement(matrix, 0) + linkIncrement(matrix, 1), 1e-12);
    EXPECT_FALSE(sampled(logPi[2]));
  }

  // anchored out of N = 1: the starved hop is bridged
  {
    std::vector<double> widomLogIncrement = noIncrements(3);
    std::vector<bool> widomHasIncrement = noFlags(3);
    widomLogIncrement[1] = -0.7;
    widomHasIncrement[1] = true;
    const std::vector<double> logPi =
        tmmcReconstructLogPi(matrix, 3, nLambda, widomLogIncrement, widomHasIncrement, 2, logZero, logZeroThreshold);
    EXPECT_NEAR(logPi[1] - logPi[0], linkIncrement(matrix, 0) + linkIncrement(matrix, 1), 1e-12);
    EXPECT_NEAR(logPi[2] - logPi[1], -0.7, 1e-12);
  }
}

TEST(TMMC_LOGPI_RECONSTRUCTION, nothing_sampled_leaves_every_state_at_the_sentinel)
{
  const std::vector<double3> matrix(4, double3(0.0, 0.0, 0.0));
  const std::vector<double> logPi =
      tmmcReconstructLogPi(matrix, 4, 1, noIncrements(4), noFlags(4), 0, logZero, logZeroThreshold);
  for (double value : logPi) EXPECT_FALSE(sampled(value));
}
