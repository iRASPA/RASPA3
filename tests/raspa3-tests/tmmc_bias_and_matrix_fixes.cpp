#include <gtest/gtest.h>

import std;

import double3;
import forcefield;
import simulationbox;
import component;
import system;
import transition_matrix;
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_parallel_tempering_swap;

// Regression tests for the transition-matrix bias of Gibbs transfers and replica swaps, the collection-matrix
// update and the Wang-Landau handover.

namespace
{

System makeTMMCBox(std::size_t numberOfMolecules, std::size_t minN, std::size_t maxN)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component co2 = Component::makeCO2(forceField, 0, true);
  co2.mc_moves_probabilities.setProbability(Move::Types::Translation, 1.0);
  co2.mc_moves_probabilities.setProbability(Move::Types::Swap, 1.0);
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {co2}, {},
                {numberOfMolecules}, 5);
  system.tmmc.doTMMC = true;
  system.tmmc.useBias = true;
  system.tmmc.useTMBias = true;
  system.tmmc.useWangLandau = false;
  system.tmmc.rejectOutOfBound = true;
  system.tmmc.minMacrostate = minN;
  system.tmmc.maxMacrostate = maxN;
  system.tmmc.initialize();
  return system;
}

}  // namespace

// Both boxes of a Gibbs pair rebuild the same distribution of the shared macrostate; the transfer must carry one
// bias, not the product of the two (which squared the flattening and drove the walk to 1 / Pi).
TEST(TMMC_BIAS_AND_MATRIX_FIXES, gibbs_transfer_carries_the_bias_of_one_box)
{
  std::vector<System> boxes;
  boxes.push_back(makeTMMCBox(3, 0, 10));
  boxes.push_back(makeTMMCBox(5, 0, 10));
  for (std::size_t n = 0; n <= 10; ++n)
  {
    // the same distribution seen from either box: bias_A(N_A) = bias_B(8 - N_A)
    boxes[0].tmmc.bias[n] = -0.7 * static_cast<double>(n);
    boxes[1].tmmc.bias[n] = n <= 8 ? -0.7 * static_cast<double>(8 - n) : 0.0;
  }
  System& first = boxes[0];
  System& second = boxes[1];
  const double single = std::exp(first.tmmc.bias[4] - first.tmmc.bias[3]);

  // first gains a molecule from second, and the reverse: both see the bias of the first box
  EXPECT_NEAR(System::gibbsTransferTMMCBiasFactor(first, second, 0, 3, 5), single, 1e-12);
  EXPECT_NEAR(System::gibbsTransferTMMCBiasFactor(second, first, 0, 4, 4), 1.0 / single, 1e-12);

  // the window of the other box still refuses a transfer out of it
  second.tmmc.minMacrostate = 5;
  second.tmmc.initialize();
  EXPECT_EQ(System::gibbsTransferTMMCBiasFactor(first, second, 0, 3, 5), 0.0);
}

// A non-finite acceptance must not poison the collection matrix: std::clamp passes a NaN through.
TEST(TMMC_BIAS_AND_MATRIX_FIXES, collection_matrix_ignores_non_finite_acceptances)
{
  TransitionMatrix matrix;
  matrix.doTMMC = true;
  matrix.minMacrostate = 0;
  matrix.maxMacrostate = 4;
  matrix.initialize();

  const double nan = std::numeric_limits<double>::quiet_NaN();
  matrix.updateMatrix(double3(0.0, nan, nan), 1);
  EXPECT_EQ(matrix.cmatrix[1].x, 0.0);
  EXPECT_EQ(matrix.cmatrix[1].y, 0.0);
  EXPECT_EQ(matrix.cmatrix[1].z, 0.0);

  // an overflowed insertion acceptance is an accepted insertion
  const double inf = std::numeric_limits<double>::max() * 10.0;
  matrix.updateMatrix(double3(0.0, 1.0 - inf, inf), 2);
  EXPECT_EQ(matrix.cmatrix[2].z, 1.0);
  EXPECT_EQ(matrix.cmatrix[2].y, 0.0);

  // an acceptance slightly above one is clamped, not dropped
  matrix.updateMatrix(double3(0.0, -0.5, 1.5), 3);
  EXPECT_EQ(matrix.cmatrix[3].z, 1.0);
  EXPECT_EQ(matrix.cmatrix[3].y, 0.0);
}

// Without a transition-matrix bias to take over, the handover freezes Wang-Landau instead of leaving it running.
TEST(TMMC_BIAS_AND_MATRIX_FIXES, handover_freezes_wang_landau_without_a_tm_bias)
{
  TransitionMatrix matrix;
  matrix.doTMMC = true;
  matrix.useBias = true;
  matrix.useTMBias = false;
  matrix.useWangLandau = true;
  matrix.minMacrostate = 0;
  matrix.maxMacrostate = 4;
  matrix.initialize();
  matrix.switchToTMBias();
  EXPECT_EQ(matrix.wangLandauFactor, 0.0);

  const std::vector<double> before = matrix.bias;
  matrix.visitWangLandau(2, 0);
  EXPECT_EQ(matrix.bias, before);
}

TEST(TMMC_BIAS_AND_MATRIX_FIXES, reweighted_loading_from_ln_pi)
{
  TransitionMatrix matrix;
  matrix.doTMMC = true;
  matrix.minMacrostate = 0;
  matrix.maxMacrostate = 2;
  matrix.initialize();
  matrix.lnpi = {0.0, std::log(2.0), 0.0};
  matrix.histogram = {1, 1, 1};
  // both links carry accepted transitions in both directions
  matrix.cmatrix = {double3(0.0, 0.5, 0.5), double3(0.5, 0.0, 0.5), double3(0.5, 0.5, 0.0)};
  const std::optional<std::pair<double, double>> moments = matrix.reweightedMoleculeCountMeanAndVariance();
  ASSERT_TRUE(moments.has_value());
  EXPECT_NEAR(moments->first, 1.0, 1e-12);
  EXPECT_NEAR(moments->second, 0.5, 1e-12);
}

// A configuration outside a window that rejects out-of-bound states is refused by the replica swap, whether or
// not that window uses a bias (it was accepted when the bias was off).
TEST(TMMC_BIAS_AND_MATRIX_FIXES, replica_swap_respects_the_window_without_a_bias)
{
  System wide = makeTMMCBox(4, 0, 10);
  System narrow = makeTMMCBox(1, 0, 2);
  wide.tmmc.useBias = false;
  narrow.tmmc.useBias = false;
  EXPECT_FALSE(MC_Moves::ParallelTemperingLogAcceptance(wide, narrow).has_value());

  narrow.tmmc.rejectOutOfBound = false;
  EXPECT_TRUE(MC_Moves::ParallelTemperingLogAcceptance(wide, narrow).has_value());
}
