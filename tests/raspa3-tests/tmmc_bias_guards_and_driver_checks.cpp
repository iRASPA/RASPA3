#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import atom;
import double3;
import forcefield;
import simulationbox;
import component;
import system;
import transition_matrix;
import parallel_tmmc;
import mc_moves_acceptance;
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_parallel_tempering_swap;

// Regression tests for the Metropolis test under a zero transition-matrix bias, the window bound of the Gibbs box
// that does not carry the bias, the unclamped log-bias of a replica swap, the weight of a configuration on an
// (N, lambda) chain, ln Pi after the statistics are dropped, and the checks of the parallel driver.

namespace
{

class ScopedCurrentPath
{
 public:
  explicit ScopedCurrentPath(const std::filesystem::path& path) : original_(std::filesystem::current_path())
  {
    std::filesystem::create_directories(path);
    std::filesystem::current_path(path);
  }

  ~ScopedCurrentPath()
  {
    std::error_code ignored;
    const std::filesystem::path used = std::filesystem::current_path();
    std::filesystem::current_path(original_);
    std::filesystem::remove_all(used, ignored);
  }

 private:
  std::filesystem::path original_;
};

System makeTMMCBox(std::size_t numberOfMolecules, std::size_t minN, std::size_t maxN,
                   Move::Types swapMove = Move::Types::Swap)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component co2 = Component::makeCO2(forceField, 0, true);
  co2.mc_moves_probabilities.setProbability(Move::Types::Translation, 1.0);
  co2.mc_moves_probabilities.setProbability(swapMove, 1.0);
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

System makeMethaneBox(std::size_t numberOfMolecules, double pressure)
{
  const ForceField forceField({{"CH4", false, 16.04246, 0.0, 0.0, 6, false}}, {{158.5, 3.72}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Swap, 1.0);
  const Component methane(forceField, "methane", 190.564, 4599200, 0.01142,
                          {Atom({0, 0, 0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21, probabilities,
                          std::nullopt, false);
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, pressure, 1.0, {}, {methane}, {},
                {numberOfMolecules}, 5);
  system.tmmc.minMacrostate = 0;
  system.tmmc.maxMacrostate = 12;
  return system;
}

std::string refusalOfParallelTMMC(System system, std::vector<double> temperatures, ParallelTMMCParameters parameters)
{
  try
  {
    ParallelTMMC simulation(std::move(system), std::move(temperatures), parameters);
  }
  catch (const std::runtime_error& error)
  {
    return error.what();
  }
  return {};
}

}  // namespace

// A zero bias factor is a refused move whatever the probability: zero times an overflowed probability is not a
// number, and under -ffast-math the comparison against it need not fail.
TEST(TMMC_BIAS_GUARDS, a_zero_bias_refuses_an_overflowed_probability)
{
  const double overflowed = std::numeric_limits<double>::infinity();
  EXPECT_FALSE(MC_Moves::acceptedWithBias(0.0, 0.0, overflowed));
  EXPECT_FALSE(MC_Moves::acceptedWithBias(0.0, 0.0, 1.0e308));
  EXPECT_FALSE(MC_Moves::acceptedWithBias(0.0, 0.0, 1.0));
  EXPECT_FALSE(MC_Moves::acceptedWithBias(0.0, 1.0, 0.0));

  // a biased move with an overflowed probability, or with a product beyond a double, is taken
  EXPECT_TRUE(MC_Moves::acceptedWithBias(0.999, std::exp(-500.0), overflowed));
  EXPECT_TRUE(MC_Moves::acceptedWithBias(0.999, std::exp(500.0), 1.0e200));

  // and the ordinary test is the product
  EXPECT_TRUE(MC_Moves::acceptedWithBias(0.24, 0.5, 0.5));
  EXPECT_FALSE(MC_Moves::acceptedWithBias(0.26, 0.5, 0.5));
  EXPECT_TRUE(MC_Moves::acceptedWithBias(0.5, std::exp(500.0), std::exp(-500.0) * 0.6));
  EXPECT_FALSE(MC_Moves::acceptedWithBias(0.5, std::exp(500.0), std::exp(-500.0) * 0.4));
}

// A box above its window is let back in by the transfers that bring it closer, as the box that carries the bias
// is. Held to the new count alone, the transfer out of it was refused while the window of the other box refused
// the reverse, and neither box moved for the whole run.
TEST(TMMC_BIAS_GUARDS, a_gibbs_box_outside_its_window_is_let_back_in)
{
  std::vector<System> boxes;
  boxes.push_back(makeTMMCBox(10, 10, 20));
  boxes.push_back(makeTMMCBox(22, 10, 20));
  System& first = boxes[0];
  System& second = boxes[1];

  // the second box (22, above [10, 20]) gives a molecule to the first: toward its window
  EXPECT_EQ(System::gibbsTransferTMMCBiasFactor(first, second, 0, 10, 22), 1.0);
  // and may not take one: away from it
  EXPECT_EQ(System::gibbsTransferTMMCBiasFactor(second, first, 0, 22, 12), 0.0);
  // inside its window a step out is refused as before
  EXPECT_EQ(System::gibbsTransferTMMCBiasFactor(first, second, 0, 12, 10), 0.0);
  EXPECT_EQ(System::gibbsTransferTMMCBiasFactor(first, second, 0, 12, 11), 1.0);
}

// Three boxes can not share the rule of the pair
TEST(TMMC_BIAS_GUARDS, gibbs_transfers_of_the_tracked_component_are_reported)
{
  System swapping = makeTMMCBox(3, 0, 10);
  EXPECT_FALSE(swapping.tmmcTracksGibbsTransfers());
  System transferring = makeTMMCBox(3, 0, 10, Move::Types::GibbsSwapCBMC);
  EXPECT_TRUE(transferring.tmmcTracksGibbsTransfers());
  transferring.tmmc.doTMMC = false;
  EXPECT_FALSE(transferring.tmmcTracksGibbsTransfers());
}

// The log-bias of a replica swap is added to the other logarithms and compared as one. Clamped per replica to
// +/- 500, a change of +700 in one replica and -685 in the other gave 0 where the swap has +15.
TEST(TMMC_BIAS_GUARDS, replica_swap_adds_the_unclamped_log_biases)
{
  std::vector<System> boxes;
  boxes.push_back(makeTMMCBox(3, 0, 10));
  boxes.push_back(makeTMMCBox(5, 0, 10));
  const std::optional<double> unbiased = MC_Moves::ParallelTemperingLogAcceptance(boxes[0], boxes[1]);
  ASSERT_TRUE(unbiased.has_value());

  boxes[0].tmmc.bias[5] = 700.0;   // replica 0 goes from 3 to 5 molecules
  boxes[1].tmmc.bias[3] = -685.0;  // replica 1 from 5 to 3
  const std::optional<double> biased = MC_Moves::ParallelTemperingLogAcceptance(boxes[0], boxes[1]);
  ASSERT_TRUE(biased.has_value());
  EXPECT_NEAR(*biased - *unbiased, 15.0, 1e-9);
}

// On an (N, lambda) chain the moves leave the bias of the GC lambda histogram out; the weight of a configuration
// must leave it out too.
TEST(TMMC_BIAS_GUARDS, the_weight_leaves_out_a_lambda_bias_the_chain_does_not_apply)
{
  System system = makeTMMCBox(3, 0, 10, Move::Types::SwapCFCMC);
  Component& component = system.components[0];
  ASSERT_TRUE(component.hasFractionalMolecule);
  ASSERT_GT(component.lambdaGC.numberOfSamplePoints, 1uz);
  std::ranges::fill(component.lambdaGC.biasFactor, 2.0);

  // a plain N chain: the histogram's bias is applied by the moves and removed by the weight
  EXPECT_FALSE(system.tmmcOwnsGCLambda(0));
  EXPECT_NEAR(system.weight(), std::exp(-2.0), 1e-12);

  system.tmmc.numberOfLambdaBins = component.lambdaGC.numberOfSamplePoints;
  system.tmmc.initialize();
  EXPECT_TRUE(system.tmmcOwnsGCLambda(0));
  EXPECT_DOUBLE_EQ(system.weight(), 1.0);
}

// ln Pi belongs to the statistics: dropped with them, while the bias stays.
TEST(TMMC_BIAS_GUARDS, dropping_the_statistics_drops_ln_pi_and_keeps_the_bias)
{
  TransitionMatrix matrix;
  matrix.doTMMC = true;
  matrix.useBias = true;
  matrix.minMacrostate = 0;
  matrix.maxMacrostate = 3;
  matrix.initialize();
  for (std::size_t n = 0; n <= 3; ++n)
  {
    matrix.cmatrix[n] = double3(n > 0 ? 0.2 : 0.0, 0.5, n < 3 ? 0.6 : 0.0);
    matrix.updateHistogram(n);
  }
  matrix.recomputeLnPiAndBias();
  const std::vector<double> bias = matrix.bias;
  ASSERT_NE(matrix.lnpi[0], matrix.lnpi[3]);

  matrix.clearStatisticsKeepBias();
  EXPECT_EQ(matrix.bias, bias);
  for (std::size_t n = 0; n <= 3; ++n)
  {
    EXPECT_EQ(matrix.lnpi[n], 0.0) << "state " << n;
    EXPECT_EQ(matrix.forward_lnpi[n], 0.0) << "state " << n;
    EXPECT_EQ(matrix.reverse_lnpi[n], 0.0) << "state " << n;
  }
}

// The walkers insert at the fugacity of the reference pressure: none given is refused, with the other numbers of
// the run, and a given ideal-gas Rosenbluth weight belongs to one temperature.
TEST(PARALLEL_TMMC_CHECKS, refuses_a_run_it_can_not_sample)
{
  ParallelTMMCParameters parameters;
  parameters.numberOfWindows = 2;

  EXPECT_EQ(refusalOfParallelTMMC(makeMethaneBox(0, 1e5), {300.0}, parameters), "");

  const std::string noPressure = refusalOfParallelTMMC(makeMethaneBox(0, 0.0), {300.0}, parameters);
  EXPECT_NE(noPressure.find("reference pressure"), std::string::npos) << noPressure;

  const std::string noTemperature = refusalOfParallelTMMC(makeMethaneBox(0, 1e5), {300.0, 0.0}, parameters);
  EXPECT_NE(noTemperature.find("temperature"), std::string::npos) << noTemperature;

  ParallelTMMCParameters noBlocks = parameters;
  noBlocks.numberOfBlocks = 0;
  const std::string blocks = refusalOfParallelTMMC(makeMethaneBox(0, 1e5), {300.0}, noBlocks);
  EXPECT_NE(blocks.find("NumberOfBlocks"), std::string::npos) << blocks;

  ParallelTMMCParameters emptyRange = parameters;
  emptyRange.reweightingPressureRange = {0.0, 1e5};
  const std::string range = refusalOfParallelTMMC(makeMethaneBox(0, 1e5), {300.0}, emptyRange);
  EXPECT_NE(range.find("reweighting pressure range"), std::string::npos) << range;

  System weighted = makeMethaneBox(0, 1e5);
  weighted.components[0].idealGasRosenbluthWeight = 0.5;
  EXPECT_EQ(refusalOfParallelTMMC(weighted, {300.0}, parameters), "");
  const std::string ladder = refusalOfParallelTMMC(weighted, {300.0, 350.0}, parameters);
  EXPECT_NE(ladder.find("IdealGasRosenbluthWeight"), std::string::npos) << ladder;
}

// ln Pi(N), the Widom increments and the coexistence densities are those of one volume: a volume move on a bulk box
// (which the input reader lets through, refusing it next to a framework only) is refused by the driver.
TEST(PARALLEL_TMMC_CHECKS, refuses_volume_moves)
{
  ParallelTMMCParameters parameters;
  parameters.numberOfWindows = 2;

  for (const Move::Types move : {Move::Types::VolumeChange, Move::Types::AnisotropicVolumeChange})
  {
    System system = makeMethaneBox(0, 1e5);
    system.mc_moves_probabilities.setProbability(move, 0.1);
    const std::string refusal = refusalOfParallelTMMC(std::move(system), {300.0}, parameters);
    EXPECT_NE(refusal.find("fixed volume"), std::string::npos) << refusal;
  }
}

// An update interval of zero (a system put together in code) is no periodic rebuild, not a remainder by zero.
TEST(TMMC_BIAS_GUARDS, an_update_interval_of_zero_never_rebuilds_the_bias)
{
  TransitionMatrix matrix;
  matrix.doTMMC = true;
  matrix.useBias = true;
  matrix.useTMBias = true;
  matrix.minMacrostate = 0;
  matrix.maxMacrostate = 3;
  matrix.updateTMEvery = 0;
  matrix.initialize();
  for (std::size_t n = 0; n <= 3; ++n)
  {
    matrix.cmatrix[n] = double3(n > 0 ? 0.2 : 0.0, 0.5, n < 3 ? 0.6 : 0.0);
    matrix.updateHistogram(n);
  }
  matrix.numberOfSteps = 10;
  matrix.adjustBias();
  EXPECT_EQ(matrix.numberOfUpdates, 0uz);
  for (std::size_t n = 0; n <= 3; ++n) EXPECT_EQ(matrix.bias[n], 0.0) << "state " << n;
}

// The windows are filled from below. A template that holds more molecules than a window starts at is brought down
// to it; grown upward only, the lower windows started above their range.
TEST(PARALLEL_TMMC_CHECKS, every_walker_starts_at_the_bottom_of_its_window)
{
  ScopedCurrentPath scratch(std::filesystem::temp_directory_path() /
                            std::format("raspa3_test_{}_parallel_tmmc_grow_down", test_support::uniqueSuffix()));

  ParallelTMMCParameters parameters;
  parameters.numberOfWindows = 2;
  parameters.randomSeed = 11;
  ParallelTMMC simulation(makeMethaneBox(9, 1e5), {300.0}, parameters);
  simulation.setup();

  for (const System& walker : simulation.systems)
  {
    EXPECT_EQ(walker.numberOfIntegerMoleculesPerComponent[0], walker.tmmc.minMacrostate);
  }
}
