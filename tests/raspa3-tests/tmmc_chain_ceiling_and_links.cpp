#include <gtest/gtest.h>

import std;

import double3;
import forcefield;
import simulationbox;
import atom;
import component;
import system;
import randomnumbers;
import transition_matrix;
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_volume;
import monte_carlo;
import simulation_schedule;

// Regression tests for the ceiling of the (N, lambda) chain, the unsupported links of ln Pi, the pair swap driven
// by the partner of the TMMC component and the volume factor of a switched-off Gibbs placeholder.

namespace
{

TransitionMatrix makeLambdaChain(std::size_t minN, std::size_t maxN, std::size_t lambdaBins)
{
  TransitionMatrix matrix;
  matrix.doTMMC = true;
  matrix.useBias = true;
  matrix.useTMBias = true;
  matrix.rejectOutOfBound = true;
  matrix.minMacrostate = minN;
  matrix.maxMacrostate = maxN;
  matrix.numberOfLambdaBins = lambdaBins;
  matrix.initialize();
  return matrix;
}

}  // namespace

// The states (maxN, k > 0) hold maxN molecules and a partly coupled one more: they are past the ceiling. At a
// ceiling set by the packing they cannot be visited, and a Wang-Landau histogram that counted them never went flat.
TEST(TMMC_CHAIN_CEILING_AND_LINKS, lambda_chain_ends_at_the_integer_ceiling)
{
  TransitionMatrix matrix = makeLambdaChain(2, 4, 5);
  EXPECT_EQ(matrix.numberOfChainStates(), 11uz);
  EXPECT_EQ(matrix.cmatrix.size(), 11uz);
  EXPECT_EQ(matrix.wangLandauHistogram.size(), 11uz);
  EXPECT_EQ(matrix.chainIndex(4, 0), 10uz);
  EXPECT_EQ(matrix.chainIndex(4, 1), TransitionMatrix::outsideChain);
  EXPECT_EQ(matrix.chainIndex(5, 0), TransitionMatrix::outsideChain);
  EXPECT_EQ(matrix.chainIndex(1, 4), TransitionMatrix::outsideChain);

  // the lambda hop up from the ceiling is refused, the one below it is not
  EXPECT_EQ(matrix.biasFactor(4, 4, 1, 0), 0.0);
  EXPECT_EQ(matrix.biasFactor(3, 3, 4, 3), 1.0);
  // the crossing (3, last) -> (4, 0) is a step, the insertion from the ceiling is not
  EXPECT_EQ(matrix.biasFactor(4, 3), 1.0);
  EXPECT_EQ(matrix.biasFactor(5, 4), 0.0);

  // a visit past the ceiling is not counted
  matrix.updateHistogram(4, 2);
  matrix.visitWangLandau(4, 2);
  EXPECT_EQ(std::accumulate(matrix.histogram.begin(), matrix.histogram.end(), 0uz), 0uz);

  // every state of the chain can be visited: Wang-Landau anneals once all of them have been
  matrix.useWangLandau = true;
  matrix.wangLandauCheckEvery = 11uz;
  const double before = matrix.wangLandauFactor;
  for (std::size_t N = 2; N <= 4; ++N)
  {
    for (std::size_t k = 0; k < 5; ++k)
    {
      if (matrix.chainIndex(N, k) != TransitionMatrix::outsideChain) matrix.visitWangLandau(N, k);
    }
  }
  EXPECT_LT(matrix.wangLandauFactor, before);
}

// A walk that starts outside the chain is let back in by steps toward it (it was refused every step, and stayed).
TEST(TMMC_CHAIN_CEILING_AND_LINKS, a_walk_outside_the_chain_can_step_back_into_it)
{
  TransitionMatrix chain = makeLambdaChain(2, 4, 5);
  // from (4, 3), past the ceiling: down in lambda is allowed, up is not
  EXPECT_EQ(chain.biasFactor(4, 4, 2, 3), 1.0);
  EXPECT_EQ(chain.biasFactor(4, 4, 4, 3), 0.0);

  TransitionMatrix plain = makeLambdaChain(2, 4, 1);
  EXPECT_EQ(plain.biasFactor(5, 6), 1.0);
  EXPECT_EQ(plain.biasFactor(7, 6), 0.0);
  EXPECT_EQ(plain.biasFactor(1, 0), 1.0);
  EXPECT_EQ(plain.biasFactor(0, 1), 0.0);
  // inside the chain nothing changed
  EXPECT_EQ(plain.biasFactor(5, 4), 0.0);
  EXPECT_EQ(plain.biasFactor(1, 2), 0.0);
}

// Across a link without accepted transitions in both directions ln Pi is carried flat (for the bias); that step
// is not an estimate, and the reweighted loading keeps to the longest run of supported links.
TEST(TMMC_CHAIN_CEILING_AND_LINKS, reweighted_loading_keeps_to_the_supported_run)
{
  TransitionMatrix matrix = makeLambdaChain(0, 5, 1);
  matrix.histogram = {3, 3, 3, 3, 3, 3};
  // links 0-1 and 1-2 supported, 2-3 has no accepted deletion, 3-4 and 4-5 supported... but 4-5 has no insertion
  matrix.cmatrix = {double3(0.0, 0.5, 0.5), double3(0.5, 0.0, 0.5), double3(0.5, 0.0, 0.5),
                    double3(0.0, 0.5, 0.5), double3(0.5, 0.5, 0.0), double3(0.0, 1.0, 0.0)};
  EXPECT_TRUE(matrix.linkSupported(0));
  EXPECT_TRUE(matrix.linkSupported(1));
  EXPECT_FALSE(matrix.linkSupported(2));
  EXPECT_TRUE(matrix.linkSupported(3));
  EXPECT_FALSE(matrix.linkSupported(4));
  EXPECT_FALSE(matrix.linkSupported(5));

  const auto visited = matrix.visitedMoleculeRange();
  const auto supported = matrix.supportedMoleculeRange();
  ASSERT_TRUE(visited.has_value());
  ASSERT_TRUE(supported.has_value());
  EXPECT_EQ(visited.value(), (std::pair<std::size_t, std::size_t>{0uz, 5uz}));
  EXPECT_EQ(supported.value(), (std::pair<std::size_t, std::size_t>{0uz, 2uz}));

  // the loading is that of N = 0..2 (equal weights from the symmetric rows), not of the plateau above it
  matrix.statisticsFileName =
      (std::filesystem::temp_directory_path() /
       std::format("raspa3_tmmc_links_{}.txt", std::chrono::steady_clock::now().time_since_epoch().count()))
          .string();
  matrix.recomputeLnPiAndBias();
  EXPECT_NEAR(matrix.lnpi[3], matrix.lnpi[2], 1e-12);
  const auto moments = matrix.reweightedMoleculeCountMeanAndVariance();
  ASSERT_TRUE(moments.has_value());
  EXPECT_NEAR(moments->first, 1.0, 1e-12);
  EXPECT_NEAR(moments->second, 2.0 / 3.0, 1e-12);

  // the statistics file numbers the runs
  std::ifstream file(matrix.statisticsFileName);
  std::string line;
  std::vector<int> segments;
  bool header = false;
  while (std::getline(file, line))
  {
    if (line.empty() || line[0] == '#') continue;
    if (!header)
    {
      EXPECT_TRUE(line.ends_with("segment"));
      header = true;
      continue;
    }
    segments.push_back(std::stoi(line.substr(line.find_last_of(' ') + 1)));
  }
  EXPECT_EQ(segments, (std::vector<int>{0, 0, 0, 1, 1, 2}));
  file.close();
  std::filesystem::remove(matrix.statisticsFileName);
}

// On an (N, lambda) chain the step N -> N + 1 is supported only when every lambda link on the way is.
TEST(TMMC_CHAIN_CEILING_AND_LINKS, supported_range_of_a_lambda_chain_needs_every_lambda_link)
{
  TransitionMatrix matrix = makeLambdaChain(0, 2, 3);
  ASSERT_EQ(matrix.numberOfChainStates(), 7uz);
  std::fill(matrix.histogram.begin(), matrix.histogram.end(), 1uz);
  std::fill(matrix.cmatrix.begin(), matrix.cmatrix.end(), double3(0.5, 0.0, 0.5));
  EXPECT_EQ(matrix.supportedMoleculeRange().value(), (std::pair<std::size_t, std::size_t>{0uz, 2uz}));
  // a starved hop (1, 1) -> (1, 2) cuts N = 2 off
  matrix.cmatrix[matrix.chainIndex(1, 1)].z = 0.0;
  EXPECT_EQ(matrix.supportedMoleculeRange().value(), (std::pair<std::size_t, std::size_t>{0uz, 1uz}));
}

// The switched-off placeholder of the serial Gibbs CFCMC swap is scaled with its box but carries no volume
// factor: the move that brings the fractional molecule into a box draws its position anew and weighs the boxes
// V_B / V_A, one factor for the box that holds it.
TEST(TMMC_CHAIN_CEILING_AND_LINKS, gibbs_placeholder_carries_no_volume_factor)
{
  auto makeBox = [](std::size_t numberOfMolecules)
  {
    const ForceField forceField({{"CH4", false, 16.04246, 0.0, 0.0, 6, false}}, {{158.5, 3.72}},
                                ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
    MCMoveProbabilities probabilities;
    probabilities.setProbability(Move::Types::Translation, 1.0);
    probabilities.setProbability(Move::Types::GibbsSwapCFCMC, 1.0);
    Component methane(forceField, "methane", 190.564, 4599200, 0.01142,
                      {Atom({0, 0, 0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21, probabilities, std::nullopt,
                      false);
    return System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {methane}, {},
                  {numberOfMolecules}, 5);
  };
  std::vector<System> boxes{makeBox(20), makeBox(10)};
  MonteCarlo mc(SimulationSchedule{1, 0, 0, 1, 1000, 10000, 5000, 5000}, boxes, 7uz, 5, false);
  mc.setup();
  ASSERT_TRUE(mc.systems[0].containsTheFractionalMolecule);
  ASSERT_FALSE(mc.systems[1].containsTheFractionalMolecule);
  ASSERT_EQ(mc.systems[0].moleculeData.size(), 21uz);
  ASSERT_EQ(mc.systems[1].moleculeData.size(), 11uz);
  EXPECT_DOUBLE_EQ(MC_Moves::numberOfScaledMolecules(mc.systems[0]), 21.0);
  EXPECT_DOUBLE_EQ(MC_Moves::numberOfScaledMolecules(mc.systems[1]), 10.0);

  // the count follows the fractional molecule
  std::swap(mc.systems[0].containsTheFractionalMolecule, mc.systems[1].containsTheFractionalMolecule);
  EXPECT_DOUBLE_EQ(MC_Moves::numberOfScaledMolecules(mc.systems[0]), 20.0);
  EXPECT_DOUBLE_EQ(MC_Moves::numberOfScaledMolecules(mc.systems[1]), 11.0);
}
