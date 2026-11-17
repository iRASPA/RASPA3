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
import int3;
import framework;
import property_block_average;
import property_lambda_probability_histogram;
import mc_moves_swap_cfcmc;

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
  ASSERT_TRUE(mc.systems[0].containsTheFractionalMolecule(0));
  ASSERT_FALSE(mc.systems[1].containsTheFractionalMolecule(0));
  ASSERT_EQ(mc.systems[0].moleculeData.size(), 21uz);
  ASSERT_EQ(mc.systems[1].moleculeData.size(), 11uz);
  EXPECT_DOUBLE_EQ(MC_Moves::numberOfScaledMolecules(mc.systems[0]), 21.0);
  EXPECT_DOUBLE_EQ(MC_Moves::numberOfScaledMolecules(mc.systems[1]), 10.0);

  // the count follows the fractional molecule
  std::swap(mc.systems[0].containsTheFractionalMoleculeOfComponent, mc.systems[1].containsTheFractionalMoleculeOfComponent);
  EXPECT_DOUBLE_EQ(MC_Moves::numberOfScaledMolecules(mc.systems[0]), 20.0);
  EXPECT_DOUBLE_EQ(MC_Moves::numberOfScaledMolecules(mc.systems[1]), 11.0);
}

// The weights of a biased run are exp(-bias) with the bias normalized to a minimum of zero: at most one. A block
// whose weights added up to less than one was divided by one (max(1, weight)), and a block was only used when its
// weight exceeded half of max(1, weight of the first): the averages were pulled to zero and the blocks dropped.
TEST(TMMC_CHAIN_CEILING_AND_LINKS, block_averages_divide_by_the_weight_they_accumulated)
{
  BlockAverage<double> average(5);
  for (std::size_t block = 0; block < 5uz; ++block)
  {
    average.addSample(block, 10.0 + static_cast<double>(block), 1e-3);
    average.addSample(block, 20.0 + static_cast<double>(block), 3e-3);
  }
  EXPECT_NEAR(average.averaged(0), 17.5, 1e-12);
  EXPECT_NEAR(average.averaged(4), 21.5, 1e-12);
  EXPECT_NEAR(average.averaged(), 19.5, 1e-12);
  for (std::size_t block = 0; block < 5uz; ++block) EXPECT_TRUE(average.filled(block));
  const auto [mean, error] = average.statistics([](double value) { return value; });
  EXPECT_NEAR(mean, 19.5, 1e-12);
  EXPECT_GT(error, 0.0);

  // the scale of the weights is immaterial, and an empty block stays out
  BlockAverage<double> scaled(5);
  for (std::size_t block = 0; block < 4uz; ++block)
  {
    scaled.addSample(block, 10.0 + static_cast<double>(block), 1e-80);
    scaled.addSample(block, 20.0 + static_cast<double>(block), 3e-80);
  }
  EXPECT_NEAR(scaled.averaged(0), 17.5, 1e-9);
  EXPECT_NEAR(scaled.averaged(), 19.0, 1e-9);
  EXPECT_FALSE(scaled.filled(4));
  EXPECT_EQ(scaled.averaged(4), 0.0);
}

// A component whose only fractional molecule belongs to a pair, a group or a conventional Gibbs CFCMC move has no
// GC lambda coordinate: Wang-Landau lowered bin 0 of its idle GC histogram on every step.
TEST(TMMC_CHAIN_CEILING_AND_LINKS, wang_landau_leaves_an_idle_gc_histogram_alone)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component swapping = Component::makeMethane(forceField, 0);
  swapping.mc_moves_probabilities.setProbability(Move::Types::SwapCFCMC, 1.0);
  Component gibbs = Component::makeMethane(forceField, 1);
  gibbs.name = "gibbs";
  gibbs.mc_moves_probabilities.setProbability(Move::Types::GibbsConventionalCFCMC, 1.0);
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {swapping, gibbs}, {},
                {0, 4}, 5);

  ASSERT_TRUE(system.components[0].hasFractionalMolecule);
  ASSERT_TRUE(system.components[1].hasFractionalMolecule);
  EXPECT_TRUE(system.usesGCLambda(0));
  EXPECT_FALSE(system.usesGCLambda(1));
  EXPECT_TRUE(system.lambdaWangLandauIsActive(0));
  EXPECT_FALSE(system.lambdaWangLandauIsActive(1));

  // the idle histogram does not enter the weight either
  system.components[1].lambdaGC.biasFactor.assign(system.components[1].lambdaGC.biasFactor.size(), 700.0);
  EXPECT_DOUBLE_EQ(system.weight(), 1.0);
}

// One shift for every histogram of every system left a shallow histogram hundreds above a deep one, its weight at
// zero. Each histogram takes its own minimum; the GC histograms of a component that the serial Gibbs swaps carry
// between the boxes keep their difference, which the box switch compares.
TEST(TMMC_CHAIN_CEILING_AND_LINKS, lambda_biases_of_several_systems_are_normalized_one_by_one)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component methane = Component::makeMethane(forceField, 0);
  methane.mc_moves_probabilities.setProbability(Move::Types::SwapCFCMC, 1.0);
  auto makeSystem = [&]()
  { return System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {methane}, {}, {4}, 5); };

  std::vector<System> systems{makeSystem(), makeSystem()};
  auto fill = [](PropertyLambdaProbabilityHistogram& histogram, double low, double high)
  {
    histogram.biasFactor.assign(histogram.biasFactor.size(), high);
    histogram.biasFactor.front() = low;
  };
  fill(systems[0].components[0].lambdaGC, -1000.0, -990.0);
  fill(systems[1].components[0].lambdaGC, -20.0, -15.0);
  fill(systems[1].components[0].lambdaGibbs, -300.0, -299.0);
  System::normalizeLambdaBiasFactors(systems);
  EXPECT_DOUBLE_EQ(systems[0].components[0].lambdaGC.biasFactor.front(), 0.0);
  EXPECT_DOUBLE_EQ(systems[0].components[0].lambdaGC.biasFactor.back(), 10.0);
  EXPECT_DOUBLE_EQ(systems[1].components[0].lambdaGC.biasFactor.front(), 0.0);
  EXPECT_DOUBLE_EQ(systems[1].components[0].lambdaGC.biasFactor.back(), 5.0);
  EXPECT_DOUBLE_EQ(systems[1].components[0].lambdaGibbs.biasFactor.front(), 0.0);

  // coupled by a serial Gibbs CFCMC swap: one shift for the two boxes
  for (System& system : systems) system.numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC[0] = 1;
  fill(systems[0].components[0].lambdaGC, -40.0, -30.0);
  fill(systems[1].components[0].lambdaGC, -20.0, -15.0);
  System::normalizeLambdaBiasFactors(systems);
  EXPECT_DOUBLE_EQ(systems[0].components[0].lambdaGC.biasFactor.front(), 0.0);
  EXPECT_DOUBLE_EQ(systems[1].components[0].lambdaGC.biasFactor.front(), 20.0);
  EXPECT_DOUBLE_EQ(systems[1].components[0].lambdaGC.biasFactor.back(), 25.0);
}

// A group that holds one molecule of the tracked component steps the matrix by one whichever member it is: the
// central molecule or a single satellite (a cation driving a 1:1 group with the adsorbate was refused, the pair
// swap of the same two accepted).
TEST(TMMC_CHAIN_CEILING_AND_LINKS, group_swap_steps_the_matrix_through_a_single_tracked_satellite)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component cation = Component::makeMethane(forceField, 0);
  cation.name = "cation";
  cation.type = Component::Type::Cation;
  Component adsorbate = Component::makeMethane(forceField, 1);
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {cation, adsorbate}, {},
                {0, 0}, 5);
  system.tmmc.doTMMC = true;
  ASSERT_EQ(system.tmmcComponent(), 1uz);

  system.components[0].groupComponentIds = {1};
  EXPECT_TRUE(system.groupSwapFitsTMMC(0));
  ASSERT_TRUE(system.groupSwapTMMCMember(0).has_value());
  EXPECT_EQ(system.groupSwapTMMCMember(0).value(), 1uz);

  // two tracked satellites move the macrostate by two
  system.components[0].groupComponentIds = {1, 1};
  EXPECT_FALSE(system.groupSwapFitsTMMC(0));
  EXPECT_FALSE(system.groupSwapTMMCMember(0).has_value());

  // a group without a tracked molecule is no step of the matrix, and fits it
  system.components[0].groupComponentIds = {0};
  EXPECT_TRUE(system.groupSwapFitsTMMC(0));
  EXPECT_FALSE(system.groupSwapTMMCMember(0).has_value());

  // the tracked central molecule, as before
  system.components[1].groupComponentIds = {0};
  EXPECT_TRUE(system.groupSwapFitsTMMC(1));
  EXPECT_EQ(system.groupSwapTMMCMember(1).value(), 1uz);
}

// The framework mass counts the cations of a fixed number. A cation that is inserted and deleted made the mass,
// and the loadings per kg of framework, depend on the count the run ended on.
TEST(TMMC_CHAIN_CEILING_AND_LINKS, framework_mass_leaves_out_swappable_cations)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  const Framework framework = Framework::makeFAU(forceField, int3(1, 1, 1));
  Component cation = Component::makeMethane(forceField, 0);
  cation.name = "cation";
  cation.type = Component::Type::Cation;
  System system(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {framework}, {cation}, {}, {3}, 5);
  ASSERT_TRUE(system.framework.has_value());
  const double bare = system.framework->mass;

  system.components[0].swappable = false;
  EXPECT_NEAR(system.frameworkMass().value(), bare + 3.0 * system.components[0].totalMass, 1e-9);
  system.components[0].swappable = true;
  EXPECT_NEAR(system.frameworkMass().value(), bare, 1e-9);
}

// A plain N chain under a lambda of several bins (a system put together in code: the input reader gives such a
// run the (N, lambda) chain): the endpoint crossings were recorded as the N -> N +/- 1 transitions, their ratio
// carrying the Wang-Landau factors of the hops in between. The move is refused and records nothing.
TEST(TMMC_CHAIN_CEILING_AND_LINKS, cfcmc_swap_is_refused_on_an_n_chain_with_several_lambda_bins)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component methane = Component::makeMethane(forceField, 0);
  methane.mc_moves_probabilities.setProbability(Move::Types::SwapCFCMC, 1.0);
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {methane}, {}, {4}, 5);
  ASSERT_GT(system.components[0].lambdaGC.numberOfSamplePoints, 1uz);
  system.tmmc.doTMMC = true;
  system.tmmc.maxMacrostate = 20;
  system.tmmc.initialize();
  ASSERT_FALSE(system.tmmc.lambdaChain());

  RandomNumber random(12);
  const std::size_t binBefore = system.components[0].lambdaGC.currentBin;
  for (std::size_t step = 0; step < 50uz; ++step)
  {
    const auto [energy, Pacc] = MC_Moves::swapMove_CFCMC(random, system, 0);
    EXPECT_FALSE(energy.has_value());
    EXPECT_EQ(Pacc.x + Pacc.y + Pacc.z, 0.0);
  }
  EXPECT_EQ(system.components[0].lambdaGC.currentBin, binBefore);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 4uz);

  // without the matrix the move runs
  system.tmmc.doTMMC = false;
  bool moved = false;
  for (std::size_t step = 0; step < 200uz && !moved; ++step)
  {
    moved = MC_Moves::swapMove_CFCMC(random, system, 0).first.has_value();
  }
  EXPECT_TRUE(moved);
}
