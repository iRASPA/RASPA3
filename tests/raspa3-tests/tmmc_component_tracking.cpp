#include <gtest/gtest.h>

import std;

import double3;
import randomnumbers;
import atom;
import forcefield;
import component;
import system;
import simulationbox;
import monte_carlo_transition_matrix;
import simulation_schedule;
import mc_moves_move_types;
import mc_moves_probabilities;

// The transition matrix is of the molecule count of one component, the adsorbate. A cation beside it is a
// component too, picked for its own moves: those moves must neither see the adsorbate's bias (a cation count
// outside the adsorbate's window would be refused as out of bound) nor write into the collection matrix.
// A group swap whose group holds more than one molecule of the tracked component moves the macrostate by
// more than one and does not fit the nearest-neighbour matrix.

namespace
{

ForceField makeForceField()
{
  return ForceField({{"CH4", false, 16.04246, 0.0, 0.0, 6, false}}, {{158.5, 3.72}},
                    ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
}

Component makeMethane(const ForceField& forceField, const MCMoveProbabilities& probabilities,
                      std::size_t lambdaBins = 21)
{
  return Component(forceField, "methane", 190.564, 4599200, 0.01142, {Atom({0, 0, 0}, 0.0, 1.0, 0, 0, 0, false, false)},
                   {}, {}, 5, lambdaBins, probabilities, std::nullopt, false);
}

Component makeCation(const ForceField& forceField, const MCMoveProbabilities& probabilities)
{
  Component cation(forceField, "cation", 190.564, 4599200, 0.01142, {Atom({0, 0, 0}, 0.0, 1.0, 0, 0, 1, false, false)},
                   {}, {}, 5, 21, probabilities, std::nullopt, false);
  cation.type = Component::Type::Cation;
  return cation;
}

MCMoveProbabilities translating()
{
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  return probabilities;
}

MCMoveProbabilities swapping()
{
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Swap, 1.0);
  return probabilities;
}

// the cation first, so that "the first component" and "the adsorbate" differ
System makeCationAdsorbateSystem(const MCMoveProbabilities& cationMoves, const MCMoveProbabilities& adsorbateMoves,
                                 std::size_t cations, std::size_t adsorbates, double pressure)
{
  const ForceField forceField = makeForceField();
  return System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, pressure, 1.0, {},
                {makeCation(forceField, cationMoves), makeMethane(forceField, adsorbateMoves)}, {},
                {cations, adsorbates}, 5);
}

void enableTMMC(System& system, std::size_t minMacrostate, std::size_t maxMacrostate,
                std::size_t numberOfLambdaBins = 1)
{
  system.tmmc.doTMMC = true;
  system.tmmc.useBias = true;
  system.tmmc.useTMBias = true;
  system.tmmc.rejectOutOfBound = true;
  system.tmmc.minMacrostate = minMacrostate;
  system.tmmc.maxMacrostate = maxMacrostate;
  system.tmmc.numberOfLambdaBins = numberOfLambdaBins;
  system.tmmc.initialize();
}

// one adsorbate alone, with the given moves
System makeAdsorbateSystem(const MCMoveProbabilities& moves, std::size_t adsorbates, double pressure)
{
  const ForceField forceField = makeForceField();
  return System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, pressure, 1.0, {},
                {makeMethane(forceField, moves)}, {}, {adsorbates}, 5);
}

bool allZero(const std::vector<double3>& matrix)
{
  return std::ranges::all_of(matrix, [](const double3& row) { return row.x == 0.0 && row.y == 0.0 && row.z == 0.0; });
}

}  // namespace

TEST(TMMC_COMPONENT_TRACKING, the_matrix_is_of_the_first_adsorbate_and_tracks_only_that_component)
{
  System system = makeCationAdsorbateSystem(translating(), translating(), 4, 2, 1e4);
  EXPECT_EQ(system.adsorbateComponentCount(), 1uz);
  EXPECT_EQ(system.tmmcComponent(), 1uz);

  // nothing is tracked while TMMC is off
  EXPECT_FALSE(system.tmmcTracks(0));
  EXPECT_FALSE(system.tmmcTracks(1));

  enableTMMC(system, 0, 10);
  EXPECT_FALSE(system.tmmcTracks(0));
  EXPECT_TRUE(system.tmmcTracks(1));
}

TEST(TMMC_COMPONENT_TRACKING, an_untracked_component_sees_no_bias_even_outside_the_window)
{
  System system = makeCationAdsorbateSystem(translating(), translating(), 9, 0, 1e4);
  enableTMMC(system, 0, 2);

  // the raw matrix refuses a 9 -> 10 transition as out of bound; a cation move is not its business
  EXPECT_DOUBLE_EQ(system.tmmc.biasFactor(10, 9), 0.0);
  EXPECT_DOUBLE_EQ(system.tmmcBiasFactor(0, 10, 9), 1.0);
  EXPECT_DOUBLE_EQ(system.tmmcBiasFactor(0, 8, 9), 1.0);

  // the tracked component gets the matrix's own answer
  EXPECT_DOUBLE_EQ(system.tmmcBiasFactor(1, 1, 0), system.tmmc.biasFactor(1, 0));
  EXPECT_DOUBLE_EQ(system.tmmcBiasFactor(1, 3, 2), 0.0);
}

TEST(TMMC_COMPONENT_TRACKING, a_group_with_tracked_satellites_does_not_fit_the_matrix)
{
  System system = makeCationAdsorbateSystem(translating(), translating(), 4, 2, 1e4);
  enableTMMC(system, 0, 10);
  Component& cation = system.components[0];
  Component& adsorbate = system.components[1];
  cation.maximumGroupDistance = 5.0;
  adsorbate.maximumGroupDistance = 5.0;

  // an adsorbate with two adsorbate satellites: N -> N +/- 3
  adsorbate.groupComponentIds = {1, 1};
  EXPECT_EQ(system.groupMoleculeCountChange(1, 1), 3uz);
  EXPECT_FALSE(system.groupSwapFitsTMMC(1));

  // an adsorbate with cation satellites only: N -> N +/- 1, the ordinary swap of the tracked count
  adsorbate.groupComponentIds = {0, 0};
  EXPECT_EQ(system.groupMoleculeCountChange(1, 1), 1uz);
  EXPECT_TRUE(system.groupSwapFitsTMMC(1));

  // a cation with an adsorbate satellite: the tracked count changes on a move the matrix does not see
  cation.groupComponentIds = {1};
  EXPECT_EQ(system.groupMoleculeCountChange(0, 1), 1uz);
  EXPECT_FALSE(system.groupSwapFitsTMMC(0));

  // a cation with cation satellites: the tracked count is untouched
  cation.groupComponentIds = {0, 0};
  EXPECT_EQ(system.groupMoleculeCountChange(0, 1), 0uz);
  EXPECT_TRUE(system.groupSwapFitsTMMC(0));

  // with TMMC off every group fits
  system.tmmc.doTMMC = false;
  adsorbate.groupComponentIds = {1, 1};
  EXPECT_TRUE(system.groupSwapFitsTMMC(1));
}

// A pair swap changes the count of its partner too: driven by a cation with the adsorbate as partner, it moves
// the tracked count on a move the matrix neither biases nor records, so it does not fit; driven by the adsorbate
// it is the ordinary step.
TEST(TMMC_COMPONENT_TRACKING, a_pair_swap_with_the_tracked_partner_does_not_fit_the_matrix)
{
  System system = makeCationAdsorbateSystem(translating(), translating(), 4, 2, 1e4);
  enableTMMC(system, 0, 10);
  system.components[0].pairComponentId = 1;
  EXPECT_FALSE(system.pairSwapFitsTMMC(0));
  system.components[0].pairComponentId = std::nullopt;
  system.components[1].pairComponentId = 0;
  EXPECT_TRUE(system.pairSwapFitsTMMC(1));
  EXPECT_TRUE(system.pairSwapFitsTMMC(0));
  system.tmmc.doTMMC = false;
  system.components[0].pairComponentId = 1;
  EXPECT_TRUE(system.pairSwapFitsTMMC(0));
}

// The CFCMC pair and group swaps walk their lambda inside one macrostate of the N chain, biased by their own
// Wang-Landau factors: only as a direct crossing, with a single bin, are they steps of the matrix.
TEST(TMMC_COMPONENT_TRACKING, cfcmc_pair_and_group_swaps_are_steps_of_the_matrix_only_as_direct_crossings)
{
  System walking = makeAdsorbateSystem(translating(), 10, 1e4);
  enableTMMC(walking, 0, 40, 1);
  ASSERT_FALSE(walking.tmmc.lambdaChain());
  ASSERT_EQ(walking.components[0].lambdaPairSwap.numberOfSamplePoints, 21uz);
  EXPECT_FALSE(walking.cfcmcSwapFitsTMMC(0, walking.components[0].lambdaPairSwap));
  EXPECT_FALSE(walking.cfcmcSwapFitsTMMC(0, walking.components[0].lambdaGroupSwap));
  walking.tmmc.doTMMC = false;
  EXPECT_TRUE(walking.cfcmcSwapFitsTMMC(0, walking.components[0].lambdaPairSwap));

  const ForceField forceField = makeForceField();
  System crossing = System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {},
                           {makeMethane(forceField, translating(), 1)}, {}, {10}, 5);
  enableTMMC(crossing, 0, 40, 1);
  EXPECT_TRUE(crossing.cfcmcSwapFitsTMMC(0, crossing.components[0].lambdaPairSwap));
  EXPECT_TRUE(crossing.cfcmcSwapFitsTMMC(0, crossing.components[0].lambdaGroupSwapCB));

  // a cation is never the matrix's business
  System mixed = makeCationAdsorbateSystem(translating(), translating(), 4, 2, 1e4);
  enableTMMC(mixed, 0, 10);
  EXPECT_TRUE(mixed.cfcmcSwapFitsTMMC(0, mixed.components[0].lambdaPairSwap));
}

// Nine cations beside an adsorbate whose window is 0--2: the cations swap at a pressure where insertions are
// accepted, so their count must move away from nine (before, every cation insertion was refused as out of the
// adsorbate's window and every deletion was biased to zero), while the collection matrix, which is of the
// adsorbate that never moves, stays empty and the histogram only ever sees the adsorbate at zero.
TEST(TMMC_COMPONENT_TRACKING, cation_swaps_are_neither_refused_by_the_window_nor_written_into_the_matrix)
{
  std::vector<System> systems{makeCationAdsorbateSystem(swapping(), translating(), 9, 0, 1e7)};
  enableTMMC(systems[0], 0, 2);
  ASSERT_EQ(systems[0].numberOfIntegerMoleculesPerComponent[0], 9uz);

  RandomNumber random(5);
  SimulationSchedule schedule{0, 0, 3, 0, 1000, 10000, 5000, 5000};
  MonteCarloTransitionMatrix tmmc(schedule, systems, random, 5);
  for (System& walker : tmmc.systems) walker.runningEnergies = walker.computeTotalEnergies();
  tmmc.simulationStage = MonteCarloTransitionMatrix::SimulationStage::Initialization;
  for (std::size_t cycle = 0; cycle < 3; ++cycle) tmmc.performCycle();

  const System& system = tmmc.systems[0];
  EXPECT_NE(system.numberOfIntegerMoleculesPerComponent[0], 9uz) << "no cation swap was accepted";
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 0uz);
  EXPECT_TRUE(allZero(system.tmmc.cmatrix)) << "cation swaps were recorded into the adsorbate's matrix";

  const std::vector<std::size_t>& histogram = system.tmmc.histogram;
  ASSERT_EQ(histogram.size(), 3uz);
  EXPECT_GT(histogram[0], 0uz);
  EXPECT_EQ(histogram[1] + histogram[2], 0uz);
}

// The integer swaps (N -> N +/- 1 with lambda unchanged) are steps of a plain N chain but not of the (N, lambda)
// chain of the CFCMC swap: there they would be recorded on the (N, k) row as lambda steps and biased by the
// endpoint crossing.
TEST(TMMC_COMPONENT_TRACKING, integer_swaps_are_steps_of_a_plain_chain_but_not_of_a_lambda_chain)
{
  System plain = makeAdsorbateSystem(swapping(), 10, 1e4);
  enableTMMC(plain, 0, 40, 1);
  EXPECT_TRUE(plain.tmmcTracksIntegerSwap(0));
  EXPECT_DOUBLE_EQ(plain.tmmcIntegerSwapBiasFactor(0, 11, 10), plain.tmmc.biasFactor(11, 10));

  System chain = makeAdsorbateSystem(swapping(), 10, 1e4);
  enableTMMC(chain, 0, 40, 21);
  ASSERT_TRUE(chain.tmmc.lambdaChain());
  EXPECT_TRUE(chain.tmmcTracks(0));
  EXPECT_FALSE(chain.tmmcTracksIntegerSwap(0));
  EXPECT_DOUBLE_EQ(chain.tmmcIntegerSwapBiasFactor(0, 11, 10), 0.0);
  EXPECT_DOUBLE_EQ(chain.tmmcIntegerSwapBiasFactor(0, 9, 10), 0.0);
  // the CFCMC endpoint crossing keeps the chain's own bias
  EXPECT_DOUBLE_EQ(chain.tmmcBiasFactor(0, 11, 10), chain.tmmc.biasFactor(11, 10));

  // a cation is never the matrix's business, chain or not
  System mixed = makeCationAdsorbateSystem(swapping(), translating(), 9, 0, 1e4);
  enableTMMC(mixed, 0, 2, 21);
  EXPECT_FALSE(mixed.tmmcTracksIntegerSwap(0));
  EXPECT_DOUBLE_EQ(mixed.tmmcIntegerSwapBiasFactor(0, 10, 9), 1.0);
}

// Ten molecules swapping by the plain insertion/deletion move while the matrix is an (N, lambda) chain: every
// such swap is refused and nothing is recorded, so the count stays at ten and the matrix stays empty (before,
// the swaps went through, biased by the endpoint crossing, and their acceptances landed on lambda-step entries).
TEST(TMMC_COMPONENT_TRACKING, integer_swaps_on_a_lambda_chain_are_refused_and_not_recorded)
{
  std::vector<System> systems{makeAdsorbateSystem(swapping(), 10, 1e7)};
  enableTMMC(systems[0], 0, 40, 21);

  RandomNumber random(9);
  SimulationSchedule schedule{0, 0, 3, 0, 1000, 10000, 5000, 5000};
  MonteCarloTransitionMatrix tmmc(schedule, systems, random, 5);
  for (System& walker : tmmc.systems) walker.runningEnergies = walker.computeTotalEnergies();
  tmmc.simulationStage = MonteCarloTransitionMatrix::SimulationStage::Initialization;
  for (std::size_t cycle = 0; cycle < 3; ++cycle) tmmc.performCycle();

  const System& system = tmmc.systems[0];
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 10uz);
  EXPECT_TRUE(allZero(system.tmmc.cmatrix));

  // the same swaps on a plain chain do move the count and do record
  std::vector<System> plainSystems{makeAdsorbateSystem(swapping(), 10, 1e7)};
  enableTMMC(plainSystems[0], 0, 40, 1);
  RandomNumber plainRandom(9);
  MonteCarloTransitionMatrix plain(schedule, plainSystems, plainRandom, 5);
  for (System& walker : plain.systems) walker.runningEnergies = walker.computeTotalEnergies();
  plain.simulationStage = MonteCarloTransitionMatrix::SimulationStage::Initialization;
  for (std::size_t cycle = 0; cycle < 3; ++cycle) plain.performCycle();
  EXPECT_NE(plain.systems[0].numberOfIntegerMoleculesPerComponent[0], 10uz);
  EXPECT_FALSE(allZero(plain.systems[0].tmmc.cmatrix));
}

// The Gibbs CFCMC transfer turns the fractional molecule whole at its lambda bin and makes a molecule of the
// other box fractional at that same bin: on an (N, lambda) chain that is an off-chain jump in both boxes, like
// an integer swap, and is refused and not recorded. On a plain chain the transfers run and record as usual.
TEST(TMMC_COMPONENT_TRACKING, gibbs_cfcmc_transfers_on_a_lambda_chain_are_refused_and_not_recorded)
{
  // the driver writes its output files relative to the working directory
  const std::filesystem::path previous = std::filesystem::current_path();
  const std::filesystem::path scratch =
      std::filesystem::temp_directory_path() /
      std::format("raspa3_tmmc_gibbs_chain_{}", std::chrono::steady_clock::now().time_since_epoch().count());
  std::filesystem::create_directories(scratch);
  std::filesystem::current_path(scratch);

  for (const std::size_t numberOfLambdaBins : {21uz, 1uz})
  {
    MCMoveProbabilities transferring;
    transferring.setProbability(Move::Types::Translation, 1.0);
    transferring.setProbability(Move::Types::GibbsSwapCFCMC, 1.0);
    std::vector<System> systems{makeAdsorbateSystem(transferring, 10, 1e4), makeAdsorbateSystem(transferring, 10, 1e4)};
    for (System& system : systems) enableTMMC(system, 0, 40, numberOfLambdaBins);
    const bool chain = numberOfLambdaBins > 1uz;
    ASSERT_EQ(systems[0].tmmc.lambdaChain(), chain);

    RandomNumber random(13);
    SimulationSchedule schedule{0, 0, 0, 4, 1000, 10000, 5000, 5000};
    MonteCarloTransitionMatrix tmmc(schedule, systems, random, 5);
    tmmc.setup();
    for (System& system : tmmc.systems) system.runningEnergies = system.computeTotalEnergies();
    tmmc.simulationStage = MonteCarloTransitionMatrix::SimulationStage::Equilibration;
    for (std::size_t cycle = 0; cycle < 4; ++cycle) tmmc.performCycle();

    const std::size_t nA = tmmc.systems[0].numberOfIntegerMoleculesPerComponent[0];
    const std::size_t nB = tmmc.systems[1].numberOfIntegerMoleculesPerComponent[0];
    const bool recorded = !allZero(tmmc.systems[0].tmmc.cmatrix) || !allZero(tmmc.systems[1].tmmc.cmatrix);
    if (chain)
    {
      // the transfers are refused and nothing is recorded (the fractional molecule, which the setup adds on
      // top of the ten integer ones, only moves in lambda)
      EXPECT_EQ(nA, 10uz) << "on the chain";
      EXPECT_EQ(nB, 10uz) << "on the chain";
      EXPECT_FALSE(recorded) << "on the chain";
    }
    else
    {
      EXPECT_TRUE(recorded) << "on the plain chain";
    }
  }

  std::filesystem::current_path(previous);
  std::filesystem::remove_all(scratch);
}
