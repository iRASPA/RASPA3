#include <gtest/gtest.h>

import std;

import int3;
import uint3;
import double3;
import archive;
import randomnumbers;
import atom;
import forcefield;
import framework;
import component;
import system;
import simulationbox;
import monte_carlo;
import monte_carlo_transition_matrix;
import transition_matrix;
import simulation_schedule;
import mc_moves_move_types;
import mc_moves_probabilities;
import cbmc_grow_context;

// Driver-level bookkeeping of the Monte Carlo drivers: the per-step sampling with a single system, the
// restart archive, and the first-run preparation of the transition-matrix driver.

namespace
{

// A single Lennard-Jones component with CFCMC swaps (so that lambda sampling is active) in an empty box;
// 'swapMove' selects the CFCMC move (the single-box swap, or the serial Gibbs swap for a pair of boxes).
System makeMethaneCFCMCSystem(std::size_t numberOfMolecules, Move::Types swapMove = Move::Types::SwapCFCMC,
                              std::size_t numberOfLambdaBins = 21uz)
{
  const ForceField forceField({{"CH4", false, 16.04246, 0.0, 0.0, 6, false}}, {{158.5, 3.72}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(swapMove, 1.0);
  Component methane(forceField, "methane", 190.564, 4599200, 0.01142,
                    {Atom({0, 0, 0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, numberOfLambdaBins, probabilities,
                    std::nullopt, false);
  methane.lambdaGC.computeDUdlambda = true;
  return System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {methane}, {},
                {numberOfMolecules}, 5);
}

SimulationSchedule makeSchedule(std::size_t production, std::size_t initialization, std::size_t equilibration)
{
  return {production, 0, initialization, equilibration, 1000, 10000, 5000, 5000};
}

}  // namespace

// With a single system the "second" system of the adjacent pair aliases the first: the lambda bookkeeping of
// that system must be sampled once per step (as it is for every system of a multi-system run), not twice.
TEST(MC_DRIVER_BOOKKEEPING, single_system_is_sampled_once_per_step)
{
  std::vector<System> systems{makeMethaneCFCMCSystem(20)};
  MonteCarlo mc(makeSchedule(1, 0, 1), systems, 7uz, 5, false);
  mc.setup();
  for (System& system : mc.systems) system.runningEnergies = system.computeTotalEnergies();

  // one component: max(number of molecules, 20) * 1 steps per cycle (the count is taken before the cycle,
  // the CFCMC swaps may change it)
  const std::size_t expectedSteps = std::max(mc.systems[0].numberOfMolecules(), 20uz);
  mc.simulationStage = MonteCarlo::SimulationStage::Equilibration;
  mc.performCycle();
  EXPECT_EQ(mc.systems[0].components[0].lambdaGC.occupancyTotal, expectedSteps);

  // the transition-matrix driver takes the same path
  std::vector<System> tmmcSystems{makeMethaneCFCMCSystem(20)};
  RandomNumber random(7);
  MonteCarloTransitionMatrix tmmc(makeSchedule(1, 0, 1), tmmcSystems, random, 5);
  for (System& system : tmmc.systems) system.runningEnergies = system.computeTotalEnergies();
  const std::size_t expectedTMMCSteps = std::max(tmmc.systems[0].numberOfMolecules(), 20uz);
  tmmc.simulationStage = MonteCarloTransitionMatrix::SimulationStage::Equilibration;
  tmmc.performCycle();
  EXPECT_EQ(tmmc.systems[0].components[0].lambdaGC.occupancyTotal, expectedTMMCSteps);
}

// Independent boxes with their own CFCMC swaps each own an active fractional molecule: the setup must leave the
// flag on in every one of them (it was switched off beyond the first box as if they shared one, which froze the
// Wang-Landau learning and the lambda sampling there). Only the serial Gibbs swap shares a fractional molecule,
// and only its later boxes hold a switched-off placeholder.
TEST(MC_DRIVER_BOOKKEEPING, independent_boxes_keep_their_fractional_molecules_active)
{
  std::vector<System> independent{makeMethaneCFCMCSystem(20), makeMethaneCFCMCSystem(10)};
  MonteCarlo mc(makeSchedule(1, 0, 1), independent, 7uz, 5, false);
  mc.setup();
  for (const System& system : mc.systems)
  {
    EXPECT_TRUE(system.containsTheFractionalMolecule(0));
    EXPECT_TRUE(system.lambdaWangLandauIsActive(0));
  }
  RandomNumber random(7);
  std::vector<System> tmmcIndependent{makeMethaneCFCMCSystem(20), makeMethaneCFCMCSystem(10)};
  MonteCarloTransitionMatrix tmmc(makeSchedule(1, 0, 1), tmmcIndependent, random, 5);
  for (const System& system : tmmc.systems) EXPECT_TRUE(system.containsTheFractionalMolecule(0));

  std::vector<System> gibbs{makeMethaneCFCMCSystem(20, Move::Types::GibbsSwapCFCMC),
                            makeMethaneCFCMCSystem(10, Move::Types::GibbsSwapCFCMC)};
  MonteCarlo serial(makeSchedule(1, 0, 1), gibbs, 7uz, 5, false);
  serial.setup();
  EXPECT_TRUE(serial.systems[0].containsTheFractionalMolecule(0));
  EXPECT_FALSE(serial.systems[1].containsTheFractionalMolecule(0));

  // the GC lambda bias reweights the box that holds the active fractional molecule, not the placeholder's box
  for (System& system : serial.systems)
  {
    system.components[0].lambdaGC.biasFactor[system.components[0].lambdaGC.currentBin] = 2.0;
  }
  EXPECT_NEAR(serial.systems[0].weight(), std::exp(-2.0), 1e-12);
  EXPECT_NEAR(serial.systems[1].weight(), 1.0, 1e-12);
}

// The (deprecated) system-level move probabilities are merged into the components before the constructor decides
// which components are swappable and which need a fractional molecule: a CFCMC swap given only on the system
// level must still allocate the fractional molecule and its dU/dlambda group.
TEST(MC_DRIVER_BOOKKEEPING, system_level_move_probabilities_drive_the_fractional_bookkeeping)
{
  const ForceField forceField({{"CH4", false, 16.04246, 0.0, 0.0, 6, false}}, {{158.5, 3.72}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  Component methane(forceField, "methane", 190.564, 4599200, 0.01142,
                    {Atom({0, 0, 0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21, MCMoveProbabilities(),
                    std::nullopt, false);
  methane.lambdaGC.computeDUdlambda = true;

  MCMoveProbabilities systemProbabilities;
  systemProbabilities.setProbability(Move::Types::Translation, 1.0);
  systemProbabilities.setProbability(Move::Types::SwapCFCMC, 1.0);

  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {methane}, {}, {4}, 5,
                systemProbabilities);

  EXPECT_GT(system.components[0].mc_moves_probabilities.getProbability(Move::Types::SwapCFCMC), 0.0);
  EXPECT_TRUE(system.components[0].swappable);
  EXPECT_TRUE(system.components[0].hasFractionalMolecule);
  EXPECT_EQ(system.numberOfGCFractionalMoleculesPerComponent_CFCMC[0], 1uz);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 1uz);
  EXPECT_NE(system.components[0].lambdaGC.dUdlambdaGroupId, 0);
  EXPECT_EQ(system.swappableComponents, std::vector<std::size_t>{0uz});
}

// The transition-matrix driver runs the lambda bookkeeping of the conventional Gibbs CFCMC move (Wang-Landau
// sampling and occupancy) like the plain Monte Carlo driver does.
TEST(MC_DRIVER_BOOKKEEPING, transition_matrix_samples_the_gibbs_lambda_histogram)
{
  const ForceField forceField({{"CH4", false, 16.04246, 0.0, 0.0, 6, false}}, {{158.5, 3.72}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::GibbsConventionalCFCMC, 1.0);
  Component methane(forceField, "methane", 190.564, 4599200, 0.01142,
                    {Atom({0, 0, 0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21, probabilities, std::nullopt,
                    false);
  methane.lambdaGibbs.computeDUdlambda = true;

  std::vector<System> systems{
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {methane}, {}, {10}, 5),
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {methane}, {}, {10}, 5)};
  ASSERT_TRUE(systems[0].usesGibbsConventionalCFCMC());

  // the driver writes its output files relative to the working directory
  const std::filesystem::path previous = std::filesystem::current_path();
  const std::filesystem::path scratch = std::filesystem::temp_directory_path() / "raspa3_tmmc_gibbs_lambda_test";
  std::filesystem::create_directories(scratch);
  std::filesystem::current_path(scratch);
  {
    RandomNumber random(5);
    MonteCarloTransitionMatrix tmmc(makeSchedule(1, 0, 1), systems, random, 5);
    tmmc.setup();
    for (System& system : tmmc.systems) system.runningEnergies = system.computeTotalEnergies();

    // one component: max(total number of molecules, 20) steps per cycle (the count is taken before the cycle)
    const std::size_t steps =
        std::max(tmmc.systems[0].numberOfMolecules() + tmmc.systems[1].numberOfMolecules(), 20uz);

    tmmc.simulationStage = MonteCarloTransitionMatrix::SimulationStage::Equilibration;
    tmmc.performCycle();

    // every step samples the occupancy of both lambda histograms in both systems of the pair
    std::size_t gibbsOccupancy{0uz};
    std::size_t gcOccupancy{0uz};
    for (const System& system : tmmc.systems)
    {
      gibbsOccupancy += system.components[0].lambdaGibbs.occupancyTotal;
      gcOccupancy += system.components[0].lambdaGC.occupancyTotal;
    }
    EXPECT_EQ(gcOccupancy, 2uz * steps);
    EXPECT_EQ(gibbsOccupancy, 2uz * steps);
  }
  std::filesystem::current_path(previous);
  std::filesystem::remove_all(scratch);
}

// A production cycle has a fixed number of steps: a cycle of max(20, N) steps with the live molecule count N
// samples the states of high N less often per step than those of low N, which biases the cycle-sampled
// averages of an open ensemble. The stages before production keep the live count.
TEST(MC_DRIVER_BOOKKEEPING, production_cycle_has_a_fixed_number_of_steps)
{
  // the average of the recorded counts, rounded; at least 20 steps per component
  ProductionCycleLength length{};
  EXPECT_FALSE(length.fixed);
  length.fix(7uz, 2uz);
  EXPECT_TRUE(length.fixed);
  EXPECT_EQ(length.stepsPerCycle, 40uz);
  length.clear();
  for (const std::size_t count : {100uz, 101uz, 103uz}) length.sample(count);
  length.fix(7uz, 2uz);
  EXPECT_EQ(length.stepsPerCycle, 2uz * 101uz);

  std::vector<System> systems{makeMethaneCFCMCSystem(30)};
  MonteCarlo mc(makeSchedule(1, 0, 1), systems, 7uz, 5, false);
  mc.setup();
  for (System& system : mc.systems) system.runningEnergies = system.computeTotalEnergies();
  const auto stepsOfACycle = [&]()
  {
    const std::size_t before = mc.systems[0].components[0].lambdaGC.occupancyTotal;
    mc.performCycle();
    return mc.systems[0].components[0].lambdaGC.occupancyTotal - before;
  };

  // before production the length follows the live molecule count
  mc.simulationStage = MonteCarlo::SimulationStage::Equilibration;
  const std::size_t liveCount = mc.systems[0].numberOfMolecules();
  EXPECT_EQ(stepsOfACycle(), std::max(liveCount, 20uz));

  // in production it is the fixed length, whatever the molecule count
  mc.simulationStage = MonteCarlo::SimulationStage::Production;
  mc.productionCycleLength.clear();
  mc.productionCycleLength.sample(57uz);
  mc.productionCycleLength.fix(mc.systems[0].numberOfMolecules(), 1uz);
  ASSERT_NE(mc.systems[0].numberOfMolecules(), 57uz);
  EXPECT_EQ(stepsOfACycle(), 57uz);
  EXPECT_EQ(stepsOfACycle(), 57uz);

  // a production stage resumed from a restart file without a stored length fixes it at its first cycle
  mc.productionCycleLength.clear();
  const std::size_t resumedCount = mc.systems[0].numberOfMolecules();
  EXPECT_EQ(stepsOfACycle(), std::max(resumedCount, 20uz));
  EXPECT_TRUE(mc.productionCycleLength.fixed);
  EXPECT_EQ(mc.productionCycleLength.stepsPerCycle, std::max(resumedCount, 20uz));

  // the transition-matrix driver takes the same path
  std::vector<System> tmmcSystems{makeMethaneCFCMCSystem(30)};
  RandomNumber random(7);
  MonteCarloTransitionMatrix tmmc(makeSchedule(1, 0, 1), tmmcSystems, random, 5);
  for (System& system : tmmc.systems) system.runningEnergies = system.computeTotalEnergies();
  tmmc.simulationStage = MonteCarloTransitionMatrix::SimulationStage::Production;
  tmmc.productionCycleLength.sample(43uz);
  tmmc.productionCycleLength.fix(tmmc.systems[0].numberOfMolecules(), 1uz);
  const std::size_t before = tmmc.systems[0].components[0].lambdaGC.occupancyTotal;
  tmmc.performCycle();
  EXPECT_EQ(tmmc.systems[0].components[0].lambdaGC.occupancyTotal - before, 43uz);
}

// The fixed length of a production cycle is set by a full run and survives a binary restart.
TEST(MC_DRIVER_BOOKKEEPING, archive_restores_the_production_cycle_length)
{
  std::vector<System> systems{makeMethaneCFCMCSystem(40)};
  MonteCarlo original(makeSchedule(2, 2, 4), systems, 11uz, 5, false);
  original.run();
  ASSERT_TRUE(original.productionCycleLength.fixed);
  // the second half of the equilibration was recorded
  EXPECT_EQ(original.productionCycleLength.numberOfSamples, 2uz);
  EXPECT_GE(original.productionCycleLength.stepsPerCycle, 20uz);

  const std::filesystem::path path =
      std::filesystem::temp_directory_path() / "raspa3_monte_carlo_cycle_length.bin";
  {
    std::ofstream stream(path, std::ios::binary);
    Archive<std::ofstream> archive(stream);
    archive << original;
  }
  MonteCarlo restored;
  {
    std::ifstream stream(path, std::ios::binary);
    Archive<std::ifstream> archive(stream);
    archive >> restored;
  }
  std::filesystem::remove(path);

  EXPECT_TRUE(restored.productionCycleLength.fixed);
  EXPECT_EQ(restored.productionCycleLength.stepsPerCycle, original.productionCycleLength.stepsPerCycle);
  EXPECT_EQ(restored.productionCycleLength.moleculeCountSum, original.productionCycleLength.moleculeCountSum);
  EXPECT_EQ(restored.productionCycleLength.numberOfSamples, original.productionCycleLength.numberOfSamples);
}

// The overall cycle number indexes the property-evolution samples; it has to survive a binary restart.
TEST(MC_DRIVER_BOOKKEEPING, archive_restores_the_absolute_cycle_number)
{
  std::vector<System> systems{makeMethaneCFCMCSystem(10)};
  MonteCarlo original(makeSchedule(4, 3, 2), systems, 11uz, 5, false);
  original.run();
  ASSERT_EQ(original.absoluteCurrentCycle, 9uz);

  const std::filesystem::path path = std::filesystem::temp_directory_path() / "raspa3_monte_carlo_bookkeeping.bin";
  {
    std::ofstream stream(path, std::ios::binary);
    Archive<std::ofstream> archive(stream);
    archive << original;
  }
  MonteCarlo restored;
  {
    std::ifstream stream(path, std::ios::binary);
    Archive<std::ifstream> archive(stream);
    archive >> restored;
  }
  std::filesystem::remove(path);

  EXPECT_EQ(restored.absoluteCurrentCycle, original.absoluteCurrentCycle);
  EXPECT_EQ(restored.currentCycle, original.currentCycle);
  EXPECT_EQ(restored.simulationStage, original.simulationStage);
}

// The transition-matrix driver prepares the systems like the plain Monte Carlo driver: the fractional molecule of a
// Gibbs CFCMC swap is active in the first box only (its copy in the other box must carry groupId 0, otherwise it
// contributes a spurious dU/dlambda), and the requested framework interpolation grids are built.
TEST(MC_DRIVER_BOOKKEEPING, transition_matrix_setup_prepares_fractional_molecules_and_grids)
{
  // charges on, so that the Ewald-real grid is built as well; only the VDW grid size is given explicitly and the
  // Coulomb grid size must then follow from its spacing (not read from the unset optional)
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  const std::size_t methaneType = forceField.findPseudoAtom("CH4").value();
  forceField.settings.gridPseudoAtomIndices = {methaneType};
  forceField.settings.numberOfVDWGridPoints = uint3(8, 8, 8);
  forceField.settings.spacingCoulombGrid = 1.5;
  forceField.settings.numberOfGridTestPoints = 10;

  Framework itq29 = Framework::makeITQ29(forceField, int3(1, 1, 1));

  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::GibbsSwapCFCMC, 1.0);
  Component methane = Component::makeMethane(forceField, 0);
  methane.mc_moves_probabilities = probabilities;
  methane.lambdaGC.computeDUdlambda = true;

  std::vector<System> systems{
      System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {itq29}, {methane}, {}, {4}, 5),
      System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {itq29}, {methane}, {}, {6}, 5)};
  ASSERT_EQ(systems[0].numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC[0], 1uz);

  // the driver writes its output files relative to the working directory
  const std::filesystem::path previous = std::filesystem::current_path();
  const std::filesystem::path scratch = std::filesystem::temp_directory_path() / "raspa3_tmmc_setup_test";
  std::filesystem::create_directories(scratch);
  std::filesystem::current_path(scratch);
  {
    RandomNumber random(3);
    MonteCarloTransitionMatrix tmmc(makeSchedule(0, 0, 0), systems, random, 5);
    tmmc.setup();

    EXPECT_TRUE(tmmc.systems[0].containsTheFractionalMolecule(0));
    EXPECT_FALSE(tmmc.systems[1].containsTheFractionalMolecule(0));
    for (std::size_t systemId : {0uz, 1uz})
    {
      System& system = tmmc.systems[systemId];
      const std::size_t fractional = system.indexOfGibbsSwapFractionalMoleculesPerComponent_CFCMC(0);
      const std::uint8_t expectedGroupId = systemId == 0 ? system.components[0].lambdaGC.dUdlambdaGroupId : 0;
      for (const Atom& atom : system.spanOfMolecule(0, fractional))
      {
        EXPECT_EQ(atom.groupId, expectedGroupId) << "system " << systemId;
      }
      EXPECT_TRUE(system.interpolationGrids[methaneType].has_value()) << "system " << systemId;
      EXPECT_TRUE(system.interpolationGrids.back().has_value()) << "system " << systemId;  // Ewald-real grid
    }
    EXPECT_EQ(tmmc.streams.size(), 2uz);
  }
  std::filesystem::current_path(previous);
  std::filesystem::remove_all(scratch);
}


// Equipped for a transition matrix over the (N, lambda) chain: the matrix biases the chain, lambda included. The
// matrix has the lambda grid of the CFCMC swap of its component, as the driver sets it for a system put together
// in code (a matrix of one lambda bin beside a swap on several is made a chain by the driver, with the arrays of
// the matrix still sized for one bin until its initialization stage); a single bin is a matrix over N alone.
System makeTMMCChainSystem(std::size_t numberOfLambdaBins)
{
  System system = makeMethaneCFCMCSystem(10, Move::Types::SwapCFCMC, numberOfLambdaBins);
  system.tmmc.doTMMC = true;
  system.tmmc.useBias = true;
  system.tmmc.useTMBias = true;
  system.tmmc.minMacrostate = 0;
  system.tmmc.maxMacrostate = 40;
  system.tmmc.numberOfLambdaBins = numberOfLambdaBins;
  system.tmmc.updateTMEvery = 100;
  system.tmmc.initialize();
  return system;
}

// On an (N, lambda) chain the transition-matrix bias flattens lambda itself, and the Wang-Landau bias of the
// lambda moves has to stay off: the endpoint jumps of the chain carried both, and were biased twice. Without
// the chain the lambda Wang-Landau bias is what flattens lambda, and it has to move.
TEST(MC_DRIVER_BOOKKEEPING, transition_matrix_on_a_lambda_chain_leaves_the_lambda_wang_landau_bias_alone)
{
  for (bool chain : {true, false})
  {
    std::vector<System> systems{makeTMMCChainSystem(chain ? 21uz : 1uz)};
    ASSERT_EQ(systems[0].tmmc.lambdaChain(), chain);
    ASSERT_TRUE(systems[0].lambdaWangLandauIsActive(0uz));

    RandomNumber random(11);
    MonteCarloTransitionMatrix tmmc(makeSchedule(0, 0, 3), systems, random, 5);
    for (System& system : tmmc.systems) system.runningEnergies = system.computeTotalEnergies();
    tmmc.simulationStage = MonteCarloTransitionMatrix::SimulationStage::Equilibration;
    for (std::size_t cycle = 0; cycle < 3; ++cycle) tmmc.performCycle();

    const std::vector<double>& biasFactor = tmmc.systems[0].components[0].lambdaGC.biasFactor;
    const bool moved = std::ranges::any_of(biasFactor, [](double value) { return value != 0.0; });
    EXPECT_EQ(moved, !chain) << (chain ? "on the chain" : "without the chain");
  }
}

// Wang-Landau on the transition matrix itself (tmmc.useWangLandau): while it owns the bias, every visit
// lowers the bias of the visited state and the collection-matrix rebuild leaves the bias alone; production
// starts by handing the bias to the transition matrix, -ln Pi from the collection matrix. Before, this driver
// never visited and never switched, so a system with the flag set ran with a bias of zero throughout.
TEST(MC_DRIVER_BOOKKEEPING, transition_matrix_wang_landau_is_driven_and_handed_over_at_production)
{
  std::vector<System> systems{makeTMMCChainSystem(1uz)};
  systems[0].tmmc.useWangLandau = true;

  const std::filesystem::path previous = std::filesystem::current_path();
  const std::filesystem::path scratch = std::filesystem::temp_directory_path() / "raspa3_tmmc_wang_landau_test";
  std::filesystem::create_directories(scratch);
  std::filesystem::current_path(scratch);
  {
    RandomNumber random(21);
    MonteCarloTransitionMatrix tmmc(makeSchedule(0, 0, 2), systems, random, 5);
    tmmc.setup();
    for (System& system : tmmc.systems) system.runningEnergies = system.computeTotalEnergies();
    tmmc.simulationStage = MonteCarloTransitionMatrix::SimulationStage::Equilibration;
    for (std::size_t cycle = 0; cycle < 2; ++cycle) tmmc.performCycle();

    const TransitionMatrix& matrix = tmmc.systems[0].tmmc;
    EXPECT_TRUE(matrix.useWangLandau);
    bool visitedStatesArePenalized = false;
    for (std::size_t i = 0; i < matrix.histogram.size(); ++i)
    {
      if (matrix.histogram[i] > 0uz)
      {
        EXPECT_LT(matrix.bias[i], 0.0) << "visited state " << i << " carries no Wang-Landau penalty";
        visitedStatesArePenalized = true;
      }
    }
    EXPECT_TRUE(visitedStatesArePenalized);

    // the bias the hand-over gives: -ln Pi of the collection matrix as it stands at the end of the equilibration
    TransitionMatrix handedOver = matrix;
    handedOver.switchToTMBias();
    ASSERT_FALSE(handedOver.useWangLandau);
    bool biasIsOfTheMatrix = false;
    for (std::size_t i = 0; i < handedOver.bias.size(); ++i)
    {
      EXPECT_DOUBLE_EQ(handedOver.bias[i], -handedOver.lnpi[i]) << "state " << i;
      if (handedOver.bias[i] != matrix.bias[i]) biasIsOfTheMatrix = true;
    }
    EXPECT_TRUE(biasIsOfTheMatrix) << "the hand-over left the Wang-Landau bias in place";

    // production starts with that bias and drops the statistics it was derived from, ln Pi with them (a ln Pi
    // left in place was written next to an empty collection matrix)
    tmmc.production();
    EXPECT_FALSE(matrix.useWangLandau);
    for (std::size_t i = 0; i < matrix.bias.size(); ++i)
    {
      EXPECT_DOUBLE_EQ(matrix.bias[i], handedOver.bias[i]) << "state " << i;
    }
  }
  std::filesystem::current_path(previous);
  std::filesystem::remove_all(scratch);
}

// The macrostate is the count of the adsorbate. A cation beside it is a component too, picked for its own
// moves; what it counts to must not be written into the adsorbate's histogram. Nine cations in a range of
// twelve, and an adsorbate at a pressure where it stays below five: the histogram above five is empty.
TEST(MC_DRIVER_BOOKKEEPING, transition_matrix_samples_the_count_of_the_adsorbate_and_not_of_a_cation)
{
  const ForceField forceField({{"CH4", false, 16.04246, 0.0, 0.0, 6, false}}, {{158.5, 3.72}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  MCMoveProbabilities swapping;
  swapping.setProbability(Move::Types::Translation, 1.0);
  swapping.setProbability(Move::Types::Swap, 1.0);
  const Component methane(forceField, "methane", 190.564, 4599200, 0.01142,
                          {Atom({0, 0, 0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21, swapping, std::nullopt,
                          false);
  MCMoveProbabilities moving;
  moving.setProbability(Move::Types::Translation, 1.0);
  Component cation(forceField, "cation", 190.564, 4599200, 0.01142, {Atom({0, 0, 0}, 0.0, 1.0, 0, 0, 1, false, false)},
                   {}, {}, 5, 21, moving, std::nullopt, false);
  cation.type = Component::Type::Cation;

  std::vector<System> systems{System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1.0, 1.0, {},
                                     {methane, cation}, {}, {0, 9}, 5)};
  System& system = systems[0];
  system.tmmc.doTMMC = true;
  system.tmmc.useBias = true;
  system.tmmc.useTMBias = true;
  system.tmmc.minMacrostate = 0;
  system.tmmc.maxMacrostate = 12;
  system.tmmc.initialize();
  ASSERT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 9uz);

  RandomNumber random(3);
  MonteCarloTransitionMatrix tmmc(makeSchedule(0, 2, 0), systems, random, 5);
  for (System& walker : tmmc.systems) walker.runningEnergies = walker.computeTotalEnergies();
  tmmc.simulationStage = MonteCarloTransitionMatrix::SimulationStage::Initialization;
  for (std::size_t cycle = 0; cycle < 2; ++cycle) tmmc.performCycle();

  const std::vector<std::size_t>& histogram = tmmc.systems[0].tmmc.histogram;
  ASSERT_EQ(histogram.size(), 13uz);
  std::size_t low = 0, high = 0;
  for (std::size_t n = 0; n < histogram.size(); ++n) (n < 5 ? low : high) += histogram[n];
  EXPECT_GT(low, 0uz);
  EXPECT_EQ(high, 0uz) << "visits above five: " << high;
  EXPECT_LT(tmmc.systems[0].numberOfIntegerMoleculesPerComponent[0], 5uz);
}

// 'containsTheFractionalMolecule' is the per-component flag of the serial Gibbs CFCMC swaps: the fractional molecule
// of a component that uses them is active in one box at a time. A co-adsorbate with a box-local CFCMC swap keeps its
// own fractional molecule in every box, biased and present throughout: its Wang-Landau updates, its lambda
// histogram and its factor in the weight of the box must not follow the flag of the other component.
TEST(MC_DRIVER_BOOKKEEPING, box_local_cfcmc_component_stays_active_in_the_box_without_the_gibbs_fractional_molecule)
{
  const ForceField forceField({{"CH4", false, 16.04246, 0.0, 0.0, 6, false}}, {{158.5, 3.72}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  const auto makeComponent = [&](const std::string& name, Move::Types swapMove)
  {
    MCMoveProbabilities probabilities;
    probabilities.setProbability(Move::Types::Translation, 1.0);
    probabilities.setProbability(swapMove, 1.0);
    return Component(forceField, name, 190.564, 4599200, 0.01142, {Atom({0, 0, 0}, 0.0, 1.0, 0, 0, 0, false, false)},
                     {}, {}, 5, 21, probabilities, std::nullopt, false);
  };
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {},
                {makeComponent("gibbs", Move::Types::GibbsSwapCFCMC), makeComponent("local", Move::Types::SwapCFCMC)},
                {}, {10, 10}, 5);
  ASSERT_GT(system.numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC[0], 0uz);
  ASSERT_EQ(system.numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC[1], 0uz);
  ASSERT_TRUE(system.usesGCLambda(0uz));
  ASSERT_TRUE(system.usesGCLambda(1uz));

  system.components[0].lambdaGC.biasFactor[system.components[0].lambdaGC.currentBin] = -3.0;
  system.components[1].lambdaGC.biasFactor[system.components[1].lambdaGC.currentBin] = -2.0;
  const double gibbsWeight = system.components[0].lambdaGC.weight();
  const double localWeight = system.components[1].lambdaGC.weight();
  ASSERT_NE(localWeight, 1.0);

  // the box that holds the Gibbs fractional molecule: both are active
  system.setContainsTheFractionalMolecules(true);
  EXPECT_TRUE(system.gcFractionalMoleculeIsActive(0uz));
  EXPECT_TRUE(system.gcFractionalMoleculeIsActive(1uz));
  EXPECT_DOUBLE_EQ(system.weight(), gibbsWeight * localWeight);

  // the other box: the Gibbs fractional molecule is a switched-off placeholder, the local one is as it was
  system.setContainsTheFractionalMolecules(false);
  EXPECT_FALSE(system.gcFractionalMoleculeIsActive(0uz));
  EXPECT_TRUE(system.gcFractionalMoleculeIsActive(1uz));
  EXPECT_FALSE(system.lambdaWangLandauIsActive(0uz));
  EXPECT_TRUE(system.lambdaWangLandauIsActive(1uz));
  EXPECT_DOUBLE_EQ(system.weight(), localWeight);
}

// The inner cut-off of the dual cut-off scheme never exceeds the full one: an automatic full cut-off is only known
// once the box is and shrinks with it, and an inner cut-off left above it reached past the cell list and the
// minimum image.
TEST(MC_DRIVER_BOOKKEEPING, inner_cut_off_of_the_dual_cut_off_scheme_is_held_at_the_full_cut_off)
{
  System system = makeMethaneCFCMCSystem(5);
  system.forceField.settings.useDualCutOff = true;

  system.forceField.settings.dualCutOff = 6.0;
  EXPECT_DOUBLE_EQ(system.makeGrowContext().cutOffMoleculeVDW, 6.0);
  EXPECT_DOUBLE_EQ(system.makeGrowContext(CBMC::CutOffMode::Full).cutOffMoleculeVDW,
                   system.forceField.cutOffMoleculeVDW);

  system.forceField.settings.dualCutOff = system.forceField.cutOffMoleculeVDW + 5.0;
  for (const CBMC::CutOffMode mode : {CBMC::CutOffMode::Growth, CBMC::CutOffMode::Inner})
  {
    const CBMC::GrowContext context = system.makeGrowContext(mode);
    EXPECT_DOUBLE_EQ(context.cutOffMoleculeVDW, system.forceField.cutOffMoleculeVDW);
    EXPECT_DOUBLE_EQ(context.cutOffFrameworkVDW, system.forceField.cutOffFrameworkVDW);
    EXPECT_DOUBLE_EQ(context.cutOffCoulomb, system.forceField.cutOffCoulomb);
  }
}

// The first of a pair of neighbouring integers is the system a Monte Carlo step is taken in: it is uniform over
// all of them. Drawn as the lower of the pair and swapped with probability one half, the first and the last were
// the first of the pair half as often as the ones in between.
TEST(MC_DRIVER_BOOKKEEPING, the_system_of_a_step_is_drawn_uniformly)
{
  RandomNumber random(11);
  const std::size_t size = 4;
  const std::size_t draws = 400000;
  std::vector<std::size_t> counts(size, 0);
  for (std::size_t i = 0; i < draws; ++i)
  {
    const auto [first, second] = random.randomPairAdjacentIntegers(size);
    ASSERT_LT(first, size);
    ASSERT_LT(second, size);
    ASSERT_EQ(first > second ? first - second : second - first, 1uz);
    ++counts[first];
  }
  for (std::size_t count : counts)
  {
    EXPECT_NEAR(static_cast<double>(count) / static_cast<double>(draws), 0.25, 0.005);
  }

  // one integer is its own pair, two are each other's neighbour
  EXPECT_EQ(random.randomPairAdjacentIntegers(1), (std::pair<std::size_t, std::size_t>(0, 0)));
  for (std::size_t i = 0; i < 100; ++i)
  {
    const auto [first, second] = random.randomPairAdjacentIntegers(2);
    EXPECT_EQ(first + second, 1uz);
  }
}
