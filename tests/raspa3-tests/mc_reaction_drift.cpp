#include <gtest/gtest.h>

import std;

import double3;
import int3;
import forcefield;
import framework;
import component;
import system;
import monte_carlo;
import simulationbox;
import running_energy;
import mc_moves_move_types;
import mc_moves_probabilities;
import reaction;
import reactions;
import mc_moves_reaction_common;
import mc_moves_reaction_cfcmc;
import mc_moves_reaction_conventional_cfcmc;
import mc_moves_reaction_conventional_cbcfcmc;
import randomnumbers;

namespace {

constexpr double kReactionDriftVDWCutoff = 50.0;
constexpr size_t kReactionDriftCycles = 20;
constexpr size_t kReactionDriftInitializationCycles = 5;
constexpr size_t kReactionDriftEquilibrationCycles = 5;
constexpr size_t kReactionDriftSingleMoveAttempts = 250;

std::filesystem::path repositoryRoot()
{
  return std::filesystem::path(__FILE__).parent_path().parent_path().parent_path();
}

ForceField makeAlkaneForceField()
{
  return ForceField({{"CH3", false, 15.04, 0.0, 0.0, 6, false},
                     {"CH2", false, 14.03, 0.0, 0.0, 6, false}},
                    {{98.0, 3.75}, {46.0, 3.95}}, ForceField::MixingRule::Lorentz_Berthelot, kReactionDriftVDWCutoff,
                     kReactionDriftVDWCutoff, 12.0, true, true, false);
}

ForceField makeZeoliteAlkaneForceField()
{
  return ForceField({{"-", false, 0.0, 0.0, 0.0, 0, false},
                     {"Si", true, 28.0855, 2.05, 0.0, 14, false},
                     {"Al", true, 26.982, 2.05, 0.0, 13, false},
                     {"O", true, 15.999, -1.025, 0.0, 8, false},
                     {"Na+", false, 12.0, 1.0, 0.0, 6, false},
                     {"Cl-", false, 15.9994, -1.0, 0.0, 8, false},
                     {"CH4", false, 16.04246, 0.0, 0.0, 6, false},
                     {"C_co2", false, 12.0, 0.6512, 0.0, 6, false},
                     {"O_co2", false, 15.9994, -0.3256, 0.0, 8, false},
                     {"Ow", false, 15.9996, 0.0, 0.0, 8, false},
                     {"Hw", false, 1.0008, 0.241, 0.0, 1, false},
                     {"Lw", false, 0.0, -0.241, 0.0, 0, false},
                     {"probe-He", false, 4.002602, 0.0, 0.0, 2, false},
                     {"probe-Ar", false, 39.948, 0.0, 0.0, 18, false},
                     {"probe-CH4", false, 16.04246, 0.0, 0.0, 6, false},
                     {"probe-N2", false, 14.00674, 0.0, 0.0, 6, false},
                     {"CH3", false, 15.04, 0.0, 0.0, 6, false},
                     {"CH2", false, 14.03, 0.0, 0.0, 6, false}},
                    {{1.0, 1.0},
                     {22.0, 2.30},
                     {22.0, 2.30},
                     {53.0, 3.30},
                     {15.0966, 2.65755},
                     {142.562, 3.51932},
                     {158.5, 3.72},
                     {29.933, 2.745},
                     {85.671, 3.017},
                     {89.633, 3.097},
                     {0.0, 1.0},
                     {0.0, 1.0},
                     {10.9, 2.64},
                     {124.070, 3.38},
                     {158.5, 3.72},
                     {91.5, 3.681},
                     {98.0, 3.75},
                     {46.0, 3.95}},
                    ForceField::MixingRule::Lorentz_Berthelot, kReactionDriftVDWCutoff, kReactionDriftVDWCutoff, 12.0,
                     true, false, true);
}

Component makeAlkaneFromExample(const ForceField &forceField, std::size_t componentId, std::string_view name,
                                const MCMoveProbabilities &probabilities)
{
  const std::filesystem::path moleculePath =
      repositoryRoot() / "examples/basic/4_mc_binary_mixture_propane_butane_in_box" / name;
  return Component(Component::Type::Adsorbate, componentId, forceField, std::string(name), moleculePath.string(), 5,
                   21, probabilities, std::nullopt, false);
}

Move::Types reactionMoveFromProbabilities(const MCMoveProbabilities &probabilities)
{
  if (probabilities.getProbability(Move::Types::ReactionCFCMC) > 0.0) return Move::Types::ReactionCFCMC;
  if (probabilities.getProbability(Move::Types::ReactionCBCFCMC) > 0.0) return Move::Types::ReactionCBCFCMC;
  if (probabilities.getProbability(Move::Types::ReactionConventionalCFCMC) > 0.0)
    return Move::Types::ReactionConventionalCFCMC;
  if (probabilities.getProbability(Move::Types::ReactionConventionalCBCFCMC) > 0.0)
    return Move::Types::ReactionConventionalCBCFCMC;
  return Move::Types::ReactionCBMC;
}

void checkEnergyDrift(System &s)
{
  constexpr double tolerance = 1e-6;
  RunningEnergy recomputedEnergies = s.computeTotalEnergies();
  RunningEnergy drift = s.runningEnergies - recomputedEnergies;

  EXPECT_NEAR(drift.potentialEnergy(), 0.0, tolerance);
  EXPECT_NEAR(drift.externalFieldVDW, 0.0, tolerance);
  EXPECT_NEAR(drift.frameworkMoleculeVDW, 0.0, tolerance);
  EXPECT_NEAR(drift.moleculeMoleculeVDW, 0.0, tolerance);
  EXPECT_NEAR(drift.externalFieldCharge, 0.0, tolerance);
  EXPECT_NEAR(drift.frameworkMoleculeCharge, 0.0, tolerance);
  EXPECT_NEAR(drift.moleculeMoleculeCharge, 0.0, tolerance);
  EXPECT_NEAR(drift.ewald_fourier, 0.0, tolerance);
  EXPECT_NEAR(drift.ewald_self, 0.0, tolerance);
  EXPECT_NEAR(drift.ewald_exclusion, 0.0, tolerance);
  EXPECT_NEAR(drift.intraVDW, 0.0, tolerance);
  EXPECT_NEAR(drift.intraCoul, 0.0, tolerance);
  EXPECT_NEAR(drift.tail, 0.0, tolerance);
  EXPECT_NEAR(drift.polarization, 0.0, tolerance);
  EXPECT_NEAR(drift.dudlambdaVDW, 0.0, tolerance);
  EXPECT_NEAR(drift.dudlambdaCharge, 0.0, tolerance);
  EXPECT_NEAR(drift.dudlambdaEwald, 0.0, tolerance);
}

System makeMultiComponentReactionSystem(const MCMoveProbabilities& systemProbabilities,
                                        std::vector<std::size_t> reactantStoichiometry,
                                        std::vector<std::size_t> productStoichiometry,
                                        std::vector<std::size_t> initialCounts)
{
  const ForceField forceField = makeAlkaneForceField();
  SimulationBox box = SimulationBox(30.0, 30.0, 30.0);

  MCMoveProbabilities componentProbabilities = MCMoveProbabilities();
  componentProbabilities.setProbability(Move::Types::Translation, 1.0);
  componentProbabilities.setProbability(Move::Types::Rotation, 1.0);
  componentProbabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);

  Component propane = makeAlkaneFromExample(forceField, 0, "propane", componentProbabilities);
  Component butane = makeAlkaneFromExample(forceField, 1, "butane", componentProbabilities);
  propane.lnPartitionFunction = 1.0;
  butane.lnPartitionFunction = 1.0;

  System system = System(forceField, box, false, 500.0, 1e4, 1.0, {}, {propane, butane}, {}, std::move(initialCounts),
                         5, systemProbabilities);

  Reaction reaction(0, std::move(reactantStoichiometry), std::move(productStoichiometry));
  reaction.reactionMove = reactionMoveFromProbabilities(systemProbabilities);
  system.reactions.list.push_back(reaction);
  return system;
}

System makePropaneButaneReactionSystem(const MCMoveProbabilities &systemProbabilities)
{
  return makeMultiComponentReactionSystem(systemProbabilities, {2, 0}, {0, 2}, {40, 20});
}

System makeCO2N2SerialReactionSystem(const MCMoveProbabilities& systemProbabilities,
                                     std::vector<std::size_t> reactantStoichiometry,
                                     std::vector<std::size_t> productStoichiometry,
                                     std::vector<std::size_t> initialCounts)
{
  const std::filesystem::path exampleDir =
      repositoryRoot() / "examples/basic/2_mc_co2_n2_in_two_independent_boxes";
  ForceField forceField = ForceField::readForceField(exampleDir.string(), "force_field.json").value();
  forceField.chargeMethod = ForceField::ChargeMethod::Ewald;
  forceField.useCharge = true;
  forceField.cutOffFrameworkVDW = kReactionDriftVDWCutoff;
  forceField.cutOffMoleculeVDW = kReactionDriftVDWCutoff;
  forceField.cutOffCoulomb = 12.0;
  SimulationBox box = SimulationBox(40.0, 40.0, 40.0);

  MCMoveProbabilities componentProbabilities = MCMoveProbabilities();
  componentProbabilities.setProbability(Move::Types::Translation, 1.0);
  componentProbabilities.setProbability(Move::Types::Rotation, 1.0);
  componentProbabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);

  Component co2(Component::Type::Adsorbate, 0, forceField, "CO2", (exampleDir / "CO2").string(), 5, 21,
                componentProbabilities, std::nullopt, false);
  Component n2(Component::Type::Adsorbate, 1, forceField, "N2", (exampleDir / "N2").string(), 5, 21,
               componentProbabilities, std::nullopt, false);
  co2.lnPartitionFunction = 55.0;
  n2.lnPartitionFunction = 53.0;

  System system = System(forceField, box, false, 300.0, 101300.0, 1.0, {}, {co2, n2}, {}, std::move(initialCounts),
                         5, systemProbabilities);

  Reaction reaction(0, std::move(reactantStoichiometry), std::move(productStoichiometry));
  reaction.reactionMove = reactionMoveFromProbabilities(systemProbabilities);
  system.reactions.list.push_back(reaction);
  return system;
}

System makeCO2N2TwoSerialReactionSystem(const MCMoveProbabilities& systemProbabilities,
                                        std::vector<std::size_t> initialCounts,
                                        std::vector<std::size_t> secondReactants = {2, 0},
                                        std::vector<std::size_t> secondProducts = {0, 2})
{
  System system = makeCO2N2SerialReactionSystem(systemProbabilities, {2, 0}, {0, 2}, std::move(initialCounts));
  Reaction reaction2(1, std::move(secondReactants), std::move(secondProducts));
  reaction2.reactionMove = Move::Types::ReactionCFCMC;
  system.reactions.list.push_back(reaction2);
  return system;
}

System makeMixedParallelSerialReactionSystem(const MCMoveProbabilities& systemProbabilities)
{
  const ForceField forceField = makeAlkaneForceField();
  SimulationBox box = SimulationBox(30.0, 30.0, 30.0);

  MCMoveProbabilities componentProbabilities = MCMoveProbabilities();
  componentProbabilities.setProbability(Move::Types::Translation, 1.0);
  componentProbabilities.setProbability(Move::Types::Rotation, 1.0);
  componentProbabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);

  Component propane = makeAlkaneFromExample(forceField, 0, "propane", componentProbabilities);
  Component butane = makeAlkaneFromExample(forceField, 1, "butane", componentProbabilities);
  propane.lnPartitionFunction = 1.0;
  butane.lnPartitionFunction = 1.0;

  System system = System(forceField, box, false, 500.0, 1e4, 1.0, {}, {propane, butane}, {},
                         std::vector<std::size_t>{40, 0}, 5, systemProbabilities);

  Reaction parallelReaction(0, {2, 0}, {0, 2});
  parallelReaction.reactionMove = Move::Types::ReactionConventionalCFCMC;
  system.reactions.list.push_back(parallelReaction);

  Reaction serialReaction(1, {2, 0}, {0, 2});
  serialReaction.reactionMove = Move::Types::ReactionCFCMC;
  system.reactions.list.push_back(serialReaction);
  return system;
}

System makeTwoParallelReactionSystem(const MCMoveProbabilities& systemProbabilities)
{
  System system = makePropaneButaneReactionSystem(systemProbabilities);

  Reaction reaction2(1, {2, 0}, {0, 2});
  reaction2.reactionMove = Move::Types::ReactionConventionalCFCMC;
  system.reactions.list.push_back(reaction2);
  return system;
}

System makeZeoliteMultiComponentReactionSystem(const MCMoveProbabilities& systemProbabilities,
                                               std::vector<std::size_t> reactantStoichiometry,
                                               std::vector<std::size_t> productStoichiometry,
                                               std::vector<std::size_t> initialCounts)
{
  const ForceField forceField = makeZeoliteAlkaneForceField();
  Framework f = Framework::makeFAU(forceField, int3(1, 1, 1));

  MCMoveProbabilities componentProbabilities = MCMoveProbabilities();
  componentProbabilities.setProbability(Move::Types::Translation, 1.0);
  componentProbabilities.setProbability(Move::Types::Rotation, 1.0);
  componentProbabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);

  Component propane = makeAlkaneFromExample(forceField, 0, "propane", componentProbabilities);
  Component butane = makeAlkaneFromExample(forceField, 1, "butane", componentProbabilities);
  propane.lnPartitionFunction = 1.0;
  butane.lnPartitionFunction = 1.0;

  System system = System(forceField, std::nullopt, false, 500.0, 1e4, 1.0, {f}, {propane, butane}, {},
                         std::move(initialCounts), 5, systemProbabilities);

  Reaction reaction(0, std::move(reactantStoichiometry), std::move(productStoichiometry));
  reaction.reactionMove = reactionMoveFromProbabilities(systemProbabilities);
  system.reactions.list.push_back(reaction);
  return system;
}

System makeZeolitePropaneButaneReactionSystem(const MCMoveProbabilities& systemProbabilities)
{
  return makeZeoliteMultiComponentReactionSystem(systemProbabilities, {2, 0}, {0, 2}, {12, 6});
}

void runReactionDriftTest(std::vector<System> systems)
{
  const size_t numberOfCycles{kReactionDriftCycles};
  const size_t numberOfInitializationCycles{kReactionDriftInitializationCycles};
  const size_t numberOfEquilibrationCycles{kReactionDriftEquilibrationCycles};
  const size_t printEvery{1000};
  const size_t writeBinaryRestartEvery{10000};
  const size_t rescaleWangLandauEvery{5000};
  const size_t optimizeMCMovesEvery{5000};
  const size_t numberOfBlocks{5};
  const bool outputToFiles{false};

  MonteCarlo mc = MonteCarlo(numberOfCycles, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery, systems, {},
                             numberOfBlocks, outputToFiles);

  mc.run();

  for (System &s : mc.systems)
  {
    checkEnergyDrift(s);
  }
}

void runReactionBoundaryDriftTest(std::vector<System> systems)
{
  const size_t numberOfCycles{kReactionDriftCycles};
  const size_t numberOfInitializationCycles{kReactionDriftInitializationCycles};
  const size_t numberOfEquilibrationCycles{kReactionDriftEquilibrationCycles};
  const size_t printEvery{1000};
  const size_t writeBinaryRestartEvery{10000};
  const size_t rescaleWangLandauEvery{5000};
  const size_t optimizeMCMovesEvery{5000};
  const size_t numberOfBlocks{5};
  const bool outputToFiles{false};

  MonteCarlo mc = MonteCarlo(numberOfCycles, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery, systems, {},
                             numberOfBlocks, outputToFiles);

  mc.run();

  for (System &s : mc.systems)
  {
    checkEnergyDrift(s);
  }
}

}  // namespace

TEST(MC_REACTION_DRIFT, reaction_cbmc)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCBMC, 1.0);

  runReactionDriftTest({makePropaneButaneReactionSystem(systemProbabilities)});
}

void setupLambdaOnlyReaction(System& system)
{
  for (Component& component : system.components)
  {
    component.mc_moves_probabilities = MCMoveProbabilities();
  }
  system.rescaleMoveProbabilities();

  Reaction& reaction = system.reactions.list[0];
  reaction.currentLambda = 0.5;
  reaction.maximumLambdaChange = 0.05;
  MC_Moves::ReactionCommon::setReactionFractionalScaling(system, reaction, reaction.currentLambda);
  system.runningEnergies = system.computeTotalEnergies();
}

TEST(MC_REACTION_DRIFT, reaction_cfcmc)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionConventionalCFCMC, 1.0);

  System system = makePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupLambdaOnlyReaction(system);
  runReactionDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_cfcmc_cbmc)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionConventionalCBCFCMC, 1.0);

  System system = makePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupLambdaOnlyReaction(system);
  runReactionDriftTest({std::move(system)});
}

void setupBoundaryCrossingReaction(System& system)
{
  for (Component& component : system.components)
  {
    component.mc_moves_probabilities = MCMoveProbabilities();
  }
  system.rescaleMoveProbabilities();

  Reaction& reaction = system.reactions.list[0];
  reaction.currentLambda = 0.85;
  reaction.maximumLambdaChange = 0.3;
  MC_Moves::ReactionCommon::setReactionFractionalScaling(system, reaction, reaction.currentLambda);
  system.runningEnergies = system.computeTotalEnergies();
}

void setupCBMCStoichiometryBoundaryReaction(System& system)
{
  system.runningEnergies = system.computeTotalEnergies();
}

TEST(MC_REACTION_DRIFT, reaction_cfcmc_boundary)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionConventionalCFCMC, 1.0);

  System system = makePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupBoundaryCrossingReaction(system);
  runReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_cfcmc_cbmc_boundary)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionConventionalCBCFCMC, 1.0);

  System system = makePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupBoundaryCrossingReaction(system);
  runReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_cbmc_boundary)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCBMC, 1.0);

  System system = makeMultiComponentReactionSystem(systemProbabilities, {2, 0}, {0, 2}, {4, 2});
  setupCBMCStoichiometryBoundaryReaction(system);
  runReactionDriftTest({std::move(system)});
}

void setupSerialLambdaOnlyReaction(System& system)
{
  for (Component& component : system.components)
  {
    component.mc_moves_probabilities = MCMoveProbabilities();
  }
  system.rescaleMoveProbabilities();

  Reaction& reaction = system.reactions.list[0];
  reaction.currentLambda = 0.25;
  reaction.maximumLambdaChange = 0.05;
  reaction.lambdaSwitchPoint = 2.0;
  MC_Moves::ReactionCommon::setSerialReactionFractionalScaling(system, reaction, reaction.currentLambda);
  system.runningEnergies = system.computeTotalEnergies();
}

void setupSerialWholeMoleculeReaction(System& system)
{
  for (Component& component : system.components)
  {
    component.mc_moves_probabilities = MCMoveProbabilities();
  }
  system.rescaleMoveProbabilities();

  Reaction& reaction = system.reactions.list[0];
  reaction.currentLambda = 0.85;
  reaction.maximumLambdaChange = 0.05;
  reaction.lambdaSwitchPoint = 0.5;
  MC_Moves::ReactionCommon::setSerialReactionFractionalScaling(system, reaction, reaction.currentLambda);
  system.runningEnergies = system.computeTotalEnergies();
}

void setupSerialBoundaryCrossingReaction(System& system, bool fromAboveSwitchPoint = false)
{
  for (Component& component : system.components)
  {
    component.mc_moves_probabilities = MCMoveProbabilities();
  }
  system.rescaleMoveProbabilities();

  Reaction& reaction = system.reactions.list[0];
  reaction.currentLambda = fromAboveSwitchPoint ? 0.55 : 0.45;
  reaction.maximumLambdaChange = 0.3;
  reaction.lambdaSwitchPoint = 0.5;
  MC_Moves::ReactionCommon::setSerialReactionFractionalScaling(system, reaction, reaction.currentLambda);
  system.runningEnergies = system.computeTotalEnergies();
}

void runSerialReactionBoundaryDriftTest(std::vector<System> systems)
{
  const size_t numberOfCycles{kReactionDriftCycles};
  const size_t numberOfInitializationCycles{kReactionDriftInitializationCycles};
  const size_t numberOfEquilibrationCycles{kReactionDriftEquilibrationCycles};
  const size_t printEvery{1000};
  const size_t writeBinaryRestartEvery{10000};
  const size_t rescaleWangLandauEvery{5000};
  const size_t optimizeMCMovesEvery{5000};
  const size_t numberOfBlocks{5};
  const bool outputToFiles{false};

  MonteCarlo mc = MonteCarlo(numberOfCycles, numberOfInitializationCycles, numberOfEquilibrationCycles, printEvery,
                             writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery, systems, {},
                             numberOfBlocks, outputToFiles);

  mc.run();

  for (System &s : mc.systems)
  {
    checkEnergyDrift(s);
  }
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_lambda)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupSerialLambdaOnlyReaction(system);
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupSerialLambdaOnlyReaction(system);
  system.reactions.list[0].lambdaSwitchPoint = 0.5;
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_cbmc)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCBCFCMC, 1.0);

  System system = makePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupSerialLambdaOnlyReaction(system);
  system.reactions.list[0].lambdaSwitchPoint = 0.5;
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_whole_molecule)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupSerialWholeMoleculeReaction(system);
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_cbmc_whole_molecule)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCBCFCMC, 1.0);

  System system = makePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupSerialWholeMoleculeReaction(system);
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_boundary)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupSerialBoundaryCrossingReaction(system);
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_cbmc_boundary)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCBCFCMC, 1.0);

  System system = makePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupSerialBoundaryCrossingReaction(system);
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_whole_molecule_boundary)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupSerialBoundaryCrossingReaction(system, true);
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_cbmc_whole_molecule_boundary)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCBCFCMC, 1.0);

  System system = makePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupSerialBoundaryCrossingReaction(system, true);
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_multi_stoichiometry_boundary)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupSerialBoundaryCrossingReaction(system);
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_multi_stoichiometry)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupSerialLambdaOnlyReaction(system);
  system.reactions.list[0].lambdaSwitchPoint = 0.5;
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_cbmc_multi_stoichiometry_boundary)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCBCFCMC, 1.0);

  System system = makePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupSerialBoundaryCrossingReaction(system);
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_cbmc_multi_stoichiometry)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCBCFCMC, 1.0);

  System system = makePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupSerialLambdaOnlyReaction(system);
  system.reactions.list[0].lambdaSwitchPoint = 0.5;
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_mixed_reactants)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupSerialLambdaOnlyReaction(system);
  system.reactions.list[0].lambdaSwitchPoint = 0.5;
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_co2_n2_ewald_no_heal_drift)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makeCO2N2SerialReactionSystem(systemProbabilities, {2, 0}, {0, 2}, {60, 20});
  system.createReactionFractionalMolecules();
  setupSerialLambdaOnlyReaction(system);
  system.reactions.list[0].lambdaSwitchPoint = 0.5;

  RandomNumber random(42);
  for (int move = 0; move < 500; ++move)
  {
    std::optional<RunningEnergy> delta = MC_Moves::reactionMove_CFCMC(random, system);
    if (delta)
    {
      system.runningEnergies += delta.value();
    }
  }
  checkEnergyDrift(system);
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_co2_n2_ewald)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makeCO2N2SerialReactionSystem(systemProbabilities, {2, 0}, {0, 2}, {60, 20});
  ASSERT_TRUE(system.forceField.useCharge);
  system.createReactionFractionalMolecules();
  setupSerialLambdaOnlyReaction(system);
  system.reactions.list[0].lambdaSwitchPoint = 0.5;
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_co2_n2_ewald_multi)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCBCFCMC, 1.0);

  System system = makeCO2N2SerialReactionSystem(systemProbabilities, {2, 0}, {0, 2}, {50, 25});
  ASSERT_TRUE(system.forceField.useCharge);
  system.createReactionFractionalMolecules();
  setupSerialLambdaOnlyReaction(system);
  system.reactions.list[0].lambdaSwitchPoint = 0.5;
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_mixed_reactants_boundary)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupSerialBoundaryCrossingReaction(system);
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_co2_n2_ewald_boundary)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makeCO2N2SerialReactionSystem(systemProbabilities, {2, 0}, {0, 2}, {60, 20});
  ASSERT_TRUE(system.forceField.useCharge);
  system.createReactionFractionalMolecules();
  setupSerialBoundaryCrossingReaction(system);
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_co2_n2_ewald_multi_boundary)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCBCFCMC, 1.0);

  System system = makeCO2N2SerialReactionSystem(systemProbabilities, {2, 0}, {0, 2}, {50, 25});
  ASSERT_TRUE(system.forceField.useCharge);
  system.createReactionFractionalMolecules();
  setupSerialBoundaryCrossingReaction(system);
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_co2_n2_ewald_mixed_boundary)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makeCO2N2SerialReactionSystem(systemProbabilities, {2, 0}, {0, 2}, {40, 30});
  ASSERT_TRUE(system.forceField.useCharge);
  system.createReactionFractionalMolecules();
  setupSerialBoundaryCrossingReaction(system);
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_co2_n2_two_reactions_init)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makeCO2N2SerialReactionSystem(systemProbabilities, {2, 0}, {0, 2}, {60, 20});
  Reaction reaction2(1, {0, 0}, {0, 0});
  reaction2.reactionMove = Move::Types::ReactionCFCMC;
  system.reactions.list.push_back(reaction2);
  EXPECT_EQ(system.reactions.list.size(), 2uz);
  system.createReactionFractionalMolecules();
  system.checkMoleculeIds();
  system.runningEnergies = system.computeTotalEnergies();
  EXPECT_TRUE(std::isfinite(system.runningEnergies.potentialEnergy()));
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_co2_n2_two_reactions_full)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makeCO2N2TwoSerialReactionSystem(systemProbabilities, {60, 20});
  system.createReactionFractionalMolecules();
  system.checkMoleculeIds();
  system.runningEnergies = system.computeTotalEnergies();
  EXPECT_TRUE(std::isfinite(system.runningEnergies.potentialEnergy()));
}

// regression test: two serial reactions with different stoichiometries ({2,0}->{0,2} and {1,1}->{2,0}).
// when reaction 0 flips to the product side its CO2 fractional molecules disappear, and the fractional
// molecule indices of reaction 1 must shift down accordingly. stale static indices made reaction 1 point
// into the integer-molecule region, so the same molecule could be treated as fractional and integer at
// once (NaN in the Ewald exclusion difference).
TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_co2_n2_two_reactions_side_flip_drift)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makeCO2N2TwoSerialReactionSystem(systemProbabilities, {60, 20}, {1, 1}, {2, 0});
  system.createReactionFractionalMolecules();

  for (Reaction& reaction : system.reactions.list)
  {
    reaction.currentLambda = 0.45;
    reaction.maximumLambdaChange = 0.3;
    reaction.lambdaSwitchPoint = 0.5;
    MC_Moves::ReactionCommon::setSerialReactionFractionalScaling(system, reaction, reaction.currentLambda);
  }
  system.runningEnergies = system.computeTotalEnergies();

  RandomNumber random(42);
  for (int move = 0; move < 1000; ++move)
  {
    std::optional<RunningEnergy> delta = MC_Moves::reactionMove_CFCMC(random, system);
    if (delta)
    {
      system.runningEnergies += delta.value();
    }
    ASSERT_TRUE(std::isfinite(system.runningEnergies.potentialEnergy()));

    // the active fractional molecule ids of every serial reaction must lie inside the fractional region
    for (const Reaction& reaction : system.reactions.list)
    {
      const std::vector<std::vector<std::size_t>>& activeIds = reaction.fractionalSideIsReactants
                                                                   ? reaction.reactantFractionalMoleculeIds
                                                                   : reaction.productFractionalMoleculeIds;
      for (std::size_t componentId = 0; componentId < activeIds.size(); ++componentId)
      {
        for (const std::size_t moleculeId : activeIds[componentId])
        {
          ASSERT_LT(moleculeId, system.numberOfFractionalMoleculesPerComponent[componentId]);
        }
      }
    }
  }
  system.checkMoleculeIds();
  checkEnergyDrift(system);
}

// serial and parallel Rx/CFC reactions can coexist in one system; each reaction declares the move
// that drives it and gets slots in the corresponding fractional-molecule region
TEST(MC_REACTION_DRIFT, reaction_mixed_parallel_serial_fractionals_init)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionConventionalCFCMC, 0.5);
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 0.5);

  System system = makeMixedParallelSerialReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  system.checkMoleculeIds();

  const Reaction& parallelReaction = system.reactions.list[0];
  const Reaction& serialReaction = system.reactions.list[1];
  EXPECT_TRUE(parallelReaction.isParallelRxCFC());
  EXPECT_TRUE(serialReaction.isSerialRxCFC());
  EXPECT_EQ(parallelReaction.reactantFractionalMoleculeIds[0].size(), 2uz);
  EXPECT_EQ(parallelReaction.productFractionalMoleculeIds[1].size(), 2uz);
  EXPECT_EQ(serialReaction.reactantFractionalMoleculeIds[0].size(), 2uz);
  EXPECT_TRUE(serialReaction.productFractionalMoleculeIds[0].empty());
  EXPECT_TRUE(serialReaction.productFractionalMoleculeIds[1].empty());
  EXPECT_EQ(system.numberOfReactionFractionalMoleculesPerComponent_CFCMC[0][0], 2uz);
  EXPECT_EQ(system.numberOfReactionFractionalMoleculesPerComponent_CFCMC[0][1], 2uz);
  EXPECT_EQ(system.numberOfReactionFractionalMoleculesPerComponent_CFCMC[1][0], 2uz);
  EXPECT_EQ(system.numberOfReactionFractionalMoleculesPerComponent_CFCMC[1][1], 2uz);

  system.runningEnergies = system.computeTotalEnergies();
  EXPECT_TRUE(std::isfinite(system.runningEnergies.potentialEnergy()));
}

TEST(MC_REACTION_DRIFT, reaction_parallel_two_reactions_init)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionConventionalCFCMC, 1.0);

  System system = makeTwoParallelReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  system.checkMoleculeIds();

  EXPECT_EQ(system.reactions.list.size(), 2uz);
  EXPECT_EQ(system.reactions.list[0].reactantFractionalMoleculeIds[0].size(), 2uz);
  EXPECT_EQ(system.reactions.list[1].reactantFractionalMoleculeIds[0].size(), 2uz);
  EXPECT_EQ(system.numberOfReactionFractionalMoleculesPerComponent_CFCMC[0][0], 2uz);
  EXPECT_EQ(system.numberOfReactionFractionalMoleculesPerComponent_CFCMC[1][0], 2uz);

  system.runningEnergies = system.computeTotalEnergies();
  EXPECT_TRUE(std::isfinite(system.runningEnergies.potentialEnergy()));
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_co2_n2_ewald_mixed)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makeCO2N2SerialReactionSystem(systemProbabilities, {2, 0}, {0, 2}, {40, 30});
  ASSERT_TRUE(system.forceField.useCharge);
  system.createReactionFractionalMolecules();
  setupSerialLambdaOnlyReaction(system);
  system.reactions.list[0].lambdaSwitchPoint = 0.5;
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_co2_n2_four_stoichiometry_single_reaction)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makeCO2N2SerialReactionSystem(systemProbabilities, {4, 0}, {0, 4}, {60, 20});
  system.createReactionFractionalMolecules();
  system.checkMoleculeIds();
  system.runningEnergies = system.computeTotalEnergies();
  EXPECT_TRUE(std::isfinite(system.runningEnergies.potentialEnergy()));
}

TEST(MC_REACTION_DRIFT, reaction_cbmc_zeolite)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCBMC, 1.0);

  runReactionDriftTest({makeZeolitePropaneButaneReactionSystem(systemProbabilities)});
}

TEST(MC_REACTION_DRIFT, reaction_cbmc_zeolite_boundary)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCBMC, 1.0);

  System system = makeZeoliteMultiComponentReactionSystem(systemProbabilities, {2, 0}, {0, 2}, {4, 2});
  setupCBMCStoichiometryBoundaryReaction(system);
  runReactionDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_cfcmc_zeolite)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionConventionalCFCMC, 1.0);

  System system = makeZeolitePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupLambdaOnlyReaction(system);
  runReactionDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_cfcmc_cbmc_zeolite)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionConventionalCBCFCMC, 1.0);

  System system = makeZeolitePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupLambdaOnlyReaction(system);
  runReactionDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_cfcmc_zeolite_boundary)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionConventionalCFCMC, 1.0);

  System system = makeZeolitePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupBoundaryCrossingReaction(system);
  runReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_cfcmc_cbmc_zeolite_boundary)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionConventionalCBCFCMC, 1.0);

  System system = makeZeolitePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupBoundaryCrossingReaction(system);
  runReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_zeolite_boundary)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makeZeolitePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupSerialBoundaryCrossingReaction(system);
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_cbmc_zeolite_boundary)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCBCFCMC, 1.0);

  System system = makeZeolitePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupSerialBoundaryCrossingReaction(system);
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_mixed_reactants_zeolite_boundary)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makeZeolitePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupSerialBoundaryCrossingReaction(system);
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_zeolite)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makeZeolitePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupSerialLambdaOnlyReaction(system);
  system.reactions.list[0].lambdaSwitchPoint = 0.5;
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_cbmc_zeolite)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCBCFCMC, 1.0);

  System system = makeZeolitePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupSerialLambdaOnlyReaction(system);
  system.reactions.list[0].lambdaSwitchPoint = 0.5;
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_mixed_reactants_zeolite)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makeZeolitePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupSerialLambdaOnlyReaction(system);
  system.reactions.list[0].lambdaSwitchPoint = 0.5;
  runSerialReactionBoundaryDriftTest({std::move(system)});
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_co2_n2_ewald_whole_molecule_single_move_drift)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makeCO2N2SerialReactionSystem(systemProbabilities, {2, 0}, {0, 2}, {60, 20});
  ASSERT_TRUE(system.forceField.useCharge);
  system.createReactionFractionalMolecules();
  setupSerialWholeMoleculeReaction(system);

  checkEnergyDrift(system);

  RandomNumber random(42);
  for (int attempt = 0; attempt < kReactionDriftSingleMoveAttempts; ++attempt)
  {
    RunningEnergy beforeRecomputed = system.computeTotalEnergies();
    std::optional<RunningEnergy> delta =
        MC_Moves::ReactionCommon::serialWholeMoleculeReactionMove(random, system, system.reactions.list[0],
                                                                  Move::Types::ReactionCFCMC, false);
    if (delta)
    {
      system.runningEnergies += delta.value();
      RunningEnergy afterRecomputed = system.computeTotalEnergies();
      RunningEnergy actualDelta = afterRecomputed - beforeRecomputed;
      EXPECT_NEAR(delta->potentialEnergy(), actualDelta.potentialEnergy(), 1e-6);
      EXPECT_NEAR(delta->ewald_exclusion, actualDelta.ewald_exclusion, 1e-6);
      checkEnergyDrift(system);
      return;
    }
  }
  FAIL() << "No accepted serial whole-molecule reaction move in " << kReactionDriftSingleMoveAttempts << " attempts";
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_co2_n2_ewald_fractional_single_move_drift)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makeCO2N2SerialReactionSystem(systemProbabilities, {2, 0}, {0, 2}, {60, 20});
  ASSERT_TRUE(system.forceField.useCharge);
  system.createReactionFractionalMolecules();
  setupSerialLambdaOnlyReaction(system);
  system.reactions.list[0].lambdaSwitchPoint = 0.5;

  checkEnergyDrift(system);

  RandomNumber random(42);
  for (int attempt = 0; attempt < kReactionDriftSingleMoveAttempts; ++attempt)
  {
    RunningEnergy beforeRecomputed = system.computeTotalEnergies();
    std::optional<RunningEnergy> delta =
        MC_Moves::ReactionCommon::serialFractionalReactionMove(random, system, system.reactions.list[0],
                                                               Move::Types::ReactionCFCMC, false);
    if (delta)
    {
      system.runningEnergies += delta.value();
      RunningEnergy afterRecomputed = system.computeTotalEnergies();
      RunningEnergy actualDelta = afterRecomputed - beforeRecomputed;
      EXPECT_NEAR(delta->potentialEnergy(), actualDelta.potentialEnergy(), 1e-6);
      EXPECT_NEAR(delta->ewald_exclusion, actualDelta.ewald_exclusion, 1e-6);
      EXPECT_NEAR(delta->ewald_fourier, actualDelta.ewald_fourier, 1e-6);
      checkEnergyDrift(system);
      return;
    }
  }
  FAIL() << "No accepted serial fractional reaction move in " << kReactionDriftSingleMoveAttempts << " attempts";
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_co2_n2_ewald_lambda_single_move_drift)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makeCO2N2SerialReactionSystem(systemProbabilities, {2, 0}, {0, 2}, {60, 20});
  ASSERT_TRUE(system.forceField.useCharge);
  system.createReactionFractionalMolecules();
  setupSerialLambdaOnlyReaction(system);
  system.reactions.list[0].lambdaSwitchPoint = 0.5;

  checkEnergyDrift(system);

  RandomNumber random(42);
  for (int attempt = 0; attempt < kReactionDriftSingleMoveAttempts; ++attempt)
  {
    RunningEnergy beforeRecomputed = system.computeTotalEnergies();
    std::optional<RunningEnergy> delta =
        MC_Moves::ReactionCommon::serialLambdaChangeMove(random, system, system.reactions.list[0],
                                                         Move::Types::ReactionCFCMC);
    if (delta)
    {
      system.runningEnergies += delta.value();
      RunningEnergy afterRecomputed = system.computeTotalEnergies();
      RunningEnergy actualDelta = afterRecomputed - beforeRecomputed;
      EXPECT_NEAR(delta->potentialEnergy(), actualDelta.potentialEnergy(), 1e-6);
      EXPECT_NEAR(delta->ewald_exclusion, actualDelta.ewald_exclusion, 1e-6);
      EXPECT_NEAR(delta->ewald_fourier, actualDelta.ewald_fourier, 1e-6);
      checkEnergyDrift(system);
      return;
    }
  }
  FAIL() << "No accepted serial lambda move in " << kReactionDriftSingleMoveAttempts << " attempts";
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_multi_stoichiometry_lambda_single_move_drift)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupSerialLambdaOnlyReaction(system);
  system.reactions.list[0].lambdaSwitchPoint = 0.5;

  checkEnergyDrift(system);

  RandomNumber random(42);
  for (int attempt = 0; attempt < kReactionDriftSingleMoveAttempts; ++attempt)
  {
    RunningEnergy beforeRecomputed = system.computeTotalEnergies();
    std::optional<RunningEnergy> delta =
        MC_Moves::ReactionCommon::serialLambdaChangeMove(random, system, system.reactions.list[0],
                                                         Move::Types::ReactionCFCMC);
    if (delta)
    {
      system.runningEnergies += delta.value();
      RunningEnergy afterRecomputed = system.computeTotalEnergies();
      RunningEnergy actualDelta = afterRecomputed - beforeRecomputed;
      EXPECT_NEAR(delta->potentialEnergy(), actualDelta.potentialEnergy(), 1e-6);
      EXPECT_NEAR(delta->tail, actualDelta.tail, 1e-6);
      checkEnergyDrift(system);
      return;
    }
  }
  FAIL() << "No accepted serial lambda move in " << kReactionDriftSingleMoveAttempts << " attempts";
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_multi_stoichiometry_single_move_drift)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupSerialLambdaOnlyReaction(system);
  system.reactions.list[0].lambdaSwitchPoint = 0.5;

  checkEnergyDrift(system);

  RandomNumber random(42);
  for (int attempt = 0; attempt < kReactionDriftSingleMoveAttempts; ++attempt)
  {
    RunningEnergy beforeRecomputed = system.computeTotalEnergies();
    std::optional<RunningEnergy> delta =
        MC_Moves::ReactionCommon::serialFractionalReactionMove(random, system, system.reactions.list[0],
                                                               Move::Types::ReactionCFCMC, false);
    if (delta)
    {
      system.runningEnergies += delta.value();
      RunningEnergy afterRecomputed = system.computeTotalEnergies();
      RunningEnergy actualDelta = afterRecomputed - beforeRecomputed;
      EXPECT_NEAR(delta->potentialEnergy(), actualDelta.potentialEnergy(), 1e-6);
      EXPECT_NEAR(delta->tail, actualDelta.tail, 1e-6);
      checkEnergyDrift(system);
      return;
    }
  }
  FAIL() << "No accepted serial fractional reaction move in " << kReactionDriftSingleMoveAttempts << " attempts";
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_fractional_single_move_drift)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupSerialLambdaOnlyReaction(system);
  system.reactions.list[0].lambdaSwitchPoint = 0.5;

  checkEnergyDrift(system);

  RandomNumber random(42);
  for (int attempt = 0; attempt < kReactionDriftSingleMoveAttempts; ++attempt)
  {
    RunningEnergy beforeRecomputed = system.computeTotalEnergies();
    std::optional<RunningEnergy> delta =
        MC_Moves::ReactionCommon::serialFractionalReactionMove(random, system, system.reactions.list[0],
                                                               Move::Types::ReactionCFCMC, false);
    if (delta)
    {
      system.runningEnergies += delta.value();
      RunningEnergy afterRecomputed = system.computeTotalEnergies();
      RunningEnergy actualDelta = afterRecomputed - beforeRecomputed;
      EXPECT_NEAR(delta->potentialEnergy(), actualDelta.potentialEnergy(), 1e-6);
      EXPECT_NEAR(delta->moleculeMoleculeVDW, actualDelta.moleculeMoleculeVDW, 1e-6);
      checkEnergyDrift(system);
      return;
    }
  }
  FAIL() << "No accepted serial fractional reaction move in " << kReactionDriftSingleMoveAttempts << " attempts";
}

TEST(MC_REACTION_DRIFT, reaction_serial_cfcmc_lambda_single_move_drift)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = makePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupSerialLambdaOnlyReaction(system);

  checkEnergyDrift(system);

  RandomNumber random(42);
  for (int attempt = 0; attempt < kReactionDriftSingleMoveAttempts; ++attempt)
  {
    RunningEnergy beforeRecomputed = system.computeTotalEnergies();
    std::optional<RunningEnergy> delta =
        MC_Moves::ReactionCommon::serialLambdaChangeMove(random, system, system.reactions.list[0],
                                                         Move::Types::ReactionCFCMC);
    if (delta)
    {
      system.runningEnergies += delta.value();
      RunningEnergy afterRecomputed = system.computeTotalEnergies();
      RunningEnergy actualDelta = afterRecomputed - beforeRecomputed;
      EXPECT_NEAR(delta->potentialEnergy(), actualDelta.potentialEnergy(), 1e-6);
      checkEnergyDrift(system);
      return;
    }
  }
  FAIL() << "No accepted serial lambda move in " << kReactionDriftSingleMoveAttempts << " attempts";
}

namespace {

// deterministic per-accept drift check for the conventional (parallel) reaction CFCMC moves,
// including the lambda boundary crossings (ForwardInsert / BackwardDelete)
void runConventionalBoundaryPerMoveDrift(bool useCBMC)
{
  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(
      useCBMC ? Move::Types::ReactionConventionalCBCFCMC : Move::Types::ReactionConventionalCFCMC, 1.0);

  System system = makePropaneButaneReactionSystem(systemProbabilities);
  system.createReactionFractionalMolecules();
  setupBoundaryCrossingReaction(system);

  checkEnergyDrift(system);

  RandomNumber random(42);
  int accepted = 0;
  for (int attempt = 0; attempt < 600; ++attempt)
  {
    const double lambdaBefore = system.reactions.list[0].currentLambda;
    RunningEnergy beforeRecomputed = system.computeTotalEnergies();
    std::optional<RunningEnergy> delta =
        useCBMC ? MC_Moves::reactionMove_ConventionalCBCFCMC(random, system)
                : MC_Moves::reactionMove_ConventionalCFCMC(random, system);
    if (delta)
    {
      ++accepted;
      SCOPED_TRACE(testing::Message() << "attempt " << attempt << " accept " << accepted << " lambda "
                                      << lambdaBefore << " -> " << system.reactions.list[0].currentLambda);
      system.runningEnergies += delta.value();
      RunningEnergy afterRecomputed = system.computeTotalEnergies();
      RunningEnergy actualDelta = afterRecomputed - beforeRecomputed;
      ASSERT_NEAR(delta->potentialEnergy(), actualDelta.potentialEnergy(), 1e-6);
      ASSERT_NEAR(delta->moleculeMoleculeVDW, actualDelta.moleculeMoleculeVDW, 1e-6);
      ASSERT_NEAR(delta->tail, actualDelta.tail, 1e-6);
      ASSERT_NEAR(delta->intraVDW, actualDelta.intraVDW, 1e-6);
    }
  }
  EXPECT_GT(accepted, 0);
  checkEnergyDrift(system);
}

}  // namespace

TEST(MC_REACTION_DRIFT, reaction_cfcmc_boundary_per_move_drift)
{
  runConventionalBoundaryPerMoveDrift(false);
}

TEST(MC_REACTION_DRIFT, reaction_cfcmc_cbmc_boundary_per_move_drift)
{
  runConventionalBoundaryPerMoveDrift(true);
}
