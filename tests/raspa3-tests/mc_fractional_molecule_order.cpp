#include <gtest/gtest.h>

import std;

import double3;
import forcefield;
import component;
import system;
import monte_carlo;
import simulationbox;
import running_energy;
import reaction;
import mc_moves_move_types;
import mc_moves_probabilities;

namespace
{

std::filesystem::path repositoryRoot()
{
  return std::filesystem::path(__FILE__).parent_path().parent_path().parent_path();
}

ForceField makeAlkaneForceField()
{
  return ForceField({{"CH3", false, 15.04, 0.0, 0.0, 6, false},
                     {"CH2", false, 14.03, 0.0, 0.0, 6, false}},
                    {{98.0, 3.75}, {46.0, 3.95}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true,
                     true, false);
}

Component makeAlkaneFromExample(const ForceField& forceField, std::size_t componentId, std::string_view name,
                                const MCMoveProbabilities& probabilities)
{
  const std::filesystem::path moleculePath =
      repositoryRoot() / "examples/basic/4_mc_binary_mixture_propane_butane_in_box" / name;
  return Component(Component::Type::Adsorbate, componentId, forceField, std::string(name), moleculePath.string(), 5,
                   21, probabilities, std::nullopt, false);
}

void expectFractionalBeforeInteger(const System& system)
{
  for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
  {
    EXPECT_EQ(system.numberOfMoleculesPerComponent[componentId],
              system.numberOfFractionalMoleculesPerComponent[componentId] +
                  system.numberOfIntegerMoleculesPerComponent[componentId]);
    for (std::size_t moleculeId = 0; moleculeId < system.numberOfFractionalMoleculesPerComponent[componentId];
         ++moleculeId)
    {
      EXPECT_TRUE(system.spanOfMolecule(componentId, moleculeId).front().isFractional);
    }
    for (std::size_t moleculeId = system.numberOfFractionalMoleculesPerComponent[componentId];
         moleculeId < system.numberOfMoleculesPerComponent[componentId]; ++moleculeId)
    {
      EXPECT_FALSE(system.spanOfMolecule(componentId, moleculeId).front().isFractional);
    }
  }
}

TEST(MC_FRACTIONAL_MOLECULE_ORDER, swap_and_serial_reaction_layout)
{
  const ForceField forceField = makeAlkaneForceField();
  SimulationBox box = SimulationBox(30.0, 30.0, 30.0);

  MCMoveProbabilities componentProbabilities = MCMoveProbabilities();
  componentProbabilities.setProbability(Move::Types::Translation, 1.0);
  componentProbabilities.setProbability(Move::Types::Rotation, 1.0);
  componentProbabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  componentProbabilities.setProbability(Move::Types::SwapCFCMC, 1.0);

  Component propane = makeAlkaneFromExample(forceField, 0, "propane", componentProbabilities);
  Component butane = makeAlkaneFromExample(forceField, 1, "butane", componentProbabilities);
  propane.lnPartitionFunction = 1.0;
  butane.lnPartitionFunction = 1.0;

  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system = System(forceField, box, false, 500.0, 1e4, 1.0, {}, {propane, butane}, {},
                         std::vector<std::size_t>{20, 10}, 5, systemProbabilities);

  Reaction reaction(0, {2, 0}, {0, 2});
  reaction.reactionMove = Move::Types::ReactionCFCMC;
  system.reactions.list.push_back(reaction);
  system.createReactionFractionalMolecules();
  system.checkMoleculeIds();

  const std::size_t gcIndex = system.indexOfFractionalMoleculeForMove(Move::Types::SwapCFCMC, 0);
  const std::size_t serialIndex0 = system.serialReactionFractionalMoleculeIndex(0, 0, 0);
  const std::size_t serialIndex1 = system.serialReactionFractionalMoleculeIndex(0, 0, 1);

  EXPECT_EQ(gcIndex, 0uz);
  EXPECT_EQ(serialIndex0, 1uz);
  EXPECT_EQ(serialIndex1, 2uz);
  EXPECT_EQ(system.reactions.list[0].reactantFractionalMoleculeIds[0][0], 1uz);
  EXPECT_EQ(system.reactions.list[0].reactantFractionalMoleculeIds[0][1], 2uz);
  expectFractionalBeforeInteger(system);
}

TEST(MC_FRACTIONAL_MOLECULE_ORDER, parallel_reaction_and_swap_layout)
{
  const ForceField forceField = makeAlkaneForceField();
  SimulationBox box = SimulationBox(30.0, 30.0, 30.0);

  MCMoveProbabilities componentProbabilities = MCMoveProbabilities();
  componentProbabilities.setProbability(Move::Types::Translation, 1.0);
  componentProbabilities.setProbability(Move::Types::Rotation, 1.0);
  componentProbabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  componentProbabilities.setProbability(Move::Types::SwapCFCMC, 1.0);

  Component propane = makeAlkaneFromExample(forceField, 0, "propane", componentProbabilities);
  Component butane = makeAlkaneFromExample(forceField, 1, "butane", componentProbabilities);
  propane.lnPartitionFunction = 1.0;
  butane.lnPartitionFunction = 1.0;

  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionConventionalCFCMC, 1.0);

  System system = System(forceField, box, false, 500.0, 1e4, 1.0, {}, {propane, butane}, {},
                         std::vector<std::size_t>{20, 10}, 5, systemProbabilities);

  Reaction reaction(0, {2, 0}, {0, 2});
  reaction.reactionMove = Move::Types::ReactionConventionalCFCMC;
  system.reactions.list.push_back(reaction);
  system.createReactionFractionalMolecules();
  system.checkMoleculeIds();

  EXPECT_EQ(system.indexOfFractionalMoleculeForMove(Move::Types::SwapCFCMC, 0), 0uz);
  EXPECT_EQ(system.parallelReactionFractionalMoleculeIndex(0, 0, false, 0), 1uz);
  EXPECT_EQ(system.parallelReactionFractionalMoleculeIndex(0, 0, false, 1), 2uz);
  EXPECT_EQ(system.parallelReactionFractionalMoleculeIndex(0, 1, true, 0), 1uz);
  EXPECT_EQ(system.parallelReactionFractionalMoleculeIndex(0, 1, true, 1), 2uz);
  EXPECT_EQ(system.reactions.list[0].reactantFractionalMoleculeIds[0][0], 1uz);
  EXPECT_EQ(system.reactions.list[0].productFractionalMoleculeIds[1][0], 1uz);
  expectFractionalBeforeInteger(system);
}

TEST(MC_FRACTIONAL_MOLECULE_ORDER, combined_moves_short_drift)
{
  const ForceField forceField = makeAlkaneForceField();
  SimulationBox box = SimulationBox(40.0, 40.0, 40.0);

  MCMoveProbabilities componentProbabilities = MCMoveProbabilities();
  componentProbabilities.setProbability(Move::Types::Translation, 1.0);
  componentProbabilities.setProbability(Move::Types::Rotation, 1.0);
  componentProbabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  componentProbabilities.setProbability(Move::Types::SwapCFCMC, 1.0);

  Component propane = makeAlkaneFromExample(forceField, 0, "propane", componentProbabilities);
  Component butane = makeAlkaneFromExample(forceField, 1, "butane", componentProbabilities);
  propane.lnPartitionFunction = 1.0;
  butane.lnPartitionFunction = 1.0;

  MCMoveProbabilities systemProbabilities = MCMoveProbabilities();
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 0.5);
  systemProbabilities.setProbability(Move::Types::ReactionConventionalCFCMC, 0.5);

  System system = System(forceField, box, false, 500.0, 1e4, 1.0, {}, {propane, butane}, {},
                         std::vector<std::size_t>{20, 10}, 5, systemProbabilities);

  Reaction serialReaction(0, {2, 0}, {0, 2});
  serialReaction.reactionMove = Move::Types::ReactionCFCMC;
  Reaction parallelReaction(1, {2, 0}, {0, 2});
  parallelReaction.reactionMove = Move::Types::ReactionConventionalCFCMC;
  system.reactions.list.push_back(serialReaction);
  system.reactions.list.push_back(parallelReaction);
  system.createReactionFractionalMolecules();
  system.checkMoleculeIds();
  system.runningEnergies = system.computeTotalEnergies();
  ASSERT_TRUE(std::isfinite(system.runningEnergies.potentialEnergy()));

  const size_t numberOfProductionCycles = 20;
  const size_t numberOfInitializationCycles = 5;
  const size_t numberOfEquilibrationCycles = 5;
  MonteCarlo mc = MonteCarlo({numberOfProductionCycles, 0, numberOfInitializationCycles, numberOfEquilibrationCycles, 100, 10000,
                              5000, 5000}, {std::move(system)}, {}, 5, false);
  mc.run();

  for (System& s : mc.systems)
  {
    s.checkMoleculeIds();
    expectFractionalBeforeInteger(s);
    const RunningEnergy recomputed = s.computeTotalEnergies();
    EXPECT_TRUE(std::isfinite(s.runningEnergies.potentialEnergy()));
    EXPECT_TRUE(std::isfinite(recomputed.potentialEnergy()));
  }
}

TEST(MC_FRACTIONAL_MOLECULE_ORDER, swap_and_gibbs_swap_layout)
{
  const ForceField forceField = makeAlkaneForceField();
  SimulationBox box = SimulationBox(30.0, 30.0, 30.0);

  MCMoveProbabilities componentProbabilities = MCMoveProbabilities();
  componentProbabilities.setProbability(Move::Types::Translation, 1.0);
  componentProbabilities.setProbability(Move::Types::Rotation, 1.0);
  componentProbabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  componentProbabilities.setProbability(Move::Types::SwapCFCMC, 1.0);
  componentProbabilities.setProbability(Move::Types::GibbsSwapCFCMC, 1.0);

  Component propane = makeAlkaneFromExample(forceField, 0, "propane", componentProbabilities);
  propane.lnPartitionFunction = 1.0;

  System system = System(forceField, box, false, 500.0, 1e4, 1.0, {}, {propane}, {}, std::vector<std::size_t>{10},
                         5, MCMoveProbabilities());

  const std::size_t gcIndex = system.indexOfFractionalMoleculeForMove(Move::Types::SwapCFCMC, 0);
  const std::size_t gibbsSwapIndex = system.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCFCMC, 0);

  EXPECT_EQ(system.numberOfGCFractionalMoleculesPerComponent_CFCMC[0], 1uz);
  EXPECT_EQ(system.numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC[0], 1uz);
  EXPECT_EQ(gcIndex, 0uz);
  EXPECT_EQ(gibbsSwapIndex, 1uz);
  EXPECT_NE(gcIndex, gibbsSwapIndex);
  expectFractionalBeforeInteger(system);
}

}  // namespace
