#include <gtest/gtest.h>

#include "../test_support.hpp"
#include "molecule_fixtures.hpp"

import std;

import atom;
import forcefield;
import component;
import system;
import simulationbox;
import running_energy;
import reaction;
import mc_moves_move_types;
import mc_moves_probabilities;

// Serial Rx/CFC keeps the fractional molecules of one side of the reaction at a time, while the layout
// reserves the larger side per component. The region after them (the Gibbs-conventional slot) has to start
// where the present molecules end, and the molecules themselves, placed decoupled one after another, must not
// be grown on top of each other.

namespace
{

ForceField makeAlkaneForceField()
{
  return ForceField({{"CH3", false, 15.04, 0.0, 0.0, 6, false}, {"CH2", false, 14.03, 0.0, 0.0, 6, false}},
                    {{98.0, 3.75}, {46.0, 3.95}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true,
                    true, false);
}

Component makeAlkane(const ForceField& forceField, std::size_t componentId, std::string_view name,
                     const MCMoveProbabilities& probabilities)
{
  TemporaryFile file(std::string(name) + ".json", molecule_fixtures::alkaneJson(name));
  Component component(Component::Type::Adsorbate, componentId, forceField, std::string(name), file.stemPath().string(),
                      5, 21, probabilities, std::nullopt, false);
  component.lnPartitionFunction = 1.0;
  return component;
}

System makeSerialReactionSystem(std::size_t reactantStoichiometry, std::size_t propanes, std::size_t butanes,
                                double boxLength, bool withGCSwaps)
{
  const ForceField forceField = makeAlkaneForceField();
  MCMoveProbabilities componentProbabilities;
  componentProbabilities.setProbability(Move::Types::Translation, 1.0);
  if (withGCSwaps) componentProbabilities.setProbability(Move::Types::SwapCFCMC, 1.0);
  MCMoveProbabilities systemProbabilities;
  systemProbabilities.setProbability(Move::Types::ReactionCFCMC, 1.0);

  System system(forceField, SimulationBox(boxLength, boxLength, boxLength), false, 500.0, 1e4, 1.0, {},
                {makeAlkane(forceField, 0, "propane", componentProbabilities),
                 makeAlkane(forceField, 1, "butane", componentProbabilities)},
                {}, std::vector<std::size_t>{propanes, butanes}, 5, systemProbabilities);
  Reaction reaction(0, {reactantStoichiometry, 0}, {0, reactantStoichiometry});
  reaction.reactionMove = Move::Types::ReactionCFCMC;
  system.reactions.list.push_back(reaction);
  system.createReactionFractionalMolecules();
  system.checkMoleculeIds();
  return system;
}

}  // namespace

TEST(SERIAL_REACTION_FRACTIONAL, the_region_after_the_serial_molecules_starts_where_the_present_ones_end)
{
  // reactants active: two propane fractional molecules present, none of butane although its capacity is two
  const System system = makeSerialReactionSystem(2, 20, 10, 30.0, true);
  ASSERT_EQ(system.numberOfSerialReactionFractionalMoleculesPerComponent_CFCMC[1], 2uz);
  EXPECT_EQ(system.numberOfPresentSerialReactionFractionalMolecules(0), 2uz);
  EXPECT_EQ(system.numberOfPresentSerialReactionFractionalMolecules(1), 0uz);

  // the GC slot, then the present serial molecules, then (empty) the Gibbs-conventional region: its start is
  // the fractional / integer boundary of the component
  for (std::size_t componentId = 0; componentId < 2; ++componentId)
  {
    EXPECT_EQ(system.indexOfGibbsConventionalFractionalMoleculesPerComponent_CFCMC(componentId),
              system.numberOfFractionalMoleculesPerComponent[componentId])
        << "component " << componentId;
  }
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 3uz);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 1uz);
}

TEST(SERIAL_REACTION_FRACTIONAL, decoupled_molecules_are_not_placed_on_top_of_each_other)
{
  // forty propane fractional molecules placed one after another at lambda = 0 in a small box: each one was
  // invisible to the overlap check of the next, so pairs landed on top of each other
  System system = makeSerialReactionSystem(40, 0, 0, 25.0, false);
  ASSERT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 40uz);

  // couple them all fully: the configuration must be free of overlaps
  for (std::size_t moleculeId = 0; moleculeId < 40; ++moleculeId)
  {
    for (Atom& atom : system.spanOfMolecule(0, moleculeId)) atom.setScaling(1.0);
  }
  const RunningEnergy energies = system.computeTotalEnergies();
  EXPECT_TRUE(std::isfinite(energies.potentialEnergy()));
  EXPECT_LT(energies.potentialEnergy(), system.forceField.energyOverlapCriteria);
}
