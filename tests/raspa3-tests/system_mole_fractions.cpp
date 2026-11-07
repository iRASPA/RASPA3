#include <gtest/gtest.h>

import std;

import atom;
import forcefield;
import component;
import system;
import simulationbox;
import mc_moves_move_types;
import mc_moves_probabilities;

// The mole fractions of the swappable components are normalized when the system is built. Given fractions are
// scaled to sum to one; with none given (all zero) the swappable components share the gas phase equally, and
// a component that is not swappable keeps whatever it has.

namespace
{

Component makeComponent(const ForceField& forceField, std::string name, bool swappable, double molFraction)
{
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  if (swappable) probabilities.setProbability(Move::Types::Swap, 1.0);
  Component component(forceField, std::move(name), 190.564, 4599200, 0.01142,
                      {Atom({0, 0, 0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21, probabilities, std::nullopt,
                      false);
  component.molFraction = molFraction;
  return component;
}

System makeSystem(std::vector<Component> components)
{
  ForceField shared({{"CH4", false, 16.04246, 0.0, 0.0, 6, false}}, {{158.5, 3.72}},
                    ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  const std::vector<std::size_t> counts(components.size(), 2uz);
  return System(shared, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, std::move(components), {}, counts,
                5);
}

ForceField makeForceField()
{
  return ForceField({{"CH4", false, 16.04246, 0.0, 0.0, 6, false}}, {{158.5, 3.72}},
                    ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
}

}  // namespace

TEST(SYSTEM_MOLE_FRACTIONS, given_fractions_are_normalized_to_one)
{
  const ForceField forceField = makeForceField();
  const System system =
      makeSystem({makeComponent(forceField, "a", true, 3.0), makeComponent(forceField, "b", true, 1.0)});
  ASSERT_TRUE(system.components[0].swappable);
  ASSERT_TRUE(system.components[1].swappable);
  EXPECT_DOUBLE_EQ(system.components[0].molFraction, 0.75);
  EXPECT_DOUBLE_EQ(system.components[1].molFraction, 0.25);
}

TEST(SYSTEM_MOLE_FRACTIONS, no_fractions_given_shares_the_gas_phase_equally)
{
  const ForceField forceField = makeForceField();
  const System system =
      makeSystem({makeComponent(forceField, "a", true, 0.0), makeComponent(forceField, "b", true, 0.0)});
  EXPECT_DOUBLE_EQ(system.components[0].molFraction, 0.5);
  EXPECT_DOUBLE_EQ(system.components[1].molFraction, 0.5);

  const System single = makeSystem({makeComponent(forceField, "a", true, 0.0)});
  EXPECT_DOUBLE_EQ(single.components[0].molFraction, 1.0);
}

TEST(SYSTEM_MOLE_FRACTIONS, a_component_that_is_not_swappable_is_left_alone)
{
  const ForceField forceField = makeForceField();
  const System system =
      makeSystem({makeComponent(forceField, "a", true, 0.0), makeComponent(forceField, "fixed", false, 0.0)});
  ASSERT_FALSE(system.components[1].swappable);
  EXPECT_DOUBLE_EQ(system.components[0].molFraction, 1.0);
  EXPECT_DOUBLE_EQ(system.components[1].molFraction, 0.0);

  // nothing swappable at all: no division by the zero count, the fractions stay as given
  const System none = makeSystem({makeComponent(forceField, "fixed", false, 0.0)});
  EXPECT_DOUBLE_EQ(none.components[0].molFraction, 0.0);
}

// A component whose count another component's move changes is swappable too: the partner of a pair swap
// configured on one ion only, and the target of an identity change. Its mole fraction is then normalized
// with the others (the pair fugacity product and the identity-change fugacity ratio use it) and it is in
// the list of swappable components (enthalpy and partial molar properties).
TEST(SYSTEM_MOLE_FRACTIONS, the_partner_of_a_pair_swap_on_one_ion_is_swappable_too)
{
  const ForceField forceField = makeForceField();
  Component cation = makeComponent(forceField, "cation", false, 1.0);
  cation.mc_moves_probabilities.setProbability(Move::Types::PairSwapCBMC, 1.0);
  cation.pairComponentId = 1;
  cation.maximumPairDistance = 5.0;
  const Component anion = makeComponent(forceField, "anion", false, 1.0);

  const System system = makeSystem({cation, anion});
  EXPECT_TRUE(system.components[0].swappable);
  EXPECT_TRUE(system.components[1].swappable);
  // the pair is one species: the salt alone has the whole gas phase, and both ions carry its fraction (the
  // pair moves use f_A f_B); the partner does not count in the normalization
  EXPECT_DOUBLE_EQ(system.components[0].molFraction, 1.0);
  EXPECT_DOUBLE_EQ(system.components[1].molFraction, 1.0);
  EXPECT_EQ(system.swappableComponents, (std::vector<std::size_t>{0, 1}));

  // beside another adsorbate the salt and that adsorbate share the gas phase, and the partner follows its driver
  Component carbonDioxide = makeComponent(forceField, "co2", true, 1.0);
  const System mixture = makeSystem({cation, anion, carbonDioxide});
  EXPECT_DOUBLE_EQ(mixture.components[0].molFraction, 0.5);
  EXPECT_DOUBLE_EQ(mixture.components[1].molFraction, 0.5);
  EXPECT_DOUBLE_EQ(mixture.components[2].molFraction, 0.5);
}

TEST(SYSTEM_MOLE_FRACTIONS, the_target_of_an_identity_change_is_swappable_too)
{
  const ForceField forceField = makeForceField();
  Component source = makeComponent(forceField, "source", false, 1.0);
  source.mc_moves_probabilities.setProbability(Move::Types::IdentityChangeCBMC, 1.0);
  source.identityChanges = {1};
  const Component target = makeComponent(forceField, "target", false, 1.0);
  const Component bystander = makeComponent(forceField, "bystander", false, 1.0);

  const System system = makeSystem({source, target, bystander});
  EXPECT_TRUE(system.components[0].swappable);
  EXPECT_TRUE(system.components[1].swappable);
  EXPECT_FALSE(system.components[2].swappable);
  EXPECT_DOUBLE_EQ(system.components[0].molFraction, 0.5);
  EXPECT_DOUBLE_EQ(system.components[1].molFraction, 0.5);
  EXPECT_DOUBLE_EQ(system.components[2].molFraction, 1.0);
  EXPECT_EQ(system.swappableComponents, (std::vector<std::size_t>{0, 1}));
}
