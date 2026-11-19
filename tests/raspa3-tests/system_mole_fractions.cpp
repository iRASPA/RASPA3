#include <gtest/gtest.h>

import std;

import atom;
import forcefield;
import component;
import system;
import simulationbox;
import mc_moves_move_types;
import mc_moves_probabilities;
import randomnumbers;
import property_lambda_probability_histogram;

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

// One anion component as the partner of two cations: it can carry the mole fraction of one salt only (it took
// that of the driver listed last), and with the CFCMC pair moves the two pairs would rescale the same fractional
// molecule at two lambdas. Refused: each pair takes a component of its own.
TEST(SYSTEM_MOLE_FRACTIONS, a_partner_shared_by_two_pair_swaps_is_refused)
{
  const ForceField forceField = makeForceField();
  for (const Move::Types move : {Move::Types::PairSwapCBMC, Move::Types::PairSwapCFCMC, Move::Types::PairSwapCBCFCMC})
  {
    Component sodium = makeComponent(forceField, "sodium", false, 1.0);
    sodium.mc_moves_probabilities.setProbability(move, 1.0);
    sodium.pairComponentId = 2;
    sodium.maximumPairDistance = 5.0;
    Component potassium = makeComponent(forceField, "potassium", false, 1.0);
    potassium.mc_moves_probabilities.setProbability(move, 1.0);
    potassium.pairComponentId = 2;
    potassium.maximumPairDistance = 5.0;
    const Component chloride = makeComponent(forceField, "chloride", false, 1.0);
    EXPECT_THROW((void)makeSystem({sodium, potassium, chloride}), std::runtime_error);
  }
}

// The CFCMC pair swap named on both ions of the pair is one pair, moved by its lower-index member: the second
// mention used to claim the fractional molecules as a second driver and was refused as two pairs sharing them.
// Named on the higher-index member alone the move never ran: that is refused, by name.
TEST(SYSTEM_MOLE_FRACTIONS, a_cfcmc_pair_swap_named_on_both_ions_is_one_pair)
{
  const ForceField forceField = makeForceField();
  for (const Move::Types move : {Move::Types::PairSwapCFCMC, Move::Types::PairSwapCBCFCMC})
  {
    Component cation = makeComponent(forceField, "cation", false, 1.0);
    cation.mc_moves_probabilities.setProbability(move, 1.0);
    cation.pairComponentId = 1;
    cation.maximumPairDistance = 5.0;
    Component anion = makeComponent(forceField, "anion", false, 1.0);
    anion.mc_moves_probabilities.setProbability(move, 1.0);
    anion.pairComponentId = 0;
    anion.maximumPairDistance = 5.0;

    EXPECT_NO_THROW((void)makeSystem({cation, anion}));
    const System system = makeSystem({cation, anion});
    EXPECT_TRUE(system.componentDrivesPairSwapLambda(0, move));
    EXPECT_FALSE(system.componentDrivesPairSwapLambda(1, move));

    Component silentCation = makeComponent(forceField, "cation", false, 1.0);
    EXPECT_THROW((void)makeSystem({silentCation, anion}), std::runtime_error);
  }
}

// The grand-canonical CFCMC swap and the Gibbs CFCMC swap both move along the one lambda of the component, each
// with a fractional molecule of its own: after either changed the bin the other's molecule kept the old scaling.
TEST(SYSTEM_MOLE_FRACTIONS, a_grand_canonical_and_a_gibbs_cfcmc_swap_on_one_component_are_refused)
{
  const ForceField forceField = makeForceField();
  Component both = makeComponent(forceField, "co2", true, 1.0);
  both.mc_moves_probabilities.setProbability(Move::Types::SwapCFCMC, 1.0);
  both.mc_moves_probabilities.setProbability(Move::Types::GibbsSwapCFCMC, 1.0);
  EXPECT_THROW((void)makeSystem({both}), std::runtime_error);

  Component gibbsOnly = makeComponent(forceField, "co2", true, 1.0);
  gibbsOnly.mc_moves_probabilities.setProbability(Move::Types::GibbsSwapCFCMC, 1.0);
  EXPECT_NO_THROW((void)makeSystem({gibbsOnly}));
}

// The pair sphere is wrapped onto the box by the insertion and looked at by minimum image in the deletion: with a
// radius above half the smallest width of the box the two count different sets of partners. Refused at the start,
// and the moves are not performed in a box that has shrunk below it.
TEST(SYSTEM_MOLE_FRACTIONS, a_pair_sphere_wider_than_the_box_is_refused)
{
  const ForceField forceField = makeForceField();
  Component cation = makeComponent(forceField, "cation", false, 1.0);
  cation.mc_moves_probabilities.setProbability(Move::Types::PairSwapCBMC, 1.0);
  cation.pairComponentId = 1;
  const Component anion = makeComponent(forceField, "anion", false, 1.0);

  cation.maximumPairDistance = 15.0;  // half the 30 Angstrom box: the largest that fits
  const System fits = makeSystem({cation, anion});
  EXPECT_TRUE(fits.pairSphereFitsBox(15.0));
  EXPECT_FALSE(fits.pairSphereFitsBox(15.1));

  cation.maximumPairDistance = 16.0;
  EXPECT_THROW((void)makeSystem({cation, anion}), std::runtime_error);
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

// The moves of the system are drawn through the tables of the components. Joined into the table of a component
// without moves of its own (fixed extra-framework ions), they made it a component "with moves": it was picked as
// often as a mobile one and then made a system move every time. They go into the tables of the components that
// have moves; into every table only when no component has one.
TEST(SYSTEM_MOLE_FRACTIONS, system_moves_are_joined_into_the_components_that_have_moves)
{
  const ForceField forceField = makeForceField();
  MCMoveProbabilities systemMoves;
  systemMoves.setProbability(Move::Types::VolumeChange, 0.1);
  const auto immobile = [&](std::string name)
  {
    return Component(forceField, std::move(name), 190.564, 4599200, 0.01142,
                     {Atom({0, 0, 0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21, MCMoveProbabilities{},
                     std::nullopt, false);
  };
  const auto volumeProbability = [](const Component& component)
  { return component.mc_moves_probabilities.getProbability(Move::Types::VolumeChange); };

  System mixed(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {},
               {makeComponent(forceField, "mobile", false, 1.0), immobile("ions")}, {}, {2, 2}, 5, systemMoves);
  EXPECT_DOUBLE_EQ(volumeProbability(mixed.components[0]), 0.1);
  EXPECT_DOUBLE_EQ(volumeProbability(mixed.components[1]), 0.0);
  RandomNumber random(5);
  for (std::size_t i = 0; i != 200; ++i) EXPECT_EQ(mixed.randomComponent(random), 0uz);

  // no component has a move of its own: the system moves are still drawn
  System fixed(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {},
               {immobile("ions"), immobile("more ions")}, {}, {2, 2}, 5, systemMoves);
  EXPECT_DOUBLE_EQ(volumeProbability(fixed.components[0]), 0.1);
  EXPECT_DOUBLE_EQ(volumeProbability(fixed.components[1]), 0.1);
}

// The Gibbs lambda histogram of a component holds a molecule only when the component has the conventional Gibbs
// CFCMC move. Its bias entered the weight of a configuration for every component with a fractional molecule of any
// kind, and the Wang-Landau sampling ran on it whenever some component of the system had the move.
TEST(SYSTEM_MOLE_FRACTIONS, the_gibbs_lambda_bias_of_a_component_without_the_move_is_not_a_weight)
{
  const ForceField forceField = makeForceField();
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::SwapCFCMC, 1.0);
  Component component(forceField, "methane", 190.564, 4599200, 0.01142,
                      {Atom({0, 0, 0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21, probabilities, std::nullopt,
                      false);
  ASSERT_FALSE(component.usesGibbsConventionalCFCMC());
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {component}, {}, {4}, 5);
  ASSERT_TRUE(system.components[0].hasFractionalMolecule);

  const double weight = system.weight();
  std::ranges::fill(system.components[0].lambdaGibbs.biasFactor, -800.0);
  EXPECT_DOUBLE_EQ(system.weight(), weight);
  EXPECT_TRUE(std::isfinite(system.weight()));
}
