#include <gtest/gtest.h>

import std;

import double3;
import atom;
import forcefield;
import component;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;
import bend_potential;
import urey_bradley_potential;
import chiral_center;
import cbmc_grow_step;
import cbmc_growth_plan;
import cbmc_rigid_tilt;
import randomnumbers;

namespace
{
// A rigid fragment {1, 2, 3} hinged on atom 1, which is a declared chiral centre with the previous bead 0, the
// fragment atom 2 and a side group 4. The bend 0-1-2 is the junction bend (spin invariant); 4-1-2 couples the
// fragment to the side group off the spin axis (spin variant).
Component makeRigidFragmentAtChiralCentre(BendType junctionBendType = BendType::Harmonic)
{
  const ForceField forceField({{"C", false, 12.0, 0.0, 0.0, 6, false}}, {{50.0, 3.5}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, false, false, false);
  ConnectivityTable connectivity(5);
  for (const auto [a, b] : std::array<std::pair<std::size_t, std::size_t>, 4>{{{0, 1}, {1, 4}, {1, 2}, {2, 3}}})
  {
    connectivity[a, b] = true;
    connectivity[b, a] = true;
  }
  Potentials::IntraMolecularPotentials potentials;
  potentials.bonds = {BondPotential({0, 1}, BondType::Harmonic, {96500.0, 1.54}),
                      BondPotential({1, 4}, BondType::Harmonic, {96500.0, 1.54})};
  potentials.bends = {BendPotential({0, 1, 2}, junctionBendType,
                                    junctionBendType == BendType::Harmonic ? std::vector<double>{62500.0, 109.5}
                                                                           : std::vector<double>{109.5}),
                      BendPotential({4, 1, 2}, BendType::Harmonic, {62500.0, 109.5})};
  ChiralCenter centre;
  centre.ids = {1, 0, 2, 4};
  centre.type = ChiralCenter::Chirality::R_Chiral;
  potentials.chiralCenters = {centre};

  Component component(
      forceField, "rigid-at-centre", 300.0, 1.0e6, 0.1,
      {Atom({-1.5, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false), Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false),
       Atom({0.8, 1.2, 0.0}, 0.0, 1.0, 0, 0, 0, false, false), Atom({1.6, 2.4, 0.3}, 0.0, 1.0, 0, 0, 0, false, false),
       Atom({0.2, -0.8, 1.2}, 0.0, 1.0, 0, 0, 0, false, false)},
      connectivity, potentials, 0, 0);
  component.buildFragmentGraph({{1, 2, 3}});
  return component;
}

std::vector<CBMC::GrowStep> planFrom(const Component& component, const std::vector<std::size_t>& placed)
{
  return CBMC::buildGrowthPlan(component.connectivityTable, component.fragmentGraph, component.intraMolecularPotentials,
                               placed);
}

bool hasBend(const std::vector<BendPotential>& bends, std::array<std::size_t, 3> ids)
{
  return std::ranges::any_of(bends, [&](const BendPotential& bend) { return bend.identifiers == ids; });
}
}  // namespace

// The spin of a rigid body hinged on a chiral centre with a kept side group decides the parity of that centre: the
// centre is classified for the spin selection of the rigid step as for a flexible attach (it was for flexible
// attaches only, so the rigid body took either parity).
TEST(CBMC_PLAN_FIXES, a_rigid_body_hinged_on_a_chiral_centre_carries_the_centre_to_the_spin)
{
  const Component component = makeRigidFragmentAtChiralCentre();
  const std::vector<CBMC::GrowStep> plan = planFrom(component, {0, 1, 4});
  ASSERT_FALSE(plan.empty());
  const CBMC::GrowStep& step = plan.front();
  ASSERT_TRUE(step.rigidBody);
  ASSERT_EQ(step.currentBead, 1uz);
  ASSERT_EQ(step.previousBead, std::optional<std::size_t>{0});

  ASSERT_EQ(step.spin.determinedChiralCenters.size(), 1uz);
  EXPECT_EQ(step.spin.determinedChiralCenters.front().ids, (std::array<std::size_t, 4>{1, 0, 2, 4}));
  EXPECT_TRUE(step.spin.parityIsHalfCircle);
  EXPECT_TRUE(step.base.determinedChiralCenters.empty());
}

// The rigid tilt samples the spin-invariant bends only; the spin-variant bend to the side group is weighted by the
// spin selection alone (the tilt summed it as well, so it was counted twice and drawn with a weight the Rosenbluth
// factor did not carry).
TEST(CBMC_PLAN_FIXES, the_rigid_tilt_leaves_the_spin_variant_bends_to_the_spin)
{
  const Component component = makeRigidFragmentAtChiralCentre();
  const std::vector<CBMC::GrowStep> plan = planFrom(component, {0, 1, 4});
  ASSERT_FALSE(plan.empty());
  const CBMC::GrowStep& step = plan.front();
  ASSERT_TRUE(step.rigidBody);
  EXPECT_TRUE(hasBend(step.rigidTilt.bends, {0, 1, 2}));
  EXPECT_FALSE(hasBend(step.rigidTilt.bends, {4, 1, 2}));
  EXPECT_TRUE(hasBend(step.spin.variantBends, {4, 1, 2}));
}

// A rigid body with a second placed atom besides its anchor cannot be regrown from the anchor alone: the plan refuses
// it (an atom of the unit itself was skipped by the check, and the step moved the unplaced atoms from the anchor while
// the second placed atom kept its old position).
TEST(CBMC_PLAN_FIXES, a_rigid_body_with_two_placed_atoms_is_refused)
{
  const Component component = makeRigidFragmentAtChiralCentre();
  EXPECT_THROW((void)planFrom(component, {0, 1, 3, 4}), std::runtime_error);
  EXPECT_NO_THROW((void)planFrom(component, {0, 1, 2, 3, 4}));
}

// A FIXED or RIGID junction bend of a rigid body is a constraint without an energy: the relaxation of the tilt only
// rolls the body about its inner bond, so the angle previous-current-inner stays on the cone it was drawn on (the
// relaxation rotated about any axis and accepted every move, which left the constrained angle uniform on the sphere).
TEST(CBMC_PLAN_FIXES, the_rigid_tilt_keeps_a_fixed_junction_bend_at_its_angle)
{
  for (const BendType type : {BendType::Fixed, BendType::Rigid})
  {
    const Component component = makeRigidFragmentAtChiralCentre(type);
    const std::vector<CBMC::GrowStep> plan = planFrom(component, {0, 1, 4});
    ASSERT_FALSE(plan.empty());
    const CBMC::GrowStep& step = plan.front();
    ASSERT_TRUE(step.rigidBody);
    ASSERT_TRUE(step.rigidTilt.junctionBend.has_value());
    const std::size_t inner = step.rigidTilt.innerBead;
    const auto innerSlot = std::ranges::find(step.nextBeads, inner);
    ASSERT_NE(innerSlot, step.nextBeads.end());

    RandomNumber random(11);
    std::vector<Atom> chain(component.atoms.begin(), component.atoms.end());
    std::size_t numberOfTilts = 0;
    for (std::size_t sample = 0; sample != 50; ++sample)
    {
      const std::optional<std::vector<Atom>> tilt =
          CBMC::generateRigidTilt(random, 10, 1.0 / (0.8314464919 * 300.0), component, chain, step);
      if (!tilt.has_value()) continue;
      ++numberOfTilts;
      const double3 toPrevious = (chain[0].position - chain[1].position).normalized();
      const double3 toInner =
          ((*tilt)[static_cast<std::size_t>(innerSlot - step.nextBeads.begin())].position - chain[1].position)
              .normalized();
      const double angle = std::acos(std::clamp(double3::dot(toPrevious, toInner), -1.0, 1.0)) * 180.0 / std::numbers::pi;
      EXPECT_NEAR(angle, 109.5, 1.0e-6);
    }
    EXPECT_GT(numberOfTilts, 0uz);
  }
}

// A FIXED Urey-Bradley distance between the ends of a FIXED bend is held by the bend only when both arms are FIXED
// bonds whose triangle has that distance: with a sampled arm the 1-3 distance varies, and with other lengths the two
// constraints contradict each other. The plan exempted it on the bend alone.
TEST(CBMC_PLAN_FIXES, a_fixed_urey_bradley_distance_needs_fixed_arms_of_its_length)
{
  const ForceField forceField({{"C", false, 12.0, 0.0, 0.0, 6, false}}, {{50.0, 3.5}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, false, false, false);
  const double angle = 109.5 * std::numbers::pi / 180.0;
  const double closing = std::sqrt(2.0 * 1.54 * 1.54 * (1.0 - std::cos(angle)));
  const auto plan = [&](BondType arms, double distance)
  {
    ConnectivityTable connectivity(3);
    for (const auto [a, b] : std::array<std::pair<std::size_t, std::size_t>, 2>{{{0, 1}, {1, 2}}})
    {
      connectivity[a, b] = true;
      connectivity[b, a] = true;
    }
    Potentials::IntraMolecularPotentials potentials;
    const std::vector<double> bondParameters =
        arms == BondType::Fixed ? std::vector<double>{1.54} : std::vector<double>{96500.0, 1.54};
    potentials.bonds = {BondPotential({0, 1}, arms, bondParameters), BondPotential({1, 2}, arms, bondParameters)};
    potentials.bends = {BendPotential({0, 1, 2}, BendType::Fixed, {109.5})};
    potentials.ureyBradleys = {UreyBradleyPotential({0, 2}, UreyBradleyType::Fixed, {distance})};
    Component component(forceField, "triangle", 300.0, 1.0e6, 0.1,
                        {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false),
                         Atom({1.54, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false),
                         Atom({2.05, 1.45, 0.0}, 0.0, 1.0, 0, 0, 0, false, false)},
                        connectivity, potentials, 0, 0);
    component.buildFragmentGraph({});
    return planFrom(component, {0});
  };
  EXPECT_NO_THROW(static_cast<void>(plan(BondType::Fixed, closing)));
  EXPECT_THROW(static_cast<void>(plan(BondType::Harmonic, closing)), std::runtime_error);
  EXPECT_THROW(static_cast<void>(plan(BondType::Fixed, closing + 0.2)), std::runtime_error);
}
