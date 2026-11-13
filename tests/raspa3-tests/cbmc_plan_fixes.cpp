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
import chiral_center;
import cbmc_grow_step;
import cbmc_growth_plan;

namespace
{
// A rigid fragment {1, 2, 3} hinged on atom 1, which is a declared chiral centre with the previous bead 0, the
// fragment atom 2 and a side group 4. The bend 0-1-2 is the junction bend (spin invariant); 4-1-2 couples the
// fragment to the side group off the spin axis (spin variant).
Component makeRigidFragmentAtChiralCentre()
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
  potentials.bends = {BendPotential({0, 1, 2}, BendType::Harmonic, {62500.0, 109.5}),
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
