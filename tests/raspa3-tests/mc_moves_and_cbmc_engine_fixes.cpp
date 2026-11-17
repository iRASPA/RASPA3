#include <gtest/gtest.h>

#include "../test_support.hpp"
#include "molecule_fixtures.hpp"

import std;

import double3;
import atom;
import molecule;
import units;
import forcefield;
import component;
import framework;
import simulationbox;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;
import bend_potential;
import torsion_potential;
import chiral_center;
import interpolation_energy_grid;
import randomnumbers;
import mc_moves_probabilities;
import property_lambda_probability_histogram;
import cbmc;
import cbmc_util;
import cbmc_grow_step;
import cbmc_growth_plan;
import system;
import running_energy;
import mc_moves;
import mc_moves_move_types;

// Regression tests of the Monte Carlo move and CBMC engine fixes: the Boltzmann factor of a first bead placed at a
// position the caller proposed, lambda steps that always move, the chirality of a rigid body hinged on a
// stereocentre, a ring closed from a placed rigid unit, and fixed bends and torsions the samplers do not hold.

namespace
{
constexpr double kTemperature = 300.0;

ForceField makeForceField()
{
  return ForceField({{"CH4", false, 16.04246, 0.0, 0.0, 6, false}, {"CH2_c", false, 14.02658, 0.0, 0.0, 6, false}},
                    {{158.5, 3.72}, {52.5, 3.91}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true,
                    true, false);
}

struct Environment
{
  ForceField forceField{makeForceField()};
  SimulationBox box{30.0, 30.0, 30.0};
  double beta{1.0 / (Units::KB * kTemperature)};
  std::optional<Framework> noFramework{};
  std::vector<std::optional<InterpolationEnergyGrid>> noGrids{forceField.pseudoAtoms.size() + 1};
  std::optional<InterpolationEnergyGrid> noExternalFieldGrid{};

  CBMC::GrowContext context(std::span<const Atom> background) const
  {
    return CBMC::GrowContext(false, forceField, box, noGrids, noExternalFieldGrid, noFramework,
                             std::span<const Atom>{}, background, beta, CBMC::CutOffMode::Full);
  }
};

std::vector<CBMC::GrowStep> planFrom(const Component& component, const std::vector<std::size_t>& placed)
{
  return CBMC::buildGrowthPlan(component.connectivityTable, component.fragmentGraph, component.intraMolecularPotentials,
                               placed);
}

ConnectivityTable connect(std::size_t numberOfBeads, std::initializer_list<std::pair<std::size_t, std::size_t>> bonds)
{
  ConnectivityTable connectivity(numberOfBeads);
  for (const auto [a, b] : bonds)
  {
    connectivity[a, b] = true;
    connectivity[b, a] = true;
  }
  return connectivity;
}

Atom carbon(double3 position) { return Atom(position, 0.0, 1.0, 0, 1, 0, false, false); }

// A three-membered rigid ring {1, 2, 3} hinged on atom 1, a stereocentre with the tail atom 0 and the two ring
// atoms: the tail lies off the plane of the ring, on the side the reference geometry gives.
Component makeRigidRingAtStereocentre(const ForceField& forceField)
{
  Potentials::IntraMolecularPotentials potentials;
  potentials.bonds = {BondPotential({0, 1}, BondType::Harmonic, {96500.0, 1.5})};
  potentials.bends = {BendPotential({0, 1, 2}, BendType::Harmonic, {62500.0, 110.4}),
                      BendPotential({0, 1, 3}, BendType::Harmonic, {62500.0, 110.4})};
  ChiralCenter centre;
  centre.ids = {1, 0, 2, 3};
  centre.type = ChiralCenter::Chirality::R_Chiral;
  potentials.chiralCenters = {centre};

  Component component(forceField, "rigid-ring-at-stereocentre", 300.0, 1.0e6, 0.1,
                      {carbon({-1.5, 0.0, 0.0}), carbon({0.0, 0.0, 0.0}), carbon({0.5, 1.2, 0.6}),
                       carbon({0.5, 1.2, -0.6})},
                      connect(4, {{0, 1}, {1, 2}, {1, 3}, {2, 3}}), potentials, 5, 21);
  component.buildFragmentGraph({{1, 2, 3}});
  return component;
}

// Isobutane-like: the branch point 1 with the previous bead 0 and the two beads 2 and 3 grown together.
Component makeBranch(const ForceField& forceField, BendType siblingBend)
{
  Potentials::IntraMolecularPotentials potentials;
  potentials.bonds = {BondPotential({0, 1}, BondType::Harmonic, {96500.0, 1.54}),
                      BondPotential({1, 2}, BondType::Harmonic, {96500.0, 1.54}),
                      BondPotential({1, 3}, BondType::Harmonic, {96500.0, 1.54})};
  potentials.bends = {BendPotential({0, 1, 2}, BendType::Fixed, {112.0}),
                      BendPotential({0, 1, 3}, BendType::Fixed, {112.0}),
                      siblingBend == BendType::Fixed ? BendPotential({2, 1, 3}, BendType::Fixed, {112.0})
                                                     : BendPotential({2, 1, 3}, BendType::Harmonic, {62500.0, 112.0})};
  return Component(forceField, "branch", 300.0, 1.0e6, 0.1,
                   {carbon({-1.5, 0.0, 0.0}), carbon({0.0, 0.0, 0.0}), carbon({0.6, 1.2, 0.7}),
                    carbon({0.6, -1.2, 0.7})},
                   connect(4, {{0, 1}, {1, 2}, {1, 3}}), potentials, 5, 21);
}
}  // namespace

// The first bead of a satellite of a pair or group move is placed at a position the move proposed. The scheme gave
// it a weight of one: the move accounts for the density of the proposal, not for the energy of the bead, and for a
// monoatomic satellite the whole of its pair energy was left out of the acceptance rule.
TEST(MC_MOVES_AND_CBMC_ENGINE, a_first_bead_at_a_proposed_position_carries_its_boltzmann_factor)
{
  Environment environment;
  TemporaryFile file("engine-fixes-methane.json", molecule_fixtures::kMethaneJson);
  Component methane(Component::Type::Adsorbate, 0, environment.forceField, "methane", file.stemPath().string(), 5, 21,
                    MCMoveProbabilities(), std::nullopt, false);
  methane.prepareGrowthPlans(environment.beta);

  // a neighbour 4 Angstrom away: inside the well of the pair potential
  Atom neighbour = methane.atoms[0];
  neighbour.position = double3(10.0, 10.0, 10.0);
  neighbour.moleculeId = 0;
  const std::vector<Atom> background{neighbour};
  const CBMC::GrowContext context = environment.context(background);
  const double3 position(14.0, 10.0, 10.0);

  RandomNumber random(7);
  const std::optional<CBMC::GrowResult> grown =
      CBMC::growNewMolecule(random, context, methane, {.componentId = 0, .moleculeId = 1},
                            {.firstBead = CBMC::FirstBeadScheme::Fixed, .firstBeadPosition = position});
  ASSERT_TRUE(grown.has_value());
  ASSERT_LT(grown->energies.potentialEnergy(), -1.0e-3);
  EXPECT_NEAR(grown->logRosenbluthWeight, -environment.beta * grown->energies.potentialEnergy(), 1.0e-12);

  const CBMC::RetraceResult retraced =
      CBMC::retraceMolecule(random, context, methane, grown->atoms, {.firstBead = CBMC::FirstBeadScheme::Fixed});
  EXPECT_NEAR(retraced.logRosenbluthWeight, grown->logRosenbluthWeight, 1.0e-12);
}

// A lambda step was the truncation of a real number of bins: nothing at all once the scaled range fell to one bin
// or less. Every proposal now moves, by one bin at least, and as often up as down.
TEST(MC_MOVES_AND_CBMC_ENGINE, a_lambda_step_always_leaves_the_bin)
{
  RandomNumber random(11);
  for (const std::size_t points : {2uz, 10uz, 41uz})
  {
    PropertyLambdaProbabilityHistogram histogram(5, points);
    histogram.setCurrentBin(points / 2);
    const int current = static_cast<int>(points / 2);
    for (const double scale : {0.01, 0.1, 0.5})
    {
      const int maximumStep = std::max(1, static_cast<int>(scale * static_cast<double>(points)));
      std::size_t up = 0;
      constexpr std::size_t numberOfDraws = 20000;
      for (std::size_t i = 0; i != numberOfDraws; ++i)
      {
        const int step = histogram.selectNewBin(random, scale) - current;
        ASSERT_NE(step, 0);
        ASSERT_LE(std::abs(step), maximumStep);
        if (step > 0) ++up;
      }
      EXPECT_NEAR(static_cast<double>(up) / numberOfDraws, 0.5, 0.02);
    }
  }
}

// A rigid body hinged on a stereocentre whose other neighbours are the previous bead and atoms of the body: the
// tilt of the body decides the parity. The centre belonged to no sampler, and the body was placed on either side.
TEST(MC_MOVES_AND_CBMC_ENGINE, a_rigid_body_hinged_on_a_stereocentre_keeps_its_chirality)
{
  Environment environment;
  Component component = makeRigidRingAtStereocentre(environment.forceField);

  const std::vector<CBMC::GrowStep> plan = planFrom(component, {0, 1});
  ASSERT_EQ(plan.size(), 1uz);
  ASSERT_TRUE(plan.front().rigidBody);
  ASSERT_EQ(plan.front().previousBead, std::optional<std::size_t>{0});
  ASSERT_EQ(plan.front().rigidTilt.chiralCenters.size(), 1uz);
  EXPECT_EQ(plan.front().rigidTilt.chiralCenters.front().ids, (std::array<std::size_t, 4>{1, 0, 2, 3}));

  component.prepareGrowthPlans(environment.beta);
  const CBMC::GrowContext context = environment.context({});
  const std::array<std::size_t, 4> ids{1, 0, 2, 3};
  const double reference = CBMC::chiralSignedVolume(ids, component.atoms);
  ASSERT_GT(std::abs(reference), 0.1);

  RandomNumber random(3);
  std::size_t numberGrown = 0;
  for (std::size_t i = 0; i != 200; ++i)
  {
    const std::optional<CBMC::GrowResult> grown =
        CBMC::growNewMolecule(random, context, component, {.componentId = 0, .moleculeId = 0});
    if (!grown.has_value()) continue;
    ++numberGrown;
    EXPECT_GT(CBMC::chiralSignedVolume(ids, grown->atoms) * reference, 0.0);
  }
  EXPECT_GT(numberGrown, 150uz);
}

// A ring closed from a placed rigid unit it contains: the flexible atoms {3, 4, 5} are bonded to the unit {0, 1, 2}
// at both ends. Seeded and spun about the anchor alone, the ring left the bond to the second atom of the unit
// stretched by Angstroms, with no weight for it. The step is held by the whole unit instead.
TEST(MC_MOVES_AND_CBMC_ENGINE, a_ring_closed_from_a_placed_rigid_unit_stays_closed_onto_it)
{
  Environment environment;
  std::string json(molecule_fixtures::kRigidUnitRingJson);
  const std::string_view from = "\"StartingBead\" : 3";
  json.replace(json.find(from), from.size(), "\"StartingBead\" : 0");
  TemporaryFile file("engine-fixes-rigid-unit-ring.json", json);
  Component molecule(Component::Type::Adsorbate, 0, environment.forceField, "ring", file.stemPath().string(), 5, 21,
                     MCMoveProbabilities(), std::nullopt, false);
  molecule.prepareGrowthPlans(environment.beta);

  const std::vector<CBMC::GrowStep> plan = molecule.growthPlan({0});
  ASSERT_EQ(plan.size(), 2uz);
  const CBMC::GrowStep& ring = plan[1];
  ASSERT_EQ(ring.kind, CBMC::GrowStep::Kind::CloseRing);
  ASSERT_TRUE(ring.ring.anchored());
  EXPECT_EQ(ring.ring.placedFrame, (std::vector<std::size_t>{0, 1, 2}));
  // nothing is left to a spin about the anchor
  EXPECT_TRUE(ring.spin.potentials.torsions.empty());
  EXPECT_TRUE(ring.spin.variantBends.empty());
  EXPECT_EQ(ring.ring.internalPotentials.bends.size(), ring.intra.bends.size());
  EXPECT_EQ(ring.ring.internalPotentials.torsions.size(), ring.intra.torsions.size());

  const CBMC::GrowContext context = environment.context({});
  RandomNumber random(5);
  std::size_t numberGrown = 0;
  for (std::size_t i = 0; i != 100; ++i)
  {
    const std::optional<CBMC::GrowResult> grown =
        CBMC::growNewMolecule(random, context, molecule, {.componentId = 0, .moleculeId = 0});
    if (!grown.has_value()) continue;
    ++numberGrown;
    for (const auto [a, b] : std::array<std::pair<std::size_t, std::size_t>, 4>{{{2, 3}, {3, 4}, {4, 5}, {5, 0}}})
    {
      const double length = (grown->atoms[a].position - grown->atoms[b].position).length();
      EXPECT_GT(length, 1.3) << "bond " << a << "-" << b;
      EXPECT_LT(length, 1.8) << "bond " << a << "-" << b;
    }
    // the unit keeps the geometry of the reference
    EXPECT_NEAR((grown->atoms[0].position - grown->atoms[2].position).length(),
                (molecule.atoms[0].position - molecule.atoms[2].position).length(), 1.0e-8);

    // the retrace of what was grown has a weight
    const CBMC::RetraceResult retraced = CBMC::retraceMolecule(random, context, molecule, grown->atoms);
    EXPECT_TRUE(std::isfinite(retraced.logRosenbluthWeight));
  }
  EXPECT_GT(numberGrown, 90uz);
}

// A FIXED bend is held where the sampler draws the angle itself: the bend previous-current-next of a grown bead.
// As the bend between two beads grown together it contributed its zero energy and was left free: refused.
TEST(MC_MOVES_AND_CBMC_ENGINE, a_fixed_bend_the_sampler_does_not_hold_is_refused)
{
  const ForceField forceField = makeForceField();
  EXPECT_NO_THROW(
      {
        const Component held = makeBranch(forceField, BendType::Harmonic);
        (void)planFrom(held, {0, 1});
      });
  EXPECT_THROW(
      {
        const Component free = makeBranch(forceField, BendType::Fixed);
        (void)planFrom(free, {0, 1});
      },
      std::runtime_error);
}

// The integrators fill the trial structure factors; the stored ones are those of the last full energy evaluation.
// The particle exchange of the grand-canonical dynamics took its Fourier energy difference against them: against
// positions the molecules had left (and against nothing at all when no full evaluation had been made). The stored
// structure factors are rebuilt before the exchange.
TEST(MC_MOVES_AND_CBMC_ENGINE, the_exchange_of_grand_canonical_dynamics_uses_the_current_structure_factors)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component co2 = Component::makeCO2(forceField, 0, true);
  co2.mc_moves_probabilities.setProbability(Move::Types::SwapCBMC, 1.0);
  co2.fugacityCoefficient = 50.0;
  System system = System(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e5, 1.0, {}, {co2}, {}, {20}, 5);
  system.runningEnergies = system.computeTotalEnergies();

  // every molecule moved as a whole, as an integration step moves it, and the energy taken from the gradients
  RandomNumber random(9);
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  for (Molecule& molecule : system.moleculeData)
  {
    const double3 shift(3.0 * (random.uniform() - 0.5), 3.0 * (random.uniform() - 0.5),
                        3.0 * (random.uniform() - 0.5));
    molecule.centerOfMassPosition += shift;
    for (std::size_t i = 0; i != molecule.numberOfAtoms; ++i) atoms[molecule.atomIndex + i].position += shift;
  }
  system.precomputeTotalGradients();

  std::size_t numberOfExchanges = 0;
  for (std::size_t attempt = 0; attempt != 2000 && numberOfExchanges < 5; ++attempt)
  {
    std::size_t exchangedMolecule{};
    if (MC_Moves::performMolecularDynamicsSwap(random, system, 0, exchangedMolecule) !=
        MC_Moves::ParticleExchangeResult::Rejected)
    {
      ++numberOfExchanges;
    }
  }
  ASSERT_GT(numberOfExchanges, 0uz);

  const double running = system.runningEnergies.potentialEnergy();
  const double recomputed = system.computeTotalEnergies().potentialEnergy();
  EXPECT_NEAR(running, recomputed, 1.0e-6 * std::max(1.0, std::abs(recomputed)));
}
