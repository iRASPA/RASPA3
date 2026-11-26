#include <gtest/gtest.h>

import std;

import atom;
import bend_potential;
import bond_potential;
import torsion_potential;
import component;
import connectivity_table;
import double3;
import forcefield;
import molecule;
import simulationbox;
import system;
import intra_molecular_potentials;
import polymer_links;
import reactive_group;
import polymer_reactions;
import property_polymer_chains;

// The chains of a polymerization per chain length: the clusters of the atom graph of the polymer table, unwrapped
// across the periodic boundaries, classified by their topology, and measured only once the force ramp of their last
// change has ended.

namespace
{
ForceField twoTypeForceField()
{
  return ForceField({{"A", false, 14.0, 0.0, 0.0, 6, false}, {"B", false, 14.0, 0.0, 0.0, 6, false}},
                    {{50.0, 3.5}, {50.0, 3.5}}, ForceField::MixingRule::Lorentz_Berthelot, 9.0, 9.0, 9.0, false, false,
                    false);
}

/// A two-bead molecule A-B (head 0, tail 1) with a harmonic bond and a reactive group 'vinyl'.
Component dimer(const ForceField& forceField)
{
  ConnectivityTable connectivityTable(2);
  connectivityTable[0, 1] = true;
  connectivityTable[1, 0] = true;
  Potentials::IntraMolecularPotentials potentials{};
  potentials.bonds.push_back(BondPotential({0, 1}, BondType::Harmonic, {96500.0, 1.54}));
  const std::vector<Atom> beads{Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false),
                                Atom({1.54, 0.0, 0.0}, 0.0, 1.0, 0, 1, 0, false, false)};
  Component component =
      Component(forceField, "dimer", 425.0, 3796000.0, 0.199, beads, connectivityTable, potentials, 5, 21);
  component.reactiveGroups.push_back(
      ReactiveGroup{.type = "vinyl", .atoms = {{"head", 0}, {"tail", 1}}, .initialState = "monomer"});
  return component;
}

/// Harmonic bond and bend templates for every type combination of A and B, and one TraPPE torsion per sequence.
void addTemplates(PolymerLinkTable& table, const ForceField& forceField)
{
  const std::array<std::string, 2> names{"A", "B"};
  for (std::size_t a = 0; a < 2; ++a)
  {
    for (std::size_t b = a; b < 2; ++b)
    {
      table.templates.bonds.push_back({.types = {names[a], names[b]}, .type = BondType::Harmonic, .parameters = {96500.0, 1.54}});
    }
  }
  for (std::size_t a = 0; a < 2; ++a)
  {
    for (std::size_t b = 0; b < 2; ++b)
    {
      for (std::size_t c = a; c < 2; ++c)
      {
        table.templates.bends.push_back(
            {.types = {names[a], names[b], names[c]}, .type = BendType::Harmonic, .parameters = {62500.0, 112.0}});
      }
    }
  }
  for (std::size_t a = 0; a < 2; ++a)
  {
    for (std::size_t b = 0; b < 2; ++b)
    {
      for (std::size_t c = 0; c < 2; ++c)
      {
        for (std::size_t d = 0; d < 2; ++d)
        {
          if (std::array<std::size_t, 4>{d, c, b, a} < std::array<std::size_t, 4>{a, b, c, d}) continue;
          table.templates.torsions.push_back({.types = {names[a], names[b], names[c], names[d]},
                                              .type = TorsionType::TraPPE,
                                              .parameters = {0.0, 355.03, -68.19, 791.32}});
        }
      }
    }
  }
  table.templates.resolve(forceField);
}

/// 'count' dimers at 'positions' (head 2m, tail 2m + 1) with initiation and propagation (bond tail-head) rules;
/// the table activated at step 1000 with a ramp of 100 steps.
System reactiveSystem(const ForceField& forceField, const SimulationBox& box, std::size_t count,
                      std::span<const double3> positions)
{
  System system = System(forceField, box, false, 300.0, 1e5, 1.0, {}, {dimer(forceField)}, {}, {count}, 5);
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  for (std::size_t i = 0; i < positions.size(); ++i) atoms[i].position = positions[i];

  PolymerLinkTable& table = system.polymerLinks;
  addTemplates(table, forceField);
  ReactionRule initiation{};
  initiation.name = "initiation";
  initiation.reactants = {{.label = "A", .type = "vinyl", .state = "monomer"}};
  initiation.products = {{.state = "radical"}};
  ReactionRule propagation{};
  propagation.name = "propagation";
  propagation.reactants = {{.label = "R", .type = "vinyl", .state = "radical"},
                           {.label = "M", .type = "vinyl", .state = "monomer"}};
  propagation.products = {{.state = "dead"}, {.state = "radical"}};
  propagation.capture = std::array<RoleRef, 2>{RoleRef{"R", "tail"}, RoleRef{"M", "head"}};
  propagation.formBonds = {{RoleRef{"R", "tail"}, RoleRef{"M", "head"}}};
  table.scheme.rules = {initiation, propagation};
  table.scheme.resolve(system.components, forceField);
  table.activate(system.components, system.numberOfMoleculesPerComponent);
  table.relaxationSteps = 100;
  table.currentStep = 1000;
  return system;
}

/// Grows the chain 0 -> 1 -> ... -> length - 1 at step 1000.
void grow(System& system, const ForceField& forceField, std::size_t length)
{
  PolymerLinkTable& table = system.polymerLinks;
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  const std::array<std::size_t, 1> first{0};
  table.apply(0, first, {}, {}, forceField, system.components, atoms, std::span<Molecule>(system.moleculeData),
              system.simulationBox, 1000);
  for (std::size_t m = 0; m + 1 < length; ++m)
  {
    const std::array<std::size_t, 2> pair{m, m + 1};
    table.apply(1, pair, {}, {}, forceField, system.components, atoms, std::span<Molecule>(system.moleculeData),
                system.simulationBox, 1000);
  }
  table.rebuild(forceField, system.components, system.numberOfMoleculesPerComponent, atoms);
}

/// Rg^2 of a set of positions with uniform weights.
double radiusOfGyrationSquared(std::span<const double3> positions)
{
  double3 centre{};
  for (const double3& r : positions) centre += r;
  centre = centre / static_cast<double>(positions.size());
  double sum = 0.0;
  for (const double3& r : positions) sum += double3::dot(r - centre, r - centre);
  return sum / static_cast<double>(positions.size());
}
}  // namespace

// Three of four dimers polymerized: one linear chain of three molecules, living (its last group is the radical),
// without the unreacted fourth dimer, measured from the positions as they are.
TEST(Polymer_chains, a_linear_chain_is_found_and_measured)
{
  const ForceField forceField = twoTypeForceField();
  const std::array<double3, 8> positions{double3(10.0, 10.0, 10.0), double3(11.3, 10.8, 10.1),
                                         double3(12.7, 10.4, 10.9), double3(14.0, 11.2, 11.0),
                                         double3(15.4, 10.7, 11.8), double3(16.7, 11.5, 11.9),
                                         double3(25.0, 25.0, 25.0), double3(26.5, 25.0, 25.0)};
  System system = reactiveSystem(forceField, SimulationBox(30.0, 30.0, 30.0), 4, positions);
  grow(system, forceField, 3);

  const std::vector<PropertyPolymerChains::Chain> chains = PropertyPolymerChains::findChains(
      system.polymerLinks, system.components, system.simulationBox, system.spanOfMoleculeAtoms());
  ASSERT_EQ(chains.size(), 1uz);
  const PropertyPolymerChains::Chain& chain = chains.front();
  EXPECT_EQ(chain.length, 3uz);
  EXPECT_EQ(chain.polymerBonds, 2uz);
  EXPECT_EQ(chain.atoms.size(), 6uz);
  EXPECT_EQ(chain.topology, PropertyPolymerChains::Topology::Linear);
  EXPECT_TRUE(chain.living);
  EXPECT_FALSE(chain.spanning);
  EXPECT_EQ(chain.lastChange, 1000u);
  EXPECT_EQ(chain.contourBonds, 5uz);
  const double3 endToEnd = positions[5] - positions[0];
  EXPECT_NEAR(std::abs(double3::dot(chain.endToEnd, endToEnd)), double3::dot(endToEnd, endToEnd), 1e-10);

  const PropertyPolymerChains::Moments values =
      PropertyPolymerChains::measure(chain, forceField, system.spanOfMoleculeAtoms(), false);
  EXPECT_NEAR(values[PropertyPolymerChains::RadiusOfGyrationSquared],
              radiusOfGyrationSquared(std::span<const double3>(positions.data(), 6)), 1e-10);
  EXPECT_NEAR(values[PropertyPolymerChains::EndToEndSquared], double3::dot(endToEnd, endToEnd), 1e-10);
  EXPECT_NEAR(values[PropertyPolymerChains::Mass], 6.0 * 14.0, 1e-10);
}

// The same chain shifted so that it crosses the periodic boundary is unwrapped: the same Rg^2 and end-to-end
// distance as the chain in one piece.
TEST(Polymer_chains, a_chain_across_the_boundary_is_unwrapped)
{
  const ForceField forceField = twoTypeForceField();
  const SimulationBox box(30.0, 30.0, 30.0);
  std::array<double3, 8> positions{double3(10.0, 10.0, 10.0), double3(11.3, 10.8, 10.1),
                                   double3(12.7, 10.4, 10.9), double3(14.0, 11.2, 11.0),
                                   double3(15.4, 10.7, 11.8), double3(16.7, 11.5, 11.9),
                                   double3(25.0, 25.0, 25.0), double3(26.5, 25.0, 25.0)};
  const double whole = radiusOfGyrationSquared(std::span<const double3>(positions.data(), 6));
  const double endToEnd = (positions[5] - positions[0]).length();
  // shift by 16.5 A along x and wrap: the boundary cuts the bond of the second molecule
  for (std::size_t i = 0; i < 6; ++i)
  {
    positions[i].x += 16.5;
    if (positions[i].x >= 30.0) positions[i].x -= 30.0;
  }
  ASSERT_LT(positions[3].x, positions[0].x);

  System system = reactiveSystem(forceField, box, 4, positions);
  grow(system, forceField, 3);
  const std::vector<PropertyPolymerChains::Chain> chains = PropertyPolymerChains::findChains(
      system.polymerLinks, system.components, system.simulationBox, system.spanOfMoleculeAtoms());
  ASSERT_EQ(chains.size(), 1uz);
  EXPECT_FALSE(chains.front().spanning);
  const PropertyPolymerChains::Moments values =
      PropertyPolymerChains::measure(chains.front(), forceField, system.spanOfMoleculeAtoms(), false);
  EXPECT_NEAR(values[PropertyPolymerChains::RadiusOfGyrationSquared], whole, 1e-10);
  EXPECT_NEAR(values[PropertyPolymerChains::EndToEnd], endToEnd, 1e-10);
}

// A chain is measured only once its last change is older than the ramp plus the lag; the key holds its length,
// topology, state and conversion bin (three of four groups reacted: conversion 0.75, bin 7 of 10).
TEST(Polymer_chains, a_chain_is_sampled_after_its_ramp)
{
  const ForceField forceField = twoTypeForceField();
  const std::array<double3, 8> positions{double3(10.0, 10.0, 10.0), double3(11.3, 10.8, 10.1),
                                         double3(12.7, 10.4, 10.9), double3(14.0, 11.2, 11.0),
                                         double3(15.4, 10.7, 11.8), double3(16.7, 11.5, 11.9),
                                         double3(25.0, 25.0, 25.0), double3(26.5, 25.0, 25.0)};
  System system = reactiveSystem(forceField, SimulationBox(30.0, 30.0, 30.0), 4, positions);
  grow(system, forceField, 3);
  PolymerLinkTable& table = system.polymerLinks;
  ASSERT_NEAR(table.conversion(), 0.75, 1e-12);

  PropertyPolymerChains property(10, 1000, PropertyPolymerChains::Sampling::Periodic, 50, 10, 2, false);
  table.currentStep = 1149;
  property.sample(table, forceField, system.components, system.simulationBox, system.spanOfMoleculeAtoms());
  EXPECT_TRUE(property.accumulators.empty());
  EXPECT_EQ(property.skippedRamping, 1u);

  table.currentStep = 1150;
  property.sample(table, forceField, system.components, system.simulationBox, system.spanOfMoleculeAtoms());
  property.sample(table, forceField, system.components, system.simulationBox, system.spanOfMoleculeAtoms());
  ASSERT_EQ(property.accumulators.size(), 1uz);
  const auto& [key, accumulator] = *property.accumulators.begin();
  EXPECT_EQ(key.length, 3uz);
  EXPECT_EQ(key.topology, PropertyPolymerChains::Topology::Linear);
  EXPECT_TRUE(key.living);
  EXPECT_EQ(key.conversionBin, 7uz);
  EXPECT_EQ(accumulator.count, 2u);
  EXPECT_NEAR(accumulator.mean(PropertyPolymerChains::NumberOfAtoms), 6.0, 1e-12);
  EXPECT_NEAR(accumulator.standardDeviation(PropertyPolymerChains::RadiusOfGyrationSquared), 0.0, 1e-8);
  EXPECT_EQ(property.numberOfSamplingEvents, 3u);

  // a minimum length above the chain's: not sampled, and not counted as skipped
  PropertyPolymerChains longOnly(10, 1000, PropertyPolymerChains::Sampling::Periodic, 0, 10, 4, false);
  longOnly.sample(table, forceField, system.components, system.simulationBox, system.spanOfMoleculeAtoms());
  EXPECT_TRUE(longOnly.accumulators.empty());
  EXPECT_EQ(longOnly.skippedRamping, 0u);
}

// The sampling before a reaction looks only at the chains of the reacting atoms: a seed in the unreacted dimer
// finds none, a seed in the chain finds it once even when both atoms of a reaction are in it.
TEST(Polymer_chains, seeds_select_the_chains_of_the_reactants)
{
  const ForceField forceField = twoTypeForceField();
  const std::array<double3, 8> positions{double3(10.0, 10.0, 10.0), double3(11.3, 10.8, 10.1),
                                         double3(12.7, 10.4, 10.9), double3(14.0, 11.2, 11.0),
                                         double3(15.4, 10.7, 11.8), double3(16.7, 11.5, 11.9),
                                         double3(25.0, 25.0, 25.0), double3(26.5, 25.0, 25.0)};
  System system = reactiveSystem(forceField, SimulationBox(30.0, 30.0, 30.0), 4, positions);
  grow(system, forceField, 3);
  const auto chainsOf = [&](std::vector<std::uint32_t> seeds)
  {
    return PropertyPolymerChains::findChains(system.polymerLinks, system.components, system.simulationBox,
                                             system.spanOfMoleculeAtoms(), seeds);
  };
  EXPECT_TRUE(chainsOf({7}).empty());
  EXPECT_EQ(chainsOf({5, 7}).size(), 1uz);
  EXPECT_EQ(chainsOf({0, 5}).size(), 1uz);
}

// A molecule with three polymer bonds is a branch point; a closed ring of molecules is cyclic; a ring that closes
// through the periodic boundary spans the box and is not measured.
TEST(Polymer_chains, topologies_and_spanning_rings)
{
  const ForceField forceField = twoTypeForceField();
  auto bondOf = [](std::uint32_t moleculeA, std::uint32_t atomA, std::uint32_t moleculeB, std::uint32_t atomB)
  { return PolymerBond{.a = {0, moleculeA, atomA}, .b = {0, moleculeB, atomB}, .reaction = 1, .formedAtStep = 0}; };

  {
    // a star: the tail of dimer 0 bonded to the heads of 1 and 2, its head to the tail of 3
    const std::array<double3, 8> positions{double3(10.0, 10.0, 10.0), double3(11.4, 10.0, 10.0),
                                           double3(12.2, 11.2, 10.0), double3(13.0, 12.4, 10.0),
                                           double3(12.2, 8.8, 10.0),  double3(13.0, 7.6, 10.0),
                                           double3(7.2, 10.0, 10.0),  double3(8.6, 10.0, 10.0)};
    System system = reactiveSystem(forceField, SimulationBox(30.0, 30.0, 30.0), 4, positions);
    PolymerLinkTable& table = system.polymerLinks;
    table.bonds = {bondOf(0, 1, 1, 0), bondOf(0, 1, 2, 0), bondOf(3, 1, 0, 0)};
    table.rebuild(forceField, system.components, system.numberOfMoleculesPerComponent, system.spanOfMoleculeAtoms());
    const std::vector<PropertyPolymerChains::Chain> chains =
        PropertyPolymerChains::findChains(table, system.components, system.simulationBox, system.spanOfMoleculeAtoms());
    ASSERT_EQ(chains.size(), 1uz);
    EXPECT_EQ(chains.front().length, 4uz);
    EXPECT_EQ(chains.front().topology, PropertyPolymerChains::Topology::Branched);
    EXPECT_FALSE(chains.front().living);
  }

  // eight dimers in a row along x, 1.4 A apart; the tail of the last bonded to the head of the first
  std::array<double3, 16> row{};
  for (std::size_t i = 0; i < row.size(); ++i) row[i] = double3(1.0 + 1.4 * static_cast<double>(i), 10.0, 10.0);
  std::vector<PolymerBond> ring{};
  for (std::uint32_t m = 0; m < 7; ++m) ring.push_back(bondOf(m, 1, m + 1, 0));
  ring.push_back(bondOf(7, 1, 0, 0));

  {
    // in a box of 22.4 A the closing bond is 1.4 A through the boundary: the ring winds around the box
    System system = reactiveSystem(forceField, SimulationBox(22.4, 30.0, 30.0), 8, row);
    PolymerLinkTable& table = system.polymerLinks;
    table.bonds = ring;
    table.rebuild(forceField, system.components, system.numberOfMoleculesPerComponent, system.spanOfMoleculeAtoms());
    const std::vector<PropertyPolymerChains::Chain> chains =
        PropertyPolymerChains::findChains(table, system.components, system.simulationBox, system.spanOfMoleculeAtoms());
    ASSERT_EQ(chains.size(), 1uz);
    EXPECT_EQ(chains.front().length, 8uz);
    EXPECT_EQ(chains.front().topology, PropertyPolymerChains::Topology::Cyclic);
    EXPECT_TRUE(chains.front().spanning);

    PropertyPolymerChains property(10, 1000, PropertyPolymerChains::Sampling::Periodic, 0, 10, 2, false);
    property.sample(table, forceField, system.components, system.simulationBox, system.spanOfMoleculeAtoms());
    EXPECT_TRUE(property.accumulators.empty());
    EXPECT_EQ(property.skippedSpanning, 1u);
  }

  {
    // the same ring folded in a large box (the closing bond made short in space): cyclic, but not spanning
    std::array<double3, 16> folded{};
    for (std::size_t i = 0; i < folded.size(); ++i)
    {
      const double angle = 2.0 * std::numbers::pi * static_cast<double>(i) / 16.0;
      folded[i] = double3(15.0 + 3.6 * std::cos(angle), 15.0 + 3.6 * std::sin(angle), 15.0);
    }
    System system = reactiveSystem(forceField, SimulationBox(30.0, 30.0, 30.0), 8, folded);
    PolymerLinkTable& table = system.polymerLinks;
    table.bonds = ring;
    table.rebuild(forceField, system.components, system.numberOfMoleculesPerComponent, system.spanOfMoleculeAtoms());
    const std::vector<PropertyPolymerChains::Chain> chains =
        PropertyPolymerChains::findChains(table, system.components, system.simulationBox, system.spanOfMoleculeAtoms());
    ASSERT_EQ(chains.size(), 1uz);
    EXPECT_EQ(chains.front().topology, PropertyPolymerChains::Topology::Cyclic);
    EXPECT_FALSE(chains.front().spanning);
    EXPECT_NEAR(PropertyPolymerChains::measure(chains.front(), forceField, system.spanOfMoleculeAtoms(), false)
                    [PropertyPolymerChains::RadiusOfGyrationSquared],
                3.6 * 3.6, 1e-10);
  }
}
