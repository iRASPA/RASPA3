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

// The polymer link table between reactions and in its network statistics: the adjacency a reaction changes is seen
// by the angular strain of the next candidate of the same scan (before the derived term lists are rebuilt), the
// order of two retypings of one atom at the same step is the order of the reactions, and the nodes of the network
// statistics are the fragments of the molecules (a molecule a reaction split is two nodes).

namespace
{
ForceField twoTypeForceField()
{
  return ForceField({{"A", false, 14.0, 0.0, 0.0, 6, false}, {"B", false, 14.0, 0.0, 0.0, 6, false}},
                    {{50.0, 3.5}, {50.0, 3.5}}, ForceField::MixingRule::Lorentz_Berthelot, 9.0, 9.0, 9.0, false, false,
                    false);
}

/// A two-bead molecule A-B (head 0, tail 1) with a harmonic bond, a reactive group 'vinyl' on it.
Component dimer(const ForceField& forceField, const std::string& name, bool reactive)
{
  ConnectivityTable connectivityTable(2);
  connectivityTable[0, 1] = true;
  connectivityTable[1, 0] = true;
  Potentials::IntraMolecularPotentials potentials{};
  potentials.bonds.push_back(BondPotential({0, 1}, BondType::Harmonic, {96500.0, 1.54}));
  const std::vector<Atom> beads{Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false),
                                Atom({1.54, 0.0, 0.0}, 0.0, 1.0, 0, 1, 0, false, false)};
  Component component =
      Component(forceField, name, 425.0, 3796000.0, 0.199, beads, connectivityTable, potentials, 5, 21);
  if (reactive)
  {
    component.reactiveGroups.push_back(
        ReactiveGroup{.type = "vinyl", .atoms = {{"head", 0}, {"tail", 1}}, .initialState = "monomer"});
  }
  return component;
}

/// A single-bead molecule A.
Component monomer(const ForceField& forceField, const std::string& name)
{
  ConnectivityTable connectivityTable(1);
  Potentials::IntraMolecularPotentials potentials{};
  const std::vector<Atom> beads{Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false)};
  return Component(forceField, name, 425.0, 3796000.0, 0.199, beads, connectivityTable, potentials, 5, 21);
}

/// Templates for every type combination of A and B: harmonic bonds and bends, a TraPPE torsion.
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

/// Three dimers in a zigzag; initiation, propagation (tail-head bond) and scission (the head-tail bond of a dead
/// group broken) rules; the table activated, the first group a radical.
System reactiveSystem(const ForceField& forceField)
{
  System system = System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {},
                         {dimer(forceField, "dimer", true)}, {}, {3}, 5);
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  const std::array<double3, 6> positions{double3(10.0, 10.0, 10.0), double3(11.3, 10.8, 10.1),
                                         double3(12.7, 10.4, 10.9), double3(14.0, 11.2, 11.0),
                                         double3(15.4, 10.7, 11.8), double3(16.7, 11.5, 11.9)};
  for (std::size_t i = 0; i < 6; ++i) atoms[i].position = positions[i];

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
  ReactionRule scission{};
  scission.name = "scission";
  scission.reactants = {{.label = "D", .type = "vinyl", .state = "dead"}};
  scission.products = {{.state = "cleaved"}};
  scission.breakBonds = {{RoleRef{"D", "head"}, RoleRef{"D", "tail"}}};
  table.scheme.rules = {initiation, propagation, scission};
  table.scheme.resolve(system.components, forceField);
  table.activate(system.components, system.numberOfMoleculesPerComponent);
  table.currentStep = 1000;
  const std::array<std::size_t, 1> first{0};
  table.apply(0, first, {}, {}, forceField, system.components, atoms, std::span<Molecule>(system.moleculeData),
              system.simulationBox, 1000);
  return system;
}
}  // namespace

TEST(Polymer_links, adjacency_follows_a_reaction_before_the_rebuild)
{
  const ForceField forceField = twoTypeForceField();
  System system = reactiveSystem(forceField);
  PolymerLinkTable& table = system.polymerLinks;
  const std::vector<Component>& components = system.components;
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  table.rebuild(forceField, components, system.numberOfMoleculesPerComponent, atoms);

  // atoms: molecule m has head 2m and tail 2m + 1
  const std::map<std::uint32_t, std::uint16_t> noRetyping{};
  const double before = table.angularStrain(3, 4, system.simulationBox, components, atoms, noRetyping);

  // propagation 0 -> 1 forms the bond tail(0)-head(1); the strain of the candidate tail(1)-head(2) then holds the
  // torsion tail(0)-head(1)-tail(1)-head(2) through the new bond
  const std::array<std::size_t, 2> pair{0, 1};
  table.apply(1, pair, {}, {}, forceField, components, atoms, std::span<Molecule>(system.moleculeData),
              system.simulationBox, 1000);
  EXPECT_TRUE(std::ranges::contains(table.neighboursOf(1, components), 2u));
  EXPECT_TRUE(std::ranges::contains(table.neighboursOf(2, components), 1u));
  const double withoutRebuild = table.angularStrain(3, 4, system.simulationBox, components, atoms, noRetyping);
  table.rebuild(forceField, components, system.numberOfMoleculesPerComponent, atoms);
  const double withRebuild = table.angularStrain(3, 4, system.simulationBox, components, atoms, noRetyping);
  EXPECT_NEAR(withoutRebuild, withRebuild, 1e-10);
  EXPECT_GT(std::abs(withRebuild - before), 1.0);

  // scission of the dead group 0 breaks head(0)-tail(0): gone from the adjacency at once
  const std::array<std::size_t, 1> dead{0};
  table.apply(2, dead, {}, {}, forceField, components, atoms, std::span<Molecule>(system.moleculeData),
              system.simulationBox, 1000);
  EXPECT_FALSE(std::ranges::contains(table.neighboursOf(0, components), 1u));
  EXPECT_EQ(table.neighboursOf(1, components), std::vector<std::uint32_t>{2u});
  const double broken = table.angularStrain(3, 4, system.simulationBox, components, atoms, noRetyping);
  table.rebuild(forceField, components, system.numberOfMoleculesPerComponent, atoms);
  EXPECT_NEAR(broken, table.angularStrain(3, 4, system.simulationBox, components, atoms, noRetyping), 1e-10);
  EXPECT_FALSE(std::ranges::contains(table.neighboursOf(0, components), 1u));
}

TEST(Polymer_links, network_nodes_are_the_fragments_of_a_split_molecule)
{
  const ForceField forceField = twoTypeForceField();
  // an initiator dimer (molecule 0 of component 0) and two monomers (component 1); each initiator atom bonded to
  // one monomer; the initiator split by a broken bond, its halves 11 A apart in a 20 A box
  System system = System(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1e5, 1.0, {},
                         {dimer(forceField, "initiator", true), monomer(forceField, "monomer")}, {}, {1, 2}, 5);
  std::span<Atom> atoms = system.spanOfMoleculeAtoms();
  atoms[0].position = double3(1.0, 10.0, 10.0);
  atoms[1].position = double3(12.0, 10.0, 10.0);
  atoms[2].position = double3(2.5, 10.0, 10.0);
  atoms[3].position = double3(13.5, 10.0, 10.0);

  PolymerLinkTable& table = system.polymerLinks;
  // activation lays out the atoms (a scheme is required, the group of the initiator takes no part here)
  ReactionRule rule{};
  rule.name = "initiation";
  rule.reactants = {{.label = "A", .type = "vinyl", .state = "monomer"}};
  rule.products = {{.state = "radical"}};
  table.scheme.rules = {rule};
  table.scheme.resolve(system.components, forceField);
  table.activate(system.components, system.numberOfMoleculesPerComponent);
  table.bonds.push_back(PolymerBond{PolymerAtomRef{0, 0, 0}, PolymerAtomRef{1, 0, 0}, 0, 0, 1.5});
  table.bonds.push_back(PolymerBond{PolymerAtomRef{0, 0, 1}, PolymerAtomRef{1, 1, 0}, 0, 0, 1.5});

  // intact initiator: one cluster of three nodes
  PolymerLinkTable::NetworkStatistics intact = table.networkStatistics(system.simulationBox, system.components, atoms);
  EXPECT_EQ(intact.moleculesInNetwork, 3uz);
  EXPECT_EQ(intact.numberOfClusters, 1uz);
  EXPECT_EQ(intact.largestCluster, 3uz);
  EXPECT_EQ(intact.cycleRank, 0uz);
  EXPECT_EQ(intact.spanningClusters, 0uz);

  // split initiator: its halves are two nodes, each in its own cluster, and no cycle runs through the molecule
  table.brokenBonds.push_back(BrokenBond{PolymerAtomRef{0, 0, 0}, PolymerAtomRef{0, 0, 1}, 1, 0});
  PolymerLinkTable::NetworkStatistics split = table.networkStatistics(system.simulationBox, system.components, atoms);
  EXPECT_EQ(split.moleculesInNetwork, 4uz);
  EXPECT_EQ(split.numberOfClusters, 2uz);
  EXPECT_EQ(split.largestCluster, 2uz);
  EXPECT_EQ(split.cycleRank, 0uz);
  EXPECT_EQ(split.spanningClusters, 0uz);
  EXPECT_FALSE(split.percolates[0]);

  // a bond between the two monomers closes the chain half-monomer-monomer-half: with the halves as separate nodes
  // it is a chain (no cycle), with the molecule as one node it would be a ring
  table.bonds.push_back(PolymerBond{PolymerAtomRef{1, 0, 0}, PolymerAtomRef{1, 1, 0}, 0, 0, 1.5});
  atoms[3].position = double3(19.5, 10.0, 10.0);
  atoms[1].position = double3(18.0, 10.0, 10.0);
  PolymerLinkTable::NetworkStatistics ring = table.networkStatistics(system.simulationBox, system.components, atoms);
  EXPECT_EQ(ring.numberOfClusters, 1uz);
  EXPECT_EQ(ring.cycleRank, 0uz);
  EXPECT_EQ(ring.spanningClusters, 0uz);
}
