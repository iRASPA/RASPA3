#include <gtest/gtest.h>

#include "../test_support.hpp"
#include "molecule_fixtures.hpp"

import std;

import double3;
import atom;
import molecule;
import forcefield;
import component;
import system;
import monte_carlo;
import simulationbox;
import running_energy;
import mc_moves_move_types;
import mc_moves_probabilities;

// The persistent per-atom tags: always a permutation of 0 ... N-1, handed out by insertion, kept dense by
// deletion (the atoms holding the top tags take the freed ones), following the atoms through the exchanges
// within a system and between replicas, and staying with the slot in an exchange between systems.

namespace
{

ForceField makeAlkaneForceField()
{
  return ForceField({{"CH3", false, 15.04, 0.0, 0.0, 6, false}, {"CH2", false, 14.03, 0.0, 0.0, 6, false}},
                    {{98.0, 3.75}, {46.0, 3.95}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true,
                    true, false);
}

Component makeAlkaneFromExample(const ForceField& forceField, std::size_t componentId, std::string_view name,
                                const MCMoveProbabilities& probabilities)
{
  TemporaryFile file(std::string(name) + ".json", molecule_fixtures::alkaneJson(name));
  return Component(Component::Type::Adsorbate, componentId, forceField, std::string(name), file.stemPath().string(), 5,
                   21, probabilities, std::nullopt, false);
}

// 20 propane (3 atoms) followed by 10 butane (4 atoms): 100 atoms, no framework
System makePropaneButaneSystem(const MCMoveProbabilities& componentProbabilities = MCMoveProbabilities())
{
  const ForceField forceField = makeAlkaneForceField();
  Component propane = makeAlkaneFromExample(forceField, 0, "propane", componentProbabilities);
  Component butane = makeAlkaneFromExample(forceField, 1, "butane", componentProbabilities);
  propane.lnPartitionFunction = 1.0;
  butane.lnPartitionFunction = 1.0;
  return System(forceField, SimulationBox(40.0, 40.0, 40.0), false, 500.0, 1e4, 1.0, {}, {propane, butane}, {},
                std::vector<std::size_t>{20, 10}, 5, MCMoveProbabilities());
}

std::vector<std::uint32_t> tagsOfMolecule(const System& system, std::size_t component, std::size_t molecule)
{
  const std::span<const Atom> atoms = system.spanOfMolecule(component, molecule);
  const std::size_t first = static_cast<std::size_t>(atoms.data() - system.atomData.data());
  return {system.atomTags.begin() + static_cast<std::ptrdiff_t>(first),
          system.atomTags.begin() + static_cast<std::ptrdiff_t>(first + atoms.size())};
}

}  // namespace

TEST(ATOM_TAGS, construction_numbers_the_atoms_in_storage_order)
{
  const System system = makePropaneButaneSystem();

  ASSERT_EQ(system.atomData.size(), 100uz);
  ASSERT_EQ(system.atomTags.size(), 100uz);
  for (std::size_t i = 0; i < system.atomTags.size(); ++i)
  {
    EXPECT_EQ(system.atomTags[i], static_cast<std::uint32_t>(i));
    EXPECT_EQ(system.atomIndexOfTag[i], static_cast<std::uint32_t>(i));
  }
  EXPECT_TRUE(system.retaggedAtoms.empty());
  EXPECT_NO_THROW(system.checkAtomTags());
}

TEST(ATOM_TAGS, deletion_hands_the_freed_tags_to_the_atoms_holding_the_top_ones)
{
  System system = makePropaneButaneSystem();

  // propane 5 occupies atoms 15, 16, 17; the last butane holds the top tags 96 ... 99
  system.deleteMolecule(0, 5, system.spanOfMolecule(0, 5));

  ASSERT_EQ(system.atomTags.size(), 97uz);
  EXPECT_NO_THROW(system.checkAtomTags());
  const std::vector<std::pair<std::uint32_t, std::uint32_t>> expectedFirst{{97, 15}, {98, 16}, {99, 17}};
  EXPECT_EQ(system.retaggedAtoms, expectedFirst);
  EXPECT_EQ(tagsOfMolecule(system, 1, 9), (std::vector<std::uint32_t>{96, 15, 16, 17}));
  // every other atom keeps its tag: the propanes before the hole and the butanes before the last one
  for (std::size_t i = 0; i < 15; ++i) EXPECT_EQ(system.atomTags[i], static_cast<std::uint32_t>(i));
  for (std::size_t i = 15; i < 93; ++i) EXPECT_EQ(system.atomTags[i], static_cast<std::uint32_t>(i + 3));

  // deleting that butane frees 96 (out of range, dropped) and 15, 16, 17, which go to the holders of 93, 94, 95
  system.deleteMolecule(1, 9, system.spanOfMolecule(1, 9));

  ASSERT_EQ(system.atomTags.size(), 93uz);
  EXPECT_NO_THROW(system.checkAtomTags());
  const std::vector<std::pair<std::uint32_t, std::uint32_t>> expectedSecond{{97, 15}, {98, 16}, {99, 17},
                                                                            {93, 15}, {94, 16}, {95, 17}};
  EXPECT_EQ(system.retaggedAtoms, expectedSecond);
  EXPECT_EQ(tagsOfMolecule(system, 1, 8), (std::vector<std::uint32_t>{92, 15, 16, 17}));

  // an insertion takes the next tags and shifts the butanes without renaming them
  const std::vector<std::uint32_t> butane0Before = tagsOfMolecule(system, 1, 0);
  std::vector<Atom> atoms(system.spanOfMolecule(0, 0).begin(), system.spanOfMolecule(0, 0).end());
  for (Atom& atom : atoms) atom.position += double3(1.0, 1.0, 1.0);
  const Molecule molecule = system.moleculeData[0];
  system.insertMolecule(0, molecule, atoms);

  ASSERT_EQ(system.atomTags.size(), 96uz);
  EXPECT_NO_THROW(system.checkAtomTags());
  EXPECT_EQ(tagsOfMolecule(system, 0, 19), (std::vector<std::uint32_t>{93, 94, 95}));
  EXPECT_EQ(tagsOfMolecule(system, 1, 0), butane0Before);
  EXPECT_EQ(system.retaggedAtoms.size(), 6uz);
}

TEST(ATOM_TAGS, swap_within_a_system_moves_the_tags_with_the_atoms)
{
  System system = makePropaneButaneSystem();
  const std::vector<std::uint32_t> tags2 = tagsOfMolecule(system, 0, 2);
  const std::vector<std::uint32_t> tags7 = tagsOfMolecule(system, 0, 7);
  const double3 position2 = system.spanOfMolecule(0, 2).front().position;
  const double3 position7 = system.spanOfMolecule(0, 7).front().position;

  system.swapMolecules(0, 2, 7);

  EXPECT_EQ(tagsOfMolecule(system, 0, 7), tags2);
  EXPECT_EQ(tagsOfMolecule(system, 0, 2), tags7);
  EXPECT_EQ(system.spanOfMolecule(0, 7).front().position, position2);
  EXPECT_EQ(system.spanOfMolecule(0, 2).front().position, position7);
  EXPECT_EQ(system.atomIndexOfTag[tags2.front()], 21u);
  EXPECT_EQ(system.atomIndexOfTag[tags7.front()], 6u);
  EXPECT_NO_THROW(system.checkAtomTags());
  EXPECT_TRUE(system.retaggedAtoms.empty());
}

TEST(ATOM_TAGS, exchange_between_systems_leaves_the_tags_with_the_slots)
{
  System systemA = makePropaneButaneSystem();
  System systemB = makePropaneButaneSystem();
  systemB.deleteMolecule(0, 0, systemB.spanOfMolecule(0, 0));
  const std::vector<std::uint32_t> tagsA = systemA.atomTags;
  const std::vector<std::uint32_t> tagsB = systemB.atomTags;
  const double3 positionA = systemA.spanOfMolecule(1, 3).front().position;
  const double3 positionB = systemB.spanOfMolecule(1, 4).front().position;

  systemA.exchangeMolecule(systemB, 1, 3, 4);

  EXPECT_EQ(systemA.spanOfMolecule(1, 3).front().position, positionB);
  EXPECT_EQ(systemB.spanOfMolecule(1, 4).front().position, positionA);
  EXPECT_EQ(systemA.atomTags, tagsA);
  EXPECT_EQ(systemB.atomTags, tagsB);
  EXPECT_NO_THROW(systemA.checkAtomTags());
  EXPECT_NO_THROW(systemB.checkAtomTags());

  // undone by the same call
  systemA.exchangeMolecule(systemB, 1, 3, 4);
  EXPECT_EQ(systemA.spanOfMolecule(1, 3).front().position, positionA);
  EXPECT_EQ(systemB.spanOfMolecule(1, 4).front().position, positionB);
  EXPECT_EQ(systemA.atomTags, tagsA);
  EXPECT_EQ(systemB.atomTags, tagsB);
}

TEST(ATOM_TAGS, replica_exchange_moves_the_tags_with_the_configuration)
{
  System systemA = makePropaneButaneSystem();
  System systemB = makePropaneButaneSystem();
  // different counts and a history in B: its tags are no longer the identity
  systemB.deleteMolecule(0, 3, systemB.spanOfMolecule(0, 3));
  systemB.deleteMolecule(1, 1, systemB.spanOfMolecule(1, 1));
  const std::vector<std::uint32_t> tagsA = systemA.atomTags;
  const std::vector<std::uint32_t> tagsB = systemB.atomTags;
  const std::vector<Atom> atomsA = systemA.atomData;
  const std::vector<Atom> atomsB = systemB.atomData;
  std::vector<std::uint32_t> identity(tagsB.size());
  std::iota(identity.begin(), identity.end(), std::uint32_t{0});
  ASSERT_NE(tagsB, identity);

  systemA.exchangeConfiguration(systemB);

  ASSERT_EQ(systemA.atomData.size(), atomsB.size());
  ASSERT_EQ(systemB.atomData.size(), atomsA.size());
  EXPECT_EQ(systemA.atomTags, tagsB);
  EXPECT_EQ(systemB.atomTags, tagsA);
  for (std::size_t i = 0; i < atomsB.size(); ++i) EXPECT_EQ(systemA.atomData[i].position, atomsB[i].position);
  EXPECT_NO_THROW(systemA.checkAtomTags());
  EXPECT_NO_THROW(systemB.checkAtomTags());
  EXPECT_EQ(systemA.numberOfMoleculesPerComponent, (std::vector<std::size_t>{19, 9}));
  EXPECT_EQ(systemB.numberOfMoleculesPerComponent, (std::vector<std::size_t>{20, 10}));
}

TEST(ATOM_TAGS, short_run_with_insertions_and_deletions_keeps_the_tags_dense)
{
  MCMoveProbabilities componentProbabilities = MCMoveProbabilities();
  componentProbabilities.setProbability(Move::Types::Translation, 1.0);
  componentProbabilities.setProbability(Move::Types::Rotation, 1.0);
  componentProbabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  componentProbabilities.setProbability(Move::Types::SwapCBMC, 1.0);
  componentProbabilities.setProbability(Move::Types::SwapCFCMC, 1.0);
  System system = makePropaneButaneSystem(componentProbabilities);
  system.runningEnergies = system.computeTotalEnergies();

  MonteCarlo mc = MonteCarlo({20, 0, 5, 5, 100, 10000, 5000, 5000}, {std::move(system)}, {}, 5, false);
  mc.run();

  for (const System& s : mc.systems)
  {
    EXPECT_NO_THROW(s.checkAtomTags());
    EXPECT_FALSE(s.retaggedAtoms.empty());
    for (const auto& [oldTag, newTag] : s.retaggedAtoms)
    {
      EXPECT_LT(newTag, oldTag);
    }
  }
}
