#include <gtest/gtest.h>

#include "../test_support.hpp"
#include "molecule_fixtures.hpp"

import std;

import archive;
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
import randomnumbers;
import mc_moves_gibbs_swap_cfcmc;

// The persistent per-atom tags: fixed for the life of an atom, handed out lowest-free-first by insertion, freed
// by deletion, following the atoms through the exchanges within a system and between replicas, and staying with
// the slot in an exchange between systems. 'atomIndexOfTag' is the inverse table, as long as the highest tag in
// use plus one, 'System::noAtom' marking a free tag.

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
  EXPECT_EQ(system.atomIndexOfTag.size(), 100uz);
  EXPECT_NO_THROW(system.checkAtomTags());
}

TEST(ATOM_TAGS, deletion_frees_the_tags_and_every_other_atom_keeps_its_own)
{
  System system = makePropaneButaneSystem();

  // propane 5 occupies atoms 15, 16, 17 and holds the tags 15, 16, 17
  system.deleteMolecule(0, 5, system.spanOfMolecule(0, 5));

  ASSERT_EQ(system.atomTags.size(), 97uz);
  EXPECT_NO_THROW(system.checkAtomTags());
  // the table keeps its length: the highest tag, 99, is still in use; the freed slots are marked free
  EXPECT_EQ(system.atomIndexOfTag.size(), 100uz);
  for (std::uint32_t tag : {15u, 16u, 17u}) EXPECT_EQ(system.atomIndexOfTag[tag], System::noAtom);
  // every surviving atom keeps its tag: the propanes before the hole and everything after it, shifted down
  for (std::size_t i = 0; i < 15; ++i) EXPECT_EQ(system.atomTags[i], static_cast<std::uint32_t>(i));
  for (std::size_t i = 15; i < 97; ++i) EXPECT_EQ(system.atomTags[i], static_cast<std::uint32_t>(i + 3));
  EXPECT_EQ(tagsOfMolecule(system, 1, 9), (std::vector<std::uint32_t>{96, 97, 98, 99}));

  // deleting the last butane frees the top tags 96 ... 99, and the table shrinks to the highest tag in use
  system.deleteMolecule(1, 9, system.spanOfMolecule(1, 9));

  ASSERT_EQ(system.atomTags.size(), 93uz);
  EXPECT_NO_THROW(system.checkAtomTags());
  EXPECT_EQ(system.atomIndexOfTag.size(), 96uz);
  EXPECT_EQ(tagsOfMolecule(system, 1, 8), (std::vector<std::uint32_t>{92, 93, 94, 95}));

  // an insertion takes the lowest free tags, 15, 16, 17, and shifts the butanes without renaming them
  const std::vector<std::uint32_t> butane0Before = tagsOfMolecule(system, 1, 0);
  std::vector<Atom> atoms(system.spanOfMolecule(0, 0).begin(), system.spanOfMolecule(0, 0).end());
  for (Atom& atom : atoms) atom.position += double3(1.0, 1.0, 1.0);
  const Molecule molecule = system.moleculeData[0];
  system.insertMolecule(0, molecule, atoms);

  ASSERT_EQ(system.atomTags.size(), 96uz);
  EXPECT_NO_THROW(system.checkAtomTags());
  EXPECT_EQ(system.atomIndexOfTag.size(), 96uz);
  EXPECT_EQ(tagsOfMolecule(system, 0, 19), (std::vector<std::uint32_t>{15, 16, 17}));
  EXPECT_EQ(tagsOfMolecule(system, 1, 0), butane0Before);
  EXPECT_EQ(system.atomIndexOfTag[15], 57u);

  // with no free tag left, the next insertion extends the table
  system.insertMolecule(0, molecule, atoms);

  ASSERT_EQ(system.atomTags.size(), 99uz);
  EXPECT_NO_THROW(system.checkAtomTags());
  EXPECT_EQ(system.atomIndexOfTag.size(), 99uz);
  EXPECT_EQ(tagsOfMolecule(system, 0, 20), (std::vector<std::uint32_t>{96, 97, 98}));
}

TEST(ATOM_TAGS, a_tag_follows_its_atom_through_a_history_of_deletions_and_insertions)
{
  System system = makePropaneButaneSystem();
  // remember every atom by its tag
  std::vector<double3> positionOfTag(system.atomData.size());
  for (std::size_t i = 0; i < system.atomData.size(); ++i)
    positionOfTag[system.atomTags[i]] = system.atomData[i].position;

  std::vector<Atom> atoms(system.spanOfMolecule(1, 2).begin(), system.spanOfMolecule(1, 2).end());
  const Molecule molecule = system.moleculeData[22];
  system.deleteMolecule(0, 3, system.spanOfMolecule(0, 3));
  system.deleteMolecule(1, 7, system.spanOfMolecule(1, 7));
  system.insertMolecule(1, molecule, atoms);
  system.deleteMolecule(0, 0, system.spanOfMolecule(0, 0));
  system.insertMolecule(1, molecule, atoms);
  system.insertMolecule(1, molecule, atoms);
  ASSERT_EQ(system.atomData.size(), 102uz);
  EXPECT_NO_THROW(system.checkAtomTags());

  // the tags 0 ... 2 (propane 0), 9 ... 11 (propane 3) and 88 ... 91 (butane 7) were freed: 10 tags, 12 atoms
  // inserted, so the table has grown by 2
  EXPECT_EQ(system.atomIndexOfTag.size(), 102uz);
  const std::set<std::uint32_t> freed{0, 1, 2, 9, 10, 11, 88, 89, 90, 91};
  for (std::uint32_t tag = 0; tag < 100; ++tag)
  {
    if (freed.contains(tag)) continue;
    const std::uint32_t index = system.atomIndexOfTag[tag];
    ASSERT_LT(index, system.atomData.size());
    EXPECT_EQ(system.atomTags[index], tag);
    EXPECT_EQ(system.atomData[index].position, positionOfTag[tag]);
  }
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
  EXPECT_EQ(systemA.atomIndexOfTag.size(), 100uz);
  EXPECT_EQ(systemB.atomIndexOfTag.size(), 100uz);
  for (std::size_t i = 0; i < atomsB.size(); ++i) EXPECT_EQ(systemA.atomData[i].position, atomsB[i].position);
  EXPECT_NO_THROW(systemA.checkAtomTags());
  EXPECT_NO_THROW(systemB.checkAtomTags());
  EXPECT_EQ(systemA.numberOfMoleculesPerComponent, (std::vector<std::size_t>{19, 9}));
  EXPECT_EQ(systemB.numberOfMoleculesPerComponent, (std::vector<std::size_t>{20, 10}));
}

TEST(ATOM_TAGS, restart_rebuilds_the_inverse_table_with_its_free_tags)
{
  System original = makePropaneButaneSystem();
  original.deleteMolecule(0, 3, original.spanOfMolecule(0, 3));
  original.deleteMolecule(1, 1, original.spanOfMolecule(1, 1));
  ASSERT_EQ(original.atomIndexOfTag.size(), 100uz);

  const std::filesystem::path path =
      std::filesystem::temp_directory_path() / std::format("raspa3_atom_tags_{}.bin", test_support::uniqueSuffix());
  {
    std::ofstream stream(path, std::ios::binary);
    Archive<std::ofstream> archive(stream);
    archive << original;
  }
  System restored = makePropaneButaneSystem();
  {
    std::ifstream stream(path, std::ios::binary);
    Archive<std::ifstream> archive(stream);
    archive >> restored;
  }
  std::filesystem::remove(path);

  EXPECT_EQ(restored.atomTags, original.atomTags);
  EXPECT_EQ(restored.atomIndexOfTag, original.atomIndexOfTag);
  EXPECT_NO_THROW(restored.checkAtomTags());
}

TEST(ATOM_TAGS, short_run_with_insertions_and_deletions_keeps_the_tags_consistent)
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
    // the table is bounded by the largest atom count the system has held, never by the number of insertions
    EXPECT_GE(s.atomIndexOfTag.size(), s.atomData.size());
    EXPECT_LE(s.atomIndexOfTag.size(), 2 * s.atomData.size() + 100);
  }
}

// The lambda interchange of the serial Gibbs CFCMC move: the fractional molecule of A becomes a whole molecule
// where it is, a whole molecule of B becomes the fractional molecule where it is, and the inactive placeholder
// moves from B to A. Only the placeholder changes box, so every other atom keeps its tag and its position.
TEST(ATOM_TAGS, gibbs_cfcmc_lambda_interchange_keeps_the_tags_of_the_atoms_that_stay)
{
  const ForceField forceField = makeAlkaneForceField();
  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::GibbsSwapCFCMC, 1.0);
  const Component propane = makeAlkaneFromExample(forceField, 0, "propane", probabilities);
  const SimulationBox box = SimulationBox(40.0, 40.0, 40.0);
  std::vector<System> systems{System(forceField, box, false, 500.0, 1e4, 1.0, {}, {propane}, {},
                                     std::vector<std::size_t>{2}, 5, MCMoveProbabilities()),
                              System(forceField, box, false, 500.0, 1e4, 1.0, {}, {propane}, {},
                                     std::vector<std::size_t>{20}, 5, MCMoveProbabilities())};
  MonteCarlo mc = MonteCarlo({0, 0, 0, 0, 1000, 10000, 5000, 5000}, systems, 42uz, 5, false);
  mc.setup();
  systems = mc.systems;
  for (System& system : systems)
  {
    system.runningEnergies = system.computeTotalEnergies();
  }

  RandomNumber random(42);
  std::size_t fractionalMoleculeSystem = 0;
  for (int attempt = 0; attempt < 5000; ++attempt)
  {
    System& systemA = systems[0].containsTheFractionalMolecule ? systems[0] : systems[1];
    System& systemB = systems[0].containsTheFractionalMolecule ? systems[1] : systems[0];
    const std::size_t countA = systemA.numberOfMoleculesPerComponent[0];
    std::map<std::uint32_t, double3> positionOfTagA;
    std::map<std::uint32_t, double3> positionOfTagB;
    for (std::size_t i = 0; i < systemA.atomTags.size(); ++i)
    {
      positionOfTagA[systemA.atomTags[i]] = systemA.atomData[i].position;
    }
    for (std::size_t i = 0; i < systemB.atomTags.size(); ++i)
    {
      positionOfTagB[systemB.atomTags[i]] = systemB.atomData[i].position;
    }
    const std::vector<std::uint32_t> placeholderTags =
        tagsOfMolecule(systemB, 0, systemB.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCFCMC, 0));

    MC_Moves::GibbsSwapMove_CFCMC(random, systemA, systemB, 0, fractionalMoleculeSystem);

    // the other branches (shuffle, lambda change) leave the molecule counts alone
    if (systemA.numberOfMoleculesPerComponent[0] != countA + 1)
    {
      continue;
    }

    ASSERT_NO_THROW(systemA.checkAtomTags());
    ASSERT_NO_THROW(systemB.checkAtomTags());
    ASSERT_EQ(systemA.atomTags.size(), positionOfTagA.size() + 3);
    ASSERT_EQ(systemB.atomTags.size(), positionOfTagB.size() - 3);
    for (const auto& [tag, position] : positionOfTagA)
    {
      const std::uint32_t index = systemA.atomIndexOfTag[tag];
      ASSERT_NE(index, System::noAtom) << "tag " << tag << " of A was freed";
      EXPECT_EQ(systemA.atomData[index].position, position) << "tag " << tag << " of A names another atom";
    }
    for (const std::uint32_t tag : placeholderTags)
    {
      EXPECT_EQ(systemB.atomIndexOfTag[tag], System::noAtom) << "placeholder tag " << tag << " is still in use";
      positionOfTagB.erase(tag);
    }
    for (const auto& [tag, position] : positionOfTagB)
    {
      const std::uint32_t index = systemB.atomIndexOfTag[tag];
      ASSERT_NE(index, System::noAtom) << "tag " << tag << " of B was freed";
      EXPECT_EQ(systemB.atomData[index].position, position) << "tag " << tag << " of B names another atom";
    }
    return;
  }
  FAIL() << "No accepted lambda interchange in 5000 attempts";
}
