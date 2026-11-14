#include <gtest/gtest.h>

#include "../test_support.hpp"
#include "molecule_fixtures.hpp"

import std;

import int3;
import archive;
import atom;
import forcefield;
import pseudo_atom;
import framework;
import component;
import system;
import simulationbox;
import write_lammps_data;
import sample_movies;
import mc_moves_move_types;
import mc_moves_probabilities;
import json;
import monte_carlo;
import molecular_dynamics;

// Binary restart round trips of the state that a restart read used to lose: the unit-cell data of a framework,
// the mixing rule of a force field, the source of a pseudo-atom, the file and frame bookkeeping of the movie and
// LAMMPS writers, and the conformation reservoirs a System rebuilds on read (they are not stored).

namespace
{

template <typename T>
void roundTrip(const T& original, T& restored, std::string_view name)
{
  const std::filesystem::path path =
      std::filesystem::temp_directory_path() / std::format("raspa3_{}_{}.bin", name, test_support::uniqueSuffix());
  {
    std::ofstream stream(path, std::ios::binary);
    Archive<std::ofstream> archive(stream);
    archive << original;
  }
  {
    std::ifstream stream(path, std::ios::binary);
    Archive<std::ifstream> archive(stream);
    archive >> restored;
  }
  std::filesystem::remove(path);
}

}  // namespace

TEST(BINARY_RESTART_ROUND_TRIP, framework_keeps_its_unit_cell_data)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  const Framework original = Framework::makeFAU(forceField, int3(1, 1, 1));
  ASSERT_GT(original.unitCellMass, 0.0);
  Framework restored;
  roundTrip(original, restored, "framework");
  EXPECT_DOUBLE_EQ(restored.unitCellMass, original.unitCellMass);
  EXPECT_EQ(restored.fractionalUnitCellAtoms.size(), original.fractionalUnitCellAtoms.size());
  EXPECT_EQ(restored.uniqueAtomTypes, original.uniqueAtomTypes);
}

TEST(BINARY_RESTART_ROUND_TRIP, force_field_keeps_its_mixing_rule_and_pseudo_atoms_their_source)
{
  ForceField original = ForceField::makeZeoliteForceField(12.0, true, false, true);
  original.mixingRule = ForceField::MixingRule::Jorgensen;
  original.pseudoAtoms.front().source = "a reference";
  ForceField restored;
  roundTrip(original, restored, "force_field");
  EXPECT_EQ(restored.mixingRule, ForceField::MixingRule::Jorgensen);
  EXPECT_EQ(restored.pseudoAtoms.front().source, "a reference");
  EXPECT_TRUE(restored.pseudoAtoms.front() == original.pseudoAtoms.front());

  // equality compares the polarizability and the oxidation state as well
  PseudoAtom other = original.pseudoAtoms.front();
  other.polarizability += 1.0;
  EXPECT_FALSE(other == original.pseudoAtoms.front());
}

TEST(BINARY_RESTART_ROUND_TRIP, movie_and_lammps_writers_keep_their_bookkeeping)
{
  SampleMovie movie(3, 25, true, std::string("run-a"));
  movie.modelNumber = 17;
  SampleMovie restoredMovie;
  roundTrip(movie, restoredMovie, "movie");
  EXPECT_EQ(restoredMovie.sampleEvery, 25uz);
  EXPECT_TRUE(restoredMovie.restrictToBox);
  EXPECT_EQ(restoredMovie.tag, std::optional<std::string>("run-a"));
  EXPECT_EQ(restoredMovie.modelNumber, 17);

  WriteLammpsData lammps(4, 50);
  lammps.modelNumber = 9;
  WriteLammpsData restoredLammps;
  roundTrip(lammps, restoredLammps, "lammps");
  EXPECT_EQ(restoredLammps.sampleEvery, 50uz);
  EXPECT_EQ(restoredLammps.systemId, 4uz);
  EXPECT_EQ(restoredLammps.modelNumber, 9);
}

TEST(BINARY_RESTART_ROUND_TRIP, system_rebuilds_its_conformation_reservoirs)
{
  const ForceField forceField({{"CH2_c", false, 14.02658, 0.0, 0.0, 6, false}}, {{52.5, 3.91}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, true, false);
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::ReinsertionCBMC, 1.0);
  TemporaryFile file("cyclohexane.json", molecule_fixtures::kCyclohexaneJson);
  const Component cyclohexane(Component::Type::Adsorbate, 0, forceField, "cyclohexane", file.stemPath().string(), 5, 21,
                              probabilities, std::nullopt, false);
  const System original(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {cyclohexane}, {}, {2},
                        5);
  const std::vector<std::vector<Atom>>& reservoir = original.components[0].conformationReservoir;
  ASSERT_FALSE(reservoir.empty());

  System restored;
  roundTrip(original, restored, "system");
  const std::vector<std::vector<Atom>>& restoredReservoir = restored.components[0].conformationReservoir;
  ASSERT_EQ(restoredReservoir.size(), reservoir.size());
  for (std::size_t atom = 0; atom < reservoir.front().size(); ++atom)
  {
    EXPECT_EQ(restoredReservoir.front()[atom].position, reservoir.front()[atom].position) << "atom " << atom;
  }
}

// The options a flexible framework regenerates its van der Waals pair list and intramolecular topology from: they
// were not stored, so a restart regenerated the pair list with the default exclusions and 1-4 scalings.
TEST(BINARY_RESTART_ROUND_TRIP, framework_keeps_its_intramolecular_generation_options)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework original = Framework::makeFAU(forceField, int3(1, 1, 1));
  original.excludeIntra12Interactions = false;
  original.excludeIntra13Interactions = false;
  original.excludeIntraBondInteractions = true;
  original.excludeIntraBendInteractions = true;
  original.intra14VanDerWaalsScaling = 0.5;
  original.intra14ChargeChargeScaling = 0.8333;
  original.frameworkDefinitionName = "fau_flexible";
  original.bondDefinitions = {FrameworkPotentialDefinition<2>{{0, 1}, "HARMONIC", {1000.0, 1.6}}};
  original.bendDefinitions = {FrameworkPotentialDefinition<3>{{1, 0, 1}, "HARMONIC", {500.0, 109.5}}};
  original.torsionDefinitions = {FrameworkPotentialDefinition<4>{{1, 0, 1, 0}, "SIX_TERM_DIHEDRAL", {1.0, 2.0}}};
  original.improperTorsionDefinitions = {};

  Framework restored;
  roundTrip(original, restored, "framework_generation_options");
  EXPECT_FALSE(restored.excludeIntra12Interactions);
  EXPECT_FALSE(restored.excludeIntra13Interactions);
  EXPECT_TRUE(restored.excludeIntraBondInteractions);
  EXPECT_TRUE(restored.excludeIntraBendInteractions);
  EXPECT_EQ(restored.intra14VanDerWaalsScaling, 0.5);
  EXPECT_EQ(restored.intra14ChargeChargeScaling, 0.8333);
  EXPECT_EQ(restored.frameworkDefinitionName, "fau_flexible");
  ASSERT_EQ(restored.bondDefinitions.size(), 1uz);
  EXPECT_EQ(restored.bondDefinitions[0].atomTypes, (std::array<std::size_t, 2>{0, 1}));
  EXPECT_EQ(restored.bondDefinitions[0].potentialType, "HARMONIC");
  EXPECT_EQ(restored.bondDefinitions[0].parameters, (std::vector<double>{1000.0, 1.6}));
  ASSERT_EQ(restored.bendDefinitions.size(), 1uz);
  EXPECT_EQ(restored.bendDefinitions[0].atomTypes, (std::array<std::size_t, 3>{1, 0, 1}));
  ASSERT_EQ(restored.torsionDefinitions.size(), 1uz);
  EXPECT_EQ(restored.torsionDefinitions[0].parameters, (std::vector<double>{1.0, 2.0}));
  EXPECT_TRUE(restored.improperTorsionDefinitions.empty());
}

// The JSON output of a run: its header (version, seed, initial conditions) is written by setup(), which a resumed run
// skips, so the final output of a restarted MonteCarlo or MolecularDynamics run had lost it.
TEST(BINARY_RESTART_ROUND_TRIP, drivers_keep_the_header_of_their_json_output)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  const Component methane = Component::makeMethane(forceField, 0);
  const std::vector<System> systems{
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {methane}, {}, {2}, 5)};

  MonteCarlo mc({0, 0, 0, 0, 1000, 10000, 5000, 5000}, systems, 42uz, 5, false);
  mc.outputJsons[0]["version"] = "3.1.0";
  mc.outputJsons[0]["seed"] = 42;
  MonteCarlo restoredMC;
  roundTrip(mc, restoredMC, "mc_json");
  ASSERT_EQ(restoredMC.outputJsons.size(), 1uz);
  EXPECT_EQ(restoredMC.outputJsons[0]["version"], "3.1.0");
  EXPECT_EQ(restoredMC.outputJsons[0]["seed"], 42);

  MolecularDynamics md({0, 0, 0, 0, 1000, 10000, 5000, 5000}, systems, 42uz, 5, false);
  md.outputJsons[0]["initialization"]["units"] = "K";
  MolecularDynamics restoredMD;
  roundTrip(md, restoredMD, "md_json");
  ASSERT_EQ(restoredMD.outputJsons.size(), 1uz);
  EXPECT_EQ(restoredMD.outputJsons[0]["initialization"]["units"], "K");
}
