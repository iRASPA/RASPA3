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
