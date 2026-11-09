#include <gtest/gtest.h>

import std;

import double3;
import int3;
import atom;
import pseudo_atom;
import vdwparameters;
import forcefield;
import component;
import molecule;
import simulationbox;
import system;
import averages;
import cif_reader;
import charge_equilibration_wilmer_snurr;

TEST(UNINITIALIZED_MEMBERS, default_constructed_atom_has_zero_charge_group_and_fractional_flag)
{
  // construct on top of non-zero bytes so a member without a default initializer shows up
  alignas(Atom) std::array<std::byte, sizeof(Atom)> storage;
  storage.fill(std::byte{0xff});
  const Atom* atom = new (storage.data()) Atom;
  EXPECT_EQ(atom->charge, 0.0);
  EXPECT_EQ(atom->groupId, 0);
  EXPECT_EQ(atom->isFractional, 0);
}

TEST(UNINITIALIZED_MEMBERS, flexible_molecule_from_initial_positions_has_its_center_of_mass)
{
  ForceField forceField({{"O", false, 16.0, 0.0, 0.0, 8, false}, {"H", false, 1.0, 0.0, 0.0, 1, false}},
                        {{80.0, 3.0}, {10.0, 2.0}}, ForceField::MixingRule::Lorentz_Berthelot, 8.0, 8.0, 8.0, false,
                        false, false);
  const std::uint16_t oxygen = static_cast<std::uint16_t>(*forceField.findPseudoAtom("O"));
  const std::uint16_t hydrogen = static_cast<std::uint16_t>(*forceField.findPseudoAtom("H"));
  Component diatomic(forceField, "diatomic", 64.0, 1.0e6, 0.02,
                     {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, oxygen, 0, false, false),
                      Atom({1.0, 0.0, 0.0}, 0.0, 1.0, 0, hydrogen, 0, false, false)},
                     {}, {}, 5, 21);
  diatomic.rigid = false;

  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {diatomic},
                {{double3(5.0, 6.0, 7.0), double3(6.7, 6.0, 7.0)}}, {0}, 5);

  ASSERT_EQ(system.moleculeData.size(), 1uz);
  const double3 expected = (16.0 * double3(5.0, 6.0, 7.0) + 1.0 * double3(6.7, 6.0, 7.0)) / 17.0;
  EXPECT_NEAR(system.moleculeData[0].centerOfMassPosition.x, expected.x, 1.0e-12);
  EXPECT_NEAR(system.moleculeData[0].centerOfMassPosition.y, expected.y, 1.0e-12);
  EXPECT_NEAR(system.moleculeData[0].centerOfMassPosition.z, expected.z, 1.0e-12);
}

TEST(UNINITIALIZED_MEMBERS, block_error_beyond_the_t_table_uses_the_t_quantile)
{
  // 95% two-sided t critical values: the table up to 20 degrees of freedom, the expansion beyond
  EXPECT_DOUBLE_EQ(studentTCriticalValue(20), 2.086);
  EXPECT_NEAR(studentTCriticalValue(21), 2.080, 1.0e-3);
  EXPECT_NEAR(studentTCriticalValue(29), 2.045, 1.0e-3);
  EXPECT_NEAR(studentTCriticalValue(99), 1.984, 1.0e-3);

  // 30 blocks, alternating around the mean: the error is t(29) * s / sqrt(30)
  std::vector<double> blocks(30);
  for (std::size_t i = 0; i < blocks.size(); ++i) blocks[i] = (i % 2 == 0) ? 1.0 : -1.0;
  const double error = blockErrorEstimate<double>(blocks, 0.0);
  EXPECT_NEAR(error, studentTCriticalValue(29) * std::sqrt(30.0 / 29.0) / std::sqrt(30.0), 1.0e-12);
}

TEST(UNINITIALIZED_MEMBERS, ddl2_cif_atom_sites_are_read)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, false);
  constexpr std::string_view cif = R"(
data_test
_cell.length_a 10.0
_cell.length_b 11.0
_cell.length_c 12.0
_cell.angle_alpha 90.0
_cell.angle_beta 90.0
_cell.angle_gamma 90.0
loop_
_atom_site.label
_atom_site.type_symbol
_atom_site.fract_x
_atom_site.fract_y
_atom_site.fract_z
_atom_site.charge
Si1 Si 0.1 0.2 0.3 1.5
O1  O  0.4 0.5 0.6 -0.75
)";

  const auto parsed = CIFReader::readCIFString(std::string(cif), forceField, CIFReader::UseChargesFrom::CIF_File);
  ASSERT_TRUE(parsed.has_value());
  const auto& [simulationBox, spaceGroupHallNumber, definedAtoms, unitCellAtoms] = parsed.value();
  EXPECT_NEAR(simulationBox.lengthB, 11.0, 1.0e-12);
  ASSERT_EQ(definedAtoms.size(), 2uz);
  EXPECT_EQ(static_cast<std::size_t>(definedAtoms[0].type), *forceField.findPseudoAtom("Si"));
  EXPECT_EQ(static_cast<std::size_t>(definedAtoms[1].type), *forceField.findPseudoAtom("O"));
  EXPECT_NEAR(definedAtoms[1].position.x, 0.4, 1.0e-12);
  EXPECT_NEAR(definedAtoms[1].position.y, 0.5, 1.0e-12);
  EXPECT_NEAR(definedAtoms[1].position.z, 0.6, 1.0e-12);
  EXPECT_NEAR(definedAtoms[0].charge, 1.5, 1.0e-12);
  EXPECT_NEAR(definedAtoms[1].charge, -0.75, 1.0e-12);
}

TEST(UNINITIALIZED_MEMBERS, automatic_coulomb_cutoff_follows_the_box_with_explicit_ewald_parameters)
{
  ForceField forceField({{"O", false, 16.0, 0.0, 0.0, 8, false}}, {{80.0, 3.0}},
                        ForceField::MixingRule::Lorentz_Berthelot, 12.0, 11.0, 12.0, false, false, true);
  forceField.cutOffFrameworkVDWAutomatic = true;
  forceField.cutOffMoleculeVDWAutomatic = true;
  forceField.cutOffCoulombAutomatic = true;
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.3;
  forceField.numberOfWaveVectors = int3(7, 8, 9);

  forceField.initializeAutomaticCutOff(SimulationBox(20.0, 30.0, 40.0));

  EXPECT_NEAR(forceField.cutOffCoulomb, 10.0, 1.0e-12);
  EXPECT_EQ(forceField.EwaldAlpha, 0.3);
  EXPECT_EQ(forceField.numberOfWaveVectors.z, 9);

  // the molecule-molecule line prints the molecule-molecule cutoff
  forceField.cutOffMoleculeVDW = 7.25;
  EXPECT_NE(forceField.printCutOffAutoStatus().find("Molecule-Molecule VDW:   7.25000"), std::string::npos)
      << forceField.printCutOffAutoStatus();
}

TEST(UNINITIALIZED_MEMBERS, charge_equilibration_rejects_unmeasured_ionization_energies_and_anions)
{
  // the table has room for ten ionization energies per element but only the measured ones are filled in: Ba has
  // two, so Ba2+ needs the (unmeasured, zero) third and used to get a negative hardness and NaN charges; Ta has
  // one, Hf four and Au two; index 0 of the table is the blank entry
  const auto forceFieldFor = [](std::string name, std::size_t atomicNumber, std::int64_t oxidationState)
  {
    ForceField forceField({{name, true, 100.0, 0.0, 0.0, atomicNumber, false}, {"O", true, 16.0, 0.0, 0.0, 8, false}},
                          {{20.0, 3.0}, {50.0, 3.0}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0,
                          false, false, true);
    forceField.pseudoAtoms[0].oxidationState = oxidationState;
    return forceField;
  };
  std::vector<Atom> atoms{Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false),
                          Atom({2.0, 0.0, 0.0}, 0.0, 1.0, 0, 1, 0, false, false)};
  const SimulationBox box(20.0, 20.0, 20.0);
  for (const auto& [name, atomicNumber, oxidationState] :
       {std::tuple<std::string, std::size_t, std::int64_t>{"Ba", 56, 2},
        {"Ta", 73, 5},
        {"Hf", 72, 4},
        {"Au", 79, 3},
        {"", 0, 0}})
  {
    EXPECT_THROW(ChargeEquilibration::computeChargeEquilibration(forceFieldFor(name, atomicNumber, oxidationState), box,
                                                                 atoms, ChargeEquilibration::Type::PeriodicEwaldSum),
                 std::runtime_error)
        << name << " " << oxidationState;
  }

  // an anion: the table holds ionization energies only (the sign used to be dropped, reading the cation's data)
  EXPECT_THROW(ChargeEquilibration::computeChargeEquilibration(forceFieldFor("O", 8, -2), box, atoms,
                                                               ChargeEquilibration::Type::PeriodicEwaldSum),
               std::runtime_error);

  // measured data (Zn has eight ionization energies; Ba+ uses its second) gives finite charges that sum to zero
  for (const auto& [name, atomicNumber, oxidationState] :
       {std::tuple<std::string, std::size_t, std::int64_t>{"Zn", 30, 2}, {"Ba", 56, 1}})
  {
    std::vector<Atom> pair = atoms;
    ChargeEquilibration::computeChargeEquilibration(forceFieldFor(name, atomicNumber, oxidationState), box, pair,
                                                    ChargeEquilibration::Type::PeriodicEwaldSum);
    EXPECT_EQ(std::to_string(pair[0].charge).find("nan"), std::string::npos) << name;
    EXPECT_GT(pair[0].charge, 0.0) << name;
    EXPECT_NEAR(pair[0].charge + pair[1].charge, 0.0, 1e-9) << name;
  }
}

TEST(UNINITIALIZED_MEMBERS, charge_equilibration_rejects_elements_beyond_the_reference_table)
{
  ForceField forceField({{"U", true, 238.0, 0.0, 0.0, 92, false}, {"O", true, 16.0, 0.0, 0.0, 8, false}},
                        {{20.0, 3.0}, {50.0, 3.0}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, false,
                        false, true);
  std::vector<Atom> atoms{Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false),
                          Atom({2.0, 0.0, 0.0}, 0.0, 1.0, 0, 1, 0, false, false)};
  EXPECT_THROW(ChargeEquilibration::computeChargeEquilibration(forceField, SimulationBox(20.0, 20.0, 20.0), atoms,
                                                               ChargeEquilibration::Type::PeriodicEwaldSum),
               std::runtime_error);
}
