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
import interactions_ewald;

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

TEST(UNINITIALIZED_MEMBERS, periodic_charge_equilibration_gives_a_small_cell_the_charges_of_its_supercell)
{
  // CsCl in a 4 A cubic cell. Caesium's hardness is small (3.4 eV), so its orbital-overlap correction reaches
  // beyond the cell: with the sites' own images and the pairs' further images left out of the sum (only the
  // minimum image was kept), the 1x1x1 cell got other charges than its 2x2x2 supercell, in which those images are
  // atoms of their own. The same crystal must give the same charges whatever cell describes it.
  const ForceField forceField({{"Cs", true, 132.9, 0.0, 0.0, 55, false}, {"Cl", true, 35.45, 0.0, 0.0, 17, false}},
                              {{20.0, 3.0}, {50.0, 3.0}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0,
                              false, false, true);
  const double a = 4.0;
  std::vector<Atom> cell{Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false),
                         Atom({0.5 * a, 0.5 * a, 0.5 * a}, 0.0, 1.0, 0, 1, 0, false, false)};
  ChargeEquilibration::computeChargeEquilibration(forceField, SimulationBox(a, a, a), cell,
                                                  ChargeEquilibration::Type::PeriodicEwaldSum);
  EXPECT_GT(cell[0].charge, 0.0);
  EXPECT_NEAR(cell[0].charge + cell[1].charge, 0.0, 1e-9);

  std::vector<Atom> supercell;
  for (int u = 0; u < 2; ++u)
    for (int v = 0; v < 2; ++v)
      for (int w = 0; w < 2; ++w)
      {
        const double3 shift(a * u, a * v, a * w);
        supercell.emplace_back(double3(0.0, 0.0, 0.0) + shift, 0.0, 1.0, 0, 0, 0, false, false);
        supercell.emplace_back(double3(0.5 * a, 0.5 * a, 0.5 * a) + shift, 0.0, 1.0, 0, 1, 0, false, false);
      }
  ChargeEquilibration::computeChargeEquilibration(forceField, SimulationBox(2.0 * a, 2.0 * a, 2.0 * a), supercell,
                                                  ChargeEquilibration::Type::PeriodicEwaldSum);
  // the two agree to the Ewald sum's convergence (6e-6 here); the minimum-image version was 0.18 e apart
  for (const Atom& atom : supercell)
  {
    EXPECT_NEAR(atom.charge, cell[static_cast<std::size_t>(atom.type)].charge, 1e-4) << atom.type;
  }
}

TEST(UNINITIALIZED_MEMBERS, periodic_charge_equilibration_is_exact_for_a_skewed_triclinic_cell)
{
  // a skewed triclinic cell (75, 110 and 80 degrees, perpendicular widths of about 5 A). The Ewald matrix picks
  // the real-space image by rounding fractional coordinates, which is the nearest image for every separation
  // below half the smallest perpendicular width (its fractional components are then all below one half), and
  // sums the wave vectors inside the physical cutoff sphere. Either going wrong gives the 1x1x1 cell another
  // potential matrix, and other charges, than its 2x2x2 supercell.
  const double degree = std::numbers::pi / 180.0;
  const SimulationBox box(5.5, 6.0, 6.5, 75.0 * degree, 110.0 * degree, 80.0 * degree);
  const SimulationBox superBox(11.0, 12.0, 13.0, 75.0 * degree, 110.0 * degree, 80.0 * degree);
  const std::array<double3, 2> fractional{double3(0.0, 0.0, 0.0), double3(0.25, 0.5, 0.625)};  // 4 A apart

  std::vector<Atom> cell;
  for (std::uint16_t type = 0; type < 2; ++type)
    cell.emplace_back(box.cell * fractional[type], 0.0, 1.0, 0, type, 0, false, false);
  std::vector<Atom> supercell;
  for (int u = 0; u < 2; ++u)
    for (int v = 0; v < 2; ++v)
      for (int w = 0; w < 2; ++w)
        for (std::uint16_t type = 0; type < 2; ++type)
          supercell.emplace_back(box.cell * (fractional[type] + double3(u, v, w)), 0.0, 1.0, 0, type, 0, false, false);

  // the Coulomb potential matrix: the potential at atom 0 from all images of atom 0 and of atom 1 is, in the
  // supercell, the sum over the eight copies of each
  std::vector<std::complex<double>> ex, ey, ez, exy;
  std::vector<double> V(4);
  Interactions::computeEwaldFourierChargeEquilibrationPotentialMatrix(ex, ey, ez, exy, box, cell, V);
  std::vector<double> W(16 * 16);
  Interactions::computeEwaldFourierChargeEquilibrationPotentialMatrix(ex, ey, ez, exy, superBox, supercell, W);
  double W00 = 0.0, W01 = 0.0;
  for (std::size_t j = 0; j < 16; ++j) (j % 2 == 0 ? W00 : W01) += W[j];
  // to 1e-6, the truncation of the Ewald sums at a tolerance of 1e-8 per term, now that the wave vectors along
  // the skewed directions are all summed (bounded by the perpendicular widths and an index sphere, they left the
  // two cells 6e-6 apart)
  EXPECT_NEAR(V[0], W00, 2e-6);
  EXPECT_NEAR(V[1], W01, 2e-6);

  // and the charges
  const ForceField forceField({{"Cs", true, 132.9, 0.0, 0.0, 55, false}, {"Cl", true, 35.45, 0.0, 0.0, 17, false}},
                              {{20.0, 3.0}, {50.0, 3.0}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0,
                              false, false, true);
  ChargeEquilibration::computeChargeEquilibration(forceField, box, cell, ChargeEquilibration::Type::PeriodicEwaldSum);
  EXPECT_GT(cell[0].charge, 0.0);
  EXPECT_NEAR(cell[0].charge + cell[1].charge, 0.0, 1e-9);
  ChargeEquilibration::computeChargeEquilibration(forceField, superBox, supercell,
                                                  ChargeEquilibration::Type::PeriodicEwaldSum);
  for (const Atom& atom : supercell)
  {
    EXPECT_NEAR(atom.charge, cell[static_cast<std::size_t>(atom.type)].charge, 1e-4) << atom.type;
  }
}

TEST(UNINITIALIZED_MEMBERS, metal_organic_framework_force_field_starts_metals_from_eqeq_charge_centres)
{
  // Wilmer and Snurr's chargecenters.dat: metals start the charge equilibration from their usual oxidation state
  // (all elements used to start from zero)
  const ForceField forceField = ForceField::makeMetalOrganicFrameworkForceField(12.0, true, false, false);
  const auto oxidationState = [&](const std::string& name)
  { return forceField.pseudoAtoms[*forceField.findPseudoAtom(name)].oxidationState; };
  EXPECT_EQ(oxidationState("Zn"), 2);
  EXPECT_EQ(oxidationState("Cu"), 2);
  EXPECT_EQ(oxidationState("Zr"), 4);
  EXPECT_EQ(oxidationState("Al"), 3);
  EXPECT_EQ(oxidationState("Mg"), 2);
  EXPECT_EQ(oxidationState("Fe"), 2);
  EXPECT_EQ(oxidationState("Ba"), 1);  // Ba has two measured ionization energies, so EQeq stops at 1
  EXPECT_EQ(oxidationState("H"), 0);
  EXPECT_EQ(oxidationState("C"), 0);
  EXPECT_EQ(oxidationState("O"), 0);
  EXPECT_EQ(oxidationState("Ta"), 0);
  EXPECT_EQ(ForceField::defaultOxidationState("Xx"), 0);

  // every default lies within the measured ionization energies, so the charge equilibration accepts it
  std::vector<Atom> atoms{Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false),
                          Atom({2.0, 0.0, 0.0}, 0.0, 1.0, 0, *forceField.findPseudoAtom("O"), 0, false, false)};
  for (const std::string& metal : {"Zn", "Cu", "Zr", "Al", "Mg", "Fe", "Ba", "Hf", "Pt", "Au"})
  {
    atoms[0].type = static_cast<std::uint16_t>(*forceField.findPseudoAtom(metal));
    EXPECT_NO_THROW(ChargeEquilibration::computeChargeEquilibration(forceField, SimulationBox(10.0, 10.0, 10.0), atoms,
                                                                    ChargeEquilibration::Type::PeriodicEwaldSum))
        << metal;
  }
}

TEST(UNINITIALIZED_MEMBERS, cif_site_label_names_the_pseudo_atom_before_the_element_does)
{
  // a force field with a generic Zn and a Zn1 with an oxidation state, as RASPA2 force fields have: a site labelled
  // Zn1 with element Zn is the Zn1 pseudo-atom (the element used to be looked up first and took the generic Zn)
  ForceField forceField({{"Zn", true, 65.38, 0.0, 0.0, 30, false},
                         {"Zn1", true, 65.38, 0.0, 0.0, 30, false},
                         {"O", true, 15.9994, 0.0, 0.0, 8, false}},
                        {{20.0, 3.0}, {20.0, 3.0}, {50.0, 3.0}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0,
                        12.0, false, false, true);
  forceField.pseudoAtoms[1].oxidationState = 2;
  constexpr std::string_view cif = R"(
data_test
_cell_length_a 10.0
_cell_length_b 10.0
_cell_length_c 10.0
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 90.0
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Zn1 Zn 0.1 0.2 0.3
Zn2 Zn 0.6 0.2 0.3
O1  O  0.4 0.5 0.6
O2  O2- 0.4 0.5 0.1
)";
  const auto parsed = CIFReader::readCIFString(std::string(cif), forceField, CIFReader::UseChargesFrom::PseudoAtoms);
  ASSERT_TRUE(parsed.has_value());
  const auto& definedAtoms = std::get<2>(parsed.value());
  ASSERT_EQ(definedAtoms.size(), 4uz);
  EXPECT_EQ(static_cast<std::size_t>(definedAtoms[0].type), *forceField.findPseudoAtom("Zn1"));
  // Zn2 is no pseudo-atom: the element is used, not Zn1 with its digits stripped to Zn either way
  EXPECT_EQ(static_cast<std::size_t>(definedAtoms[1].type), *forceField.findPseudoAtom("Zn"));
  // a label that is only a name falls back to the element, also when the element carries a charge suffix
  EXPECT_EQ(static_cast<std::size_t>(definedAtoms[2].type), *forceField.findPseudoAtom("O"));
  EXPECT_EQ(static_cast<std::size_t>(definedAtoms[3].type), *forceField.findPseudoAtom("O"));
}

TEST(UNINITIALIZED_MEMBERS, charge_equilibration_of_no_atoms_does_nothing)
{
  // an empty framework used to write the total-charge row past an empty vector and loop from 1 towards SIZE_MAX
  const ForceField forceField = ForceField::makeMetalOrganicFrameworkForceField(12.0, true, false, false);
  std::vector<Atom> none;
  for (ChargeEquilibration::Type type :
       {ChargeEquilibration::Type::PeriodicEwaldSum, ChargeEquilibration::Type::NonPeriodic,
        ChargeEquilibration::Type::PeriodicDirectSum})
  {
    EXPECT_NO_THROW(
        ChargeEquilibration::computeChargeEquilibration(forceField, SimulationBox(10.0, 10.0, 10.0), none, type));
    EXPECT_TRUE(none.empty());
  }
}

TEST(UNINITIALIZED_MEMBERS, cif_charge_equilibration_gives_a_redundant_site_the_charge_of_the_site_it_duplicates)
{
  // P -1 with the oxygen listed twice, as O1 at x and O2 at -x: the orbits coincide, the expansion keeps the later
  // atom's copies, and the earlier atom averaged its charge over no atoms (0/0, NaN; -ffast-math hides a NaN from
  // std::isnan, hence the text check)
  const ForceField forceField = ForceField::makeMetalOrganicFrameworkForceField(12.0, true, false, false);
  constexpr std::string_view cif = R"(
data_test
_cell_length_a 10.0
_cell_length_b 10.0
_cell_length_c 10.0
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 90.0
_symmetry_space_group_name_H-M 'P -1'
_symmetry_Int_Tables_number 2
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
O1 O 0.1 0.2 0.3
O2 O 0.9 0.8 0.7
Zn1 Zn 0.5 0.5 0.0
)";
  const auto parsed =
      CIFReader::readCIFString(std::string(cif), forceField, CIFReader::UseChargesFrom::ChargeEquilibration);
  ASSERT_TRUE(parsed.has_value());
  const auto& [simulationBox, spaceGroupHallNumber, definedAtoms, unitCellAtoms] = parsed.value();
  ASSERT_EQ(definedAtoms.size(), 3uz);
  ASSERT_EQ(unitCellAtoms.size(), 3uz);  // two oxygens and the zinc on the inversion centre
  for (const Atom& atom : definedAtoms)
  {
    EXPECT_EQ(std::to_string(atom.charge).find("nan"), std::string::npos) << atom.moleculeId;
  }
  EXPECT_LT(definedAtoms[0].charge, 0.0);
  EXPECT_NEAR(definedAtoms[0].charge, definedAtoms[1].charge, 1e-12);
  EXPECT_GT(definedAtoms[2].charge, 0.0);
  double total = 0.0;
  for (const Atom& atom : unitCellAtoms)
  {
    EXPECT_EQ(std::to_string(atom.charge).find("nan"), std::string::npos);
    total += atom.charge;
  }
  EXPECT_NEAR(total, 0.0, 1e-9);
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

TEST(UNINITIALIZED_MEMBERS, direct_sum_charge_equilibration_does_not_change_when_the_structure_is_translated)
{
  // the direct sum adds the images of a pair in a window of cells centred on the pair's separation. Centred on the
  // raw separation, the window moved by a whole cell when the structure was translated so that the pair's raw
  // separation changed from -0.4a to +0.6a (it is the same minimum image, 0.4a); the truncated lattice sum is not
  // shift-invariant, so the same crystal got other charges after a translation. The window is now centred on the
  // minimum image. The matrix is also filled from its upper triangle, which the symmetry of the result checks.
  const ForceField forceField({{"Cs", true, 132.9, 0.0, 0.0, 55, false}, {"Cl", true, 35.45, 0.0, 0.0, 17, false}},
                              {{20.0, 3.0}, {50.0, 3.0}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0,
                              false, false, true);
  const double a = 4.0;
  const SimulationBox box(a, a, a);

  std::vector<Atom> cell{Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false),
                         Atom({0.4 * a, 0.4 * a, 0.4 * a}, 0.0, 1.0, 0, 1, 0, false, false)};
  ChargeEquilibration::computeChargeEquilibration(forceField, box, cell, ChargeEquilibration::Type::PeriodicDirectSum);
  EXPECT_GT(cell[0].charge, 0.0);
  EXPECT_NEAR(cell[0].charge + cell[1].charge, 0.0, 1e-9);

  // the same crystal translated by 0.7a along every axis and wrapped into the cell
  std::vector<Atom> translated{Atom({0.7 * a, 0.7 * a, 0.7 * a}, 0.0, 1.0, 0, 0, 0, false, false),
                               Atom({0.1 * a, 0.1 * a, 0.1 * a}, 0.0, 1.0, 0, 1, 0, false, false)};
  ChargeEquilibration::computeChargeEquilibration(forceField, box, translated,
                                                  ChargeEquilibration::Type::PeriodicDirectSum);
  EXPECT_NEAR(translated[0].charge, cell[0].charge, 1e-9);
  EXPECT_NEAR(translated[1].charge, cell[1].charge, 1e-9);

  // and with the two atoms listed the other way round
  std::vector<Atom> reordered{Atom({0.1 * a, 0.1 * a, 0.1 * a}, 0.0, 1.0, 0, 1, 0, false, false),
                              Atom({0.7 * a, 0.7 * a, 0.7 * a}, 0.0, 1.0, 0, 0, 0, false, false)};
  ChargeEquilibration::computeChargeEquilibration(forceField, box, reordered,
                                                  ChargeEquilibration::Type::PeriodicDirectSum);
  EXPECT_NEAR(reordered[0].charge, cell[1].charge, 1e-9);
  EXPECT_NEAR(reordered[1].charge, cell[0].charge, 1e-9);
}
