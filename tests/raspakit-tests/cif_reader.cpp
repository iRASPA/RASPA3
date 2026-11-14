#include <gtest/gtest.h>

import std;

import atom;
import cif_reader;
import forcefield;
import framework;

namespace
{

// IZA-style: tetrahedral sites labelled T1, T2 with _atom_site_type_symbol Si.
constexpr std::string_view kIzaTSitesCif = R"(
data_test
_cell_length_a 10.0
_cell_length_b 10.0
_cell_length_c 10.0
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 90.0
_symmetry_space_group_name_H-M 'P 1'
_symmetry_Int_Tables_number 1
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
T1 Si 0.1 0.2 0.3
O1 O  0.4 0.5 0.6
)";

}  // namespace

TEST(cif_reader, iza_tetrahedral_site_labels_map_to_the_type_symbol)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, false);
  const auto cif = CIFReader::readCIFString(std::string(kIzaTSitesCif), forceField,
                                              CIFReader::UseChargesFrom::CIF_File);
  ASSERT_TRUE(cif.has_value());

  auto [simulationBox, spaceGroupHallNumber, definedAtoms, unitCellAtoms] = cif.value();
  ASSERT_EQ(definedAtoms.size(), 2);

  std::optional<std::size_t> silicon = forceField.findPseudoAtom("Si");
  std::optional<std::size_t> oxygen = forceField.findPseudoAtom("O");
  ASSERT_TRUE(silicon.has_value());
  ASSERT_TRUE(oxygen.has_value());

  EXPECT_EQ(static_cast<std::size_t>(definedAtoms[0].type), silicon.value());
  EXPECT_EQ(static_cast<std::size_t>(definedAtoms[1].type), oxygen.value());

  Framework framework(forceField, "test", simulationBox, spaceGroupHallNumber, definedAtoms, unitCellAtoms,
                      {1, 1, 1});
  EXPECT_EQ(framework.unitCellAtoms.size(), 2);
}

TEST(cif_reader, site_label_is_used_when_the_type_symbol_is_not_a_pseudo_atom)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, false);
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
Si1 Xx 0.0 0.0 0.0
)";

  const auto parsed = CIFReader::readCIFString(std::string(cif), forceField, CIFReader::UseChargesFrom::CIF_File);
  ASSERT_TRUE(parsed.has_value());
  auto [simulationBox, spaceGroupHallNumber, definedAtoms, unitCellAtoms] = parsed.value();
  ASSERT_EQ(definedAtoms.size(), 1);
  EXPECT_EQ(static_cast<std::size_t>(definedAtoms[0].type), forceField.findPseudoAtom("Si").value());
}

TEST(cif_reader, unknown_type_symbol_is_invalid_force_field)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, false);
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
Zz1 Zz 0.0 0.0 0.0
)";

  const auto parsed = CIFReader::readCIFString(std::string(cif), forceField, CIFReader::UseChargesFrom::CIF_File);
  ASSERT_FALSE(parsed.has_value());
  EXPECT_EQ(parsed.error(), CIFReader::ParseError::invalidForceField);
}

TEST(cif_reader, type_symbol_oxidation_states_in_either_notation_and_case_strip_to_the_element)
{
  // a force field with the elements and one pseudo-atom named after a full type symbol
  const ForceField forceField({{"Na", true, 22.99, 1.0, 0.0, 11, false},
                               {"Cl", true, 35.45, -1.0, 0.0, 17, false},
                               {"Zn", true, 65.38, 2.0, 0.0, 30, false},
                               {"Fe", true, 55.85, 3.0, 0.0, 26, false},
                               {"O", true, 15.999, -2.0, 0.0, 8, false},
                               {"O2-", true, 15.999, -2.0, 0.0, 8, false}},
                              {{10.0, 3.0}, {10.0, 3.0}, {10.0, 3.0}, {10.0, 3.0}, {10.0, 3.0}, {10.0, 3.0}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, true, false);
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
Na1 Na+  0.0 0.0 0.0
Cl1 Cl-  0.1 0.0 0.0
Zn1 Zn+2 0.2 0.0 0.0
Fe1 Fe3+ 0.3 0.0 0.0
X1  O-2  0.4 0.0 0.0
X2  O2-  0.5 0.0 0.0
X3  ZN   0.6 0.0 0.0
X4  cl   0.7 0.0 0.0
X5  O2   0.8 0.0 0.0
)";

  const auto parsed = CIFReader::readCIFString(std::string(cif), forceField, CIFReader::UseChargesFrom::PseudoAtoms);
  ASSERT_TRUE(parsed.has_value());
  auto [simulationBox, spaceGroupHallNumber, definedAtoms, unitCellAtoms] = parsed.value();
  ASSERT_EQ(definedAtoms.size(), 9);

  const std::vector<std::string> expected{"Na", "Cl", "Zn", "Fe", "O", "O2-", "Zn", "Cl", "O"};
  for (std::size_t i = 0; i < expected.size(); ++i)
  {
    EXPECT_EQ(static_cast<std::size_t>(definedAtoms[i].type), forceField.findPseudoAtom(expected[i]).value())
        << "site " << i << " should be " << expected[i];
  }
}

// A file that names its space group 'unknown', lists its operations under the current IUCr tag with a comment among
// the loop tags, quotes operations that contain spaces, and quotes an atom label. The name 'unknown' matched the
// placeholder entry of the space-group table, which has no operations and expanded every atom to none; the comment
// ended the loop header; and the words of a quoted value were read as values of their own.
TEST(cif_reader, unknown_space_group_quoted_operations_and_modern_tags)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, false);
  constexpr std::string_view cif = R"(
data_test
_cell_length_a 10.0
_cell_length_b 10.0
_cell_length_c 10.0
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 90.0
_symmetry_space_group_name_H-M 'unknown'
loop_
_space_group_symop_id
# the operations of P -1, quoted, with spaces
_space_group_symop_operation_xyz
1 'x, y, z'
2 '-x, -y, -z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
'Si 1' Si 0.1 0.2 0.3
O1 O  0.4 0.5 0.6
)";
  const auto read = CIFReader::readCIFString(std::string(cif), forceField, CIFReader::UseChargesFrom::CIF_File);
  ASSERT_TRUE(read.has_value());
  auto [simulationBox, spaceGroupHallNumber, definedAtoms, unitCellAtoms] = read.value();
  ASSERT_EQ(definedAtoms.size(), 2uz);
  EXPECT_NEAR(definedAtoms[0].position.x, 0.1, 1e-12);
  EXPECT_NEAR(definedAtoms[1].position.z, 0.6, 1e-12);
  // the inversion doubles both atoms
  EXPECT_EQ(unitCellAtoms.size(), 4uz);
}

// The space group by its Hall symbol under the IUCr core names, '_space_group_name_Hall' (CIF1) and
// '_space_group.name_Hall' (CIF2): only the deprecated '_symmetry_space_group_name_Hall' was read, and a file with the
// standard name fell back to P 1 (no symmetry operations applied).
TEST(cif_reader, hall_symbol_under_the_iucr_names)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, false);
  std::vector<std::size_t> numbers;
  for (const std::string tag : {"_symmetry_space_group_name_Hall", "_space_group_name_Hall", "_space_group.name_Hall"})
  {
    const std::string cif = std::format(R"(
data_test
_cell_length_a 10.0
_cell_length_b 10.0
_cell_length_c 10.0
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 90.0
{} '-P 4 2 3'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si1 Si 0.1 0.2 0.3
)",
                                        tag);
    const auto parsed = CIFReader::readCIFString(cif, forceField, CIFReader::UseChargesFrom::CIF_File);
    ASSERT_TRUE(parsed.has_value()) << tag;
    numbers.push_back(std::get<1>(parsed.value()));
  }
  EXPECT_NE(numbers[0], 1uz);
  EXPECT_EQ(numbers[1], numbers[0]);
  EXPECT_EQ(numbers[2], numbers[0]);
}

// A file with site labels but no _atom_site_type_symbol (common in CSD, ICSD and RASPA2 files) reads its atoms from
// the labels; every row used to be skipped, leaving a framework without atoms. A data-block keyword after the loop
// ends it instead of being read as the next row.
TEST(cif_reader, site_labels_without_type_symbols_are_read)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, false);
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
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si1 0.1 0.2 0.3
O1  0.4 0.5 0.6
data_second
)";

  const auto parsed = CIFReader::readCIFString(std::string(cif), forceField, CIFReader::UseChargesFrom::CIF_File);
  ASSERT_TRUE(parsed.has_value());
  auto [simulationBox, spaceGroupHallNumber, definedAtoms, unitCellAtoms] = parsed.value();
  ASSERT_EQ(definedAtoms.size(), 2);
  EXPECT_EQ(static_cast<std::size_t>(definedAtoms[0].type), forceField.findPseudoAtom("Si").value());
  EXPECT_EQ(static_cast<std::size_t>(definedAtoms[1].type), forceField.findPseudoAtom("O").value());
}
