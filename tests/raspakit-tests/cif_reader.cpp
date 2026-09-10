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
