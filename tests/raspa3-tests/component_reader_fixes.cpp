#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import forcefield;
import component;
import bond_bend_potential;
import mc_moves_probabilities;

// The molecule-file reader: the indices it is handed have to lie inside the molecule, and the three-body
// bond-bend takes three of them.

namespace
{

ForceField makeAlkaneForceField()
{
  return ForceField({{"CH3", false, 15.04, 0.0, 0.0, 6, false}, {"CH2", false, 14.03, 0.0, 0.0, 6, false}},
                    {{98.0, 3.75}, {46.0, 3.95}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true,
                    true, false);
}

// propane with a flexible bond, with 'extra' spliced into the top-level object
std::string propaneJson(std::string_view extra)
{
  return std::format(R"({{
  "CriticalTemperature" : 369.825,
  "CriticalPressure" : 4247660.0,
  "AcentricFactor" : 0.1524,
  "pseudoAtoms" : [["CH3", [0.0, 0.0, 0.0]], ["CH2", [1.54, 0.0, 0.0]], ["CH3", [2.0, 1.4, 0.0]]],
  "Connectivity" : [[0, 1], [1, 2]],
  "Bonds" : [[["CH3", "CH2"], "HARMONIC", [96500.0, 1.54]]],
  "Bends" : [[["CH3", "CH2", "CH3"], "HARMONIC", [62500.0, 114]]]{}
}})",
                     extra);
}

Component readPropane(const ForceField& forceField, std::string_view extra)
{
  TemporaryFile file("propane.json", propaneJson(extra));
  return Component(Component::Type::Adsorbate, 0, forceField, "propane", file.stemPath().string(), 5, 21,
                   MCMoveProbabilities(), std::nullopt, false);
}

}  // namespace

// A bond-bend (stretch-bend) is a three-body term: three identifiers, and a molecule of three atoms can have one.
TEST(COMPONENT_READER_FIXES, bond_bend_takes_three_identifiers)
{
  const ForceField forceField = makeAlkaneForceField();
  const Component propane =
      readPropane(forceField, R"(, "BondBends" : [[[0, 1, 2], "CVFF", [114.0, 800.0, 1.54, 800.0, 1.54]]])");
  ASSERT_EQ(propane.intraMolecularPotentials.bondBends.size(), 1uz);
  EXPECT_EQ(propane.intraMolecularPotentials.bondBends[0].identifiers, (std::array<std::size_t, 3>{0, 1, 2}));

  EXPECT_THROW(static_cast<void>(readPropane(
                   forceField, R"(, "BondBends" : [[[0, 1, 2, 0], "CVFF", [114.0, 800.0, 1.54, 800.0, 1.54]]])")),
               std::runtime_error);
}

// Explicit internal potentials and the CBMC starting bead index atoms of the molecule, counted from 0; an index past
// the last atom read past the molecule's atoms.
TEST(COMPONENT_READER_FIXES, indices_outside_the_molecule_are_refused)
{
  const ForceField forceField = makeAlkaneForceField();
  EXPECT_THROW(static_cast<void>(readPropane(
                   forceField, R"(, "BondBends" : [[[1, 2, 3], "CVFF", [114.0, 800.0, 1.54, 800.0, 1.54]]])")),
               std::runtime_error);
  EXPECT_THROW(static_cast<void>(readPropane(forceField, R"(, "StartingBead" : 3)")), std::runtime_error);
  EXPECT_THROW(static_cast<void>(readPropane(forceField, R"(, "StartingBead" : -1)")), std::runtime_error);
  EXPECT_EQ(readPropane(forceField, R"(, "StartingBead" : 2)").startingBead, 2uz);
}
