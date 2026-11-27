#include <gtest/gtest.h>

import std;

import units;

TEST(units, Test_units_link_from_foundationkit_alone)
{
  // the implementation of the module lived in raspakit, so a target that links only foundationkit could not
  // resolve the non-constexpr members
  EXPECT_DOUBLE_EQ(Units::KB, Units::BoltzmannConstant / Units::EnergyConversionFactor);
  EXPECT_DOUBLE_EQ(Units::KelvinToEnergy * Units::EnergyToKelvin, 1.0);
  Units::setUnits(Units::System::RASPA);
  EXPECT_EQ(Units::unitOfLengthString, "m");
}

// In reduced units every energy is in units of epsilon: a force constant given in kcal/mol (MM3, CVFF, AMBER and
// LAMMPS bonded terms), kJ/mol or eV is taken as it is, like one in Kelvin. The factors used to be zero, which
// zeroed those force constants and printed them as 0 / 0.
TEST(units, reduced_units_take_every_energy_unit_as_the_reduced_energy)
{
  Units::setUnits(Units::System::ReducedUnits);
  EXPECT_EQ(Units::KelvinToEnergy, 1.0);
  EXPECT_EQ(Units::KCalPerMolToEnergy, 1.0);
  EXPECT_EQ(Units::EnergyToKCalPerMol, 1.0);
  EXPECT_EQ(Units::EnergyToKJPerMol, 1.0);
  EXPECT_EQ(Units::EnergyToEV, 1.0);
  Units::setUnits(Units::System::RASPA);
  EXPECT_DOUBLE_EQ(Units::KCalPerMolToEnergy * Units::EnergyToKCalPerMol, 1.0);
}
