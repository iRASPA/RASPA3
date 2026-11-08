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
