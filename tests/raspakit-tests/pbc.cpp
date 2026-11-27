#include <gtest/gtest.h>

import std;

import double3;
import double3x3;
import simulationbox;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;

TEST(pbc, triclinic)
{
  // double3x3 m{ double3{25.0, 2.5e-6, 0.0}, double3{0.0, 25.0, 0.0}, double3{0.0, 0.0, 25.0} };
  // SimulationBox rect_box = SimulationBox(m, SimulationBox::Type::Rectangular);
  // SimulationBox triclinic_box = SimulationBox(m, SimulationBox::Type::Triclinic);
  // double3 dr = double3{6.8015008307321985, -8.0937084818712712, -13.749755827850278};
  // double3 res = rect_box.applyPeriodicBoundaryConditions(dr);
  // double3 s = triclinic_box.inverseCell * dr;
  // double3 r = triclinic_box.cell * s;
  //
  // double3 test = triclinic_box.applyPeriodicBoundaryConditions(dr);
}

// An angle 1e-7 rad off a right angle has a tilt in the cell matrix (3e-6 A for 30 A), which the minimum image of a
// rectangular box leaves out: such a box is triclinic. It was taken for rectangular within 1e-5 rad; a right angle
// given in degrees still makes a rectangular box.
TEST(pbc, a_cell_a_little_off_right_angles_is_triclinic)
{
  const double rightAngle = 0.5 * std::numbers::pi;
  const SimulationBox tilted(30.0, 30.0, 30.0, rightAngle + 1.0e-7, rightAngle, rightAngle);
  EXPECT_EQ(tilted.type, SimulationBox::Type::Triclinic);
  EXPECT_NE(tilted.cell.cy, 0.0);
  const double degrees = 90.0 * std::numbers::pi / 180.0;
  EXPECT_EQ(SimulationBox(30.0, 30.0, 30.0, degrees, degrees, degrees).type, SimulationBox::Type::Rectangular);
}
