#include <gtest/gtest.h>

import <vector>;
import <tuple>;
import <algorithm>;

import double3;
import double3x3;
import simulationbox;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;

TEST(pdb, triclinic)
{
  double3x3 m{ double3{25.0, 2.5e-6, 0.0}, double3{0.0, 25.0, 0.0}, double3{0.0, 0.0, 25.0} };
  SimulationBox rect_box = SimulationBox(m, SimulationBox::Type::Rectangular);
  SimulationBox triclinic_box = SimulationBox(m, SimulationBox::Type::Triclinic);
  double3 dr = double3{6.8015008307321985, -8.0937084818712712, -13.749755827850278};
  double3 res = rect_box.applyPeriodicBoundaryConditions(dr);
  double3 s = triclinic_box.inverseUnitCell * dr;
  double3 r = triclinic_box.unitCell * s;
  
  double3 test = triclinic_box.applyPeriodicBoundaryConditions(dr);
}
