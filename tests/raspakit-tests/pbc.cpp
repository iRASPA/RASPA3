#ifdef USE_LEGACY_HEADERS
#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <tuple>
#include <vector>
#endif

#ifdef USE_STD_IMPORT
#include <gtest/gtest.h>
import std;
#endif

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
