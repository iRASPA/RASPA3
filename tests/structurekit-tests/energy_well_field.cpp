#include <gtest/gtest.h>

import std;

import uint3;
import double3;
import unit_cell;
import crystal;
import pair_interactions;
import blocking_spheres;
import energy_shared_linear_probe;
import energy_shared_well_field;
import energy_well_field;

// The processor well field for a single site. Each pair is held to the ceiling, and so must their sum be:
// a voxel on top of several atoms is at the ceiling, not several times it, which is what the molecular route
// and the device both answer.

namespace
{

constexpr double edge = 10.0;
constexpr double ceiling = 1.0e7;
const uint3 gridSize{5, 5, 5};

PairInteractions probeAndAtom()
{
  PairInteractions interactions;
  interactions.numberOfTypes = 2;
  interactions.names = {"X", "O"};
  interactions.charges = {0.0, 0.0};
  interactions.parameters = {PairParameters{3.0, 100.0, 0.0}, PairParameters{3.0, 100.0, 0.0},
                             PairParameters{3.0, 100.0, 0.0}, PairParameters{3.0, 100.0, 0.0}};
  interactions.cutOffVDW = 4.0;
  interactions.cutOffCoulomb = 4.0;
  return interactions;
}

// Three atoms a tenth of an Ångström from the origin, which is the first grid point.
Crystal threeAtomsOnTheOrigin()
{
  Crystal framework;
  framework.name = "crowded";
  framework.unitCell = UnitCell(edge, edge, edge);
  framework.mass = 48.0;
  for (double3 position : {double3(0.1, 0.0, 0.0), double3(0.0, 0.1, 0.0), double3(0.0, 0.0, 0.1)})
  {
    framework.atoms.push_back(CrystalAtom{position, 0.0, 1});
    framework.fractionalPositions.push_back(position / edge);
  }
  return framework;
}

}  // namespace


TEST(energy_well_field, a_voxel_on_several_atoms_is_at_the_ceiling_and_not_above_it)
{
  PairInteractions interactions = probeAndAtom();
  Crystal framework = threeAtomsOnTheOrigin();
  LinearProbe probe = LinearProbe::singleSite(interactions, "X").value();

  WellField field = WellFieldCPU::compute(interactions, framework, probe, gridSize, 1, 0.0, {}, 1000.0, ceiling);
  ASSERT_EQ(field.numberOfVoxels(), 125u);

  EXPECT_FLOAT_EQ(field.energy[0], static_cast<float>(ceiling));
  for (std::size_t voxel = 0; voxel < field.numberOfVoxels(); ++voxel)
  {
    EXPECT_LE(field.energy[voxel], static_cast<float>(ceiling)) << "voxel " << voxel;
  }
}
