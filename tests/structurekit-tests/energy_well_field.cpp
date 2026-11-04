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
import energy_opencl_well_field;
import opencl;

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


namespace
{

// One atom and its probe in a cell smaller than the cutoff reaches, so that several images of the atom count at every
// point, near faces as much as anywhere, with the atom given at a fractional position of the caller's choosing.
Crystal oneAtomAt(const double3 &fractional, bool carryFractional)
{
  Crystal framework;
  framework.name = "moved";
  framework.unitCell = UnitCell(edge, edge, edge);
  framework.mass = 16.0;
  framework.atoms.push_back(CrystalAtom{fractional * edge, 0.0, 1});
  if (carryFractional) framework.fractionalPositions.push_back(fractional);
  return framework;
}

}  // namespace


// A lattice translation moves nothing. An atom given outside the home cell is the same crystal as one given inside it,
// so the field is the same, point for point; the images looked at have to be counted from where the atom really is.
TEST(energy_well_field, an_atom_given_outside_the_cell_makes_the_same_field)
{
  PairInteractions interactions = probeAndAtom();
  LinearProbe probe = LinearProbe::singleSite(interactions, "X").value();
  const double3 home(0.75, 0.3, 0.6);

  WellField expected =
      WellFieldCPU::compute(interactions, oneAtomAt(home, true), probe, gridSize, 1, 0.0, {}, 1000.0, ceiling);

  for (const double3 &shift : {double3(-1.0, 0.0, 0.0), double3(1.0, -1.0, 2.0), double3(-3.0, 3.0, -5.0)})
  {
    for (bool carryFractional : {true, false})
    {
      WellField moved = WellFieldCPU::compute(interactions, oneAtomAt(home + shift, carryFractional), probe, gridSize, 1,
                                              0.0, {}, 1000.0, ceiling);
      ASSERT_EQ(moved.numberOfVoxels(), expected.numberOfVoxels());
      std::size_t different = 0;
      float worst = 0.0f;
      for (std::size_t voxel = 0; voxel < moved.numberOfVoxels(); ++voxel)
      {
        const float difference = std::abs(moved.energy[voxel] - expected.energy[voxel]);
        worst = std::max(worst, difference);
        if (difference > 1.0e-3f * std::max(std::abs(expected.energy[voxel]), 1.0f)) ++different;
      }
      EXPECT_EQ(different, 0uz) << "shift " << shift.x << " " << shift.y << " " << shift.z << " carried " << carryFractional
                                << ": worst by " << worst;
    }
  }
}

// The device takes the same crystal to the same field, wherever in the lattice the atom was given.
TEST(energy_well_field, the_processor_and_the_device_agree_about_an_atom_given_outside_the_cell)
{
  OpenCL::initialize();
  if (!OpenCL::clContext.has_value()) GTEST_SKIP() << "no OpenCL device";

  PairInteractions interactions = probeAndAtom();
  LinearProbe probe = LinearProbe::singleSite(interactions, "X").value();
  const double3 home(0.75, 0.3, 0.6);

  for (const double3 &shift : {double3(0.0, 0.0, 0.0), double3(1.0, -1.0, 2.0), double3(-3.0, 3.0, -5.0)})
  {
    Crystal framework = oneAtomAt(home + shift, true);
    WellField cpu = WellFieldCPU::compute(interactions, framework, probe, gridSize, 1, 0.0, {}, 1000.0, ceiling);
    WellField gpu = WellFieldOpenCL::compute(interactions, framework, probe, gridSize, 1, 0.0, {}, 1000.0, ceiling);
    ASSERT_EQ(gpu.numberOfVoxels(), cpu.numberOfVoxels());
    std::size_t different = 0;
    for (std::size_t voxel = 0; voxel < cpu.numberOfVoxels(); ++voxel)
    {
      if (std::abs(gpu.energy[voxel] - cpu.energy[voxel]) > 1.0e-3f * std::max(std::abs(cpu.energy[voxel]), 1.0f)) ++different;
    }
    EXPECT_EQ(different, 0uz) << "shift " << shift.x << " " << shift.y << " " << shift.z;
  }
}
