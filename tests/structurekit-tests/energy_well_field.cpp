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
import energy_shared_blocking_mask;
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

namespace
{

// A cell with two of its axes at sixty degrees, an atom well away from everything, and one blocking sphere at the
// origin wide enough to reach past the middle of the cell's short diagonal. The point at fractions (0.45, 0.45)
// is 7.8 from the copy of the sphere nearest it in fractions and 5.1 from the one a cell back along a, which is
// inside the sphere.
const UnitCell skewedCell(10.0, 10.0, 10.0, 0.5 * std::numbers::pi, 0.5 * std::numbers::pi, std::numbers::pi / 3.0);
const uint3 skewedGrid{20, 20, 4};
const BlockingSphere wideSphere{double3(0.0, 0.0, 0.0), 5.5};

Crystal oneAtomInASkewedCell()
{
  Crystal framework;
  framework.name = "skewed";
  framework.unitCell = skewedCell;
  framework.mass = 16.0;
  // off the grid, so that no voxel is at the cutoff from it to the last digit and counted by one and not the other
  const double3 fractional(0.52, 0.47, 0.55);
  framework.atoms.push_back(CrystalAtom{skewedCell.cell * fractional, 0.0, 1});
  framework.fractionalPositions.push_back(fractional);
  return framework;
}

// The distance to the sphere's surface by every copy of it for two cells around, one at a time.
double distanceBySearch(const double3 &fractional)
{
  double nearest = 1.0e10;
  for (int i = -2; i <= 2; ++i)
    for (int j = -2; j <= 2; ++j)
      for (int k = -2; k <= 2; ++k)
      {
        const double3 dr = skewedCell.cell * (fractional - wideSphere.centerFractional + double3(i, j, k));
        nearest = std::min(nearest, dr.length() - wideSphere.radius);
      }
  return nearest;
}

}  // namespace

// How far a point is from a blocking sphere is by the nearest copy of the sphere, and in a skewed cell that is not
// always the copy nearest in fractions of the cell. Taking that one left points inside a wide sphere outside it.
TEST(energy_well_field, a_point_is_inside_a_blocking_sphere_by_whichever_copy_of_it_is_nearest)
{
  const std::vector<BlockingSphere> spheres{wideSphere};
  EXPECT_LT(blockingSphereDistance(double3(0.45, 0.45, 0.0), skewedCell, spheres), -0.4);

  std::size_t inside = 0;
  for (std::size_t iy = 0; iy < skewedGrid.y; ++iy)
  {
    for (std::size_t ix = 0; ix < skewedGrid.x; ++ix)
    {
      const double3 fractional(static_cast<double>(ix) / static_cast<double>(skewedGrid.x),
                               static_cast<double>(iy) / static_cast<double>(skewedGrid.y), 0.25);
      const double distance = blockingSphereDistance(fractional, skewedCell, spheres);
      EXPECT_NEAR(distance, distanceBySearch(fractional), 1.0e-9)
          << "at fractions " << fractional.x << ", " << fractional.y;
      if (distance < 0.0) ++inside;
    }
  }
  EXPECT_GT(inside, 100uz);
}

// And the device finds the same copy: the field it makes of a skewed cell with a wide sphere in it is the
// processor's, blocked at the same voxels by the same amount.
TEST(energy_well_field, the_processor_and_the_device_block_the_same_voxels_of_a_skewed_cell)
{
  OpenCL::initialize();
  if (!OpenCL::clContext.has_value()) GTEST_SKIP() << "no OpenCL device";

  PairInteractions interactions = probeAndAtom();
  LinearProbe probe = LinearProbe::singleSite(interactions, "X").value();
  const std::vector<BlockingSphere> spheres{wideSphere};
  Crystal framework = oneAtomInASkewedCell();

  WellField cpu =
      WellFieldCPU::compute(interactions, framework, probe, skewedGrid, 1, 0.0, spheres, 1000.0, ceiling);
  WellField gpu =
      WellFieldOpenCL::compute(interactions, framework, probe, skewedGrid, 1, 0.0, spheres, 1000.0, ceiling);
  ASSERT_EQ(gpu.numberOfVoxels(), cpu.numberOfVoxels());

  // the voxel at fractions (0.45, 0.45, 0), which the copy nearest in fractions leaves outside the sphere
  const std::size_t told = 9 * skewedGrid.x + 9;
  const float expected = static_cast<float>(-1000.0 * distanceBySearch(double3(0.45, 0.45, 0.0)));
  EXPECT_NEAR(cpu.energy[told], expected, 1.0e-3f * expected);
  EXPECT_NEAR(gpu.energy[told], expected, 1.0e-3f * expected);

  std::size_t different = 0;
  for (std::size_t voxel = 0; voxel < cpu.numberOfVoxels(); ++voxel)
  {
    if (std::abs(gpu.energy[voxel] - cpu.energy[voxel]) > 1.0e-3f * std::max(std::abs(cpu.energy[voxel]), 1.0f))
      ++different;
  }
  EXPECT_EQ(different, 0uz);
}
