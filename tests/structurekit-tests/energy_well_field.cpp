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
import energy_probe_energy_grid;
import energy_shared_probe_energy_grid;
import energy_shared_molecular_energy_grid;
import energy_shared_molecular_void_fraction;
import units;
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

// A mesh vertex inside a blocking sphere has the blocked ramp's energy, as the field has there, on the processor and
// on the device. The refinement used to keep the bare framework energy, so a sphere cap in an attractive blocked
// cage came out below the cap threshold and its triangles were not trimmed.
TEST(energy_well_field, a_vertex_in_a_blocking_sphere_is_refined_to_the_blocked_energy)
{
  PairInteractions interactions = probeAndAtom();
  const std::vector<BlockingSphere> spheres{wideSphere};
  Crystal framework = oneAtomInASkewedCell();
  NeighbourhoodParameters parameters{};
  parameters.probe = LinearProbe::singleSite(interactions, "X").value();
  parameters.blockedEnergyPerAngstrom = 1000.0;
  parameters.ceiling = ceiling;

  const double3 corner(0.45, 0.45, 0.0);
  const double expected = -1000.0 * distanceBySearch(corner);
  ASSERT_GT(expected, 0.0);

  std::vector<double3> corners{corner};
  std::vector<double> energies;
  WellFieldCPU::refineVertices(corners, energies, interactions, framework, parameters, spheres, -50.0);
  ASSERT_EQ(energies.size(), 1uz);
  EXPECT_NEAR(energies[0], expected, 1.0e-4 * expected);  // the vertex is quantized to 2^-20 of the cell

  OpenCL::initialize();
  if (!OpenCL::clContext.has_value()) GTEST_SKIP() << "no OpenCL device";
  corners = {corner};
  energies.clear();
  WellFieldOpenCL::refineVertices(corners, energies, interactions, framework, parameters, spheres, -50.0);
  ASSERT_EQ(energies.size(), 1uz);
  EXPECT_NEAR(energies[0], expected, 1.0e-3 * expected);
}


namespace
{

// Two sites of the probe's own kind, 'half' either side of the centre.
LinearProbe dumbbell(double half, bool headTailSymmetric)
{
  LinearProbe probe;
  probe.name = "dumbbell";
  probe.sites = {LinearProbe::Site{0, -half, 0.0, "X"}, LinearProbe::Site{0, half, 0.0, "X"}};
  probe.headTailSymmetric = headTailSymmetric;
  return probe;
}

PairInteractions smallProbeAndAtom()
{
  PairInteractions interactions = probeAndAtom();
  interactions.parameters = {PairParameters{1.0, 100.0, 0.0}, PairParameters{1.0, 100.0, 0.0},
                             PairParameters{1.0, 100.0, 0.0}, PairParameters{1.0, 100.0, 0.0}};
  return interactions;
}

}  // namespace


// An orientation with a site on top of an atom is an overlap, and is the ceiling: the pair is held there, and the
// attraction the other site feels three Ångström off is not left under it (it was, by half a unit here).
TEST(energy_well_field, an_orientation_with_a_site_inside_an_atom_is_at_the_ceiling)
{
  PairInteractions interactions = smallProbeAndAtom();
  Crystal framework = oneAtomAt(double3(0.0, 0.0, 0.15), true);
  const double lowCeiling = 1.0e4;

  // one orientation, along z: the site at +1.5 Å is on the atom, the one at -1.5 Å is 3 Å from it
  WellField field =
      WellFieldCPU::compute(interactions, framework, dumbbell(1.5, true), gridSize, 1, 0.0, {}, 1000.0, lowCeiling);
  EXPECT_EQ(field.energy[0], static_cast<float>(lowCeiling));

  OpenCL::initialize();
  if (!OpenCL::clContext.has_value()) return;
  WellField device =
      WellFieldOpenCL::compute(interactions, framework, dumbbell(1.5, true), gridSize, 1, 0.0, {}, 1000.0, lowCeiling);
  EXPECT_EQ(device.energy[0], static_cast<float>(lowCeiling));
}

// A framework without atoms has a flat field, from either builder, and not a complaint about the probe.
TEST(energy_well_field, a_framework_without_atoms_has_a_flat_field)
{
  PairInteractions interactions = probeAndAtom();
  Crystal framework;
  framework.name = "empty";
  framework.unitCell = UnitCell(edge, edge, edge);
  LinearProbe probe = LinearProbe::singleSite(interactions, "X").value();

  WellField field;
  ASSERT_NO_THROW(field = WellFieldCPU::compute(interactions, framework, probe, gridSize, 1, 0.0, {}, 1000.0, ceiling));
  ASSERT_EQ(field.numberOfVoxels(), 125u);
  for (float energy : field.energy) EXPECT_EQ(energy, 0.0f);

  OpenCL::initialize();
  if (!OpenCL::clContext.has_value()) return;
  ASSERT_NO_THROW(field =
                      WellFieldOpenCL::compute(interactions, framework, probe, gridSize, 1, 0.0, {}, 1000.0, ceiling));
  for (float energy : field.energy) EXPECT_EQ(energy, 0.0f);
}

// The device kernel holds 128 orientations and hands a field of more to the processor; the vertices of that field
// are refined there too, and were not: the refusal went uncaught and ended the run.
TEST(energy_well_field, vertices_are_refined_for_more_orientations_than_the_device_holds)
{
  OpenCL::initialize();
  if (!OpenCL::clContext.has_value()) GTEST_SKIP() << "no OpenCL device";

  PairInteractions interactions = probeAndAtom();
  Crystal framework = oneAtomAt(double3(0.5, 0.5, 0.5), true);
  NeighbourhoodParameters parameters{};
  parameters.probe = dumbbell(0.5, false);
  parameters.numberOfOrientations = 600;
  parameters.thermalEnergy = 1.0;
  parameters.blockedEnergyPerAngstrom = 1000.0;
  parameters.ceiling = ceiling;

  const std::vector<double3> start{double3(0.5, 0.5, 0.9), double3(0.1, 0.5, 0.5)};
  std::vector<double3> onProcessor = start;
  std::vector<double> processorEnergies;
  WellFieldCPU::refineVertices(onProcessor, processorEnergies, interactions, framework, parameters, {}, -50.0);

  std::vector<double3> onDevice = start;
  std::vector<double> deviceEnergies;
  ASSERT_NO_THROW(
      WellFieldOpenCL::refineVertices(onDevice, deviceEnergies, interactions, framework, parameters, {}, -50.0));
  ASSERT_EQ(deviceEnergies.size(), processorEnergies.size());
  for (std::size_t i = 0; i < deviceEnergies.size(); ++i)
  {
    EXPECT_EQ(deviceEnergies[i], processorEnergies[i]) << "vertex " << i;
  }
}

// A voxel on an atom's nucleus is that atom's. Its neighbour 3.5 Å off pulls on the point, as every atom beyond
// sigma does, and was given it for being the only one that pulls.
TEST(energy_well_field, a_voxel_inside_an_atom_is_labelled_with_that_atom)
{
  PairInteractions interactions = probeAndAtom();
  Crystal framework = oneAtomAt(double3(0.0, 0.0, 0.0), true);
  framework.atoms.push_back(CrystalAtom{double3(3.5, 0.0, 0.0), 0.0, 1});
  framework.fractionalPositions.push_back(double3(0.35, 0.0, 0.0));

  ProbeEnergyGrid grid = ProbeEnergyGridCPU::compute(interactions, framework, "X", gridSize);
  ASSERT_EQ(grid.strongestAtom.size(), 125uz);
  EXPECT_EQ(grid.strongestAtom[0], 0);
}

// A site deeper than 709 kT is more than exp holds in a double. The chemical potential is its logarithm and is
// a number all the same; the average is held at the largest double rather than at infinity.
TEST(energy_well_field, a_boltzmann_average_past_the_range_of_a_double_keeps_its_chemical_potential)
{
  const double temperature = 77.0;
  const double kT = Units::KB * temperature;

  MolecularEnergyGrid grid;
  grid.gridSize = uint3(2, 1, 1);
  grid.freeEnergy = {static_cast<float>(-800.0 * kT), 0.0f};
  grid.minimumEnergy = grid.freeEnergy;

  Crystal framework = oneAtomAt(double3(0.0, 0.0, 0.0), true);
  MolecularVoidFraction result = MolecularVoidFraction::fromGrid(grid, framework, temperature);

  // ln((e^800 + 1) / 2) = 800 - ln 2
  EXPECT_NEAR(result.excessChemicalPotential / kT, -(800.0 - std::log(2.0)), 1.0e-3);
  EXPECT_EQ(result.boltzmannAverage, std::numeric_limits<double>::max());
  EXPECT_FALSE(result.readsAsFraction);
}

// The other end: every site more than 700 kT up (a closed or repulsive framework). The logarithm is a number there
// as well, and so is the excess chemical potential; it was left at the ideal-gas zero.
TEST(energy_well_field, a_repulsive_grid_keeps_its_large_excess_chemical_potential)
{
  const double temperature = 77.0;
  const double kT = Units::KB * temperature;

  MolecularEnergyGrid grid;
  grid.gridSize = uint3(2, 1, 1);
  grid.freeEnergy = {static_cast<float>(800.0 * kT), static_cast<float>(900.0 * kT)};
  grid.minimumEnergy = grid.freeEnergy;

  Crystal framework = oneAtomAt(double3(0.0, 0.0, 0.0), true);
  MolecularVoidFraction result = MolecularVoidFraction::fromGrid(grid, framework, temperature);

  // -ln((e^-800 + e^-900) / 2) = 800 + ln 2 (to within e^-100)
  EXPECT_NEAR(result.excessChemicalPotential / kT, 800.0 + std::log(2.0), 1.0e-3);
}
