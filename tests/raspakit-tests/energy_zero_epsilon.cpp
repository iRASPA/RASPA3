#include <gtest/gtest.h>

import std;

import uint3;
import double3;
import unit_cell;
import crystal;
import pair_interactions;
import blocking_spheres;
import energy_shared_linear_probe;
import energy_shared_probe_energy_grid;
import energy_probe_energy_grid;
import energy_opencl_probe_energy_grid;
import energy_shared_molecular_energy_grid;
import energy_opencl_molecular_energy_grid;
import energy_shared_well_field;
import energy_well_field;
import energy_opencl_well_field;
import opencl;

// Pairs with no strength on the device. A framework atom with no dispersion of its own mixes to a pair of
// zero strength and nonzero size with every probe, and the pair's energy is nought at every separation. On the
// device the ratio of size to separation is a float, and on a grid point lying on the atom's centre, where the
// separation is held at its floor, the sixth power of that ratio overflows: nought times an overflow is not
// nought but not a number, and the ceiling then stood in for it. The CPU works in double, where nothing
// overflows, and is the reference.

namespace
{

constexpr double edge = 10.0;
const uint3 gridSize{8, 8, 8};

// Two types: the probe, with a strength of its own, and an atom with none. The mixed pair keeps a size.
PairInteractions probeAndStrengthlessAtom()
{
  PairInteractions interactions;
  interactions.numberOfTypes = 2;
  interactions.names = {"X", "H"};
  interactions.charges = {0.0, 0.0};
  interactions.parameters = {PairParameters{3.0, 100.0, 0.0}, PairParameters{2.75, 0.0, 0.0},
                             PairParameters{2.75, 0.0, 0.0}, PairParameters{2.5, 0.0, 0.0}};
  interactions.cutOffVDW = 4.0;
  interactions.cutOffCoulomb = 4.0;
  return interactions;
}

// One strengthless atom on the origin, which is the first point of the grid, and one atom with strength in
// the middle of the cell, so that the probe has some dispersion somewhere and is not refused outright.
Crystal atomOnTheOrigin()
{
  Crystal framework;
  framework.name = "strengthless";
  framework.unitCell = UnitCell(edge, edge, edge);
  framework.mass = 1.0;
  framework.atoms.push_back(CrystalAtom{double3(0.0, 0.0, 0.0), 0.0, 1});
  framework.fractionalPositions.push_back(double3(0.0, 0.0, 0.0));
  framework.atoms.push_back(CrystalAtom{double3(0.5 * edge, 0.5 * edge, 0.5 * edge), 0.0, 0});
  framework.fractionalPositions.push_back(double3(0.5, 0.5, 0.5));
  return framework;
}

bool noDevice()
{
  OpenCL::initialize();
  return !OpenCL::clContext.has_value();
}

}  // namespace


TEST(energy_zero_epsilon, the_probe_grid_feels_nothing_from_a_strengthless_atom_even_on_its_centre)
{
  if (noDevice()) GTEST_SKIP() << "no OpenCL device";
  PairInteractions interactions = probeAndStrengthlessAtom();
  Crystal framework = atomOnTheOrigin();

  ProbeEnergyGrid cpu = ProbeEnergyGridCPU::compute(interactions, framework, "X", gridSize);
  ProbeEnergyGrid gpu = ProbeEnergyGridOpenCL::compute(interactions, framework, "X", gridSize);

  ASSERT_EQ(gpu.energy.size(), cpu.energy.size());
  EXPECT_EQ(cpu.energy[0], 0.0f);
  EXPECT_EQ(gpu.energy[0], 0.0f) << "the point on the atom's centre";
  // Elsewhere the two agree to float precision, the device working in float where the CPU works in double.
  for (std::size_t voxel = 0; voxel < cpu.energy.size(); ++voxel)
  {
    EXPECT_NEAR(gpu.energy[voxel], cpu.energy[voxel], 1.0e-5 * std::max(1.0f, std::abs(cpu.energy[voxel])))
        << "voxel " << voxel;
  }
}


TEST(energy_zero_epsilon, the_molecular_grid_feels_nothing_from_a_strengthless_atom_even_on_its_centre)
{
  if (noDevice()) GTEST_SKIP() << "no OpenCL device";
  PairInteractions interactions = probeAndStrengthlessAtom();
  Crystal framework = atomOnTheOrigin();
  LinearProbe probe = LinearProbe::singleSite(interactions, "X").value();

  MolecularEnergyGrid grid = MolecularEnergyGridOpenCL::compute(interactions, framework, probe, gridSize, 1, 300.0);

  ASSERT_EQ(grid.numberOfVoxels(), static_cast<std::size_t>(gridSize.x) * gridSize.y * gridSize.z);
  // The origin is a cutoff and more from the atom with strength, so the only thing felt there is nothing.
  EXPECT_NEAR(grid.freeEnergy[0], 0.0f, 1.0e-4) << "the point on the strengthless atom's centre";
}


TEST(energy_zero_epsilon, the_well_field_feels_nothing_from_a_strengthless_atom_even_on_its_centre)
{
  if (noDevice()) GTEST_SKIP() << "no OpenCL device";
  PairInteractions interactions = probeAndStrengthlessAtom();
  Crystal framework = atomOnTheOrigin();
  LinearProbe probe = LinearProbe::singleSite(interactions, "X").value();
  const double ceiling = 1.0e7;

  WellField cpu = WellFieldCPU::compute(interactions, framework, probe, gridSize, 1, 0.0, {}, 1000.0, ceiling);
  WellField gpu = WellFieldOpenCL::compute(interactions, framework, probe, gridSize, 1, 0.0, {}, 1000.0, ceiling);

  ASSERT_EQ(gpu.numberOfVoxels(), cpu.numberOfVoxels());
  for (std::size_t voxel = 0; voxel < cpu.numberOfVoxels(); ++voxel)
  {
    EXPECT_NEAR(gpu.energy[voxel], cpu.energy[voxel], 1.0e-5 * std::max(1.0, std::abs(double(cpu.energy[voxel]))))
        << "voxel " << voxel;
  }
  EXPECT_NEAR(gpu.energy[0], 0.0, 1.0e-4) << "the point on the atom's centre";
}
