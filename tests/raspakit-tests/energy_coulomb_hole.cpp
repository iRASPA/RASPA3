#include <gtest/gtest.h>

import std;

import uint3;
import double3;
import unit_cell;
import crystal;
import units;
import pair_interactions;
import energy_shared_linear_probe;
import energy_shared_electrostatic_potential_grid;
import energy_shared_molecular_energy_grid;
import energy_molecular_energy_grid;
import energy_opencl_molecular_energy_grid;
import opencl;

// A charged site on top of an oppositely charged atom. The dispersion of the pair is held at the ceiling,
// but its Coulomb term, taken on its own, is a hole under the ceiling deeper than any well: at the floor the
// separation is held at, it is minus a hundred million Kelvin and more, and the Boltzmann weight of the point
// is infinite. The pair is an overlap, so its Coulomb term goes with its dispersion and the orientation is
// at the ceiling. A site with no dispersion at all is given a core of one Ångström, as in the well field.

namespace
{

constexpr double edge = 20.0;
const uint3 gridSize{8, 8, 8};
constexpr double alpha = 0.1;

// The probe with dispersion, a probe without, and the atom.
PairInteractions twoProbesAndAnAtom()
{
  PairInteractions interactions;
  interactions.numberOfTypes = 3;
  interactions.names = {"Q", "M", "O"};
  interactions.charges = {1.0, 1.0, -1.0};
  interactions.parameters = {PairParameters{3.0, 100.0, 0.0}, PairParameters{0.0, 0.0, 0.0},
                             PairParameters{3.0, 100.0, 0.0}, PairParameters{0.0, 0.0, 0.0},
                             PairParameters{0.0, 0.0, 0.0},   PairParameters{0.0, 0.0, 0.0},
                             PairParameters{3.0, 100.0, 0.0}, PairParameters{0.0, 0.0, 0.0},
                             PairParameters{3.0, 100.0, 0.0}};
  interactions.cutOffVDW = 8.0;
  interactions.cutOffCoulomb = 8.0;
  return interactions;
}

// One atom on the origin, which is the first grid point.
Crystal atomOnTheOrigin()
{
  Crystal framework;
  framework.name = "hole";
  framework.unitCell = UnitCell(edge, edge, edge);
  framework.mass = 16.0;
  framework.atoms.push_back(CrystalAtom{double3(0.0, 0.0, 0.0), -1.0, 2});
  framework.fractionalPositions.push_back(double3(0.0, 0.0, 0.0));
  return framework;
}

LinearProbe chargedSite(std::size_t type, const std::string &name)
{
  LinearProbe probe;
  probe.name = name;
  probe.headTailSymmetric = true;
  probe.sites.push_back(LinearProbe::Site{type, 0.0, 1.0, name});
  return probe;
}

ElectrostaticPotentialGrid nearHalfOnly(const Crystal &framework)
{
  ElectrostaticPotentialGrid potential;
  potential.gridSize = gridSize;
  potential.unitCell = framework.unitCell;
  potential.alpha = alpha;
  potential.cutOff = 8.0;
  potential.smoothPotential.assign(static_cast<std::size_t>(gridSize.x) * gridSize.y * gridSize.z, 0.0f);
  return potential;
}

bool noDevice()
{
  OpenCL::initialize();
  return !OpenCL::clContext.has_value();
}

}  // namespace


TEST(energy_coulomb_hole, a_charged_site_on_an_opposite_atom_is_at_the_ceiling)
{
  PairInteractions interactions = twoProbesAndAnAtom();
  Crystal framework = atomOnTheOrigin();
  LinearProbe probe = chargedSite(0, "Q");
  ElectrostaticPotentialGrid potential = nearHalfOnly(framework);

  MolecularEnergyGrid grid =
      MolecularEnergyGridCPU::compute(interactions, framework, probe, gridSize, 1, 300.0, &potential);
  ASSERT_TRUE(grid.chargesIncluded);
  EXPECT_FLOAT_EQ(grid.minimumEnergy[0], static_cast<float>(grid.ceiling));
  EXPECT_FLOAT_EQ(grid.freeEnergy[0], static_cast<float>(grid.ceiling));

  // One grid step away the pair is not overlapping, and its Coulomb term is kept in full.
  MolecularEnergyGrid without = MolecularEnergyGridCPU::compute(interactions, framework, probe, gridSize, 1, 300.0);
  const double r = edge / static_cast<double>(gridSize.x);
  const double expected = -Units::CoulombicConversionFactor * std::erfc(alpha * r) / r;
  const double felt = static_cast<double>(grid.freeEnergy[1]) - static_cast<double>(without.freeEnergy[1]);
  EXPECT_NEAR(felt, expected, 1.0e-4 * std::abs(expected));
}


TEST(energy_coulomb_hole, a_site_without_dispersion_on_an_opposite_atom_is_at_the_ceiling)
{
  PairInteractions interactions = twoProbesAndAnAtom();
  Crystal framework = atomOnTheOrigin();
  LinearProbe probe = chargedSite(1, "M");
  ElectrostaticPotentialGrid potential = nearHalfOnly(framework);

  MolecularEnergyGrid grid =
      MolecularEnergyGridCPU::compute(interactions, framework, probe, gridSize, 1, 300.0, &potential);
  ASSERT_TRUE(grid.chargesIncluded);
  EXPECT_FLOAT_EQ(grid.freeEnergy[0], static_cast<float>(grid.ceiling));

  // Beyond the core the bare charge keeps its full Coulomb term.
  const double r = edge / static_cast<double>(gridSize.x);
  const double expected = -Units::CoulombicConversionFactor * std::erfc(alpha * r) / r;
  EXPECT_NEAR(static_cast<double>(grid.freeEnergy[1]), expected, 1.0e-4 * std::abs(expected));
}


TEST(energy_coulomb_hole, the_device_answers_the_same)
{
  if (noDevice()) GTEST_SKIP() << "no OpenCL device";
  PairInteractions interactions = twoProbesAndAnAtom();
  Crystal framework = atomOnTheOrigin();
  ElectrostaticPotentialGrid potential = nearHalfOnly(framework);

  for (const LinearProbe &probe : {chargedSite(0, "Q"), chargedSite(1, "M")})
  {
    MolecularEnergyGrid cpu =
        MolecularEnergyGridCPU::compute(interactions, framework, probe, gridSize, 1, 300.0, &potential);
    MolecularEnergyGrid gpu =
        MolecularEnergyGridOpenCL::compute(interactions, framework, probe, gridSize, 1, 300.0, &potential);
    ASSERT_EQ(gpu.freeEnergy.size(), cpu.freeEnergy.size());
    EXPECT_FLOAT_EQ(gpu.freeEnergy[0], static_cast<float>(cpu.ceiling)) << probe.name;
    for (std::size_t voxel = 0; voxel < cpu.freeEnergy.size(); ++voxel)
    {
      const double scale = std::max(1.0, std::abs(static_cast<double>(cpu.freeEnergy[voxel])));
      EXPECT_NEAR(gpu.freeEnergy[voxel], cpu.freeEnergy[voxel], 1.0e-4 * scale) << probe.name << " voxel " << voxel;
    }
  }
}
