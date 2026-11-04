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

// The molecular grid walks the periodic images of every framework atom and throws an image out, before its
// separation is built, when a lower bound on that separation already exceeds the cutoff. There are two cutoffs
// in play, the dispersion one and the electrostatic one, and the walk has to reach as far as the longer of
// them: an image that is beyond the dispersion cutoff but within the electrostatic one still owes the near
// half of its Ewald sum. The tests below place the probe where every atom is in exactly that band.

namespace
{

constexpr double edge = 30.0;
const uint3 gridSize{30, 30, 30};
constexpr double alpha = 0.1;
constexpr double probeCharge = 1.0;

// A dispersion cutoff well short of the electrostatic one.
PairInteractions shortDispersionLongCoulomb()
{
  PairInteractions interactions;
  interactions.numberOfTypes = 2;
  interactions.names = {"O", "Q"};
  interactions.charges = {0.0, probeCharge};
  interactions.parameters = {PairParameters{3.0, 100.0, 0.0}, PairParameters{3.0, 100.0, 0.0},
                             PairParameters{3.0, 100.0, 0.0}, PairParameters{3.0, 100.0, 0.0}};
  interactions.cutOffVDW = 4.0;
  interactions.cutOffCoulomb = 12.0;
  return interactions;
}

// Two charged atoms on the x axis: one on the origin, one halfway across the cell. Their nearest images are
// the only ones within either cutoff.
Crystal twoChargedAtoms()
{
  Crystal framework;
  framework.name = "two charges";
  framework.unitCell = UnitCell(edge, edge, edge);
  framework.mass = 2.0;
  framework.atoms.push_back(CrystalAtom{double3(0.0, 0.0, 0.0), -1.0, 0});
  framework.fractionalPositions.push_back(double3(0.0, 0.0, 0.0));
  framework.atoms.push_back(CrystalAtom{double3(0.5 * edge, 0.0, 0.0), 0.4, 0});
  framework.fractionalPositions.push_back(double3(0.5, 0.0, 0.0));
  return framework;
}

LinearProbe chargedSite()
{
  LinearProbe probe;
  probe.name = "Q";
  probe.headTailSymmetric = true;
  probe.sites.push_back(LinearProbe::Site{1, 0.0, probeCharge, "Q"});
  return probe;
}

// A potential with no far half at all, so that whatever the grid feels from the charges is the near half
// walked pair by pair.
ElectrostaticPotentialGrid nearHalfOnly(const Crystal &framework)
{
  ElectrostaticPotentialGrid potential;
  potential.gridSize = gridSize;
  potential.unitCell = framework.unitCell;
  potential.alpha = alpha;
  potential.cutOff = 12.0;
  potential.smoothPotential.assign(static_cast<std::size_t>(gridSize.x) * gridSize.y * gridSize.z, 0.0f);
  return potential;
}

// The grid point 8 Å along x: 8 Å from the atom on the origin and 7 Å from the other, both beyond the
// dispersion cutoff and within the electrostatic one.
constexpr std::size_t probeVoxel = 8;

}  // namespace


TEST(energy_molecular_cutoff, the_near_half_of_the_ewald_sum_reaches_past_the_dispersion_cutoff)
{
  PairInteractions interactions = shortDispersionLongCoulomb();
  Crystal framework = twoChargedAtoms();
  LinearProbe probe = chargedSite();
  ElectrostaticPotentialGrid potential = nearHalfOnly(framework);

  MolecularEnergyGrid without = MolecularEnergyGridCPU::compute(interactions, framework, probe, gridSize, 1, 300.0);
  MolecularEnergyGrid with =
      MolecularEnergyGridCPU::compute(interactions, framework, probe, gridSize, 1, 300.0, &potential);
  ASSERT_TRUE(without.chargesIgnored);
  ASSERT_TRUE(with.chargesIncluded);

  // The dispersion is the same in both and cancels; what is left is the two screened pair terms.
  const double expected = Units::CoulombicConversionFactor * probeCharge *
                          (-1.0 * std::erfc(alpha * 8.0) / 8.0 + 0.4 * std::erfc(alpha * 7.0) / 7.0);
  const double felt =
      static_cast<double>(with.freeEnergy[probeVoxel]) - static_cast<double>(without.freeEnergy[probeVoxel]);
  EXPECT_NEAR(felt, expected, 1.0e-4 * std::abs(expected));
}


TEST(energy_molecular_cutoff, an_atom_pulling_only_through_its_charge_still_holds_the_point)
{
  PairInteractions interactions = shortDispersionLongCoulomb();
  Crystal framework = twoChargedAtoms();
  LinearProbe probe = chargedSite();
  ElectrostaticPotentialGrid potential = nearHalfOnly(framework);

  std::vector<std::int32_t> labels =
      MolecularEnergyGridCPU::strongestAtoms(interactions, framework, probe, gridSize, 1, 300.0, &potential);
  ASSERT_EQ(labels.size(), static_cast<std::size_t>(gridSize.x) * gridSize.y * gridSize.z);

  // The atom on the origin attracts the probe and the other repels it; both act only through their charges
  // here, and the attracting one is the stronger pull.
  EXPECT_EQ(labels[probeVoxel], 0);
}
