#include <gtest/gtest.h>

import std;

import int3;
import uint3;
import double3;
import double3x3;
import units;
import interpolation_energy_grid;
import atom;
import pseudo_atom;
import vdwparameters;
import forcefield;
import framework;
import component;
import system;
import simulationbox;
import running_energy;
import energy_status;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;

TEST(static_energy, Test_2_CO2_in_ITQ_29_1x1x1)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework f = Framework::makeITQ29(forceField, int3(1, 1, 1));
  Component c = Component::makeCO2(forceField, 0, true);

  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {c}, {}, {2}, 5);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  atomData[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  atomData[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  atomData[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  atomData[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  atomData[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  atomData[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);

  RunningEnergy energy =
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomData) +
      Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                   system.framework, frameworkAtoms, atomData);

  EXPECT_NEAR(energy.frameworkMoleculeVDW * Units::EnergyToKelvin, -1545.62921755, 1e-6);
  EXPECT_NEAR(energy.frameworkMoleculeCharge * Units::EnergyToKelvin, -592.13188606, 1e-6);
  EXPECT_NEAR(energy.moleculeMoleculeVDW * Units::EnergyToKelvin, -242.42459776, 1e-6);
  EXPECT_NEAR(energy.moleculeMoleculeCharge * Units::EnergyToKelvin, 154.11883595, 1e-6);
}

TEST(static_energy, Test_2_CO2_in_MFI_2x2x2_shifted)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component c = Component::makeCO2(forceField, 0, true);
  Framework f = Framework::makeMFI(forceField, int3(2, 2, 2));

  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {c}, {}, {2}, 5);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  atomData[0].position = double3(10.011, 4.97475 + 2.0, 1.149);
  atomData[1].position = double3(10.011, 4.97475 + 2.0, 0.0);
  atomData[2].position = double3(10.011, 4.97475 + 2.0, -1.149);
  atomData[3].position = double3(10.011, 4.97475 - 2.0, 1.149);
  atomData[4].position = double3(10.011, 4.97475 - 2.0, 0.0);
  atomData[5].position = double3(10.011, 4.97475 - 2.0, -1.149);

  RunningEnergy energy =
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomData) +
      Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                   system.framework, frameworkAtoms, atomData);

  EXPECT_NEAR(energy.frameworkMoleculeVDW * Units::EnergyToKelvin, -2525.36580663, 1e-6);
  EXPECT_NEAR(energy.frameworkMoleculeCharge * Units::EnergyToKelvin, 2167.45591472, 1e-6);
  EXPECT_NEAR(energy.moleculeMoleculeVDW * Units::EnergyToKelvin, -242.42459776, 1e-6);
  EXPECT_NEAR(energy.moleculeMoleculeCharge * Units::EnergyToKelvin, 154.11883595, 1e-6);
}

TEST(static_energy, Test_2_CO2_in_MFI_2x2x2_truncated)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, false, true, true);

  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);
  Framework f = Framework::makeMFI(forceField, int3(2, 2, 2));
  Component c = Component::makeCO2(forceField, 0, true);
  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {c}, {}, {2}, 5);

  std::span<Atom> moleculeAtomPositions = system.spanOfMoleculeAtoms();
  moleculeAtomPositions[0].position = double3(10.011, 4.97475 + 2.0, 1.149);
  moleculeAtomPositions[1].position = double3(10.011, 4.97475 + 2.0, 0.0);
  moleculeAtomPositions[2].position = double3(10.011, 4.97475 + 2.0, -1.149);
  moleculeAtomPositions[3].position = double3(10.011, 4.97475 - 2.0, 1.149);
  moleculeAtomPositions[4].position = double3(10.011, 4.97475 - 2.0, 0.0);
  moleculeAtomPositions[5].position = double3(10.011, 4.97475 - 2.0, -1.149);

  RunningEnergy energy = system.computeTotalEnergies();

  EXPECT_NEAR(energy.frameworkMoleculeVDW * Units::EnergyToKelvin, -2657.36121975, 1e-6);
  EXPECT_NEAR(energy.frameworkMoleculeCharge * Units::EnergyToKelvin, 1971.00612979, 1e-6);
  EXPECT_NEAR(energy.moleculeMoleculeVDW * Units::EnergyToKelvin, -242.94298709, 1e-6);
  EXPECT_NEAR(energy.moleculeMoleculeCharge * Units::EnergyToKelvin, 162.41877650, 1e-6);
  EXPECT_NEAR(energy.tail * Units::EnergyToKelvin, -127.81601515, 1e-6);
  // EXPECT_NEAR(energy.tail * Units::EnergyToKelvin, -127.72803736223419, 1e-6);
  EXPECT_NEAR((energy.ewald_fourier + energy.ewald_self + energy.ewald_exclusion) * Units::EnergyToKelvin,
              -1197.23909965, 1e-6);
}

// The Brick-CFCMC-style aggregated tail correction must reproduce the per-atom O(N^2) reference exactly, for
// both the tail energy and every per-group dU/dlambda accumulator, at all fractional-molecule lambda values.
TEST(static_energy, tail_correction_aggregated_matches_per_atom_reference)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, false, true, true);
  Framework f = Framework::makeMFI(forceField, int3(2, 2, 2));
  Component c = Component::makeCO2(forceField, 0, true);
  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {c}, {}, {4}, 5);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  ASSERT_GE(atomData.size(), static_cast<std::size_t>(6));

  // Tail corrections are position-independent, but give the atoms distinct positions anyway.
  for (std::size_t i = 0; i < atomData.size(); ++i)
  {
    atomData[i].position =
        double3(1.0 + 0.31 * static_cast<double>(i), 2.0 + 0.13 * static_cast<double>(i),
                3.0 + 0.19 * static_cast<double>(i));
  }

  // Tag the first molecule as fractional in dU/dlambda group 1, the second as fractional in group 2, and leave
  // the remaining molecules as integer molecules; sweep lambda over the full [0, 1] range.
  for (double lambda : {0.0, 0.15, 0.37, 0.5, 0.66, 0.83, 1.0})
  {
    for (std::size_t i = 0; i < 3; ++i) atomData[i].setScalingToFractional(lambda, std::uint8_t{1});
    for (std::size_t i = 3; i < 6; ++i) atomData[i].setScalingToFractional(1.0 - lambda, std::uint8_t{2});

    RunningEnergy reference =
        Interactions::computeInterMolecularTailEnergyReference(system.forceField, system.simulationBox, atomData);
    RunningEnergy aggregated =
        Interactions::computeInterMolecularTailEnergy(system.forceField, system.simulationBox, atomData);

    EXPECT_NEAR(aggregated.tail, reference.tail, 1e-10);
    for (std::size_t g = 0; g < maximumNumberOfDUDlambdaGroups; ++g)
    {
      EXPECT_NEAR(aggregated.dudlambdaVDW[g], reference.dudlambdaVDW[g], 1e-10);
    }
  }
}

// The aggregated difference (used on the CFCMC move hot path) must equal the per-atom tail-energy difference for
// a fractional-molecule lambda change, for both the tail energy and every per-group dU/dlambda accumulator.
TEST(static_energy, tail_correction_aggregated_difference_matches_per_atom_reference)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, false, true, true);
  Framework f = Framework::makeMFI(forceField, int3(2, 2, 2));
  Component c = Component::makeCO2(forceField, 0, true);
  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {c}, {}, {4}, 5);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  ASSERT_GE(atomData.size(), static_cast<std::size_t>(6));

  for (std::size_t i = 0; i < atomData.size(); ++i)
  {
    atomData[i].position =
        double3(1.0 + 0.31 * static_cast<double>(i), 2.0 + 0.13 * static_cast<double>(i),
                3.0 + 0.19 * static_cast<double>(i));
  }

  // Current state: first molecule fractional (group 1) at lambda = 0.4, everything else integer.
  for (std::size_t i = 0; i < 3; ++i) atomData[i].setScalingToFractional(0.4, std::uint8_t{1});

  system.computeTailCorrectionCounts();

  // Trial state: the fractional molecule changes lambda to 0.7 (positions unchanged).
  std::vector<Atom> oldFractional(atomData.begin(), atomData.begin() + 3);
  std::vector<Atom> newFractional(oldFractional.begin(), oldFractional.end());
  for (Atom& atom : newFractional) atom.setScalingToFractional(0.7, std::uint8_t{1});

  RunningEnergy reference = Interactions::computeInterMolecularTailEnergyDifference(
      system.forceField, system.simulationBox, atomData, newFractional, oldFractional);
  RunningEnergy aggregated = Interactions::computeInterMolecularTailEnergyDifferenceAggregated(
      system.forceField, system.simulationBox, system.effectiveNumberOfPseudoAtomsVDW,
      system.fractionalPseudoAtomCountsPerGroup, newFractional, oldFractional);

  EXPECT_NEAR(aggregated.tail, reference.tail, 1e-10);
  for (std::size_t g = 0; g < maximumNumberOfDUDlambdaGroups; ++g)
  {
    EXPECT_NEAR(aggregated.dudlambdaVDW[g], reference.dudlambdaVDW[g], 1e-10);
  }
}

// The virial needs the pair displacements, which the framework grids do not hold, so the pressure of a system with
// grids must equal the pressure of the same system without them, and report the same tail correction.
TEST(static_energy, molecular_pressure_is_independent_of_the_framework_grids)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, false, true, true);
  Component co2 = Component::makeCO2(forceField, 0, true);
  forceField.settings.gridPseudoAtomIndices = {static_cast<std::size_t>(co2.atoms[0].type),
                                               static_cast<std::size_t>(co2.atoms[1].type)};
  forceField.settings.numberOfVDWGridPoints = uint3(32, 32, 32);
  forceField.settings.numberOfCoulombGridPoints = uint3(32, 32, 32);
  forceField.settings.numberOfGridTestPoints = 10;
  Framework f = Framework::makeITQ29(forceField, int3(2, 2, 2));

  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {co2}, {}, {1}, 5);
  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  atomData[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  atomData[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  atomData[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  std::ostringstream log;
  system.createFrameworkInterpolationGrids(log);
  system.rebuildConfigurationDerivedState();
  ASSERT_TRUE(system.interpolationGrids[static_cast<std::size_t>(co2.atoms[0].type)].has_value());

  const std::pair<EnergyStatus, double3x3> withGrids = system.computeMolecularPressure();

  system.interpolationGrids.assign(system.interpolationGrids.size(), std::nullopt);
  const std::pair<EnergyStatus, double3x3> withoutGrids = system.computeMolecularPressure();

  const double3x3& a = withGrids.second;
  const double3x3& b = withoutGrids.second;
  EXPECT_NEAR(a.ax, b.ax, 1e-8);
  EXPECT_NEAR(a.ay, b.ay, 1e-8);
  EXPECT_NEAR(a.az, b.az, 1e-8);
  EXPECT_NEAR(a.bx, b.bx, 1e-8);
  EXPECT_NEAR(a.by, b.by, 1e-8);
  EXPECT_NEAR(a.bz, b.bz, 1e-8);
  EXPECT_NEAR(a.cx, b.cx, 1e-8);
  EXPECT_NEAR(a.cy, b.cy, 1e-8);
  EXPECT_NEAR(a.cz, b.cz, 1e-8);

  const double tail = withoutGrids.first.frameworkMoleculeEnergy.VanDerWaalsTailCorrection.energy;
  EXPECT_NE(tail, 0.0);
  EXPECT_NEAR(withGrids.first.frameworkMoleculeEnergy.VanDerWaalsTailCorrection.energy, tail, 1e-8);
}

// The framework grids tabulate the unscaled (lambda = 1) energies. A fractional molecule that is switched off
// for deletion (scaling 0, not fractional: Atom::setScalingOff, as in the pair/group CFCMC swaps) must not be
// read from the grid as if it were a fully coupled integer molecule.
TEST(static_energy, switched_off_molecule_is_not_read_from_the_framework_grid)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component co2 = Component::makeCO2(forceField, 0, true);
  forceField.settings.gridPseudoAtomIndices = {static_cast<std::size_t>(co2.atoms[0].type),
                                               static_cast<std::size_t>(co2.atoms[1].type)};
  forceField.settings.numberOfVDWGridPoints = uint3(32, 32, 32);
  forceField.settings.numberOfCoulombGridPoints = uint3(32, 32, 32);
  forceField.settings.numberOfGridTestPoints = 10;
  Framework f = Framework::makeITQ29(forceField, int3(2, 2, 2));

  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {co2}, {}, {1}, 5);
  std::ostringstream log;
  system.createFrameworkInterpolationGrids(log);
  ASSERT_TRUE(system.interpolationGrids[static_cast<std::size_t>(co2.atoms[0].type)].has_value());

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  atomData[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  atomData[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  atomData[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  const std::vector<std::optional<InterpolationEnergyGrid>> noGrids(system.interpolationGrids.size());

  // the grids are live: an integer molecule reads them (and agrees with the pair sum to interpolation accuracy)
  RunningEnergy integerFromGrid = Interactions::computeFrameworkMoleculeEnergy(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
      system.spanOfFrameworkAtoms(), atomData);
  RunningEnergy integerFromPairs = Interactions::computeFrameworkMoleculeEnergy(
      system.forceField, system.simulationBox, noGrids, system.framework, system.spanOfFrameworkAtoms(), atomData);
  EXPECT_NE(integerFromGrid.frameworkMoleculeVDW, integerFromPairs.frameworkMoleculeVDW);
  EXPECT_NEAR(integerFromGrid.frameworkMoleculeVDW, integerFromPairs.frameworkMoleculeVDW,
              0.05 * std::abs(integerFromPairs.frameworkMoleculeVDW));

  // deletion of a fractional molecule: lambda = 0.6 -> switched off
  std::vector<Atom> oldFractional(atomData.begin(), atomData.end());
  for (Atom& atom : oldFractional) atom.setScalingToFractional(0.6, std::uint8_t{1});
  std::vector<Atom> switchedOff(oldFractional.begin(), oldFractional.end());
  for (Atom& atom : switchedOff) atom.setScalingOff();

  std::optional<RunningEnergy> withGrids = Interactions::computeFrameworkMoleculeEnergyDifference(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
      system.frameworkCellList(), system.spanOfFrameworkAtoms(), switchedOff, oldFractional);
  std::optional<RunningEnergy> withoutGrids = Interactions::computeFrameworkMoleculeEnergyDifference(
      system.forceField, system.simulationBox, noGrids, system.framework, system.frameworkCellList(),
      system.spanOfFrameworkAtoms(), switchedOff, oldFractional);
  ASSERT_TRUE(withGrids.has_value());
  ASSERT_TRUE(withoutGrids.has_value());
  EXPECT_NEAR(withGrids->frameworkMoleculeVDW, withoutGrids->frameworkMoleculeVDW, 1e-10);
  EXPECT_NEAR(withGrids->frameworkMoleculeCharge, withoutGrids->frameworkMoleculeCharge, 1e-10);

  // the switched-off molecule has no framework energy: the difference is minus the fractional energy
  RunningEnergy fractionalEnergy = Interactions::computeFrameworkMoleculeEnergy(
      system.forceField, system.simulationBox, noGrids, system.framework, system.spanOfFrameworkAtoms(),
      oldFractional);
  EXPECT_NEAR(withGrids->frameworkMoleculeVDW, -fractionalEnergy.frameworkMoleculeVDW, 1e-10);
  EXPECT_NEAR(withGrids->frameworkMoleculeCharge, -fractionalEnergy.frameworkMoleculeCharge, 1e-10);

  // the brute-force variant of the difference agrees as well
  std::optional<RunningEnergy> bruteForce = Interactions::computeFrameworkMoleculeEnergyDifference(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
      system.spanOfFrameworkAtoms(), switchedOff, oldFractional);
  ASSERT_TRUE(bruteForce.has_value());
  EXPECT_NEAR(bruteForce->frameworkMoleculeVDW, -fractionalEnergy.frameworkMoleculeVDW, 1e-10);
}
