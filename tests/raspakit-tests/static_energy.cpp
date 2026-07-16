#include <gtest/gtest.h>

import std;

import int3;
import double3;
import units;
import atom;
import pseudo_atom;
import vdwparameters;
import forcefield;
import framework;
import component;
import system;
import simulationbox;
import running_energy;
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
