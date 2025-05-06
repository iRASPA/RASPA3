#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <cstddef>
#include <span>
#include <vector>

import int3;
import double3;
import factory;
import units;
import atom;
import pseudo_atom;
import vdwparameters;
import forcefield;
import framework;
import component;
import system;
import simulationbox;
import energy_factor;
import running_energy;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;

TEST(static_energy, Test_2_CO2_in_ITQ_29_1x1x1)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Framework f = TestFactories::makeITQ29(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeCO2(forceField, 0, true);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {2}, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  atomPositions[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  atomPositions[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  atomPositions[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  atomPositions[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  atomPositions[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  atomPositions[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);

  RunningEnergy energy =
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions) +
      Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                   system.framework, frameworkAtoms, atomPositions);

  EXPECT_NEAR(energy.frameworkMoleculeVDW * Units::EnergyToKelvin, -1545.62921755, 1e-6);
  EXPECT_NEAR(energy.frameworkMoleculeCharge * Units::EnergyToKelvin, -592.13188606, 1e-6);
  EXPECT_NEAR(energy.moleculeMoleculeVDW * Units::EnergyToKelvin, -242.42459776, 1e-6);
  EXPECT_NEAR(energy.moleculeMoleculeCharge * Units::EnergyToKelvin, 154.11883595, 1e-6);
}

TEST(static_energy, Test_2_CO2_in_MFI_2x2x2_shifted)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Component c = TestFactories::makeCO2(forceField, 0, true);
  Framework f = TestFactories::makeMFI_Si(forceField, int3(2, 2, 2));

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {2}, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  atomPositions[0].position = double3(10.011, 4.97475 + 2.0, 1.149);
  atomPositions[1].position = double3(10.011, 4.97475 + 2.0, 0.0);
  atomPositions[2].position = double3(10.011, 4.97475 + 2.0, -1.149);
  atomPositions[3].position = double3(10.011, 4.97475 - 2.0, 1.149);
  atomPositions[4].position = double3(10.011, 4.97475 - 2.0, 0.0);
  atomPositions[5].position = double3(10.011, 4.97475 - 2.0, -1.149);

  RunningEnergy energy =
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions) +
      Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                   system.framework, frameworkAtoms, atomPositions);

  EXPECT_NEAR(energy.frameworkMoleculeVDW * Units::EnergyToKelvin, -2525.36580663, 1e-6);
  EXPECT_NEAR(energy.frameworkMoleculeCharge * Units::EnergyToKelvin, 2167.45591472, 1e-6);
  EXPECT_NEAR(energy.moleculeMoleculeVDW * Units::EnergyToKelvin, -242.42459776, 1e-6);
  EXPECT_NEAR(energy.moleculeMoleculeCharge * Units::EnergyToKelvin, 154.11883595, 1e-6);
}

TEST(static_energy, Test_2_CO2_in_MFI_2x2x2_truncated)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, false, true, true);

  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);
  Framework f = TestFactories::makeMFI_Si(forceField, int3(2, 2, 2));
  Component c = TestFactories::makeCO2(forceField, 0, true);
  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {2}, 5);

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
