#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <cstddef>
#include <span>
#include <vector>

import int3;
import double3;
import double3x3;
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
import gradient_factor;
import energy_status;
import units;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;

TEST(Ewald, Test_2_CO2_in_Box_10_10_10)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);

  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

  Component c = TestFactories::makeCO2(forceField, 0, true);

  System system = System(0, forceField, SimulationBox(10.0, 10.0, 10.0), 300.0, 1e4, 1.0, {}, {c}, {2}, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  atomPositions[0].position = double3(-1.0, 0.0, 1.149);
  atomPositions[1].position = double3(-1.0, 0.0, 0.0);
  atomPositions[2].position = double3(-1.0, 0.0, -1.149);
  atomPositions[3].position = double3(1.0, 0.0, 1.149);
  atomPositions[4].position = double3(1.0, 0.0, 0.0);
  atomPositions[5].position = double3(1.0, 0.0, -1.149);

  Interactions::computeEwaldFourierEnergySingleIon(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                   system.forceField, system.simulationBox, double3(0.0, 0.0, 0.0),
                                                   1.0);
  system.precomputeTotalRigidEnergy();
  RunningEnergy energy = Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
      system.spanOfMoleculeAtoms());

  EXPECT_NEAR(energy.ewaldFourier() * Units::EnergyToKelvin, 90.54613836, 1e-6);
}

TEST(Ewald, Test_1_Na_1_Cl_in_Box_10_10_10_Gradient)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);

  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

  Component na = TestFactories::makeIon(forceField, 0, "Na", 6, 0.0);
  Component cl = TestFactories::makeIon(forceField, 1, "Cl", 7, 0.0);

  System system = System(0, forceField, SimulationBox(10.0, 10.0, 10.0), 300.0, 1e4, 1.0, {}, {na, cl}, {1, 1}, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::vector<Atom> atomPositions = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end());

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  RunningEnergy energy;

  system.precomputeTotalRigidEnergy();
  std::pair<EnergyStatus, double3x3> strainDerivative = Interactions::computeEwaldFourierEnergyStrainDerivative(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, system.framework, system.components,
      system.numberOfMoleculesPerComponent, system.spanOfMoleculeAtoms(), system.CoulombicFourierEnergySingleIon,
      system.netChargeFramework, system.netChargePerComponent);

  double delta = 1e-5;
  double tolerance = 1e-5;
  double3 gradient;
  for (size_t i = 0; i < atomPositions.size(); ++i)
  {
    RunningEnergy x1, x2, y1, y2, z1, z2;

    // finite difference x
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x + 0.5 * delta;
    x2 = Interactions::computeEwaldFourierEnergy(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                 system.fixedFrameworkStoredEik, system.storedEik, system.forceField,
                                                 system.simulationBox, system.components,
                                                 system.numberOfMoleculesPerComponent, atomPositions);

    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x - 0.5 * delta;
    x1 = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        std::span<const Atom>(atomPositions));
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x;

    // finite difference y
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y + 0.5 * delta;
    y2 = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        std::span<const Atom>(atomPositions));

    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y - 0.5 * delta;
    y1 = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        std::span<const Atom>(atomPositions));
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y;

    // finite difference z
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z + 0.5 * delta;
    z2 = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        std::span<const Atom>(atomPositions));

    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z - 0.5 * delta;
    z1 = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        std::span<const Atom>(atomPositions));
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z;

    gradient.x = (x2.ewaldFourier() - x1.ewaldFourier()) / delta;
    gradient.y = (y2.ewaldFourier() - y1.ewaldFourier()) / delta;
    gradient.z = (z2.ewaldFourier() - z1.ewaldFourier()) / delta;

    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.x, gradient.x, tolerance) << "Wrong x-gradient";
    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.y, gradient.y, tolerance) << "Wrong y-gradient";
    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.z, gradient.z, tolerance) << "Wrong z-gradient";
  }
}

TEST(Ewald, Test_2_CO2_in_ITQ_29_1x1x1)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);

  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

  Framework f = TestFactories::makeITQ29(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeCO2(forceField, 0, true);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {2}, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  atomPositions[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  atomPositions[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  atomPositions[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  atomPositions[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  atomPositions[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  atomPositions[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);

  system.precomputeTotalRigidEnergy();
  system.CoulombicFourierEnergySingleIon = Interactions::computeEwaldFourierEnergySingleIon(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.forceField, system.simulationBox,
      double3(0.0, 0.0, 0.0), 1.0);
  RunningEnergy energy = Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
      system.spanOfMoleculeAtoms());

  EXPECT_NEAR(system.CoulombicFourierEnergySingleIon * Units::EnergyToKelvin, 17464.2371790130, 1e-6);
  EXPECT_NEAR(energy.ewaldFourier() * Units::EnergyToKelvin, -702.65863478, 1e-6);
}

TEST(Ewald, Test_2_CO2_in_ITQ_29_2x2x2)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);

  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

  Framework f = TestFactories::makeITQ29(forceField, int3(2, 2, 2));
  Component c = TestFactories::makeCO2(forceField, 0, true);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {2}, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  atomPositions[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  atomPositions[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  atomPositions[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  atomPositions[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  atomPositions[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  atomPositions[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);

  system.precomputeTotalRigidEnergy();
  system.CoulombicFourierEnergySingleIon = Interactions::computeEwaldFourierEnergySingleIon(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.forceField, system.simulationBox,
      double3(0.0, 0.0, 0.0), 1.0);
  RunningEnergy energy = Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
      system.spanOfMoleculeAtoms());

  EXPECT_NEAR(system.CoulombicFourierEnergySingleIon * Units::EnergyToKelvin, 9673.9032373025, 1e-6);
  EXPECT_NEAR(energy.ewaldFourier() * Units::EnergyToKelvin, -721.64644486, 1e-6);
}

TEST(Ewald, Test_2_CO2_in_MFI_1x1x1)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);

  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

  Framework f = TestFactories::makeMFI_Si(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeCO2(forceField, 0, true);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {2}, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  atomPositions[0].position = double3(10.011, 4.97475 + 2.0, 1.149);
  atomPositions[1].position = double3(10.011, 4.97475 + 2.0, 0.0);
  atomPositions[2].position = double3(10.011, 4.97475 + 2.0, -1.149);
  atomPositions[3].position = double3(10.011, 4.97475 - 2.0, 1.149);
  atomPositions[4].position = double3(10.011, 4.97475 - 2.0, 0.0);
  atomPositions[5].position = double3(10.011, 4.97475 - 2.0, -1.149);

  system.precomputeTotalRigidEnergy();
  system.CoulombicFourierEnergySingleIon = Interactions::computeEwaldFourierEnergySingleIon(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.forceField, system.simulationBox,
      double3(0.0, 0.0, 0.0), 1.0);
  RunningEnergy energy = Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
      system.spanOfMoleculeAtoms());

  EXPECT_NEAR(system.CoulombicFourierEnergySingleIon * Units::EnergyToKelvin, 12028.1731827280, 1e-6);
  EXPECT_NEAR(energy.ewaldFourier() * Units::EnergyToKelvin, -1191.77790165, 1e-6);
}

TEST(Ewald, Test_2_CO2_in_MFI_2x2x2)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);

  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

  Framework f = TestFactories::makeMFI_Si(forceField, int3(2, 2, 2));
  Component c = TestFactories::makeCO2(forceField, 0, true);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {2}, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  atomPositions[0].position = double3(10.011, 4.97475 + 2.0, 1.149);
  atomPositions[1].position = double3(10.011, 4.97475 + 2.0, 0.0);
  atomPositions[2].position = double3(10.011, 4.97475 + 2.0, -1.149);
  atomPositions[3].position = double3(10.011, 4.97475 - 2.0, 1.149);
  atomPositions[4].position = double3(10.011, 4.97475 - 2.0, 0.0);
  atomPositions[5].position = double3(10.011, 4.97475 - 2.0, -1.149);

  system.precomputeTotalRigidEnergy();
  system.CoulombicFourierEnergySingleIon = Interactions::computeEwaldFourierEnergySingleIon(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.forceField, system.simulationBox,
      double3(0.0, 0.0, 0.0), 1.0);
  RunningEnergy energy = Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
      system.spanOfMoleculeAtoms());

  EXPECT_NEAR(system.CoulombicFourierEnergySingleIon * Units::EnergyToKelvin, 6309.7866899037, 1e-6);
  EXPECT_NEAR(energy.ewaldFourier() * Units::EnergyToKelvin, -1197.23909965, 1e-6);
}

TEST(Ewald, Test_20_Na_Cl_in_Box_25x25x25)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);

  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

  Component na = TestFactories::makeIon(forceField, 0, "Na", 6, 0.0);
  Component cl = TestFactories::makeIon(forceField, 1, "Cl", 7, 0.0);

  System system = System(0, forceField, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, 1.0, {}, {na, cl}, {20, 20}, 5);

  // std::fill(system.forceField.data.begin(), system.forceField.data.end(), VDWParameters(0.0, 1.0));

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::vector<Atom> atomPositions = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end());

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  for (size_t i = 0; i < 20; ++i)
  {
    atomPositions[i].charge = 1.0;
    system.atomPositions[i].charge = 1.0;
  }
  for (size_t i = 0; i < 20; ++i)
  {
    atomPositions[i + 20].charge = -1.0;
    system.atomPositions[i + 20].charge = -1.0;
  }

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  RunningEnergy energy, rigidenergy;

  system.precomputeTotalRigidEnergy();
  [[maybe_unused]] RunningEnergy factor = Interactions::computeEwaldFourierGradient(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.fixedFrameworkStoredEik,
      system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent, atomPositions);

  double delta = 1e-4;
  double tolerance = 1e-4;
  double3 gradient;
  for (size_t i = 0; i < atomPositions.size(); ++i)
  {
    RunningEnergy forward1_x, forward2_x, backward1_x, backward2_x;
    RunningEnergy forward1_y, forward2_y, backward1_y, backward2_y;
    RunningEnergy forward1_z, forward2_z, backward1_z, backward2_z;

    // finite difference x
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x + 1.0 * delta;
    forward2_x = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        system.spanOfMoleculeAtoms());
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x + 0.5 * delta;
    forward1_x = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        system.spanOfMoleculeAtoms());
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x - 0.5 * delta;
    backward1_x = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        system.spanOfMoleculeAtoms());
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x - 1.0 * delta;
    backward2_x = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        system.spanOfMoleculeAtoms());
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x;

    // finite difference y
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y + 1.0 * delta;
    forward2_y = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        system.spanOfMoleculeAtoms());
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y + 0.5 * delta;
    forward1_y = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        system.spanOfMoleculeAtoms());
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y - 0.5 * delta;
    backward1_y = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        system.spanOfMoleculeAtoms());
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y - 1.0 * delta;
    backward2_y = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        system.spanOfMoleculeAtoms());
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y;

    // finite difference z
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z + 1.0 * delta;
    forward2_z = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        system.spanOfMoleculeAtoms());
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z + 0.5 * delta;
    forward1_z = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        system.spanOfMoleculeAtoms());
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z - 0.5 * delta;
    backward1_z = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        system.spanOfMoleculeAtoms());
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z - 1.0 * delta;
    backward2_z = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        system.spanOfMoleculeAtoms());
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z;

    gradient.x = (-forward2_x.ewaldFourier() + 8.0 * forward1_x.ewaldFourier() - 8.0 * backward1_x.ewaldFourier() +
                  backward2_x.ewaldFourier()) /
                 (6.0 * delta);
    gradient.y = (-forward2_y.ewaldFourier() + 8.0 * forward1_y.ewaldFourier() - 8.0 * backward1_y.ewaldFourier() +
                  backward2_y.ewaldFourier()) /
                 (6.0 * delta);
    gradient.z = (-forward2_z.ewaldFourier() + 8.0 * forward1_z.ewaldFourier() - 8.0 * backward1_z.ewaldFourier() +
                  backward2_z.ewaldFourier()) /
                 (6.0 * delta);

    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.x, gradient.x, tolerance) << "Wrong x-gradient";
    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.y, gradient.y, tolerance) << "Wrong y-gradient";
    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.z, gradient.z, tolerance) << "Wrong z-gradient";
  }
}
