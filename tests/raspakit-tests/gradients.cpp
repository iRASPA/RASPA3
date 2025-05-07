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
import gradient_factor;
import running_energy;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;
import energy_status;
import integrators;
import integrators_compute;
import integrators_update;
import interpolation_energy_grid;

TEST(gradients, Test_2_CO2_in_ITQ_29_2x2x2_inter)
{
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);
  Framework f = TestFactories::makeITQ29(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeCO2(forceField, 0, false);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {2}, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  atomPositions[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  atomPositions[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  atomPositions[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  atomPositions[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  atomPositions[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  atomPositions[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);

  RunningEnergy factorInterMolecular = Interactions::computeInterMolecularGradient(
      system.forceField, system.simulationBox, system.spanOfMoleculeAtoms());

  EXPECT_NEAR(factorInterMolecular.moleculeMoleculeVDW * Units::EnergyToKelvin, -242.36960932, 1e-6);
  EXPECT_NEAR(factorInterMolecular.moleculeMoleculeCharge * Units::EnergyToKelvin, 0.00000000, 1e-6);

  EXPECT_NEAR(atomPositions[0].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[0].gradient.y, 90.951955774481, 1e-6);
  EXPECT_NEAR(atomPositions[0].gradient.z, 17.938271420558, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.y, 52.398564956433, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.z, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.y, 90.951955774481, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.z, -17.938271420558, 1e-6);

  EXPECT_NEAR(atomPositions[3].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[3].gradient.y, -90.951955774481, 1e-6);
  EXPECT_NEAR(atomPositions[3].gradient.z, 17.938271420558, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.y, -52.398564956433, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.z, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.y, -90.951955774481, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.z, -17.938271420558, 1e-6);

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  std::pair<EnergyStatus, double3x3> pressureInfo = Interactions::computeInterMolecularEnergyStrainDerivative(
      system.forceField, system.components, system.simulationBox, system.spanOfMoleculeAtoms());

  EXPECT_NEAR(atomPositions[0].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[0].gradient.y, 90.951955774481, 1e-6);
  EXPECT_NEAR(atomPositions[0].gradient.z, 17.938271420558, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.y, 52.398564956433, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.z, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.y, 90.951955774481, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.z, -17.938271420558, 1e-6);

  EXPECT_NEAR(atomPositions[3].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[3].gradient.y, -90.951955774481, 1e-6);
  EXPECT_NEAR(atomPositions[3].gradient.z, 17.938271420558, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.y, -52.398564956433, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.z, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.y, -90.951955774481, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.z, -17.938271420558, 1e-6);
}

TEST(gradients, Test_2_CO2_in_ITQ_29_2x2x2_framework_molecule)
{
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);
  Framework f = TestFactories::makeITQ29(forceField, int3(2, 2, 2));
  Component c = TestFactories::makeCO2(forceField, 0, false);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {2}, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  atomPositions[0].position = double3(5.93355, 7.93355, 7.08255);
  atomPositions[1].position = double3(5.93355, 7.93355, 5.93355);
  atomPositions[2].position = double3(5.93355, 7.93355, 4.78455);
  atomPositions[3].position = double3(5.93355, 3.93355, 7.08255);
  atomPositions[4].position = double3(5.93355, 3.93355, 5.93355);
  atomPositions[5].position = double3(5.93355, 3.93355, 4.78455);

  RunningEnergy factorFrameworkMolecular = Interactions::computeFrameworkMoleculeGradient(
      system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
      system.interpolationGrids);

  EXPECT_NEAR(factorFrameworkMolecular.frameworkMoleculeVDW * Units::EnergyToKelvin, -1932.15586114, 1e-6);
  EXPECT_NEAR(factorFrameworkMolecular.frameworkMoleculeCharge * Units::EnergyToKelvin, 0.00000000, 1e-6);

  EXPECT_NEAR(atomPositions[0].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[0].gradient.y, -131.516544539514, 1e-6);
  EXPECT_NEAR(atomPositions[0].gradient.z, -90.888952502176, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.y, -48.847937726172, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.z, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.y, -131.516544539514, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.z, 90.888952502176, 1e-6);

  EXPECT_NEAR(atomPositions[3].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[3].gradient.y, 131.516544539514, 1e-6);
  EXPECT_NEAR(atomPositions[3].gradient.z, -90.888952502176, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.y, 48.847937726172, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.z, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.y, 131.516544539514, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.z, 90.888952502176, 1e-6);

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  std::pair<EnergyStatus, double3x3> pressureInfo = Interactions::computeFrameworkMoleculeEnergyStrainDerivative(
      system.forceField, system.framework, system.interpolationGrids,
      system.components, system.simulationBox, system.spanOfFrameworkAtoms(),
      system.spanOfMoleculeAtoms());

  EXPECT_NEAR(atomPositions[0].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[0].gradient.y, -131.516544539514, 1e-6);
  EXPECT_NEAR(atomPositions[0].gradient.z, -90.888952502176, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.y, -48.847937726172, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.z, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.y, -131.516544539514, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.z, 90.888952502176, 1e-6);

  EXPECT_NEAR(atomPositions[3].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[3].gradient.y, 131.516544539514, 1e-6);
  EXPECT_NEAR(atomPositions[3].gradient.z, -90.888952502176, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.y, 48.847937726172, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.z, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.y, 131.516544539514, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.z, 90.888952502176, 1e-6);
}

TEST(gradients, Test_2_CO2_in_ITQ_29_2x2x2_NonEwald)
{
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);

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

  RunningEnergy factorFrameworkMolecular = Interactions::computeFrameworkMoleculeGradient(
      system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
      system.interpolationGrids);
  RunningEnergy factorInterMolecular = Interactions::computeInterMolecularGradient(
      system.forceField, system.simulationBox, system.spanOfMoleculeAtoms());

  EXPECT_NEAR(factorFrameworkMolecular.frameworkMoleculeVDW * Units::EnergyToKelvin, -1932.15586114, 1e-6);
  EXPECT_NEAR(factorFrameworkMolecular.frameworkMoleculeCharge * Units::EnergyToKelvin, 554.41444763, 1e-6);
  EXPECT_NEAR(factorInterMolecular.moleculeMoleculeVDW * Units::EnergyToKelvin, -242.36960932, 1e-6);
  EXPECT_NEAR(factorInterMolecular.moleculeMoleculeCharge * Units::EnergyToKelvin, 162.41877650, 1e-6);

  EXPECT_NEAR(atomPositions[0].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[0].gradient.y, 466.705848885190, 1e-6);
  EXPECT_NEAR(atomPositions[0].gradient.z, 223.065047349758, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.y, -1259.011389396075, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.z, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.y, 466.705848885189, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.z, -223.065047349755, 1e-6);

  EXPECT_NEAR(atomPositions[3].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[3].gradient.y, -466.705848885190, 1e-6);
  EXPECT_NEAR(atomPositions[3].gradient.z, 223.065047349758, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.y, 1259.011389396075, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.z, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.y, -466.705848885189, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.z, -223.065047349755, 1e-6);

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  std::pair<EnergyStatus, double3x3> pressureInfo1 = Interactions::computeFrameworkMoleculeEnergyStrainDerivative(
      system.forceField, system.framework, system.interpolationGrids,
      system.components, system.simulationBox, system.spanOfFrameworkAtoms(),
      system.spanOfMoleculeAtoms());

  std::pair<EnergyStatus, double3x3> pressureInfo2 = Interactions::computeInterMolecularEnergyStrainDerivative(
      system.forceField, system.components, system.simulationBox, system.spanOfMoleculeAtoms());

  EXPECT_NEAR(atomPositions[0].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[0].gradient.y, 466.705848885190, 1e-6);
  EXPECT_NEAR(atomPositions[0].gradient.z, 223.065047349758, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.y, -1259.011389396075, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.z, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.y, 466.705848885189, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.z, -223.065047349755, 1e-6);

  EXPECT_NEAR(atomPositions[3].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[3].gradient.y, -466.705848885190, 1e-6);
  EXPECT_NEAR(atomPositions[3].gradient.z, 223.065047349758, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.y, 1259.011389396075, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.z, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.x, 0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.y, -466.705848885189, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.z, -223.065047349755, 1e-6);
}

TEST(gradients, Test_2_CO2_in_ITQ_29_2x2x2_Ewald)
{
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);

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
  RunningEnergy factorEwald = Interactions::computeEwaldFourierGradient(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.fixedFrameworkStoredEik,
      system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
      system.spanOfMoleculeAtoms());

  EXPECT_NEAR(
      (factorEwald.ewald_fourier + factorEwald.ewald_self + factorEwald.ewald_exclusion) * Units::EnergyToKelvin,
      -759.67572774 + 38.02930863, 1e-4);

  EXPECT_NEAR(atomPositions[0].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[0].gradient.y, -362.766142495638, 1e-4);
  EXPECT_NEAR(atomPositions[0].gradient.z, -241.771707768611, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.y, 684.800882899876, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.z, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.y, -362.766142495638, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.z, 241.771707768619, 1e-4);

  EXPECT_NEAR(atomPositions[3].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[3].gradient.y, 362.766142495638, 1e-4);
  EXPECT_NEAR(atomPositions[3].gradient.z, -241.771707768611, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.y, -684.800882899876, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.z, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.y, 362.766142495638, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.z, 241.771707768619, 1e-4);

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  std::pair<EnergyStatus, double3x3> pressureInfo = Interactions::computeEwaldFourierEnergyStrainDerivative(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, system.framework, system.components,
      system.numberOfMoleculesPerComponent, system.spanOfMoleculeAtoms(), system.CoulombicFourierEnergySingleIon,
      system.netChargeFramework, system.netChargePerComponent);

  EXPECT_NEAR(atomPositions[0].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[0].gradient.y, -362.766142495638, 1e-4);
  EXPECT_NEAR(atomPositions[0].gradient.z, -241.771707768611, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.y, 684.800882899876, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.z, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.y, -362.766142495638, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.z, 241.771707768619, 1e-4);

  EXPECT_NEAR(atomPositions[3].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[3].gradient.y, 362.766142495638, 1e-4);
  EXPECT_NEAR(atomPositions[3].gradient.z, -241.771707768611, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.y, -684.800882899876, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.z, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.y, 362.766142495638, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.z, 241.771707768619, 1e-4);
}

TEST(gradients, Test_2_CO2_in_ITQ_29_2x2x2_Total)
{
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);

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

  [[maybe_unused]] RunningEnergy factorFrameworkMolecular = Interactions::computeFrameworkMoleculeGradient(
      system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
      system.interpolationGrids);
  [[maybe_unused]] RunningEnergy factorInterMolecular = Interactions::computeInterMolecularGradient(
      system.forceField, system.simulationBox, system.spanOfMoleculeAtoms());

  system.precomputeTotalRigidEnergy();
  [[maybe_unused]] RunningEnergy factorEwald = Interactions::computeEwaldFourierGradient(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.fixedFrameworkStoredEik,
      system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
      system.spanOfMoleculeAtoms());

  EXPECT_NEAR(atomPositions[0].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[0].gradient.y, 103.939706389550, 1e-4);
  EXPECT_NEAR(atomPositions[0].gradient.z, -18.706660418853, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.y, -574.210506496196, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.z, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.y, 103.939706389549, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.z, 18.706660418865, 1e-4);

  EXPECT_NEAR(atomPositions[3].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[3].gradient.y, -103.939706389550, 1e-4);
  EXPECT_NEAR(atomPositions[3].gradient.z, -18.706660418853, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.y, 574.210506496196, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.z, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.y, -103.939706389549, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.z, 18.706660418865, 1e-4);

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }
  [[maybe_unused]] RunningEnergy gradientEnergy = Integrators::updateGradients(
      system.spanOfMoleculeAtoms(), system.spanOfFrameworkAtoms(), system.forceField, system.simulationBox,
      system.components, system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik,
      system.fixedFrameworkStoredEik, system.interpolationGrids, system.numberOfMoleculesPerComponent);

  // EXPECT_NEAR(gradientEnergy.total()  * Units::EnergyToKelvin, -2179.338665434245, 1e-4);

  EXPECT_NEAR(atomPositions[0].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[0].gradient.y, 103.939706389550, 1e-4);
  EXPECT_NEAR(atomPositions[0].gradient.z, -18.706660418853, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.y, -574.210506496196, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.z, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.y, 103.939706389549, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.z, 18.706660418865, 1e-4);

  EXPECT_NEAR(atomPositions[3].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[3].gradient.y, -103.939706389550, 1e-4);
  EXPECT_NEAR(atomPositions[3].gradient.z, -18.706660418853, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.y, 574.210506496196, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.z, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.x, 0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.y, -103.939706389549, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.z, 18.706660418865, 1e-4);
}

/*
TEST(gradients, Test_2_CO2_in_ITQ_29_1x1x1_numerical)
{
  ForceField forceField = ForceField(
    { PseudoAtom("Si",    28.0855,   2.05,   0.0, 14, false),
      PseudoAtom("O",     15.999,   -1.025,  0.0,  8, false),
      PseudoAtom("CH4",   16.04246,  0.0,    0.0,  6, false),
      PseudoAtom("C_co2", 12.0,      0.6512, 0.0,  6, false),
      PseudoAtom("O_co2", 15.9994,  -0.3256, 0.0,  8, false),
    },
    { VDWParameters(22.0, 2.30),
      VDWParameters(53.0, 3.3),
      VDWParameters(158.5, 3.72),
      VDWParameters(29.933, 2.745),
      VDWParameters(85.671, 3.017)
    },
    ForceField::MixingRule::Lorentz_Berthelot,
    12.0,
    true,
    false);
  Framework f = Framework(0, forceField, "ITQ-29", SimulationBox(11.8671, 11.8671, 11.8671),
    517,
    { // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId,
uint8_t groupId Atom(double3(0.3683, 0.1847, 0),       2.05,  1.0, 0, 0, 0, 0), Atom(double3(0.5,    0.2179, 0),
-1.025, 1.0, 0, 1, 0, 0), Atom(double3(0.2939, 0.2939, 0),      -1.025, 1.0, 0, 1, 0, 0), Atom(double3(0.3429, 0.1098,
0.1098), -1.025, 1.0, 0, 1, 0, 0)
    },
    int3(1, 1, 1));
  Component c = Component(0, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
    { // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId,
uint8_t groupId
      //Atom(double3(0.0, 0.0,  1.149), 0.0, 1.0, 0, 4, 1, 0),
      //Atom(double3(0.0, 0.0,  0.0  ), 0.0, 1.0, 0, 3, 1, 0),
      //Atom(double3(0.0, 0.0, -1.149), 0.0, 1.0, 0, 4, 1, 0)
      Atom(double3(0.0, 0.0,  1.149), -0.3256, 1.0, 0, 4, 1, 0),
      Atom(double3(0.0, 0.0,  0.0  ),  0.6512, 1.0, 0, 3, 1, 0),
      Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)
    }, 5, 21);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, { f }, { c }, { 2 }, 5);

  //std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  //std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  //spanOfMoleculeAtoms[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  //spanOfMoleculeAtoms[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  //spanOfMoleculeAtoms[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  //spanOfMoleculeAtoms[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  //spanOfMoleculeAtoms[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  //spanOfMoleculeAtoms[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);

  //[[maybe_unused]] ForceFactor factorFrameworkMolecular =
  //  Interactions::computeFrameworkMoleculeGradient(system.forceField, system.simulationBox,
system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms());
  //[[maybe_unused]] ForceFactor factorInterMolecular =
  //  Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox,
system.spanOfMoleculeAtoms());

  //// make copy
  //std::vector<Atom> atomPositions = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end());

  //for (Atom& atom : atomPositions)
  //{
  //  atom.gradient = double3(0.0, 0.0, 0.0);
  //}

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  spanOfMoleculeAtoms[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);
  std::vector<Atom> atomPositions = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end());

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  //[[maybe_unused]] ForceFactor factorFrameworkMolecular =
  //  Interactions::computeFrameworkMoleculeGradient(system.forceField, system.simulationBox,
system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms());
  [[maybe_unused]] ForceFactor factorInterMolecular =
    Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox, system.spanOfMoleculeAtoms());


  double delta = 1e-3;
  double tolerance = 1e-5;
  double3 gradient;
  for (size_t i = 0; i < atomPositions.size(); ++i)
  {
    RunningEnergy x1, x2, y1, y2, z1, z2;

    // finite difference x
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x + 0.5 * delta;
    x2 = Interactions::computeFrameworkMoleculeEnergy(system.forceField,system.simulationBox, frameworkAtoms,
atomPositions) + Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);

    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x - 0.5 * delta;
    x1 = Interactions::computeFrameworkMoleculeEnergy(system.forceField,system.simulationBox, frameworkAtoms,
atomPositions) + Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x;

    // finite difference y
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y + 0.5 * delta;
    y2 = Interactions::computeFrameworkMoleculeEnergy(system.forceField,system.simulationBox, frameworkAtoms,
atomPositions) + Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);

    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y - 0.5 * delta;
    y1 = Interactions::computeFrameworkMoleculeEnergy(system.forceField,system.simulationBox, frameworkAtoms,
atomPositions) + Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y;

    // finite difference z
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z + 0.5 * delta;
    z2 = Interactions::computeFrameworkMoleculeEnergy(system.forceField,system.simulationBox, frameworkAtoms,
atomPositions) + Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);

    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z - 0.5 * delta;
    z1 = Interactions::computeFrameworkMoleculeEnergy(system.forceField,system.simulationBox, frameworkAtoms,
atomPositions) + Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z;

    gradient.x = (x2.moleculeMoleculeVDW + x2.moleculeMoleculeCharge + x2.frameworkMoleculeVDW  +
x2.frameworkMoleculeCharge
                - (x1.moleculeMoleculeVDW + x1.moleculeMoleculeCharge + x1.frameworkMoleculeVDW +
x1.frameworkMoleculeCharge)) / delta; gradient.y = (y2.moleculeMoleculeVDW + y2.moleculeMoleculeCharge +
y2.frameworkMoleculeVDW + y2.frameworkMoleculeCharge
                - (y1.moleculeMoleculeVDW + y1.moleculeMoleculeCharge + y1.frameworkMoleculeVDW +
y1.frameworkMoleculeCharge)) / delta; gradient.z = (z2.moleculeMoleculeVDW + z2.moleculeMoleculeCharge +
z2.frameworkMoleculeVDW + z2.frameworkMoleculeCharge
                - (z1.moleculeMoleculeVDW + z1.moleculeMoleculeCharge + z1.frameworkMoleculeVDW +
z1.frameworkMoleculeCharge)) / delta;

    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.x, gradient.x, tolerance);
    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.y, gradient.y, tolerance);
    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.z, gradient.z, tolerance);
  }
}
*/

TEST(gradients, Test_CO2_in_ITQ_29_1x1x1)
{
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);

  Framework f = TestFactories::makeITQ29(forceField, int3(2, 2, 2));
  Component c = TestFactories::makeCO2(forceField, 0, true);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {2}, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();

  spanOfMoleculeAtoms[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  std::vector<Atom> atomPositions = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end());

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  [[maybe_unused]] RunningEnergy factor = Interactions::computeInterMolecularGradient(
      system.forceField, system.simulationBox, system.spanOfMoleculeAtoms());
  [[maybe_unused]] RunningEnergy factor2 = Interactions::computeFrameworkMoleculeGradient(
      system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
      system.interpolationGrids);

  double delta = 1e-5;
  double tolerance = 1e-4;
  double3 gradient;
  for (size_t i = 0; i < atomPositions.size(); ++i)
  {
    RunningEnergy x1, x2, y1, y2, z1, z2;

    // finite difference x
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x + 0.5 * delta;
    x2 =
        Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions) +
        Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                     system.framework, frameworkAtoms, atomPositions);

    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x - 0.5 * delta;
    x1 =
        Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions) +
        Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                     system.framework, frameworkAtoms, atomPositions);
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x;

    // finite difference y
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y + 0.5 * delta;
    y2 =
        Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions) +
        Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                     system.framework, frameworkAtoms, atomPositions);

    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y - 0.5 * delta;
    y1 =
        Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions) +
        Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                     system.framework, frameworkAtoms, atomPositions);

    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y;

    // finite difference z
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z + 0.5 * delta;
    z2 =
        Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions) +
        Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                     system.framework, frameworkAtoms, atomPositions);

    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z - 0.5 * delta;
    z1 =
        Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions) +
        Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                     system.framework, frameworkAtoms, atomPositions);

    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z;

    gradient.x =
        (x2.moleculeMoleculeVDW + x2.moleculeMoleculeCharge + x2.frameworkMoleculeVDW + x2.frameworkMoleculeCharge -
         (x1.moleculeMoleculeVDW + x1.moleculeMoleculeCharge + x1.frameworkMoleculeVDW + x1.frameworkMoleculeCharge)) /
        delta;
    gradient.y =
        (y2.moleculeMoleculeVDW + y2.moleculeMoleculeCharge + y2.frameworkMoleculeVDW + y2.frameworkMoleculeCharge -
         (y1.moleculeMoleculeVDW + y1.moleculeMoleculeCharge + y1.frameworkMoleculeVDW + y1.frameworkMoleculeCharge)) /
        delta;
    gradient.z =
        (z2.moleculeMoleculeVDW + z2.moleculeMoleculeCharge + z2.frameworkMoleculeVDW + z2.frameworkMoleculeCharge -
         (z1.moleculeMoleculeVDW + z1.moleculeMoleculeCharge + z1.frameworkMoleculeVDW + z1.frameworkMoleculeCharge)) /
        delta;

    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.x, gradient.x, tolerance);
    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.y, gradient.y, tolerance);
    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.z, gradient.z, tolerance);
  }
}

TEST(gradients, Test_CH4_in_Box_25x25x25)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, false, false, false);
  Component c = TestFactories::makeMethane(forceField, 0);

  System system = System(0, forceField, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, 1.0, {}, {c}, {50}, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  std::vector<Atom> atomPositions = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end());

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  [[maybe_unused]] RunningEnergy factorFrameworkMolecular = Interactions::computeFrameworkMoleculeGradient(
      system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
      system.interpolationGrids);
  [[maybe_unused]] RunningEnergy factorInterMolecular = Interactions::computeInterMolecularGradient(
      system.forceField, system.simulationBox, system.spanOfMoleculeAtoms());

  double delta = 1e-5;
  double tolerance = 1e-6;
  double3 gradient;
  for (size_t i = 0; i < atomPositions.size(); ++i)
  {
    RunningEnergy x1, x2, y1, y2, z1, z2;

    // finite difference x
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x + 0.5 * delta;
    x2 =
        Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                     system.framework, frameworkAtoms, atomPositions) +
        Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);

    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x - 0.5 * delta;
    x1 =
        Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                     system.framework, frameworkAtoms, atomPositions) +
        Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x;

    // finite difference y
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y + 0.5 * delta;
    y2 =
        Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                     system.framework, frameworkAtoms, atomPositions) +
        Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);

    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y - 0.5 * delta;
    y1 =
        Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                     system.framework, frameworkAtoms, atomPositions) +
        Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y;

    // finite difference z
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z + 0.5 * delta;
    z2 =
        Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                     system.framework, frameworkAtoms, atomPositions) +
        Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);

    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z - 0.5 * delta;
    z1 =
        Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                     system.framework, frameworkAtoms, atomPositions) +
        Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z;

    gradient.x =
        (x2.moleculeMoleculeVDW + x2.moleculeMoleculeCharge + x2.frameworkMoleculeVDW + x2.frameworkMoleculeCharge -
         (x1.moleculeMoleculeVDW + x1.moleculeMoleculeCharge + x1.frameworkMoleculeVDW + x1.frameworkMoleculeCharge)) /
        delta;
    gradient.y =
        (y2.moleculeMoleculeVDW + y2.moleculeMoleculeCharge + y2.frameworkMoleculeVDW + y2.frameworkMoleculeCharge -
         (y1.moleculeMoleculeVDW + y1.moleculeMoleculeCharge + y1.frameworkMoleculeVDW + y1.frameworkMoleculeCharge)) /
        delta;
    gradient.z =
        (z2.moleculeMoleculeVDW + z2.moleculeMoleculeCharge + z2.frameworkMoleculeVDW + z2.frameworkMoleculeCharge -
         (z1.moleculeMoleculeVDW + z1.moleculeMoleculeCharge + z1.frameworkMoleculeVDW + z1.frameworkMoleculeCharge)) /
        delta;

    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.x, gradient.x, tolerance);
    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.y, gradient.y, tolerance);
    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.z, gradient.z, tolerance);
  }
}

TEST(gradients, Test_CO2_in_MFI_2x2x2)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Framework f = TestFactories::makeMFI_Si(forceField, int3(2, 2, 2));
  Component c = TestFactories::makeCO2(forceField, 0, true);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {10}, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  std::vector<Atom> atomPositions = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end());

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  [[maybe_unused]] RunningEnergy factor = Interactions::computeInterMolecularGradient(
      system.forceField, system.simulationBox, system.spanOfMoleculeAtoms());
  [[maybe_unused]] RunningEnergy factor2 = Interactions::computeFrameworkMoleculeGradient(
      system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
      system.interpolationGrids);

  double delta = 1e-5;
  double tolerance = 1e-4;
  double3 gradient;
  for (size_t i = 0; i < atomPositions.size(); ++i)
  {
    RunningEnergy x1, x2, y1, y2, z1, z2;

    // finite difference x
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x + 0.5 * delta;
    x2 =
        Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions) +
        Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                     system.framework, frameworkAtoms, atomPositions);

    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x - 0.5 * delta;
    x1 =
        Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions) +
        Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                     system.framework, frameworkAtoms, atomPositions);
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x;

    // finite difference y
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y + 0.5 * delta;
    y2 =
        Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions) +
        Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                     system.framework, frameworkAtoms, atomPositions);

    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y - 0.5 * delta;
    y1 =
        Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions) +
        Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                     system.framework, frameworkAtoms, atomPositions);
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y;

    // finite difference z
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z + 0.5 * delta;
    z2 =
        Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions) +
        Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                     system.framework, frameworkAtoms, atomPositions);

    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z - 0.5 * delta;
    z1 =
        Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions) +
        Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                     system.framework, frameworkAtoms, atomPositions);
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z;

    gradient.x =
        (x2.moleculeMoleculeVDW + x2.moleculeMoleculeCharge + x2.frameworkMoleculeVDW + x2.frameworkMoleculeCharge -
         (x1.moleculeMoleculeVDW + x1.moleculeMoleculeCharge + x1.frameworkMoleculeVDW + x1.frameworkMoleculeCharge)) /
        delta;
    gradient.y =
        (y2.moleculeMoleculeVDW + y2.moleculeMoleculeCharge + y2.frameworkMoleculeVDW + y2.frameworkMoleculeCharge -
         (y1.moleculeMoleculeVDW + y1.moleculeMoleculeCharge + y1.frameworkMoleculeVDW + y1.frameworkMoleculeCharge)) /
        delta;
    gradient.z =
        (z2.moleculeMoleculeVDW + z2.moleculeMoleculeCharge + z2.frameworkMoleculeVDW + z2.frameworkMoleculeCharge -
         (z1.moleculeMoleculeVDW + z1.moleculeMoleculeCharge + z1.frameworkMoleculeVDW + z1.frameworkMoleculeCharge)) /
        delta;

    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.x, gradient.x, tolerance);
    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.y, gradient.y, tolerance);
    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.z, gradient.z, tolerance);
  }
}

TEST(gradients, Test_20_Na_Cl_in_Box_25x25x25)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);

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

  [[maybe_unused]] RunningEnergy factor = Interactions::computeInterMolecularGradient(
      system.forceField, system.simulationBox, system.spanOfMoleculeAtoms());

  double delta = 1e-5;
  double tolerance = 1e-4;
  double3 gradient;
  for (size_t i = 0; i < atomPositions.size(); ++i)
  {
    RunningEnergy x1, x2, y1, y2, z1, z2;

    // finite difference x
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x + 0.5 * delta;
    x2 = Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);

    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x - 0.5 * delta;
    x1 = Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x;

    // finite difference y
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y + 0.5 * delta;
    y2 = Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);

    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y - 0.5 * delta;
    y1 = Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y;

    // finite difference z
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z + 0.5 * delta;
    z2 = Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);

    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z - 0.5 * delta;
    z1 = Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z;

    gradient.x =
        (x2.moleculeMoleculeVDW + x2.moleculeMoleculeCharge - (x1.moleculeMoleculeVDW + x1.moleculeMoleculeCharge)) /
        delta;
    gradient.y =
        (y2.moleculeMoleculeVDW + y2.moleculeMoleculeCharge - (y1.moleculeMoleculeVDW + y1.moleculeMoleculeCharge)) /
        delta;
    gradient.z =
        (z2.moleculeMoleculeVDW + z2.moleculeMoleculeCharge - (z1.moleculeMoleculeVDW + z1.moleculeMoleculeCharge)) /
        delta;

    EXPECT_NEAR(system.atomPositions[i].gradient.x, gradient.x, tolerance);
    EXPECT_NEAR(system.atomPositions[i].gradient.y, gradient.y, tolerance);
    EXPECT_NEAR(system.atomPositions[i].gradient.z, gradient.z, tolerance);
  }
}
