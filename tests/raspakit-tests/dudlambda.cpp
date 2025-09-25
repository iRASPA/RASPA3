#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <cstddef>
#include <span>
#include <tuple>
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
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;

TEST(dudlambda, Test_20_Na_Cl_in_Box_25x25x25_VDW)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Component na = TestFactories::makeIon(forceField, 0, "Na", 6, 0.0);
  Component cl = TestFactories::makeIon(forceField, 1, "Cl", 7, 0.0);
  System system =
      System(0, forceField, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, 1.0, {}, {na, cl}, {}, {20, 20}, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::vector<Atom> atomData = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end());

  for (Atom& atom : atomData)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  for (size_t i = 0; i < 20; ++i)
  {
    atomData[i].charge = 0.0;
    system.atomData[i].charge = 0.0;
  }
  for (size_t i = 0; i < 20; ++i)
  {
    atomData[i + 20].charge = 0.0;
    system.atomData[i + 20].charge = 0.0;
  }

  system.atomData[8].groupId = 1;
  system.atomData[2].groupId = 1;
  // system.atomData[12].groupId = 1;
  // system.atomData[4].groupId = 1;
  // system.atomData[14].groupId = 1;

  system.atomData[8].scalingVDW = 0.15;
  system.atomData[2].scalingVDW = 0.25;
  system.atomData[12].scalingVDW = 0.34;
  system.atomData[4].scalingVDW = 0.16;
  system.atomData[14].scalingVDW = 0.27;
  RunningEnergy factor = Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox,
                                                                     system.spanOfMoleculeAtoms());

  double delta = 1e-6;
  double tolerance = 1e-4;

  system.atomData[8].scalingVDW = 0.15 + 0.5 * delta;
  system.atomData[2].scalingVDW = 0.25 + 0.5 * delta;
  // system.atomData[12].scalingVDW = 0.34 + 0.5 * delta;
  // system.atomData[4].scalingVDW = 0.16 + 0.5 * delta;
  // system.atomData[14].scalingVDW = 0.27 + 0.5 * delta;
  RunningEnergy energyForward = Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox,
                                                                            system.spanOfMoleculeAtoms());

  system.atomData[8].scalingVDW = 0.15 - 0.5 * delta;
  system.atomData[2].scalingVDW = 0.25 - 0.5 * delta;
  // system.atomData[12].scalingVDW = 0.34 - 0.5 * delta;
  // system.atomData[4].scalingVDW = 0.16 - 0.5 * delta;
  // system.atomData[14].scalingVDW = 0.27 - 0.5 * delta;
  RunningEnergy energyBackward = Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox,
                                                                             system.spanOfMoleculeAtoms());

  double dUdlambda = ((energyForward.moleculeMoleculeVDW + energyForward.moleculeMoleculeCharge) -
                      (energyBackward.moleculeMoleculeVDW + energyBackward.moleculeMoleculeCharge)) /
                     delta;

  EXPECT_NEAR(factor.dudlambdaVDW + factor.dudlambdaCharge, dUdlambda, tolerance)
      << " ratio: " << (factor.dudlambdaVDW + factor.dudlambdaCharge) / dUdlambda << " "
      << dUdlambda / (factor.dudlambdaVDW + factor.dudlambdaCharge);
}

TEST(dudlambda, Test_20_Na_Cl_in_Box_25x25x25_Coulomb)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Component na = TestFactories::makeIon(forceField, 0, "Na", 6, 0.0);
  Component cl = TestFactories::makeIon(forceField, 1, "Cl", 7, 0.0);
  System system =
      System(0, forceField, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, 1.0, {}, {na, cl}, {}, {20, 20}, 5);

  std::fill(system.forceField.data.begin(), system.forceField.data.end(), VDWParameters(0.0, 1.0));

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::vector<Atom> atomData = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end());

  for (Atom& atom : atomData)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  for (size_t i = 0; i < 20; ++i)
  {
    atomData[i].charge = 1.0;
    system.atomData[i].charge = 1.0;
  }
  for (size_t i = 0; i < 20; ++i)
  {
    atomData[i + 20].charge = -1.0;
    system.atomData[i + 20].charge = -1.0;
  }

  system.atomData[8].groupId = 1;
  system.atomData[2].groupId = 1;
  // system.atomData[12].groupId = 1;
  // system.atomData[4].groupId = 1;
  // system.atomData[14].groupId = 1;

  system.atomData[8].scalingCoulomb = 0.45;
  system.atomData[2].scalingCoulomb = 0.5;
  system.atomData[12].scalingCoulomb = 0.4;
  system.atomData[4].scalingCoulomb = 0.6;
  system.atomData[14].scalingCoulomb = 0.7;
  RunningEnergy factor = Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox,
                                                                     system.spanOfMoleculeAtoms());

  double delta = 1e-6;
  double tolerance = 1e-4;

  system.atomData[8].scalingCoulomb = 0.45 + 0.5 * delta;
  system.atomData[2].scalingCoulomb = 0.5 + 0.5 * delta;
  // system.atomData[12].scalingCoulomb = 0.4 + 0.5 * delta;
  // system.atomData[4].scalingCoulomb = 0.6 + 0.5 * delta;
  // system.atomData[14].scalingCoulomb = 0.7 + 0.5 * delta;
  RunningEnergy energyForward = Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox,
                                                                            system.spanOfMoleculeAtoms());

  system.atomData[8].scalingCoulomb = 0.45 - 0.5 * delta;
  system.atomData[2].scalingCoulomb = 0.5 - 0.5 * delta;
  // system.atomData[12].scalingCoulomb = 0.4 - 0.5 * delta;
  // system.atomData[4].scalingCoulomb = 0.6 - 0.5 * delta;
  // system.atomData[14].scalingCoulomb = 0.7 - 0.5 * delta;
  RunningEnergy energyBackward = Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox,
                                                                             system.spanOfMoleculeAtoms());

  double dUdlambda = ((energyForward.moleculeMoleculeVDW + energyForward.moleculeMoleculeCharge) -
                      (energyBackward.moleculeMoleculeVDW + energyBackward.moleculeMoleculeCharge)) /
                     delta;

  EXPECT_NEAR(factor.dudlambdaVDW + factor.dudlambdaCharge, dUdlambda, tolerance)
      << " ratio: " << (factor.dudlambdaVDW + factor.dudlambdaCharge) / dUdlambda << " "
      << dUdlambda / (factor.dudlambdaVDW + factor.dudlambdaCharge);
}

TEST(dudlambda, Test_20_Na_Cl_in_Box_25x25x25_Fourier)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Component na = TestFactories::makeIon(forceField, 0, "Na", 6, 0.0);
  Component cl = TestFactories::makeIon(forceField, 1, "Cl", 7, 0.0);
  System system =
      System(0, forceField, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, 1.0, {}, {na, cl}, {}, {20, 20}, 5);

  std::fill(system.forceField.data.begin(), system.forceField.data.end(), VDWParameters(0.0, 1.0));

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::vector<Atom> atomData = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end());

  for (Atom& atom : atomData)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  for (size_t i = 0; i < 20; ++i)
  {
    atomData[i].charge = 1.0;
    system.atomData[i].charge = 1.0;
  }
  for (size_t i = 0; i < 20; ++i)
  {
    atomData[i + 20].charge = -1.0;
    system.atomData[i + 20].charge = -1.0;
  }

  system.atomData[8].groupId = 1;
  system.atomData[2].groupId = 1;
  // system.atomData[12].groupId = 1;
  // system.atomData[4].groupId = 1;
  // system.atomData[14].groupId = 1;

  system.atomData[8].scalingCoulomb = 0.45;
  system.atomData[2].scalingCoulomb = 0.5;
  system.atomData[12].scalingCoulomb = 0.4;
  system.atomData[4].scalingCoulomb = 0.6;
  system.atomData[14].scalingCoulomb = 0.7;
  RunningEnergy factor = Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox,
                                                                     system.spanOfMoleculeAtoms());

  double delta = 1e-5;
  double tolerance = 1e-3;

  system.atomData[8].scalingCoulomb = 0.45 + 0.5 * delta;
  system.atomData[2].scalingCoulomb = 0.5 + 0.5 * delta;
  // system.atomData[12].scalingCoulomb = 0.4 + 0.5 * delta;
  // system.atomData[4].scalingCoulomb = 0.6 + 0.5 * delta;
  // system.atomData[14].scalingCoulomb = 0.7 + 0.5 * delta;
  RunningEnergy energyForward = Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox,
                                                                            system.spanOfMoleculeAtoms());

  system.atomData[8].scalingCoulomb = 0.45 - 0.5 * delta;
  system.atomData[2].scalingCoulomb = 0.5 - 0.5 * delta;
  // system.atomData[12].scalingCoulomb = 0.4 - 0.5 * delta;
  // system.atomData[4].scalingCoulomb = 0.6 - 0.5 * delta;
  // system.atomData[14].scalingCoulomb = 0.7 - 0.5 * delta;
  RunningEnergy energyBackward = Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox,
                                                                             system.spanOfMoleculeAtoms());

  double dUdlambda = ((energyForward.moleculeMoleculeVDW + energyForward.moleculeMoleculeCharge) -
                      (energyBackward.moleculeMoleculeVDW + energyBackward.moleculeMoleculeCharge)) /
                     delta;

  EXPECT_NEAR(factor.dudlambdaVDW + factor.dudlambdaCharge, dUdlambda, tolerance)
      << " ratio: " << (factor.dudlambdaVDW + factor.dudlambdaCharge) / dUdlambda << " "
      << dUdlambda / (factor.dudlambdaVDW + factor.dudlambdaCharge);
}

TEST(dudlambda, Test_20_CO2_in_Box_25x25x25_Fourier)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Component c = TestFactories::makeCO2(forceField, 0, true);
  System system = System(0, forceField, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, 1.0, {}, {c}, {}, {20}, 5);

  std::fill(system.forceField.data.begin(), system.forceField.data.end(), VDWParameters(0.0, 1.0));

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::vector<Atom> atomData = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end());

  for (Atom& atom : atomData)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  system.atomData[8].groupId = 1;
  system.atomData[2].groupId = 1;
  // system.atomData[12].groupId = 1;
  // system.atomData[4].groupId = 1;
  // system.atomData[14].groupId = 1;

  system.atomData[8].scalingCoulomb = 0.45;
  system.atomData[2].scalingCoulomb = 0.5;
  system.atomData[12].scalingCoulomb = 0.4;
  system.atomData[4].scalingCoulomb = 0.6;
  system.atomData[14].scalingCoulomb = 0.7;
  RunningEnergy factor = Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox,
                                                                     system.spanOfMoleculeAtoms());

  double delta = 1e-5;
  double tolerance = 1e-4;

  system.atomData[8].scalingCoulomb = 0.45 + 0.5 * delta;
  system.atomData[2].scalingCoulomb = 0.5 + 0.5 * delta;
  // system.atomData[12].scalingCoulomb = 0.4 + 0.5 * delta;
  // system.atomData[4].scalingCoulomb = 0.6 + 0.5 * delta;
  // system.atomData[14].scalingCoulomb = 0.7 + 0.5 * delta;
  RunningEnergy energyForward = Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox,
                                                                            system.spanOfMoleculeAtoms());

  system.atomData[8].scalingCoulomb = 0.45 - 0.5 * delta;
  system.atomData[2].scalingCoulomb = 0.5 - 0.5 * delta;
  // system.atomData[12].scalingCoulomb = 0.4 - 0.5 * delta;
  // system.atomData[4].scalingCoulomb = 0.6 - 0.5 * delta;
  // system.atomData[14].scalingCoulomb = 0.7 - 0.5 * delta;
  RunningEnergy energyBackward = Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox,
                                                                             system.spanOfMoleculeAtoms());

  double dUdlambda = ((energyForward.moleculeMoleculeVDW + energyForward.moleculeMoleculeCharge) -
                      (energyBackward.moleculeMoleculeVDW + energyBackward.moleculeMoleculeCharge)) /
                     delta;

  EXPECT_NEAR(factor.dudlambdaVDW + factor.dudlambdaCharge, dUdlambda, tolerance)
      << " ratio: " << (factor.dudlambdaVDW + factor.dudlambdaCharge) / dUdlambda << " "
      << dUdlambda / (factor.dudlambdaVDW + factor.dudlambdaCharge);
}

TEST(dudlambda, Test_2_CO2_in_MFI_2x2x2_VDW)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Component c = TestFactories::makeCO2(forceField, 0, false);
  Framework f = TestFactories::makeMFI_Si(forceField, int3(2, 2, 2));
  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {}, {2}, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::vector<Atom> atomData = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end());

  atomData[0].position = double3(10.011, 4.97475 + 2.0, 1.149);
  atomData[1].position = double3(10.011, 4.97475 + 2.0, 0.0);
  atomData[2].position = double3(10.011, 4.97475 + 2.0, -1.149);
  atomData[3].position = double3(10.011, 4.97475 - 2.0, 1.149);
  atomData[4].position = double3(10.011, 4.97475 - 2.0, 0.0);
  atomData[5].position = double3(10.011, 4.97475 - 2.0, -1.149);

  // atomData[8].groupId = 1;
  atomData[2].groupId = 1;
  // atomData[12].groupId = 1;
  atomData[4].groupId = 1;
  // atomData[14].groupId = 1;

  // atomData[8].scalingVDW = 0.45;
  atomData[2].scalingVDW = 0.5;
  // atomData[12].scalingVDW = 0.4;
  atomData[4].scalingVDW = 0.6;
  // atomData[14].scalingVDW = 0.7;
  RunningEnergy factor = Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomData);

  double delta = 1e-5;
  double tolerance = 1e-3;

  // atomData[8].scalingVDW = 0.45 + 0.5 * delta;
  atomData[2].scalingVDW = 0.5 + 0.5 * delta;
  // atomData[12].scalingVDW = 0.4 + 0.5 * delta;
  atomData[4].scalingVDW = 0.6 + 0.5 * delta;
  // atomData[14].scalingVDW = 0.7 + 0.5 * delta;
  RunningEnergy energyForward =
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomData);

  // atomData[8].scalingVDW = 0.45 - 0.5 * delta;
  atomData[2].scalingVDW = 0.5 - 0.5 * delta;
  // atomData[12].scalingVDW = 0.4 - 0.5 * delta;
  atomData[4].scalingVDW = 0.6 - 0.5 * delta;
  // atomData[14].scalingVDW = 0.7 - 0.5 * delta;
  RunningEnergy energyBackward =
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomData);

  double dUdlambda = (energyForward.potentialEnergy() - energyBackward.potentialEnergy()) / delta;

  EXPECT_NEAR(factor.dudlambdaVDW, dUdlambda, tolerance)
      << " ratio: " << factor.dudlambdaVDW / dUdlambda << " " << dUdlambda / factor.dudlambdaVDW;
}

TEST(dudlambda, Test_2_CO2_in_MFI_2x2x2_Coulomb)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Component c = TestFactories::makeCO2(forceField, 0, true);
  Framework f = TestFactories::makeMFI_Si(forceField, int3(2, 2, 2));
  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {}, {2}, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::vector<Atom> atomData = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end());

  atomData[0].position = double3(10.011, 4.97475 + 2.0, 1.149);
  atomData[1].position = double3(10.011, 4.97475 + 2.0, 0.0);
  atomData[2].position = double3(10.011, 4.97475 + 2.0, -1.149);
  atomData[3].position = double3(10.011, 4.97475 - 2.0, 1.149);
  atomData[4].position = double3(10.011, 4.97475 - 2.0, 0.0);
  atomData[5].position = double3(10.011, 4.97475 - 2.0, -1.149);

  // atomData[8].groupId = 1;
  atomData[2].groupId = 1;
  // atomData[12].groupId = 1;
  atomData[4].groupId = 1;
  // atomData[14].groupId = 1;

  // atomData[8].scalingCoulomb = 0.45;
  atomData[2].scalingCoulomb = 0.5;
  // atomData[12].scalingCoulomb = 0.4;
  atomData[4].scalingCoulomb = 0.6;
  // atomData[14].scalingCoulomb = 0.7;
  RunningEnergy factor = Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomData);

  double delta = 1e-5;
  double tolerance = 1e-3;

  // atomData[8].scalingCoulomb = 0.45 + 0.5 * delta;
  atomData[2].scalingCoulomb = 0.5 + 0.5 * delta;
  // atomData[12].scalingCoulomb = 0.4 + 0.5 * delta;
  atomData[4].scalingCoulomb = 0.6 + 0.5 * delta;
  // atomData[14].scalingCoulomb = 0.7 + 0.5 * delta;
  RunningEnergy energyForward =
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomData);

  // atomData[8].scalingCoulomb = 0.45 - 0.5 * delta;
  atomData[2].scalingCoulomb = 0.5 - 0.5 * delta;
  // atomData[12].scalingCoulomb = 0.4 - 0.5 * delta;
  atomData[4].scalingCoulomb = 0.6 - 0.5 * delta;
  // atomData[14].scalingCoulomb = 0.7 - 0.5 * delta;
  RunningEnergy energyBackward =
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomData);

  double dUdlambda = (energyForward.potentialEnergy() - energyBackward.potentialEnergy()) / delta;

  EXPECT_NEAR(factor.dudlambdaCharge, dUdlambda, tolerance)
      << " ratio: " << factor.dudlambdaCharge / dUdlambda << " " << dUdlambda / factor.dudlambdaCharge;
}

TEST(dudlambda, Test_2_CO2_in_MFI_2x2x2_Ewald)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Component c = TestFactories::makeCO2(forceField, 0, true);
  Framework f = TestFactories::makeMFI_Si(forceField, int3(2, 2, 2));
  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {}, {2}, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::vector<Atom> atomData = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end());

  atomData[0].position = double3(10.011, 4.97475 + 2.0, 1.149);
  atomData[1].position = double3(10.011, 4.97475 + 2.0, 0.0);
  atomData[2].position = double3(10.011, 4.97475 + 2.0, -1.149);
  atomData[3].position = double3(10.011, 4.97475 - 2.0, 1.149);
  atomData[4].position = double3(10.011, 4.97475 - 2.0, 0.0);
  atomData[5].position = double3(10.011, 4.97475 - 2.0, -1.149);

  // atomData[8].groupId = 1;
  atomData[2].groupId = 1;
  // atomData[12].groupId = 1;
  atomData[4].groupId = 1;
  // atomData[14].groupId = 1;

  // atomData[8].scalingCoulomb = 0.45;
  atomData[2].scalingCoulomb = 0.5;
  // atomData[12].scalingCoulomb = 0.4;
  atomData[4].scalingCoulomb = 0.6;
  // atomData[14].scalingCoulomb = 0.7;

  system.precomputeTotalRigidEnergy();
  RunningEnergy factor = Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent, atomData);

  double delta = 1e-5;
  double tolerance = 1e-3;

  // atomData[8].scalingCoulomb = 0.45 + 0.5 * delta;
  atomData[2].scalingCoulomb = 0.5 + 0.5 * delta;
  // atomData[12].scalingCoulomb = 0.4 + 0.5 * delta;
  atomData[4].scalingCoulomb = 0.6 + 0.5 * delta;
  // atomData[14].scalingCoulomb = 0.7 + 0.5 * delta;
  RunningEnergy energyForward = Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent, atomData);

  // atomData[8].scalingCoulomb = 0.45 - 0.5 * delta;
  atomData[2].scalingCoulomb = 0.5 - 0.5 * delta;
  // atomData[12].scalingCoulomb = 0.4 - 0.5 * delta;
  atomData[4].scalingCoulomb = 0.6 - 0.5 * delta;
  // atomData[14].scalingCoulomb = 0.7 - 0.5 * delta;
  RunningEnergy energyBackward = Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent, atomData);

  double dUdlambda = (energyForward.potentialEnergy() - energyBackward.potentialEnergy()) / delta;

  EXPECT_NEAR(factor.dudlambdaEwald, dUdlambda, tolerance)
      << " ratio: " << factor.dudlambdaEwald / dUdlambda << " " << dUdlambda / factor.dudlambdaEwald;
}
