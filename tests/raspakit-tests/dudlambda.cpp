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
import force_factor;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;

TEST(dudlambda, Test_20_Na_Cl_in_Box_25x25x25_VDW)
{
  ForceField forceField = ForceField(
      {
          PseudoAtom("Si", true, 28.0855, 2.05, 0.0, 14, false),
          PseudoAtom("O", true, 15.999, -1.025, 0.0, 8, false),
          PseudoAtom("CH4", false, 16.04246, 0.0, 0.0, 6, false),
          PseudoAtom("Na+", false, 12.0, 0.0, 0.0, 6, false),
          PseudoAtom("Cl-", false, 15.9994, 0.0, 0.0, 8, false),
      },
      {VDWParameters(22.0, 2.30), VDWParameters(53.0, 3.3), VDWParameters(158.5, 3.72), VDWParameters(15.0966, 2.65755),
       VDWParameters(142.562, 3.51932)},
      ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, true);
  Component na = Component(0, forceField, "Na", 304.1282, 7377300.0, 0.22394,
                           {
                               // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                               // uint8_t componentId, uint8_t groupId
                               Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 3, 0, 0),
                           },
                           5, 21);
  Component cl = Component(1, forceField, "Cl", 304.1282, 7377300.0, 0.22394,
                           {
                               // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                               // uint8_t componentId, uint8_t groupId
                               Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 4, 1, 0),
                           },
                           5, 21);

  System system = System(0, forceField, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, 1.0, {}, {na, cl}, {20, 20}, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::vector<Atom> atomPositions = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end());

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  for (size_t i = 0; i < 20; ++i)
  {
    atomPositions[i].charge = 0.0;
    system.atomPositions[i].charge = 0.0;
  }
  for (size_t i = 0; i < 20; ++i)
  {
    atomPositions[i + 20].charge = 0.0;
    system.atomPositions[i + 20].charge = 0.0;
  }

  system.atomPositions[8].groupId = 1;
  system.atomPositions[2].groupId = 1;
  // system.atomPositions[12].groupId = 1;
  // system.atomPositions[4].groupId = 1;
  // system.atomPositions[14].groupId = 1;

  system.atomPositions[8].scalingVDW = 0.15;
  system.atomPositions[2].scalingVDW = 0.25;
  system.atomPositions[12].scalingVDW = 0.34;
  system.atomPositions[4].scalingVDW = 0.16;
  system.atomPositions[14].scalingVDW = 0.27;
  RunningEnergy factor = Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox,
                                                                     system.spanOfMoleculeAtoms());

  double delta = 1e-6;
  double tolerance = 1e-4;

  system.atomPositions[8].scalingVDW = 0.15 + 0.5 * delta;
  system.atomPositions[2].scalingVDW = 0.25 + 0.5 * delta;
  // system.atomPositions[12].scalingVDW = 0.34 + 0.5 * delta;
  // system.atomPositions[4].scalingVDW = 0.16 + 0.5 * delta;
  // system.atomPositions[14].scalingVDW = 0.27 + 0.5 * delta;
  RunningEnergy energyForward = Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox,
                                                                            system.spanOfMoleculeAtoms());

  system.atomPositions[8].scalingVDW = 0.15 - 0.5 * delta;
  system.atomPositions[2].scalingVDW = 0.25 - 0.5 * delta;
  // system.atomPositions[12].scalingVDW = 0.34 - 0.5 * delta;
  // system.atomPositions[4].scalingVDW = 0.16 - 0.5 * delta;
  // system.atomPositions[14].scalingVDW = 0.27 - 0.5 * delta;
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
  ForceField forceField = ForceField(
      {
          PseudoAtom("Si", true, 28.0855, 2.05, 0.0, 14, false),
          PseudoAtom("O", true, 15.999, -1.025, 0.0, 8, false),
          PseudoAtom("CH4", false, 16.04246, 0.0, 0.0, 6, false),
          PseudoAtom("Na+", false, 12.0, 0.0, 0.0, 6, false),
          PseudoAtom("Cl-", false, 15.9994, 0.0, 0.0, 8, false),
      },
      {VDWParameters(22.0, 2.30), VDWParameters(53.0, 3.3), VDWParameters(158.5, 3.72), VDWParameters(15.0966, 2.65755),
       VDWParameters(142.562, 3.51932)},
      ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false);
  Component na = Component(0, forceField, "Na", 304.1282, 7377300.0, 0.22394,
                           {
                               // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                               // uint8_t componentId, uint8_t groupId
                               Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 3, 0, 0),
                           },
                           5, 21);
  Component cl = Component(1, forceField, "Cl", 304.1282, 7377300.0, 0.22394,
                           {
                               // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                               // uint8_t componentId, uint8_t groupId
                               Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 4, 1, 0),
                           },
                           5, 21);

  System system = System(0, forceField, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, 1.0, {}, {na, cl}, {20, 20}, 5);

  std::fill(system.forceField.data.begin(), system.forceField.data.end(), VDWParameters(0.0, 1.0));

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

  system.atomPositions[8].groupId = 1;
  system.atomPositions[2].groupId = 1;
  // system.atomPositions[12].groupId = 1;
  // system.atomPositions[4].groupId = 1;
  // system.atomPositions[14].groupId = 1;

  system.atomPositions[8].scalingCoulomb = 0.45;
  system.atomPositions[2].scalingCoulomb = 0.5;
  system.atomPositions[12].scalingCoulomb = 0.4;
  system.atomPositions[4].scalingCoulomb = 0.6;
  system.atomPositions[14].scalingCoulomb = 0.7;
  RunningEnergy factor = Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox,
                                                                     system.spanOfMoleculeAtoms());

  double delta = 1e-6;
  double tolerance = 1e-4;

  system.atomPositions[8].scalingCoulomb = 0.45 + 0.5 * delta;
  system.atomPositions[2].scalingCoulomb = 0.5 + 0.5 * delta;
  // system.atomPositions[12].scalingCoulomb = 0.4 + 0.5 * delta;
  // system.atomPositions[4].scalingCoulomb = 0.6 + 0.5 * delta;
  // system.atomPositions[14].scalingCoulomb = 0.7 + 0.5 * delta;
  RunningEnergy energyForward = Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox,
                                                                            system.spanOfMoleculeAtoms());

  system.atomPositions[8].scalingCoulomb = 0.45 - 0.5 * delta;
  system.atomPositions[2].scalingCoulomb = 0.5 - 0.5 * delta;
  // system.atomPositions[12].scalingCoulomb = 0.4 - 0.5 * delta;
  // system.atomPositions[4].scalingCoulomb = 0.6 - 0.5 * delta;
  // system.atomPositions[14].scalingCoulomb = 0.7 - 0.5 * delta;
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
  ForceField forceField = ForceField(
      {
          PseudoAtom("Si", true, 28.0855, 2.05, 0.0, 14, false),
          PseudoAtom("O", true, 15.999, -1.025, 0.0, 8, false),
          PseudoAtom("CH4", false, 16.04246, 0.0, 0.0, 6, false),
          PseudoAtom("Na+", false, 12.0, 0.0, 0.0, 6, false),
          PseudoAtom("Cl-", false, 15.9994, 0.0, 0.0, 8, false),
      },
      {VDWParameters(22.0, 2.30), VDWParameters(53.0, 3.3), VDWParameters(158.5, 3.72), VDWParameters(15.0966, 2.65755),
       VDWParameters(142.562, 3.51932)},
      ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false);
  Component na = Component(0, forceField, "Na", 304.1282, 7377300.0, 0.22394,
                           {
                               // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                               // uint8_t componentId, uint8_t groupId
                               Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 3, 0, 0),
                           },
                           5, 21);
  Component cl = Component(1, forceField, "Cl", 304.1282, 7377300.0, 0.22394,
                           {
                               // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                               // uint8_t componentId, uint8_t groupId
                               Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 4, 1, 0),
                           },
                           5, 21);

  System system = System(0, forceField, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, 1.0, {}, {na, cl}, {20, 20}, 5);

  std::fill(system.forceField.data.begin(), system.forceField.data.end(), VDWParameters(0.0, 1.0));

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

  system.atomPositions[8].groupId = 1;
  system.atomPositions[2].groupId = 1;
  // system.atomPositions[12].groupId = 1;
  // system.atomPositions[4].groupId = 1;
  // system.atomPositions[14].groupId = 1;

  system.atomPositions[8].scalingCoulomb = 0.45;
  system.atomPositions[2].scalingCoulomb = 0.5;
  system.atomPositions[12].scalingCoulomb = 0.4;
  system.atomPositions[4].scalingCoulomb = 0.6;
  system.atomPositions[14].scalingCoulomb = 0.7;
  RunningEnergy factor = Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox,
                                                                     system.spanOfMoleculeAtoms());

  double delta = 1e-5;
  double tolerance = 1e-3;

  system.atomPositions[8].scalingCoulomb = 0.45 + 0.5 * delta;
  system.atomPositions[2].scalingCoulomb = 0.5 + 0.5 * delta;
  // system.atomPositions[12].scalingCoulomb = 0.4 + 0.5 * delta;
  // system.atomPositions[4].scalingCoulomb = 0.6 + 0.5 * delta;
  // system.atomPositions[14].scalingCoulomb = 0.7 + 0.5 * delta;
  RunningEnergy energyForward = Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox,
                                                                            system.spanOfMoleculeAtoms());

  system.atomPositions[8].scalingCoulomb = 0.45 - 0.5 * delta;
  system.atomPositions[2].scalingCoulomb = 0.5 - 0.5 * delta;
  // system.atomPositions[12].scalingCoulomb = 0.4 - 0.5 * delta;
  // system.atomPositions[4].scalingCoulomb = 0.6 - 0.5 * delta;
  // system.atomPositions[14].scalingCoulomb = 0.7 - 0.5 * delta;
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
  ForceField forceField = ForceField(
      {
          PseudoAtom("Si", true, 28.0855, 2.05, 0.0, 14, false),
          PseudoAtom("O", true, 15.999, -1.025, 0.0, 8, false),
          PseudoAtom("CH4", false, 16.04246, 0.0, 0.0, 6, false),
          PseudoAtom("C_co2", false, 12.0, 0.6512, 0.0, 6, false),
          PseudoAtom("O_co2", false, 15.9994, -0.3256, 0.0, 8, false),
      },
      {VDWParameters(22.0, 2.30), VDWParameters(53.0, 3.3), VDWParameters(158.5, 3.72), VDWParameters(29.933, 2.745),
       VDWParameters(85.671, 3.017)},
      ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false);

  Component c = Component(
      0, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
      {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId,
       // uint8_t groupId
       Atom(double3(0.0, 0.0, 1.149), -0.3256, 1.0, 0, 4, 1, 0), Atom(double3(0.0, 0.0, 0.0), 0.6512, 1.0, 0, 3, 1, 0),
       Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)},
      5, 21);

  System system = System(0, forceField, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, 1.0, {}, {c}, {20}, 5);

  std::fill(system.forceField.data.begin(), system.forceField.data.end(), VDWParameters(0.0, 1.0));

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::vector<Atom> atomPositions = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end());

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  system.atomPositions[8].groupId = 1;
  system.atomPositions[2].groupId = 1;
  // system.atomPositions[12].groupId = 1;
  // system.atomPositions[4].groupId = 1;
  // system.atomPositions[14].groupId = 1;

  system.atomPositions[8].scalingCoulomb = 0.45;
  system.atomPositions[2].scalingCoulomb = 0.5;
  system.atomPositions[12].scalingCoulomb = 0.4;
  system.atomPositions[4].scalingCoulomb = 0.6;
  system.atomPositions[14].scalingCoulomb = 0.7;
  RunningEnergy factor = Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox,
                                                                     system.spanOfMoleculeAtoms());

  double delta = 1e-5;
  double tolerance = 1e-4;

  system.atomPositions[8].scalingCoulomb = 0.45 + 0.5 * delta;
  system.atomPositions[2].scalingCoulomb = 0.5 + 0.5 * delta;
  // system.atomPositions[12].scalingCoulomb = 0.4 + 0.5 * delta;
  // system.atomPositions[4].scalingCoulomb = 0.6 + 0.5 * delta;
  // system.atomPositions[14].scalingCoulomb = 0.7 + 0.5 * delta;
  RunningEnergy energyForward = Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox,
                                                                            system.spanOfMoleculeAtoms());

  system.atomPositions[8].scalingCoulomb = 0.45 - 0.5 * delta;
  system.atomPositions[2].scalingCoulomb = 0.5 - 0.5 * delta;
  // system.atomPositions[12].scalingCoulomb = 0.4 - 0.5 * delta;
  // system.atomPositions[4].scalingCoulomb = 0.6 - 0.5 * delta;
  // system.atomPositions[14].scalingCoulomb = 0.7 - 0.5 * delta;
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
  ForceField forceField = ForceField(
      {
          PseudoAtom("Si", true, 28.0855, 2.05, 0.0, 14, false),
          PseudoAtom("O", true, 15.999, -1.025, 0.0, 8, false),
          PseudoAtom("CH4", false, 16.04246, 0.0, 0.0, 6, false),
          PseudoAtom("C_co2", false, 12.0, 0.6512, 0.0, 6, false),
          PseudoAtom("O_co2", false, 15.9994, -0.3256, 0.0, 8, false),
      },
      {VDWParameters(22.0, 2.30), VDWParameters(53.0, 3.3), VDWParameters(158.5, 3.72), VDWParameters(29.933, 2.745),
       VDWParameters(85.671, 3.017)},
      ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false);

  Framework f = Framework(0, forceField, "MFI_SI", SimulationBox(20.022, 19.899, 13.383), 292,
                          {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                           // uint8_t componentId, uint8_t groupId
                           Atom(double3(0.42238, 0.0565, -0.33598), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.30716, 0.02772, -0.1893), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.27911, 0.06127, 0.0312), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.12215, 0.06298, 0.0267), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.07128, 0.02722, -0.18551), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.18641, 0.05896, -0.32818), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.42265, -0.1725, -0.32718), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.30778, -0.13016, -0.18548), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.27554, -0.17279, 0.03109), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.12058, -0.1731, 0.02979), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.07044, -0.13037, -0.182), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.18706, -0.17327, -0.31933), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.3726, 0.0534, -0.2442), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.3084, 0.0587, -0.0789), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.2007, 0.0592, 0.0289), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.0969, 0.0611, -0.0856), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.1149, 0.0541, -0.2763), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.2435, 0.0553, -0.246), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.3742, -0.1561, -0.2372), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.3085, -0.1552, -0.0728), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.198, -0.1554, 0.0288), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.091, -0.1614, -0.0777), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.1169, -0.1578, -0.2694), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.2448, -0.1594, -0.2422), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.3047, -0.051, -0.1866), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.0768, -0.0519, -0.1769), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.4161, 0.1276, -0.3896), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.4086, -0.0017, -0.4136), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.402, -0.1314, -0.4239), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.1886, 0.1298, -0.3836), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.194, 0.0007, -0.4082), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.1951, -0.1291, -0.419), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(-0.0037, 0.0502, -0.208), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(-0.004, -0.1528, -0.2078), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.4192, -0.25, -0.354), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.1884, -0.25, -0.3538), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.2883, -0.25, 0.0579), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.1085, -0.25, 0.0611), -1.025, 1.0, 0, 1, 0, 0)},
                          int3(2, 2, 2));
  Component c = Component(
      0, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
      {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId,
       // uint8_t groupId
       Atom(double3(0.0, 0.0, 1.149), 0.0, 1.0, 0, 4, 1, 0), Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 3, 1, 0),
       Atom(double3(0.0, 0.0, -1.149), 0.0, 1.0, 0, 4, 1, 0)},
      5, 21);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {2}, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::vector<Atom> atomPositions = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end());

  atomPositions[0].position = double3(10.011, 4.97475 + 2.0, 1.149);
  atomPositions[1].position = double3(10.011, 4.97475 + 2.0, 0.0);
  atomPositions[2].position = double3(10.011, 4.97475 + 2.0, -1.149);
  atomPositions[3].position = double3(10.011, 4.97475 - 2.0, 1.149);
  atomPositions[4].position = double3(10.011, 4.97475 - 2.0, 0.0);
  atomPositions[5].position = double3(10.011, 4.97475 - 2.0, -1.149);

  // atomPositions[8].groupId = 1;
  atomPositions[2].groupId = 1;
  // atomPositions[12].groupId = 1;
  atomPositions[4].groupId = 1;
  // atomPositions[14].groupId = 1;

  // atomPositions[8].scalingVDW = 0.45;
  atomPositions[2].scalingVDW = 0.5;
  // atomPositions[12].scalingVDW = 0.4;
  atomPositions[4].scalingVDW = 0.6;
  // atomPositions[14].scalingVDW = 0.7;
  RunningEnergy factor =
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);

  double delta = 1e-5;
  double tolerance = 1e-3;

  // atomPositions[8].scalingVDW = 0.45 + 0.5 * delta;
  atomPositions[2].scalingVDW = 0.5 + 0.5 * delta;
  // atomPositions[12].scalingVDW = 0.4 + 0.5 * delta;
  atomPositions[4].scalingVDW = 0.6 + 0.5 * delta;
  // atomPositions[14].scalingVDW = 0.7 + 0.5 * delta;
  RunningEnergy energyForward =
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);

  // atomPositions[8].scalingVDW = 0.45 - 0.5 * delta;
  atomPositions[2].scalingVDW = 0.5 - 0.5 * delta;
  // atomPositions[12].scalingVDW = 0.4 - 0.5 * delta;
  atomPositions[4].scalingVDW = 0.6 - 0.5 * delta;
  // atomPositions[14].scalingVDW = 0.7 - 0.5 * delta;
  RunningEnergy energyBackward =
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);

  double dUdlambda = (energyForward.potentialEnergy() - energyBackward.potentialEnergy()) / delta;

  EXPECT_NEAR(factor.dudlambdaVDW, dUdlambda, tolerance)
      << " ratio: " << factor.dudlambdaVDW / dUdlambda << " " << dUdlambda / factor.dudlambdaVDW;
}

TEST(dudlambda, Test_2_CO2_in_MFI_2x2x2_Coulomb)
{
  ForceField forceField = ForceField(
      {
          PseudoAtom("Si", true, 28.0855, 2.05, 0.0, 14, false),
          PseudoAtom("O", true, 15.999, -1.025, 0.0, 8, false),
          PseudoAtom("CH4", false, 16.04246, 0.0, 0.0, 6, false),
          PseudoAtom("C_co2", false, 12.0, 0.6512, 0.0, 6, false),
          PseudoAtom("O_co2", false, 15.9994, -0.3256, 0.0, 8, false),
      },
      {VDWParameters(22.0, 2.30), VDWParameters(53.0, 3.3), VDWParameters(158.5, 3.72), VDWParameters(29.933, 2.745),
       VDWParameters(85.671, 3.017)},
      ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false);
  Framework f = Framework(0, forceField, "MFI_SI", SimulationBox(20.022, 19.899, 13.383), 292,
                          {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                           // uint8_t componentId, uint8_t groupId
                           Atom(double3(0.42238, 0.0565, -0.33598), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.30716, 0.02772, -0.1893), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.27911, 0.06127, 0.0312), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.12215, 0.06298, 0.0267), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.07128, 0.02722, -0.18551), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.18641, 0.05896, -0.32818), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.42265, -0.1725, -0.32718), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.30778, -0.13016, -0.18548), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.27554, -0.17279, 0.03109), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.12058, -0.1731, 0.02979), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.07044, -0.13037, -0.182), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.18706, -0.17327, -0.31933), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.3726, 0.0534, -0.2442), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.3084, 0.0587, -0.0789), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.2007, 0.0592, 0.0289), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.0969, 0.0611, -0.0856), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.1149, 0.0541, -0.2763), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.2435, 0.0553, -0.246), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.3742, -0.1561, -0.2372), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.3085, -0.1552, -0.0728), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.198, -0.1554, 0.0288), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.091, -0.1614, -0.0777), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.1169, -0.1578, -0.2694), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.2448, -0.1594, -0.2422), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.3047, -0.051, -0.1866), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.0768, -0.0519, -0.1769), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.4161, 0.1276, -0.3896), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.4086, -0.0017, -0.4136), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.402, -0.1314, -0.4239), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.1886, 0.1298, -0.3836), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.194, 0.0007, -0.4082), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.1951, -0.1291, -0.419), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(-0.0037, 0.0502, -0.208), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(-0.004, -0.1528, -0.2078), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.4192, -0.25, -0.354), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.1884, -0.25, -0.3538), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.2883, -0.25, 0.0579), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.1085, -0.25, 0.0611), -1.025, 1.0, 0, 1, 0, 0)},
                          int3(2, 2, 2));
  Component c = Component(
      0, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
      {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId,
       // uint8_t groupId
       Atom(double3(0.0, 0.0, 1.149), -0.3256, 1.0, 0, 4, 1, 0), Atom(double3(0.0, 0.0, 0.0), 0.6512, 1.0, 0, 3, 1, 0),
       Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)},
      5, 21);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {2}, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::vector<Atom> atomPositions = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end());

  atomPositions[0].position = double3(10.011, 4.97475 + 2.0, 1.149);
  atomPositions[1].position = double3(10.011, 4.97475 + 2.0, 0.0);
  atomPositions[2].position = double3(10.011, 4.97475 + 2.0, -1.149);
  atomPositions[3].position = double3(10.011, 4.97475 - 2.0, 1.149);
  atomPositions[4].position = double3(10.011, 4.97475 - 2.0, 0.0);
  atomPositions[5].position = double3(10.011, 4.97475 - 2.0, -1.149);

  // atomPositions[8].groupId = 1;
  atomPositions[2].groupId = 1;
  // atomPositions[12].groupId = 1;
  atomPositions[4].groupId = 1;
  // atomPositions[14].groupId = 1;

  // atomPositions[8].scalingCoulomb = 0.45;
  atomPositions[2].scalingCoulomb = 0.5;
  // atomPositions[12].scalingCoulomb = 0.4;
  atomPositions[4].scalingCoulomb = 0.6;
  // atomPositions[14].scalingCoulomb = 0.7;
  RunningEnergy factor =
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);

  double delta = 1e-5;
  double tolerance = 1e-3;

  // atomPositions[8].scalingCoulomb = 0.45 + 0.5 * delta;
  atomPositions[2].scalingCoulomb = 0.5 + 0.5 * delta;
  // atomPositions[12].scalingCoulomb = 0.4 + 0.5 * delta;
  atomPositions[4].scalingCoulomb = 0.6 + 0.5 * delta;
  // atomPositions[14].scalingCoulomb = 0.7 + 0.5 * delta;
  RunningEnergy energyForward =
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);

  // atomPositions[8].scalingCoulomb = 0.45 - 0.5 * delta;
  atomPositions[2].scalingCoulomb = 0.5 - 0.5 * delta;
  // atomPositions[12].scalingCoulomb = 0.4 - 0.5 * delta;
  atomPositions[4].scalingCoulomb = 0.6 - 0.5 * delta;
  // atomPositions[14].scalingCoulomb = 0.7 - 0.5 * delta;
  RunningEnergy energyBackward =
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions);

  double dUdlambda = (energyForward.potentialEnergy() - energyBackward.potentialEnergy()) / delta;

  EXPECT_NEAR(factor.dudlambdaCharge, dUdlambda, tolerance)
      << " ratio: " << factor.dudlambdaCharge / dUdlambda << " " << dUdlambda / factor.dudlambdaCharge;
}

TEST(dudlambda, Test_2_CO2_in_MFI_2x2x2_Ewald)
{
  ForceField forceField = ForceField(
      {
          PseudoAtom("Si", true, 28.0855, 2.05, 0.0, 14, false),
          PseudoAtom("O", true, 15.999, -1.025, 0.0, 8, false),
          PseudoAtom("CH4", false, 16.04246, 0.0, 0.0, 6, false),
          PseudoAtom("C_co2", false, 12.0, 0.6512, 0.0, 6, false),
          PseudoAtom("O_co2", false, 15.9994, -0.3256, 0.0, 8, false),
      },
      {VDWParameters(22.0, 2.30), VDWParameters(53.0, 3.3), VDWParameters(158.5, 3.72), VDWParameters(29.933, 2.745),
       VDWParameters(85.671, 3.017)},
      ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, true);
  Framework f = Framework(0, forceField, "MFI_SI", SimulationBox(20.022, 19.899, 13.383), 292,
                          {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                           // uint8_t componentId, uint8_t groupId
                           Atom(double3(0.42238, 0.0565, -0.33598), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.30716, 0.02772, -0.1893), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.27911, 0.06127, 0.0312), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.12215, 0.06298, 0.0267), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.07128, 0.02722, -0.18551), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.18641, 0.05896, -0.32818), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.42265, -0.1725, -0.32718), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.30778, -0.13016, -0.18548), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.27554, -0.17279, 0.03109), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.12058, -0.1731, 0.02979), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.07044, -0.13037, -0.182), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.18706, -0.17327, -0.31933), 2.05, 1.0, 0, 0, 0, 0),
                           Atom(double3(0.3726, 0.0534, -0.2442), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.3084, 0.0587, -0.0789), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.2007, 0.0592, 0.0289), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.0969, 0.0611, -0.0856), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.1149, 0.0541, -0.2763), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.2435, 0.0553, -0.246), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.3742, -0.1561, -0.2372), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.3085, -0.1552, -0.0728), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.198, -0.1554, 0.0288), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.091, -0.1614, -0.0777), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.1169, -0.1578, -0.2694), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.2448, -0.1594, -0.2422), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.3047, -0.051, -0.1866), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.0768, -0.0519, -0.1769), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.4161, 0.1276, -0.3896), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.4086, -0.0017, -0.4136), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.402, -0.1314, -0.4239), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.1886, 0.1298, -0.3836), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.194, 0.0007, -0.4082), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.1951, -0.1291, -0.419), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(-0.0037, 0.0502, -0.208), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(-0.004, -0.1528, -0.2078), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.4192, -0.25, -0.354), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.1884, -0.25, -0.3538), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.2883, -0.25, 0.0579), -1.025, 1.0, 0, 1, 0, 0),
                           Atom(double3(0.1085, -0.25, 0.0611), -1.025, 1.0, 0, 1, 0, 0)},
                          int3(2, 2, 2));
  Component c = Component(
      0, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
      {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId,
       // uint8_t groupId
       Atom(double3(0.0, 0.0, 1.149), -0.3256, 1.0, 0, 4, 1, 0), Atom(double3(0.0, 0.0, 0.0), 0.6512, 1.0, 0, 3, 1, 0),
       Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)},
      5, 21);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {2}, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::vector<Atom> atomPositions = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end());

  atomPositions[0].position = double3(10.011, 4.97475 + 2.0, 1.149);
  atomPositions[1].position = double3(10.011, 4.97475 + 2.0, 0.0);
  atomPositions[2].position = double3(10.011, 4.97475 + 2.0, -1.149);
  atomPositions[3].position = double3(10.011, 4.97475 - 2.0, 1.149);
  atomPositions[4].position = double3(10.011, 4.97475 - 2.0, 0.0);
  atomPositions[5].position = double3(10.011, 4.97475 - 2.0, -1.149);

  // atomPositions[8].groupId = 1;
  atomPositions[2].groupId = 1;
  // atomPositions[12].groupId = 1;
  atomPositions[4].groupId = 1;
  // atomPositions[14].groupId = 1;

  // atomPositions[8].scalingCoulomb = 0.45;
  atomPositions[2].scalingCoulomb = 0.5;
  // atomPositions[12].scalingCoulomb = 0.4;
  atomPositions[4].scalingCoulomb = 0.6;
  // atomPositions[14].scalingCoulomb = 0.7;

  system.precomputeTotalRigidEnergy();
  RunningEnergy factor = Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent, atomPositions);

  double delta = 1e-5;
  double tolerance = 1e-3;

  // atomPositions[8].scalingCoulomb = 0.45 + 0.5 * delta;
  atomPositions[2].scalingCoulomb = 0.5 + 0.5 * delta;
  // atomPositions[12].scalingCoulomb = 0.4 + 0.5 * delta;
  atomPositions[4].scalingCoulomb = 0.6 + 0.5 * delta;
  // atomPositions[14].scalingCoulomb = 0.7 + 0.5 * delta;
  RunningEnergy energyForward = Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent, atomPositions);

  // atomPositions[8].scalingCoulomb = 0.45 - 0.5 * delta;
  atomPositions[2].scalingCoulomb = 0.5 - 0.5 * delta;
  // atomPositions[12].scalingCoulomb = 0.4 - 0.5 * delta;
  atomPositions[4].scalingCoulomb = 0.6 - 0.5 * delta;
  // atomPositions[14].scalingCoulomb = 0.7 - 0.5 * delta;
  RunningEnergy energyBackward = Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent, atomPositions);

  double dUdlambda = (energyForward.potentialEnergy() - energyBackward.potentialEnergy()) / delta;

  EXPECT_NEAR(factor.dudlambdaEwald, dUdlambda, tolerance)
      << " ratio: " << factor.dudlambdaEwald / dUdlambda << " " << dUdlambda / factor.dudlambdaEwald;
}
