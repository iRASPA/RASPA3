#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <cstddef>
#include <span>
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
import gradient_factor;
import energy_status;
import units;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;

TEST(Ewald, Test_2_CO2_in_Box_10_10_10)
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

  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

  Component c = Component(
      0, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
      {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId,
       // uint8_t groupId
       Atom(double3(0.0, 0.0, 1.149), -0.3256, 1.0, 0, 4, 0, 0), Atom(double3(0.0, 0.0, 0.0), 0.6512, 1.0, 0, 3, 0, 0),
       Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 0, 0)},
      5, 21);

  System system = System(0, forceField, SimulationBox(10.0, 10.0, 10.0), 300.0, 1e4, 1.0, {}, {c}, {2}, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  atomPositions[0].position = double3(-1.0, 0.0, 1.149);
  atomPositions[1].position = double3(-1.0, 0.0, 0.0);
  atomPositions[2].position = double3(-1.0, 0.0, -1.149);
  atomPositions[3].position = double3(1.0, 0.0, 1.149);
  atomPositions[4].position = double3(1.0, 0.0, 0.0);
  atomPositions[5].position = double3(1.0, 0.0, -1.149);

  Interactions::computeEwaldFourierEnergySingleIon(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                   *system.forceField, *system.simulationBox, double3(0.0, 0.0, 0.0),
                                                   1.0);
  system.precomputeTotalRigidEnergy();
  RunningEnergy energy = Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      *system.forceField, *system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
      system.spanOfMoleculeAtoms());

  EXPECT_NEAR(energy.ewaldFourier() * Units::EnergyToKelvin, 90.54613836, 1e-6);
}

TEST(Ewald, Test_1_Na_1_Cl_in_Box_10_10_10_Gradient)
{
  ForceField forceField = ForceField(
      {
          PseudoAtom("Si", true, 28.0855, 2.05, 0.0, 14, false),
          PseudoAtom("O", true, 15.999, -1.025, 0.0, 8, false),
          PseudoAtom("CH4", false, 16.04246, 0.0, 0.0, 6, false),
          PseudoAtom("Na+", false, 12.0, 1.0, 0.0, 6, false),
          PseudoAtom("Cl-", false, 15.9994, -1.0, 0.0, 8, false),
      },
      {VDWParameters(22.0, 2.30), VDWParameters(53.0, 3.3), VDWParameters(158.5, 3.72), VDWParameters(75.0, 1.0),
       VDWParameters(75.0, 1.0)},
      ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, true);

  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

  Component na = Component(0, forceField, "Na", 304.1282, 7377300.0, 0.22394,
                           {
                               Atom(double3(0.0, 0.0, 0.0), 1.0, 1.0, 0, 3, 0, 0),
                           },
                           5, 21);
  Component cl = Component(1, forceField, "Cl", 304.1282, 7377300.0, 0.22394,
                           {
                               Atom(double3(0.0, 0.0, 0.0), -1.0, 1.0, 0, 4, 1, 0),
                           },
                           5, 21);

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
      *system.forceField, *system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
      system.spanOfMoleculeAtoms(), system.CoulombicFourierEnergySingleIon, system.netChargeFramework,
      system.netChargePerComponent);

  double delta = 1e-5;
  double tolerance = 1e-5;
  double3 gradient;
  for (size_t i = 0; i < atomPositions.size(); ++i)
  {
    RunningEnergy x1, x2, y1, y2, z1, z2;

    // finite difference x
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x + 0.5 * delta;
    x2 = Interactions::computeEwaldFourierEnergy(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                 system.fixedFrameworkStoredEik, system.storedEik, *system.forceField,
                                                 *system.simulationBox, system.components,
                                                 system.numberOfMoleculesPerComponent, atomPositions);

    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x - 0.5 * delta;
    x1 = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        *system.forceField, *system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        std::span<const Atom>(atomPositions));
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x;

    // finite difference y
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y + 0.5 * delta;
    y2 = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        *system.forceField, *system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        std::span<const Atom>(atomPositions));

    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y - 0.5 * delta;
    y1 = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        *system.forceField, *system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        std::span<const Atom>(atomPositions));
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y;

    // finite difference z
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z + 0.5 * delta;
    z2 = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        *system.forceField, *system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        std::span<const Atom>(atomPositions));

    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z - 0.5 * delta;
    z1 = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        *system.forceField, *system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
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

  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

  Framework f = Framework(
      0, forceField, "ITQ-29", SimulationBox(11.8671, 11.8671, 11.8671), 517,
      {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId,
       // uint8_t groupId
       Atom(double3(0.3683, 0.1847, 0), 2.05, 1.0, 0, 0, 0, 0), Atom(double3(0.5, 0.2179, 0), -1.025, 1.0, 0, 1, 0, 0),
       Atom(double3(0.2939, 0.2939, 0), -1.025, 1.0, 0, 1, 0, 0),
       Atom(double3(0.3429, 0.1098, 0.1098), -1.025, 1.0, 0, 1, 0, 0)},
      int3(1, 1, 1));
  Component c = Component(
      0, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
      {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId,
       // uint8_t groupId
       Atom(double3(0.0, 0.0, 1.149), -0.3256, 1.0, 0, 4, 1, 0), Atom(double3(0.0, 0.0, 0.0), 0.6512, 1.0, 0, 3, 1, 0),
       Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)},
      5, 21);

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
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, *system.forceField, *system.simulationBox,
      double3(0.0, 0.0, 0.0), 1.0);
  RunningEnergy energy = Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      *system.forceField, *system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
      system.spanOfMoleculeAtoms());

  EXPECT_NEAR(system.CoulombicFourierEnergySingleIon * Units::EnergyToKelvin, 17464.2371790130, 1e-6);
  EXPECT_NEAR(energy.ewaldFourier() * Units::EnergyToKelvin, -702.65863478, 1e-6);
}

TEST(Ewald, Test_2_CO2_in_ITQ_29_2x2x2)
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

  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

  Framework f = Framework(
      0, forceField, "ITQ-29", SimulationBox(11.8671, 11.8671, 11.8671), 517,
      {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId,
       // uint8_t groupId
       Atom(double3(0.3683, 0.1847, 0), 2.05, 1.0, 0, 0, 0, 0), Atom(double3(0.5, 0.2179, 0), -1.025, 1.0, 0, 1, 0, 0),
       Atom(double3(0.2939, 0.2939, 0), -1.025, 1.0, 0, 1, 0, 0),
       Atom(double3(0.3429, 0.1098, 0.1098), -1.025, 1.0, 0, 1, 0, 0)},
      int3(2, 2, 2));
  Component c = Component(
      0, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
      {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId,
       // uint8_t groupId
       Atom(double3(0.0, 0.0, 1.149), -0.3256, 1.0, 0, 4, 1, 0), Atom(double3(0.0, 0.0, 0.0), 0.6512, 1.0, 0, 3, 1, 0),
       Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)},
      5, 21);

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
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, *system.forceField, *system.simulationBox,
      double3(0.0, 0.0, 0.0), 1.0);
  RunningEnergy energy = Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      *system.forceField, *system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
      system.spanOfMoleculeAtoms());

  EXPECT_NEAR(system.CoulombicFourierEnergySingleIon * Units::EnergyToKelvin, 9673.9032373025, 1e-6);
  EXPECT_NEAR(energy.ewaldFourier() * Units::EnergyToKelvin, -721.64644486, 1e-6);
}

TEST(Ewald, Test_2_CO2_in_MFI_1x1x1)
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

  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

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
                          int3(1, 1, 1));
  Component c = Component(
      0, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
      {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId,
       // uint8_t groupId
       Atom(double3(0.0, 0.0, 1.149), -0.3256, 1.0, 0, 4, 1, 0), Atom(double3(0.0, 0.0, 0.0), 0.6512, 1.0, 0, 3, 1, 0),
       Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)},
      5, 21);

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
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, *system.forceField, *system.simulationBox,
      double3(0.0, 0.0, 0.0), 1.0);
  RunningEnergy energy = Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      *system.forceField, *system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
      system.spanOfMoleculeAtoms());

  EXPECT_NEAR(system.CoulombicFourierEnergySingleIon * Units::EnergyToKelvin, 12028.1731827280, 1e-6);
  EXPECT_NEAR(energy.ewaldFourier() * Units::EnergyToKelvin, -1191.77790165, 1e-6);
}

TEST(Ewald, Test_2_CO2_in_MFI_2x2x2)
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

  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

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

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  atomPositions[0].position = double3(10.011, 4.97475 + 2.0, 1.149);
  atomPositions[1].position = double3(10.011, 4.97475 + 2.0, 0.0);
  atomPositions[2].position = double3(10.011, 4.97475 + 2.0, -1.149);
  atomPositions[3].position = double3(10.011, 4.97475 - 2.0, 1.149);
  atomPositions[4].position = double3(10.011, 4.97475 - 2.0, 0.0);
  atomPositions[5].position = double3(10.011, 4.97475 - 2.0, -1.149);

  system.precomputeTotalRigidEnergy();
  system.CoulombicFourierEnergySingleIon = Interactions::computeEwaldFourierEnergySingleIon(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, *system.forceField, *system.simulationBox,
      double3(0.0, 0.0, 0.0), 1.0);
  RunningEnergy energy = Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      *system.forceField, *system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
      system.spanOfMoleculeAtoms());

  EXPECT_NEAR(system.CoulombicFourierEnergySingleIon * Units::EnergyToKelvin, 6309.7866899037, 1e-6);
  EXPECT_NEAR(energy.ewaldFourier() * Units::EnergyToKelvin, -1197.23909965, 1e-6);
}

TEST(Ewald, Test_20_Na_Cl_in_Box_25x25x25)
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

  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

  Component na = Component(0, forceField, "Na", 304.1282, 7377300.0, 0.22394,
                           {
                               // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                               // uint8_t componentId, uint8_t groupId
                               Atom(double3(0.0, 0.0, 0.0), 1.0, 1.0, 0, 3, 0, 0),
                           },
                           5, 21);
  Component cl = Component(1, forceField, "Cl", 304.1282, 7377300.0, 0.22394,
                           {
                               // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                               // uint8_t componentId, uint8_t groupId
                               Atom(double3(0.0, 0.0, 0.0), -1.0, 1.0, 0, 4, 1, 0),
                           },
                           5, 21);

  System system = System(0, forceField, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, 1.0, {}, {na, cl}, {20, 20}, 5);

  // std::fill(system.forceField->data.begin(), system.forceField->data.end(), VDWParameters(0.0, 1.0));

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
      *system.forceField, *system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
      atomPositions);

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
        *system.forceField, *system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        system.spanOfMoleculeAtoms());
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x + 0.5 * delta;
    forward1_x = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        *system.forceField, *system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        system.spanOfMoleculeAtoms());
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x - 0.5 * delta;
    backward1_x = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        *system.forceField, *system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        system.spanOfMoleculeAtoms());
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x - 1.0 * delta;
    backward2_x = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        *system.forceField, *system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        system.spanOfMoleculeAtoms());
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x;

    // finite difference y
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y + 1.0 * delta;
    forward2_y = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        *system.forceField, *system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        system.spanOfMoleculeAtoms());
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y + 0.5 * delta;
    forward1_y = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        *system.forceField, *system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        system.spanOfMoleculeAtoms());
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y - 0.5 * delta;
    backward1_y = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        *system.forceField, *system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        system.spanOfMoleculeAtoms());
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y - 1.0 * delta;
    backward2_y = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        *system.forceField, *system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        system.spanOfMoleculeAtoms());
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y;

    // finite difference z
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z + 1.0 * delta;
    forward2_z = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        *system.forceField, *system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        system.spanOfMoleculeAtoms());
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z + 0.5 * delta;
    forward1_z = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        *system.forceField, *system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        system.spanOfMoleculeAtoms());
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z - 0.5 * delta;
    backward1_z = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        *system.forceField, *system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
        system.spanOfMoleculeAtoms());
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z - 1.0 * delta;
    backward2_z = Interactions::computeEwaldFourierEnergy(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        *system.forceField, *system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
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
