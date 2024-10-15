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
import energy_status;
import running_energy;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;
import energy_status;

TEST(energy_decomposition, CO2_Methane_in_Box)
{
  ForceField forceField = ForceField(
      {PseudoAtom("Si", true, 28.0855, 2.05, 0.0, 14, false), PseudoAtom("O", true, 15.999, -1.025, 0.0, 8, false),
       PseudoAtom("CH4", false, 16.04246, 0.0, 0.0, 6, false), PseudoAtom("C_co2", false, 12.0, 0.6512, 0.0, 6, false),
       PseudoAtom("O_co2", false, 15.9994, -0.3256, 0.0, 8, false)},
      {VDWParameters(22.0, 2.30), VDWParameters(53.0, 3.3), VDWParameters(158.5, 3.72), VDWParameters(29.933, 2.745),
       VDWParameters(85.671, 3.017)},
      ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);

  Component methane = Component(0, forceField, "methane", 190.564, 45599200, 0.01142,
                                {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 2, 0, 0)}, 5, 21);

  Component co2 = Component(
      1, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
      {// double3 position, double charge, double lambda, uint16_t type, uint8_t componentId, uint32_t moleculeId
       Atom(double3(0.0, 0.0, 1.149), -0.3256, 1.0, 0, 4, 1, 0), Atom(double3(0.0, 0.0, 0.0), 0.6512, 1.0, 0, 3, 1, 0),
       Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)},
      5, 21);

  System system =
      System(0, forceField, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, 1.0, {}, {methane, co2}, {15, 30}, 5);

  RunningEnergy energy = system.computeTotalEnergies();

  std::pair<EnergyStatus, double3x3> strainDerivative = Interactions::computeInterMolecularEnergyStrainDerivative(
      system.forceField, system.components, system.simulationBox, system.atomPositions);
  strainDerivative.first.sumTotal();

  EXPECT_NEAR(energy.moleculeMoleculeVDW + energy.moleculeMoleculeCharge, strainDerivative.first.totalEnergy.energy,
              1e-6);
}

TEST(energy_decomposition, CO2_Methane_in_Box_Ewald)
{
  ForceField forceField = ForceField(
      {PseudoAtom("Si", true, 28.0855, 2.05, 0.0, 14, false), PseudoAtom("O", true, 15.999, -1.025, 0.0, 8, false),
       PseudoAtom("CH4", false, 16.04246, 0.0, 0.0, 6, false), PseudoAtom("C_co2", false, 12.0, 0.6512, 0.0, 6, false),
       PseudoAtom("O_co2", false, 15.9994, -0.3256, 0.0, 8, false)},
      {VDWParameters(22.0, 2.30), VDWParameters(53.0, 3.3), VDWParameters(158.5, 3.72), VDWParameters(29.933, 2.745),
       VDWParameters(85.671, 3.017)},
      ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, true);

  Component methane = Component(0, forceField, "methane", 190.564, 45599200, 0.01142,
                                {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 2, 0, 0)}, 5, 21);

  Component co2 = Component(
      1, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
      {// double3 position, double charge, double lambda, uint16_t type, uint8_t componentId, uint32_t moleculeId
       Atom(double3(0.0, 0.0, 1.149), -0.3256, 1.0, 0, 4, 1, 0), Atom(double3(0.0, 0.0, 0.0), 0.6512, 1.0, 0, 3, 1, 0),
       Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)},
      5, 21);

  Component co2_2 = Component(
      2, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
      {// double3 position, double charge, double lambda, uint16_t type, uint8_t componentId, uint32_t moleculeId
       Atom(double3(0.0, 0.0, 1.149), -0.3256, 1.0, 0, 4, 1, 0), Atom(double3(0.0, 0.0, 0.0), 0.6512, 1.0, 0, 3, 1, 0),
       Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)},
      5, 21);

  System system =
      System(0, forceField, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, 1.0, {}, {co2, co2_2}, {15, 30}, 5);

  system.forceField.EwaldAlpha = 0.25;
  system.forceField.numberOfWaveVectors = int3(8, 8, 8);

  RunningEnergy energy = system.computeTotalEnergies();

  Interactions::computeEwaldFourierEnergySingleIon(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                   system.forceField, system.simulationBox, double3(0.0, 0.0, 0.0),
                                                   1.0);

  system.precomputeTotalRigidEnergy();
  std::pair<EnergyStatus, double3x3> strainDerivative = Interactions::computeEwaldFourierEnergyStrainDerivative(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, system.frameworkComponents, system.components,
      system.numberOfMoleculesPerComponent, system.spanOfMoleculeAtoms(), system.CoulombicFourierEnergySingleIon,
      system.netChargeFramework, system.netChargePerComponent);

  strainDerivative.first.sumTotal();

  EXPECT_NEAR(energy.ewald_fourier + energy.ewald_self + energy.ewald_exclusion,
              strainDerivative.first.totalEnergy.energy, 1e-6);
}

inline std::pair<EnergyStatus, double3x3> pair_acc(const std::pair<EnergyStatus, double3x3> &lhs,
                                                   const std::pair<EnergyStatus, double3x3> &rhs)
{
  return std::make_pair(lhs.first + rhs.first, lhs.second + rhs.second);
}

TEST(energy_decomposition, CO2_Methane_in_Framework)
{
  ForceField forceField = ForceField(
      {PseudoAtom("Si", true, 28.0855, 2.05, 0.0, 14, false), PseudoAtom("O", true, 15.999, -1.025, 0.0, 8, false),
       PseudoAtom("CH4", false, 16.04246, 0.0, 0.0, 6, false), PseudoAtom("C_co2", false, 12.0, 0.6512, 0.0, 6, false),
       PseudoAtom("O_co2", false, 15.9994, -0.3256, 0.0, 8, false)},
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

  Component methane = Component(0, forceField, "methane", 190.564, 45599200, 0.01142,
                                {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 2, 0, 0)}, 5, 21);

  Component co2 = Component(
      1, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
      {// double3 position, double charge, double lambda, uint16_t type, uint8_t componentId, uint32_t moleculeId
       Atom(double3(0.0, 0.0, 1.149), -0.3256, 1.0, 0, 4, 1, 0), Atom(double3(0.0, 0.0, 0.0), 0.6512, 1.0, 0, 3, 1, 0),
       Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)},
      5, 21);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {methane, co2}, {10, 15}, 5);

  system.precomputeTotalRigidEnergy();

  RunningEnergy energy = system.computeTotalEnergies();

  RunningEnergy energyForces = system.computeTotalGradients();

  std::pair<EnergyStatus, double3x3> strainDerivative = system.computeMolecularPressure();

  EXPECT_NEAR(energy.potentialEnergy(), strainDerivative.first.totalEnergy.energy, 1e-6);
  EXPECT_NEAR(energy.potentialEnergy(), energyForces.potentialEnergy(), 1e-6);
}
