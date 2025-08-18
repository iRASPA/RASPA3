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
import energy_status;
import running_energy;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;
import energy_status;
import integrators;
import integrators_update;
import integrators_compute;
import interpolation_energy_grid;

TEST(energy_decomposition, CO2_Methane_in_Box)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);

  Component methane = TestFactories::makeMethane(forceField, 0);
  Component co2 = TestFactories::makeCO2(forceField, 1, true);

  System system =
      System(0, forceField, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, 1.0, {}, {methane, co2}, {}, {15, 30}, 5);

  RunningEnergy energy = system.computeTotalEnergies();

  std::pair<EnergyStatus, double3x3> strainDerivative = Interactions::computeInterMolecularEnergyStrainDerivative(
      system.forceField, system.components, system.simulationBox, system.atomData);
  strainDerivative.first.sumTotal();

  EXPECT_NEAR(energy.moleculeMoleculeVDW + energy.moleculeMoleculeCharge, strainDerivative.first.totalEnergy.energy,
              1e-6);
}

TEST(energy_decomposition, CO2_Methane_in_Box_Ewald)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);

  Component methane = TestFactories::makeMethane(forceField, 0);
  Component co2 = TestFactories::makeCO2(forceField, 1, true);
  Component co2_2 = TestFactories::makeCO2(forceField, 2, true);

  System system =
      System(0, forceField, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, 1.0, {}, {co2, co2_2}, {}, {15, 30}, 5);

  system.forceField.EwaldAlpha = 0.25;
  system.forceField.numberOfWaveVectors = int3(8, 8, 8);

  RunningEnergy energy = system.computeTotalEnergies();

  Interactions::computeEwaldFourierEnergySingleIon(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                   system.forceField, system.simulationBox, double3(0.0, 0.0, 0.0),
                                                   1.0);

  system.precomputeTotalRigidEnergy();
  std::pair<EnergyStatus, double3x3> strainDerivative = Interactions::computeEwaldFourierEnergyStrainDerivative(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, system.framework, system.components,
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
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Framework f = TestFactories::makeMFI_Si(forceField, int3(2, 2, 2));
  Component methane = TestFactories::makeMethane(forceField, 0);
  Component co2 = TestFactories::makeCO2(forceField, 1, true);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {methane, co2}, {}, {10, 15}, 5);

  system.precomputeTotalRigidEnergy();

  RunningEnergy energy = system.computeTotalEnergies();

  RunningEnergy energyForces = Integrators::updateGradients(
      system.spanOfMoleculeAtoms(), system.spanOfFrameworkAtoms(), system.forceField, system.simulationBox,
      system.components, system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik,
      system.fixedFrameworkStoredEik, system.interpolationGrids, system.numberOfMoleculesPerComponent);

  std::pair<EnergyStatus, double3x3> strainDerivative = system.computeMolecularPressure();

  EXPECT_NEAR(energy.potentialEnergy(), strainDerivative.first.totalEnergy.energy, 1e-6);
  EXPECT_NEAR(energy.potentialEnergy(), energyForces.potentialEnergy(), 1e-6);
}
