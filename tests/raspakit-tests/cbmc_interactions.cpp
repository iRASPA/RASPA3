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
import cbmc_interactions;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;
import cbmc_interactions;
import cbmc_interactions_external_field;
import cbmc_interactions_framework_molecule;
import cbmc_interactions_intermolecular;

TEST(cbmc_interactions, framework_molecule_1)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, false);
  Framework f = TestFactories::makeITQ29(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeMethane(forceField, 0);
  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {1}, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  atomPositions[0].position = double3(0.5 * 11.8671, 0.5 * 11.8671, 0.5 * 11.8671);

  RunningEnergy energy =
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions) +
      Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                   system.framework, frameworkAtoms, atomPositions);

  EXPECT_NEAR(energy.frameworkMoleculeVDW * Units::EnergyToKelvin, -337.77056357, 1e-6);

  std::optional<RunningEnergy> frameworkMoleculeEnergy =
      CBMC::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                           system.framework, frameworkAtoms, 12.0, 12.0, atomPositions, -1);

  EXPECT_NEAR(frameworkMoleculeEnergy->frameworkMoleculeVDW * Units::EnergyToKelvin, -337.77056357, 1e-6);
}

TEST(cbmc_interactions, framework_molecule_2)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, false);
  Framework f = TestFactories::makeITQ29(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeMethane(forceField, 0);
  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {2}, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  atomPositions[0].position = double3(0.5 * 11.8671, 0.5 * 11.8671, 0.5 * 11.8671);
  atomPositions[1].position = double3(0.6 * 11.8671, 0.7 * 11.8671, 0.65 * 11.8671);

  RunningEnergy energy =
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions) +
      Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                   system.framework, frameworkAtoms, atomPositions);

  EXPECT_NEAR(energy.frameworkMoleculeVDW * Units::EnergyToKelvin, -1599.10322574, 1e-6);
  EXPECT_NEAR(energy.moleculeMoleculeVDW * Units::EnergyToKelvin, 2352.42793591, 1e-6);

  std::optional<RunningEnergy> frameworkMoleculeEnergy =
      CBMC::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                           system.framework, frameworkAtoms, 12.0, 12.0, atomPositions, -1);

  std::optional<RunningEnergy> interMoleculeEnergy1 = CBMC::computeInterMolecularEnergy(
      system.forceField, system.simulationBox, atomPositions, 12.0, 12.0, {atomPositions.begin(), 1}, -1);

  std::optional<RunningEnergy> interMoleculeEnergy2 = CBMC::computeInterMolecularEnergy(
      system.forceField, system.simulationBox, atomPositions, 12.0, 12.0, {atomPositions.begin() + 1, 1}, -1);

  std::optional<RunningEnergy> interMoleculeEnergy3 = CBMC::computeInterMolecularEnergy(
      system.forceField, system.simulationBox, atomPositions, 12.0, 12.0, {atomPositions.begin(), 2}, -1);

  EXPECT_NEAR(frameworkMoleculeEnergy->frameworkMoleculeVDW * Units::EnergyToKelvin, -1599.10322574, 1e-6);
  EXPECT_NEAR(interMoleculeEnergy1->moleculeMoleculeVDW * Units::EnergyToKelvin, 2352.42793591, 1e-6);
  EXPECT_NEAR(interMoleculeEnergy2->moleculeMoleculeVDW * Units::EnergyToKelvin, 2352.42793591, 1e-6);
  EXPECT_NEAR(interMoleculeEnergy3->moleculeMoleculeVDW * Units::EnergyToKelvin / 2.0, 2352.42793591, 1e-6);
}
