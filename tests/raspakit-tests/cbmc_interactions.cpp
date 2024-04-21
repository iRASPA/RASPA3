#include <gtest/gtest.h>

#include <vector>
#include <tuple>
#include <algorithm>
#include <span>
#include <complex>

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
    ForceField forceField = ForceField(
    { PseudoAtom("Si", 28.0855, 2.05, 14, false),
      PseudoAtom("O", 15.999, -1.025, 8, false),
      PseudoAtom("CH4", 16.04246, 0.0, 6, false)
    },
    { VDWParameters(22.0, 2.30),
      VDWParameters(53.0, 3.3),
      VDWParameters(158.5, 3.72) },
    ForceField::MixingRule::Lorentz_Berthelot,
    12.0,
    true,
    false);
  Framework f = Framework(0, forceField, "ITQ-29", SimulationBox(11.8671, 11.8671, 11.8671),
    517,
    { // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
      Atom(double3(0.3683, 0.1847, 0),       2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.5,    0.2179, 0),      -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.2939, 0.2939, 0),      -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.3429, 0.1098, 0.1098), -1.025, 1.0, 0, 1, 0, 0)
    },
    int3(1, 1, 1));
  Component c = Component(0, forceField, "methane", 190.564, 45599200, 0.01142,
    { // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
      Atom(double3(0.0, 0.0,  0.0),    0.0, 1.0, 0, 2, 1, 0) 
    }, 5, 21);

  System system = System(0, std::nullopt, 300.0, 1e4, forceField, { f }, { c }, { 1 }, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  atomPositions[0].position = double3(0.5 * 11.8671, 0.5 * 11.8671, 0.5 * 11.8671);

  RunningEnergy energy;
  Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions, energy);
  Interactions::computeFrameworkMoleculeEnergy(system.forceField,system.simulationBox, frameworkAtoms, atomPositions, energy);

  EXPECT_NEAR(energy.frameworkMoleculeVDW * Units::EnergyToKelvin, -337.77056357, 1e-6);

  std::optional<RunningEnergy> frameworkMoleculeEnergy = 
    CBMC::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, frameworkAtoms, 12.0, 12.0, atomPositions, -1);

  EXPECT_NEAR(frameworkMoleculeEnergy->frameworkMoleculeVDW * Units::EnergyToKelvin, -337.77056357, 1e-6);
}

TEST(cbmc_interactions, framework_molecule_2)
{
    ForceField forceField = ForceField(
    { PseudoAtom("Si", 28.0855, 2.05, 14, false),
      PseudoAtom("O", 15.999, -1.025, 8, false),
      PseudoAtom("CH4", 16.04246, 0.0, 6, false)
    },
    { VDWParameters(22.0, 2.30),
      VDWParameters(53.0, 3.3),
      VDWParameters(158.5, 3.72) },
    ForceField::MixingRule::Lorentz_Berthelot,
    12.0,
    true,
    false);
  Framework f = Framework(0, forceField, "ITQ-29", SimulationBox(11.8671, 11.8671, 11.8671),
    517,
    { // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
      Atom(double3(0.3683, 0.1847, 0),       2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.5,    0.2179, 0),      -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.2939, 0.2939, 0),      -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.3429, 0.1098, 0.1098), -1.025, 1.0, 0, 1, 0, 0)
    },
    int3(1, 1, 1));
  Component c = Component(0, forceField, "methane", 190.564, 45599200, 0.01142,
    { // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
      Atom(double3(0.0, 0.0,  0.0),    0.0, 1.0, 0, 2, 1, 0) 
    }, 5, 21);

  System system = System(0, std::nullopt, 300.0, 1e4, forceField, { f }, { c }, { 2 }, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  atomPositions[0].position = double3(0.5 * 11.8671, 0.5 * 11.8671, 0.5 * 11.8671);
  atomPositions[1].position = double3(0.6 * 11.8671, 0.7 * 11.8671, 0.65 * 11.8671);

  RunningEnergy energy;
  Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions, energy);
  Interactions::computeFrameworkMoleculeEnergy(system.forceField,system.simulationBox, frameworkAtoms, atomPositions, energy);

  EXPECT_NEAR(energy.frameworkMoleculeVDW * Units::EnergyToKelvin, -1599.10322574, 1e-6);
  EXPECT_NEAR(energy.moleculeMoleculeVDW * Units::EnergyToKelvin, 2352.42793591, 1e-6);


  std::optional<RunningEnergy> frameworkMoleculeEnergy = 
    CBMC::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, frameworkAtoms, 12.0, 12.0, atomPositions, -1);

  std::optional<RunningEnergy> interMoleculeEnergy1 = 
    CBMC::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions, 12.0, 12.0, {atomPositions.begin(), 1}, -1);

  std::optional<RunningEnergy> interMoleculeEnergy2 = 
    CBMC::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions, 12.0, 12.0, {atomPositions.begin()+1, 1}, -1);

  std::optional<RunningEnergy> interMoleculeEnergy3 = 
    CBMC::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions, 12.0, 12.0, {atomPositions.begin(), 2}, -1);

  EXPECT_NEAR(frameworkMoleculeEnergy->frameworkMoleculeVDW * Units::EnergyToKelvin, -1599.10322574, 1e-6);
  EXPECT_NEAR(interMoleculeEnergy1->moleculeMoleculeVDW * Units::EnergyToKelvin, 2352.42793591, 1e-6);
  EXPECT_NEAR(interMoleculeEnergy2->moleculeMoleculeVDW * Units::EnergyToKelvin, 2352.42793591, 1e-6);
  EXPECT_NEAR(interMoleculeEnergy3->moleculeMoleculeVDW * Units::EnergyToKelvin / 2.0, 2352.42793591, 1e-6);
}
