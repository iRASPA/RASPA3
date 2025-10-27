#ifdef USE_LEGACY_HEADERS
#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <cstddef>
#include <span>
#include <vector>
#endif

#ifdef USE_STD_IMPORT
#include <gtest/gtest.h>
import std;
#endif

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
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;

TEST(spacegroup, TestLennardJonesVDWTwoMethanes)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, false, false, false);
  Component c = TestFactories::makeMethane(forceField, 0);

  System system = System(0, forceField, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, 1.0, {}, {c}, {}, {2}, 5);

  system.atomData[0].position = double3(0.0, 0.0, 0.0);
  system.atomData[1].position = double3(0.0, 0.0, 3.72 * std::pow(2.0, 1.0 / 6.0));

  RunningEnergy energy =
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, system.atomData);

  EXPECT_NEAR(energy.moleculeMoleculeVDW * Units::EnergyToKelvin, -158.5, 1e-12);
}

TEST(spacegroup, TestLennardJonesVDWMethaneInITQ_29_P1)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, false);
  Component c = TestFactories::makeMethane(forceField, 0);

  Framework f =
      Framework(0, forceField, "ITQ-29", SimulationBox(11.8671, 11.8671, 11.8671), 1,
                {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                 // uint8_t componentId, uint8_t groupId
                 Atom(double3(0.368300000000, 0.184700000000, 0.000000000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.500000000000, 0.217900000000, 0.000000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.293900000000, 0.293900000000, 0.000000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.342900000000, 0.109800000000, 0.109800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.815300000000, 0.368300000000, 0.000000000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.782100000000, 0.500000000000, 0.000000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.706100000000, 0.293900000000, 0.000000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.890200000000, 0.342900000000, 0.109800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.631700000000, 0.815300000000, 0.000000000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.500000000000, 0.782100000000, 0.000000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.706100000000, 0.706100000000, 0.000000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.657100000000, 0.890200000000, 0.109800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.184700000000, 0.631700000000, 0.000000000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.217900000000, 0.500000000000, 0.000000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.293900000000, 0.706100000000, 0.000000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.109800000000, 0.657100000000, 0.109800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.368300000000, 0.815300000000, 0.000000000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.342900000000, 0.890200000000, 0.890200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.184700000000, 0.368300000000, 0.000000000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.109800000000, 0.342900000000, 0.890200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.631700000000, 0.184700000000, 0.000000000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.657100000000, 0.109800000000, 0.890200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.815300000000, 0.631700000000, 0.000000000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.890200000000, 0.657100000000, 0.890200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.000000000000, 0.368300000000, 0.184700000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.000000000000, 0.500000000000, 0.217900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.000000000000, 0.293900000000, 0.293900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.109800000000, 0.342900000000, 0.109800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.631700000000, 0.000000000000, 0.184700000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.500000000000, 0.000000000000, 0.217900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.706100000000, 0.000000000000, 0.293900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.657100000000, 0.109800000000, 0.109800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.000000000000, 0.631700000000, 0.184700000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.000000000000, 0.706100000000, 0.293900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.890200000000, 0.657100000000, 0.109800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.368300000000, 0.000000000000, 0.184700000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.293900000000, 0.000000000000, 0.293900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.342900000000, 0.890200000000, 0.109800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.000000000000, 0.631700000000, 0.815300000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.000000000000, 0.500000000000, 0.782100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.000000000000, 0.706100000000, 0.706100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.109800000000, 0.657100000000, 0.890200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.368300000000, 0.000000000000, 0.815300000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.500000000000, 0.000000000000, 0.782100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.293900000000, 0.000000000000, 0.706100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.342900000000, 0.109800000000, 0.890200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.000000000000, 0.368300000000, 0.815300000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.000000000000, 0.293900000000, 0.706100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.890200000000, 0.342900000000, 0.890200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.631700000000, 0.000000000000, 0.815300000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.706100000000, 0.000000000000, 0.706100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.657100000000, 0.890200000000, 0.890200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.184700000000, 0.000000000000, 0.368300000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.217900000000, 0.000000000000, 0.500000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.109800000000, 0.109800000000, 0.342900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.184700000000, 0.000000000000, 0.631700000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.109800000000, 0.890200000000, 0.657100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.000000000000, 0.184700000000, 0.631700000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.000000000000, 0.217900000000, 0.500000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.109800000000, 0.109800000000, 0.657100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.815300000000, 0.000000000000, 0.631700000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.782100000000, 0.000000000000, 0.500000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.890200000000, 0.109800000000, 0.657100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.000000000000, 0.815300000000, 0.631700000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.000000000000, 0.782100000000, 0.500000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.890200000000, 0.890200000000, 0.657100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.815300000000, 0.000000000000, 0.368300000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.890200000000, 0.890200000000, 0.342900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.000000000000, 0.815300000000, 0.368300000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.109800000000, 0.890200000000, 0.342900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.000000000000, 0.184700000000, 0.368300000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.890200000000, 0.109800000000, 0.342900000000), -1.025, 1.0, 0, 1, 0, false, false)},
                int3(1, 1, 1));

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {}, {1}, 5);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  atomData[0].position = double3(0.5 * 11.8671, 0.5 * 11.8671, 0.5 * 11.8671);

  RunningEnergy energy = Interactions::computeFrameworkMoleculeEnergy(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework, frameworkAtoms, atomData);

  EXPECT_NEAR(energy.frameworkMoleculeVDW * Units::EnergyToKelvin, -337.77056357, 1e-6);
}

TEST(spacegroup, TestLennardJonesVDWMethaneInITQ_29_2x2x2_P1)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, false);
  Component c = TestFactories::makeMethane(forceField, 0);

  Framework f =
      Framework(0, forceField, "ITQ-29", SimulationBox(11.8671, 11.8671, 11.8671), 1,
                {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                 // uint8_t componentId, uint8_t groupId
                 Atom(double3(0.368300000000, 0.184700000000, 0.000000000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.500000000000, 0.217900000000, 0.000000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.293900000000, 0.293900000000, 0.000000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.342900000000, 0.109800000000, 0.109800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.815300000000, 0.368300000000, 0.000000000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.782100000000, 0.500000000000, 0.000000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.706100000000, 0.293900000000, 0.000000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.890200000000, 0.342900000000, 0.109800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.631700000000, 0.815300000000, 0.000000000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.500000000000, 0.782100000000, 0.000000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.706100000000, 0.706100000000, 0.000000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.657100000000, 0.890200000000, 0.109800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.184700000000, 0.631700000000, 0.000000000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.217900000000, 0.500000000000, 0.000000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.293900000000, 0.706100000000, 0.000000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.109800000000, 0.657100000000, 0.109800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.368300000000, 0.815300000000, 0.000000000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.342900000000, 0.890200000000, 0.890200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.184700000000, 0.368300000000, 0.000000000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.109800000000, 0.342900000000, 0.890200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.631700000000, 0.184700000000, 0.000000000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.657100000000, 0.109800000000, 0.890200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.815300000000, 0.631700000000, 0.000000000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.890200000000, 0.657100000000, 0.890200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.000000000000, 0.368300000000, 0.184700000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.000000000000, 0.500000000000, 0.217900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.000000000000, 0.293900000000, 0.293900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.109800000000, 0.342900000000, 0.109800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.631700000000, 0.000000000000, 0.184700000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.500000000000, 0.000000000000, 0.217900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.706100000000, 0.000000000000, 0.293900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.657100000000, 0.109800000000, 0.109800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.000000000000, 0.631700000000, 0.184700000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.000000000000, 0.706100000000, 0.293900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.890200000000, 0.657100000000, 0.109800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.368300000000, 0.000000000000, 0.184700000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.293900000000, 0.000000000000, 0.293900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.342900000000, 0.890200000000, 0.109800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.000000000000, 0.631700000000, 0.815300000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.000000000000, 0.500000000000, 0.782100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.000000000000, 0.706100000000, 0.706100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.109800000000, 0.657100000000, 0.890200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.368300000000, 0.000000000000, 0.815300000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.500000000000, 0.000000000000, 0.782100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.293900000000, 0.000000000000, 0.706100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.342900000000, 0.109800000000, 0.890200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.000000000000, 0.368300000000, 0.815300000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.000000000000, 0.293900000000, 0.706100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.890200000000, 0.342900000000, 0.890200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.631700000000, 0.000000000000, 0.815300000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.706100000000, 0.000000000000, 0.706100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.657100000000, 0.890200000000, 0.890200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.184700000000, 0.000000000000, 0.368300000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.217900000000, 0.000000000000, 0.500000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.109800000000, 0.109800000000, 0.342900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.184700000000, 0.000000000000, 0.631700000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.109800000000, 0.890200000000, 0.657100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.000000000000, 0.184700000000, 0.631700000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.000000000000, 0.217900000000, 0.500000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.109800000000, 0.109800000000, 0.657100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.815300000000, 0.000000000000, 0.631700000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.782100000000, 0.000000000000, 0.500000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.890200000000, 0.109800000000, 0.657100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.000000000000, 0.815300000000, 0.631700000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.000000000000, 0.782100000000, 0.500000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.890200000000, 0.890200000000, 0.657100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.815300000000, 0.000000000000, 0.368300000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.890200000000, 0.890200000000, 0.342900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.000000000000, 0.815300000000, 0.368300000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.109800000000, 0.890200000000, 0.342900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.000000000000, 0.184700000000, 0.368300000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.890200000000, 0.109800000000, 0.342900000000), -1.025, 1.0, 0, 1, 0, false, false)},
                int3(2, 2, 2));
  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {}, {1}, 5);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  atomData[0].position = double3(0.5 * 11.8671, 0.5 * 11.8671, 0.5 * 11.8671);

  RunningEnergy energy = Interactions::computeFrameworkMoleculeEnergy(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework, frameworkAtoms, atomData);

  EXPECT_NEAR(energy.frameworkMoleculeVDW * Units::EnergyToKelvin, -602.89568378, 1e-6);
}

TEST(spacegroup, TestLennardJonesVDWMethaneInITQ_29)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, false);
  Framework f = TestFactories::makeITQ29(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeMethane(forceField, 0);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {}, {1}, 5);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  atomData[0].position = double3(0.5 * 11.8671, 0.5 * 11.8671, 0.5 * 11.8671);

  RunningEnergy energy = Interactions::computeFrameworkMoleculeEnergy(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework, frameworkAtoms, atomData);

  EXPECT_NEAR(energy.frameworkMoleculeVDW * Units::EnergyToKelvin, -337.77056357, 1e-6);
}

TEST(spacegroup, TestLennardJonesVDWMethaneInMFI)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, false);
  Framework f = TestFactories::makeMFI_Si(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeMethane(forceField, 0);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {}, {1}, 5);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  atomData[0].position = double3(10.011, 4.097475, 0.0);

  RunningEnergy energy = Interactions::computeFrameworkMoleculeEnergy(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework, frameworkAtoms, atomData);

  EXPECT_NEAR(energy.frameworkMoleculeVDW * Units::EnergyToKelvin, -1784.82292180, 1e-6);
}

TEST(spacegroup, TestLennardJonesVDWMethaneInMFI2x2x2)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, false);
  Framework f = TestFactories::makeMFI_Si(forceField, int3(2, 2, 2));
  Component c = TestFactories::makeMethane(forceField, 0);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {}, {1}, 5);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  atomData[0].position = double3(10.011, 4.097475, 0.0);

  RunningEnergy energy = Interactions::computeFrameworkMoleculeEnergy(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework, frameworkAtoms, atomData);

  EXPECT_NEAR(energy.frameworkMoleculeVDW * Units::EnergyToKelvin, -1828.89015075, 1e-6);
}

TEST(spacegroup, TestLennardJonesVDWMethaneInMFI_P1)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, false);
  Component c = TestFactories::makeMethane(forceField, 0);

  Framework f =
      Framework(0, forceField, "MFI_SI", SimulationBox(20.022, 19.899, 13.383), 1,
                {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                 // uint8_t componentId, uint8_t groupId
                 Atom(double3(0.422380000000, 0.056500000000, 0.664020000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.307160000000, 0.027720000000, 0.810700000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.279110000000, 0.061270000000, 0.031200000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.122150000000, 0.062980000000, 0.026700000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.071280000000, 0.027220000000, 0.814490000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.186410000000, 0.058960000000, 0.671820000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.422650000000, 0.827500000000, 0.672820000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.307780000000, 0.869840000000, 0.814520000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.275540000000, 0.827210000000, 0.031090000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.120580000000, 0.826900000000, 0.029790000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.070440000000, 0.869630000000, 0.818000000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.187060000000, 0.826730000000, 0.680670000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.372600000000, 0.053400000000, 0.755800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.308400000000, 0.058700000000, 0.921100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.200700000000, 0.059200000000, 0.028900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.096900000000, 0.061100000000, 0.914400000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.114900000000, 0.054100000000, 0.723700000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.243500000000, 0.055300000000, 0.754000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.374200000000, 0.843900000000, 0.762800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.308500000000, 0.844800000000, 0.927200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.198000000000, 0.844600000000, 0.028800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.091000000000, 0.838600000000, 0.922300000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.116900000000, 0.842200000000, 0.730600000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.244800000000, 0.840600000000, 0.757800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.304700000000, 0.949000000000, 0.813400000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.076800000000, 0.948100000000, 0.823100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.416100000000, 0.127600000000, 0.610400000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.408600000000, 0.998300000000, 0.586400000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.402000000000, 0.868600000000, 0.576100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.188600000000, 0.129800000000, 0.616400000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.194000000000, 0.000700000000, 0.591800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.195100000000, 0.870900000000, 0.581000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.996300000000, 0.050200000000, 0.792000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.996000000000, 0.847200000000, 0.792200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.419200000000, 0.750000000000, 0.646000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.188400000000, 0.750000000000, 0.646200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.288300000000, 0.750000000000, 0.057900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.108500000000, 0.750000000000, 0.061100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.077620000000, 0.943500000000, 0.164020000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.192840000000, 0.972280000000, 0.310700000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.220890000000, 0.938730000000, 0.531200000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.377850000000, 0.937020000000, 0.526700000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.428720000000, 0.972780000000, 0.314490000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.313590000000, 0.941040000000, 0.171820000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.077350000000, 0.172500000000, 0.172820000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.192220000000, 0.130160000000, 0.314520000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.224460000000, 0.172790000000, 0.531090000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.379420000000, 0.173100000000, 0.529790000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.429560000000, 0.130370000000, 0.318000000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.312940000000, 0.173270000000, 0.180670000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.127400000000, 0.946600000000, 0.255800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.191600000000, 0.941300000000, 0.421100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.299300000000, 0.940800000000, 0.528900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.403100000000, 0.938900000000, 0.414400000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.385100000000, 0.945900000000, 0.223700000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.256500000000, 0.944700000000, 0.254000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.125800000000, 0.156100000000, 0.262800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.191500000000, 0.155200000000, 0.427200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.302000000000, 0.155400000000, 0.528800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.409000000000, 0.161400000000, 0.422300000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.383100000000, 0.157800000000, 0.230600000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.255200000000, 0.159400000000, 0.257800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.195300000000, 0.051000000000, 0.313400000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.423200000000, 0.051900000000, 0.323100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.083900000000, 0.872400000000, 0.110400000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.091400000000, 0.001700000000, 0.086400000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.098000000000, 0.131400000000, 0.076100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.311400000000, 0.870200000000, 0.116400000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.306000000000, 0.999300000000, 0.091800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.304900000000, 0.129100000000, 0.081000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.503700000000, 0.949800000000, 0.292000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.504000000000, 0.152800000000, 0.292200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.080800000000, 0.250000000000, 0.146000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.311600000000, 0.250000000000, 0.146200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.211700000000, 0.250000000000, 0.557900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.391500000000, 0.250000000000, 0.561100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.922380000000, 0.443500000000, 0.835980000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.807160000000, 0.472280000000, 0.689300000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.779110000000, 0.438730000000, 0.468800000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.622150000000, 0.437020000000, 0.473300000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.571280000000, 0.472780000000, 0.685510000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.686410000000, 0.441040000000, 0.828180000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.922650000000, 0.672500000000, 0.827180000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.807780000000, 0.630160000000, 0.685480000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.775540000000, 0.672790000000, 0.468910000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.620580000000, 0.673100000000, 0.470210000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.570440000000, 0.630370000000, 0.682000000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.687060000000, 0.673270000000, 0.819330000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.872600000000, 0.446600000000, 0.744200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.808400000000, 0.441300000000, 0.578900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.700700000000, 0.440800000000, 0.471100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.596900000000, 0.438900000000, 0.585600000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.614900000000, 0.445900000000, 0.776300000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.743500000000, 0.444700000000, 0.746000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.874200000000, 0.656100000000, 0.737200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.808500000000, 0.655200000000, 0.572800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.698000000000, 0.655400000000, 0.471200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.591000000000, 0.661400000000, 0.577700000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.616900000000, 0.657800000000, 0.769400000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.744800000000, 0.659400000000, 0.742200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.804700000000, 0.551000000000, 0.686600000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.576800000000, 0.551900000000, 0.676900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.916100000000, 0.372400000000, 0.889600000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.908600000000, 0.501700000000, 0.913600000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.902000000000, 0.631400000000, 0.923900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.688600000000, 0.370200000000, 0.883600000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.694000000000, 0.499300000000, 0.908200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.695100000000, 0.629100000000, 0.919000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.496300000000, 0.449800000000, 0.708000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.496000000000, 0.652800000000, 0.707800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.919200000000, 0.750000000000, 0.854000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.688400000000, 0.750000000000, 0.853800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.788300000000, 0.750000000000, 0.442100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.608500000000, 0.750000000000, 0.438900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.577620000000, 0.556500000000, 0.335980000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.692840000000, 0.527720000000, 0.189300000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.720890000000, 0.561270000000, 0.968800000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.877850000000, 0.562980000000, 0.973300000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.928720000000, 0.527220000000, 0.185510000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.813590000000, 0.558960000000, 0.328180000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.577350000000, 0.327500000000, 0.327180000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.692220000000, 0.369840000000, 0.185480000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.724460000000, 0.327210000000, 0.968910000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.879420000000, 0.326900000000, 0.970210000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.929560000000, 0.369630000000, 0.182000000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.812940000000, 0.326730000000, 0.319330000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.627400000000, 0.553400000000, 0.244200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.691600000000, 0.558700000000, 0.078900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.799300000000, 0.559200000000, 0.971100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.903100000000, 0.561100000000, 0.085600000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.885100000000, 0.554100000000, 0.276300000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.756500000000, 0.555300000000, 0.246000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.625800000000, 0.343900000000, 0.237200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.691500000000, 0.344800000000, 0.072800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.802000000000, 0.344600000000, 0.971200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.909000000000, 0.338600000000, 0.077700000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.883100000000, 0.342200000000, 0.269400000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.755200000000, 0.340600000000, 0.242200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.695300000000, 0.449000000000, 0.186600000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.923200000000, 0.448100000000, 0.176900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.583900000000, 0.627600000000, 0.389600000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.591400000000, 0.498300000000, 0.413600000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.598000000000, 0.368600000000, 0.423900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.811400000000, 0.629800000000, 0.383600000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.806000000000, 0.500700000000, 0.408200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.804900000000, 0.370900000000, 0.419000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.003700000000, 0.550200000000, 0.208000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.004000000000, 0.347200000000, 0.207800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.580800000000, 0.250000000000, 0.354000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.811600000000, 0.250000000000, 0.353800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.711700000000, 0.250000000000, 0.942100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.891500000000, 0.250000000000, 0.938900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.577620000000, 0.943500000000, 0.335980000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.692840000000, 0.972280000000, 0.189300000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.720890000000, 0.938730000000, 0.968800000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.877850000000, 0.937020000000, 0.973300000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.928720000000, 0.972780000000, 0.185510000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.813590000000, 0.941040000000, 0.328180000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.577350000000, 0.172500000000, 0.327180000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.692220000000, 0.130160000000, 0.185480000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.724460000000, 0.172790000000, 0.968910000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.879420000000, 0.173100000000, 0.970210000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.929560000000, 0.130370000000, 0.182000000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.812940000000, 0.173270000000, 0.319330000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.627400000000, 0.946600000000, 0.244200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.691600000000, 0.941300000000, 0.078900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.799300000000, 0.940800000000, 0.971100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.903100000000, 0.938900000000, 0.085600000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.885100000000, 0.945900000000, 0.276300000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.756500000000, 0.944700000000, 0.246000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.625800000000, 0.156100000000, 0.237200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.691500000000, 0.155200000000, 0.072800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.802000000000, 0.155400000000, 0.971200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.909000000000, 0.161400000000, 0.077700000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.883100000000, 0.157800000000, 0.269400000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.755200000000, 0.159400000000, 0.242200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.695300000000, 0.051000000000, 0.186600000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.923200000000, 0.051900000000, 0.176900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.583900000000, 0.872400000000, 0.389600000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.591400000000, 0.001700000000, 0.413600000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.598000000000, 0.131400000000, 0.423900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.811400000000, 0.870200000000, 0.383600000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.806000000000, 0.999300000000, 0.408200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.804900000000, 0.129100000000, 0.419000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.003700000000, 0.949800000000, 0.208000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.004000000000, 0.152800000000, 0.207800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.922380000000, 0.056500000000, 0.835980000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.807160000000, 0.027720000000, 0.689300000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.779110000000, 0.061270000000, 0.468800000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.622150000000, 0.062980000000, 0.473300000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.571280000000, 0.027220000000, 0.685510000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.686410000000, 0.058960000000, 0.828180000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.922650000000, 0.827500000000, 0.827180000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.807780000000, 0.869840000000, 0.685480000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.775540000000, 0.827210000000, 0.468910000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.620580000000, 0.826900000000, 0.470210000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.570440000000, 0.869630000000, 0.682000000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.687060000000, 0.826730000000, 0.819330000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.872600000000, 0.053400000000, 0.744200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.808400000000, 0.058700000000, 0.578900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.700700000000, 0.059200000000, 0.471100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.596900000000, 0.061100000000, 0.585600000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.614900000000, 0.054100000000, 0.776300000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.743500000000, 0.055300000000, 0.746000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.874200000000, 0.843900000000, 0.737200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.808500000000, 0.844800000000, 0.572800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.698000000000, 0.844600000000, 0.471200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.591000000000, 0.838600000000, 0.577700000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.616900000000, 0.842200000000, 0.769400000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.744800000000, 0.840600000000, 0.742200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.804700000000, 0.949000000000, 0.686600000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.576800000000, 0.948100000000, 0.676900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.916100000000, 0.127600000000, 0.889600000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.908600000000, 0.998300000000, 0.913600000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.902000000000, 0.868600000000, 0.923900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.688600000000, 0.129800000000, 0.883600000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.694000000000, 0.000700000000, 0.908200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.695100000000, 0.870900000000, 0.919000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.496300000000, 0.050200000000, 0.708000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.496000000000, 0.847200000000, 0.707800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.077620000000, 0.556500000000, 0.164020000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.192840000000, 0.527720000000, 0.310700000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.220890000000, 0.561270000000, 0.531200000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.377850000000, 0.562980000000, 0.526700000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.428720000000, 0.527220000000, 0.314490000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.313590000000, 0.558960000000, 0.171820000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.077350000000, 0.327500000000, 0.172820000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.192220000000, 0.369840000000, 0.314520000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.224460000000, 0.327210000000, 0.531090000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.379420000000, 0.326900000000, 0.529790000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.429560000000, 0.369630000000, 0.318000000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.312940000000, 0.326730000000, 0.180670000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.127400000000, 0.553400000000, 0.255800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.191600000000, 0.558700000000, 0.421100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.299300000000, 0.559200000000, 0.528900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.403100000000, 0.561100000000, 0.414400000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.385100000000, 0.554100000000, 0.223700000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.256500000000, 0.555300000000, 0.254000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.125800000000, 0.343900000000, 0.262800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.191500000000, 0.344800000000, 0.427200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.302000000000, 0.344600000000, 0.528800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.409000000000, 0.338600000000, 0.422300000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.383100000000, 0.342200000000, 0.230600000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.255200000000, 0.340600000000, 0.257800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.195300000000, 0.449000000000, 0.313400000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.423200000000, 0.448100000000, 0.323100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.083900000000, 0.627600000000, 0.110400000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.091400000000, 0.498300000000, 0.086400000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.098000000000, 0.368600000000, 0.076100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.311400000000, 0.629800000000, 0.116400000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.306000000000, 0.500700000000, 0.091800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.304900000000, 0.370900000000, 0.081000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.503700000000, 0.550200000000, 0.292000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.504000000000, 0.347200000000, 0.292200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.422380000000, 0.443500000000, 0.664020000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.307160000000, 0.472280000000, 0.810700000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.279110000000, 0.438730000000, 0.031200000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.122150000000, 0.437020000000, 0.026700000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.071280000000, 0.472780000000, 0.814490000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.186410000000, 0.441040000000, 0.671820000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.422650000000, 0.672500000000, 0.672820000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.307780000000, 0.630160000000, 0.814520000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.275540000000, 0.672790000000, 0.031090000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.120580000000, 0.673100000000, 0.029790000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.070440000000, 0.630370000000, 0.818000000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.187060000000, 0.673270000000, 0.680670000000), 2.05, 1.0, 0, 0, 0, false, false),
                 Atom(double3(0.372600000000, 0.446600000000, 0.755800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.308400000000, 0.441300000000, 0.921100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.200700000000, 0.440800000000, 0.028900000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.096900000000, 0.438900000000, 0.914400000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.114900000000, 0.445900000000, 0.723700000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.243500000000, 0.444700000000, 0.754000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.374200000000, 0.656100000000, 0.762800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.308500000000, 0.655200000000, 0.927200000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.198000000000, 0.655400000000, 0.028800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.091000000000, 0.661400000000, 0.922300000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.116900000000, 0.657800000000, 0.730600000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.244800000000, 0.659400000000, 0.757800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.304700000000, 0.551000000000, 0.813400000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.076800000000, 0.551900000000, 0.823100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.416100000000, 0.372400000000, 0.610400000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.408600000000, 0.501700000000, 0.586400000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.402000000000, 0.631400000000, 0.576100000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.188600000000, 0.370200000000, 0.616400000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.194000000000, 0.499300000000, 0.591800000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.195100000000, 0.629100000000, 0.581000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.996300000000, 0.449800000000, 0.792000000000), -1.025, 1.0, 0, 1, 0, false, false),
                 Atom(double3(0.996000000000, 0.652800000000, 0.792200000000), -1.025, 1.0, 0, 1, 0, false, false)},
                int3(1, 1, 1));

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {}, {1}, 5);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  atomData[0].position = double3(10.011, 4.097475, 0.0);

  RunningEnergy energy = Interactions::computeFrameworkMoleculeEnergy(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework, frameworkAtoms, atomData);

  EXPECT_NEAR(energy.frameworkMoleculeVDW * Units::EnergyToKelvin, -1784.82292180, 1e-6);
}
