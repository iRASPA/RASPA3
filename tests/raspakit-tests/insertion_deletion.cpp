#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <cstddef>
#include <ranges>
#include <span>
#include <tuple>
#include <vector>

import int3;
import double3;
import double3x3;
import factory;
import units;
import molecule;
import atom;
import pseudo_atom;
import vdwparameters;
import forcefield;
import framework;
import component;
import system;
import simulationbox;

TEST(insertion_deletion, methane_number_of_molecules_per_component)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Component c = TestFactories::makeMethane(forceField, 0);
  System system = System(0, forceField, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, 1.0, {}, {c}, {20}, 5);

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 20uz);
  EXPECT_EQ(system.numberOfPseudoAtoms[0][2], 20uz);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();

  for (size_t i = 0; i != atomPositions.size(); ++i)
  {
    EXPECT_NEAR(atomPositions[i].charge, 0.0, 1e-6);
    EXPECT_NEAR(atomPositions[i].scalingVDW, 1.0, 1e-6);
    EXPECT_NEAR(atomPositions[i].scalingCoulomb, 1.0, 1e-6);
    EXPECT_EQ(atomPositions[i].moleculeId, i);
    EXPECT_EQ(atomPositions[i].type, 2);
    EXPECT_EQ(atomPositions[i].componentId, 0);
    EXPECT_EQ(atomPositions[i].groupId, 0);
  }
}

TEST(insertion_deletion, CO2_number_of_molecules_per_component)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Component c = TestFactories::makeCO2(forceField, 0, true);
  System system = System(0, forceField, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, 1.0, std::nullopt, {c}, {3}, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();

  for (size_t i = 0; i != atomPositions.size(); ++i)
  {
    EXPECT_NEAR(atomPositions[i].scalingVDW, 1.0, 1e-6);
    EXPECT_NEAR(atomPositions[i].scalingCoulomb, 1.0, 1e-6);
    EXPECT_EQ(atomPositions[i].componentId, 0);
    EXPECT_EQ(atomPositions[i].groupId, 0);
  }

  EXPECT_EQ(atomPositions[0].moleculeId, 0);
  EXPECT_EQ(atomPositions[1].moleculeId, 0);
  EXPECT_EQ(atomPositions[2].moleculeId, 0);
  EXPECT_EQ(atomPositions[3].moleculeId, 1);
  EXPECT_EQ(atomPositions[4].moleculeId, 1);
  EXPECT_EQ(atomPositions[5].moleculeId, 1);
  EXPECT_EQ(atomPositions[6].moleculeId, 2);
  EXPECT_EQ(atomPositions[7].moleculeId, 2);
  EXPECT_EQ(atomPositions[8].moleculeId, 2);

  EXPECT_EQ(atomPositions[0].type, 4);
  EXPECT_EQ(atomPositions[1].type, 3);
  EXPECT_EQ(atomPositions[2].type, 4);
  EXPECT_EQ(atomPositions[3].type, 4);
  EXPECT_EQ(atomPositions[4].type, 3);
  EXPECT_EQ(atomPositions[5].type, 4);
  EXPECT_EQ(atomPositions[6].type, 4);
  EXPECT_EQ(atomPositions[7].type, 3);
  EXPECT_EQ(atomPositions[8].type, 4);
}

TEST(insertion_deletion, CO2_Methane_number_of_molecules_per_component)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Framework f = TestFactories::makeMFI_Si(forceField, int3(2, 2, 2));
  Component methane = TestFactories::makeMethane(forceField, 0);
  Component co2 = TestFactories::makeCO2(forceField, 1, true);
  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {methane, co2}, {5, 3}, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();

  // component 0: Methane
  EXPECT_EQ(atomPositions[0].componentId, 0);
  EXPECT_EQ(atomPositions[1].componentId, 0);
  EXPECT_EQ(atomPositions[2].componentId, 0);
  EXPECT_EQ(atomPositions[3].componentId, 0);
  EXPECT_EQ(atomPositions[4].componentId, 0);
  // component 1: CO2
  EXPECT_EQ(atomPositions[5].componentId, 1);
  EXPECT_EQ(atomPositions[6].componentId, 1);
  EXPECT_EQ(atomPositions[7].componentId, 1);
  EXPECT_EQ(atomPositions[8].componentId, 1);
  EXPECT_EQ(atomPositions[9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);
  EXPECT_EQ(atomPositions[13].componentId, 1);

  // component 0: Methane
  EXPECT_EQ(atomPositions[0].moleculeId, 0);
  EXPECT_EQ(atomPositions[1].moleculeId, 1);
  EXPECT_EQ(atomPositions[2].moleculeId, 2);
  EXPECT_EQ(atomPositions[3].moleculeId, 3);
  EXPECT_EQ(atomPositions[4].moleculeId, 4);
  // component 1: CO2
  EXPECT_EQ(atomPositions[5].moleculeId, 0);
  EXPECT_EQ(atomPositions[6].moleculeId, 0);
  EXPECT_EQ(atomPositions[7].moleculeId, 0);
  EXPECT_EQ(atomPositions[8].moleculeId, 1);
  EXPECT_EQ(atomPositions[9].moleculeId, 1);
  EXPECT_EQ(atomPositions[10].moleculeId, 1);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 2);
  EXPECT_EQ(atomPositions[13].moleculeId, 2);

  // component 0: Methane
  EXPECT_EQ(atomPositions[0].type, 2);
  EXPECT_EQ(atomPositions[1].type, 2);
  EXPECT_EQ(atomPositions[2].type, 2);
  EXPECT_EQ(atomPositions[3].type, 2);
  EXPECT_EQ(atomPositions[4].type, 2);
  // component 1: CO2
  EXPECT_EQ(atomPositions[5].type, 4);
  EXPECT_EQ(atomPositions[6].type, 3);
  EXPECT_EQ(atomPositions[7].type, 4);
  EXPECT_EQ(atomPositions[8].type, 4);
  EXPECT_EQ(atomPositions[9].type, 3);
  EXPECT_EQ(atomPositions[10].type, 4);
  EXPECT_EQ(atomPositions[11].type, 4);
  EXPECT_EQ(atomPositions[12].type, 3);
  EXPECT_EQ(atomPositions[13].type, 4);

  // component 0: Methane
  EXPECT_NEAR(atomPositions[0].charge, 0.0, 1e-6);
  EXPECT_NEAR(atomPositions[1].charge, 0.0, 1e-6);
  EXPECT_NEAR(atomPositions[2].charge, 0.0, 1e-6);
  EXPECT_NEAR(atomPositions[3].charge, 0.0, 1e-6);
  EXPECT_NEAR(atomPositions[4].charge, 0.0, 1e-6);
  // component 1: CO2
  EXPECT_NEAR(atomPositions[5].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomPositions[6].charge, 0.6512, 1e-6);
  EXPECT_NEAR(atomPositions[7].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomPositions[8].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomPositions[9].charge, 0.6512, 1e-6);
  EXPECT_NEAR(atomPositions[10].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomPositions[11].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomPositions[12].charge, 0.6512, 1e-6);
  EXPECT_NEAR(atomPositions[13].charge, -0.3256, 1e-6);
}

TEST(insertion_deletion, Dynamic_CO2_Methane_number_of_molecules_per_component)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Framework f = TestFactories::makeMFI_Si(forceField, int3(2, 2, 2));
  Component methane = TestFactories::makeMethane(forceField, 0);
  Component co2 = TestFactories::makeCO2(forceField, 1, true);
  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {methane, co2}, {5, 3}, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 3);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 3);

  // insert new CO2
  system.insertMolecule(
      1, Molecule(),
      {Atom(double3(0.0, 0.0, 1.149), -0.3256, 1.0, 0, 4, 1, 0), Atom(double3(0.0, 0.0, 0.0), 0.6512, 1.0, 0, 3, 1, 0),
       Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)});

  // refresh span
  atomPositions = system.spanOfMoleculeAtoms();

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 0);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 4);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 4);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 0);

  EXPECT_EQ(atomPositions[0].moleculeId, 0);
  EXPECT_EQ(atomPositions[1].moleculeId, 1);
  EXPECT_EQ(atomPositions[2].moleculeId, 2);
  EXPECT_EQ(atomPositions[3].moleculeId, 3);
  EXPECT_EQ(atomPositions[4].moleculeId, 4);
  EXPECT_EQ(atomPositions[5].moleculeId, 0);
  EXPECT_EQ(atomPositions[6].moleculeId, 0);
  EXPECT_EQ(atomPositions[7].moleculeId, 0);
  EXPECT_EQ(atomPositions[8].moleculeId, 1);
  EXPECT_EQ(atomPositions[9].moleculeId, 1);
  EXPECT_EQ(atomPositions[10].moleculeId, 1);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 2);
  EXPECT_EQ(atomPositions[13].moleculeId, 2);
  EXPECT_EQ(atomPositions[14].moleculeId, 3);
  EXPECT_EQ(atomPositions[15].moleculeId, 3);
  EXPECT_EQ(atomPositions[16].moleculeId, 3);

  EXPECT_EQ(atomPositions[0].componentId, 0);
  EXPECT_EQ(atomPositions[1].componentId, 0);
  EXPECT_EQ(atomPositions[2].componentId, 0);
  EXPECT_EQ(atomPositions[3].componentId, 0);
  EXPECT_EQ(atomPositions[4].componentId, 0);
  EXPECT_EQ(atomPositions[5].componentId, 1);
  EXPECT_EQ(atomPositions[6].componentId, 1);
  EXPECT_EQ(atomPositions[7].componentId, 1);
  EXPECT_EQ(atomPositions[8].componentId, 1);
  EXPECT_EQ(atomPositions[9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);
  EXPECT_EQ(atomPositions[13].componentId, 1);
  EXPECT_EQ(atomPositions[14].componentId, 1);
  EXPECT_EQ(atomPositions[15].componentId, 1);
  EXPECT_EQ(atomPositions[16].componentId, 1);

  // delete Methane 2
  std::vector<Atom>::iterator iterator_methane1 = system.iteratorForMolecule(0, 2);
  system.deleteMolecule(0, 2, {iterator_methane1, 1});

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 4);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 4);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 0);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 4);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 4);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 0);

  // refresh span
  atomPositions = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomPositions[0].moleculeId, 0);
  EXPECT_EQ(atomPositions[1].moleculeId, 1);
  EXPECT_EQ(atomPositions[2].moleculeId, 2);
  EXPECT_EQ(atomPositions[3].moleculeId, 3);
  EXPECT_EQ(atomPositions[4].moleculeId, 0);
  EXPECT_EQ(atomPositions[5].moleculeId, 0);
  EXPECT_EQ(atomPositions[6].moleculeId, 0);
  EXPECT_EQ(atomPositions[7].moleculeId, 1);
  EXPECT_EQ(atomPositions[8].moleculeId, 1);
  EXPECT_EQ(atomPositions[9].moleculeId, 1);
  EXPECT_EQ(atomPositions[10].moleculeId, 2);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 2);
  EXPECT_EQ(atomPositions[13].moleculeId, 3);
  EXPECT_EQ(atomPositions[14].moleculeId, 3);
  EXPECT_EQ(atomPositions[15].moleculeId, 3);

  EXPECT_EQ(atomPositions[0].componentId, 0);
  EXPECT_EQ(atomPositions[1].componentId, 0);
  EXPECT_EQ(atomPositions[2].componentId, 0);
  EXPECT_EQ(atomPositions[3].componentId, 0);
  EXPECT_EQ(atomPositions[4].componentId, 1);
  EXPECT_EQ(atomPositions[5].componentId, 1);
  EXPECT_EQ(atomPositions[6].componentId, 1);
  EXPECT_EQ(atomPositions[7].componentId, 1);
  EXPECT_EQ(atomPositions[8].componentId, 1);
  EXPECT_EQ(atomPositions[9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);
  EXPECT_EQ(atomPositions[13].componentId, 1);
  EXPECT_EQ(atomPositions[14].componentId, 1);
  EXPECT_EQ(atomPositions[15].componentId, 1);

  // delete Methane 2
  std::vector<Atom>::iterator iterator_methane2 = system.iteratorForMolecule(0, 2);
  system.deleteMolecule(0, 2, {iterator_methane2, 1});

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 3);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 3);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 0);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 4);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 4);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 0);

  // refresh span
  atomPositions = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomPositions[0].moleculeId, 0);
  EXPECT_EQ(atomPositions[1].moleculeId, 1);
  EXPECT_EQ(atomPositions[2].moleculeId, 2);
  EXPECT_EQ(atomPositions[3].moleculeId, 0);
  EXPECT_EQ(atomPositions[4].moleculeId, 0);
  EXPECT_EQ(atomPositions[5].moleculeId, 0);
  EXPECT_EQ(atomPositions[6].moleculeId, 1);
  EXPECT_EQ(atomPositions[7].moleculeId, 1);
  EXPECT_EQ(atomPositions[8].moleculeId, 1);
  EXPECT_EQ(atomPositions[9].moleculeId, 2);
  EXPECT_EQ(atomPositions[10].moleculeId, 2);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 3);
  EXPECT_EQ(atomPositions[13].moleculeId, 3);
  EXPECT_EQ(atomPositions[14].moleculeId, 3);

  EXPECT_EQ(atomPositions[0].componentId, 0);
  EXPECT_EQ(atomPositions[1].componentId, 0);
  EXPECT_EQ(atomPositions[2].componentId, 0);
  EXPECT_EQ(atomPositions[3].componentId, 1);
  EXPECT_EQ(atomPositions[4].componentId, 1);
  EXPECT_EQ(atomPositions[5].componentId, 1);
  EXPECT_EQ(atomPositions[6].componentId, 1);
  EXPECT_EQ(atomPositions[7].componentId, 1);
  EXPECT_EQ(atomPositions[8].componentId, 1);
  EXPECT_EQ(atomPositions[9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);
  EXPECT_EQ(atomPositions[13].componentId, 1);
  EXPECT_EQ(atomPositions[14].componentId, 1);

  // insert new CO2
  system.insertMolecule(
      1, Molecule(),
      {Atom(double3(0.0, 0.0, 1.149), -0.3256, 1.0, 0, 4, 1, 0), Atom(double3(0.0, 0.0, 0.0), 0.6512, 1.0, 0, 3, 1, 0),
       Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)});

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 3);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 3);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 0);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 5);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 5);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 0);

  // refresh span
  atomPositions = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomPositions[0].moleculeId, 0);
  EXPECT_EQ(atomPositions[1].moleculeId, 1);
  EXPECT_EQ(atomPositions[2].moleculeId, 2);
  EXPECT_EQ(atomPositions[3].moleculeId, 0);
  EXPECT_EQ(atomPositions[4].moleculeId, 0);
  EXPECT_EQ(atomPositions[5].moleculeId, 0);
  EXPECT_EQ(atomPositions[6].moleculeId, 1);
  EXPECT_EQ(atomPositions[7].moleculeId, 1);
  EXPECT_EQ(atomPositions[8].moleculeId, 1);
  EXPECT_EQ(atomPositions[9].moleculeId, 2);
  EXPECT_EQ(atomPositions[10].moleculeId, 2);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 3);
  EXPECT_EQ(atomPositions[13].moleculeId, 3);
  EXPECT_EQ(atomPositions[14].moleculeId, 3);
  EXPECT_EQ(atomPositions[15].moleculeId, 4);
  EXPECT_EQ(atomPositions[16].moleculeId, 4);
  EXPECT_EQ(atomPositions[17].moleculeId, 4);

  EXPECT_EQ(atomPositions[0].componentId, 0);
  EXPECT_EQ(atomPositions[1].componentId, 0);
  EXPECT_EQ(atomPositions[2].componentId, 0);
  EXPECT_EQ(atomPositions[3].componentId, 1);
  EXPECT_EQ(atomPositions[4].componentId, 1);
  EXPECT_EQ(atomPositions[5].componentId, 1);
  EXPECT_EQ(atomPositions[6].componentId, 1);
  EXPECT_EQ(atomPositions[7].componentId, 1);
  EXPECT_EQ(atomPositions[8].componentId, 1);
  EXPECT_EQ(atomPositions[9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);
  EXPECT_EQ(atomPositions[13].componentId, 1);
  EXPECT_EQ(atomPositions[14].componentId, 1);
  EXPECT_EQ(atomPositions[15].componentId, 1);
  EXPECT_EQ(atomPositions[16].componentId, 1);
  EXPECT_EQ(atomPositions[17].componentId, 1);

  // insert new Methane
  system.insertMolecule(0, Molecule(), {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 2, 0, 0)});

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 4);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 4);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 0);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 5);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 5);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 0);

  // refresh span
  atomPositions = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomPositions[0].moleculeId, 0);
  EXPECT_EQ(atomPositions[1].moleculeId, 1);
  EXPECT_EQ(atomPositions[2].moleculeId, 2);
  EXPECT_EQ(atomPositions[3].moleculeId, 3);
  EXPECT_EQ(atomPositions[4].moleculeId, 0);
  EXPECT_EQ(atomPositions[5].moleculeId, 0);
  EXPECT_EQ(atomPositions[6].moleculeId, 0);
  EXPECT_EQ(atomPositions[7].moleculeId, 1);
  EXPECT_EQ(atomPositions[8].moleculeId, 1);
  EXPECT_EQ(atomPositions[9].moleculeId, 1);
  EXPECT_EQ(atomPositions[10].moleculeId, 2);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 2);
  EXPECT_EQ(atomPositions[13].moleculeId, 3);
  EXPECT_EQ(atomPositions[14].moleculeId, 3);
  EXPECT_EQ(atomPositions[15].moleculeId, 3);
  EXPECT_EQ(atomPositions[16].moleculeId, 4);
  EXPECT_EQ(atomPositions[17].moleculeId, 4);
  EXPECT_EQ(atomPositions[18].moleculeId, 4);

  EXPECT_EQ(atomPositions[0].componentId, 0);
  EXPECT_EQ(atomPositions[1].componentId, 0);
  EXPECT_EQ(atomPositions[2].componentId, 0);
  EXPECT_EQ(atomPositions[3].componentId, 0);
  EXPECT_EQ(atomPositions[4].componentId, 1);
  EXPECT_EQ(atomPositions[5].componentId, 1);
  EXPECT_EQ(atomPositions[6].componentId, 1);
  EXPECT_EQ(atomPositions[7].componentId, 1);
  EXPECT_EQ(atomPositions[8].componentId, 1);
  EXPECT_EQ(atomPositions[9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);
  EXPECT_EQ(atomPositions[13].componentId, 1);
  EXPECT_EQ(atomPositions[14].componentId, 1);
  EXPECT_EQ(atomPositions[15].componentId, 1);
  EXPECT_EQ(atomPositions[16].componentId, 1);
  EXPECT_EQ(atomPositions[17].componentId, 1);
  EXPECT_EQ(atomPositions[18].componentId, 1);

  // delete CO2 1
  std::vector<Atom>::iterator iterator_CO2 = system.iteratorForMolecule(1, 1);
  system.deleteMolecule(1, 1, {iterator_CO2, 3});

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 4);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 4);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 0);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 4);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 4);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 0);

  // refresh span
  atomPositions = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomPositions[0].moleculeId, 0);
  EXPECT_EQ(atomPositions[1].moleculeId, 1);
  EXPECT_EQ(atomPositions[2].moleculeId, 2);
  EXPECT_EQ(atomPositions[3].moleculeId, 3);
  EXPECT_EQ(atomPositions[4].moleculeId, 0);
  EXPECT_EQ(atomPositions[5].moleculeId, 0);
  EXPECT_EQ(atomPositions[6].moleculeId, 0);
  EXPECT_EQ(atomPositions[7].moleculeId, 1);
  EXPECT_EQ(atomPositions[8].moleculeId, 1);
  EXPECT_EQ(atomPositions[9].moleculeId, 1);
  EXPECT_EQ(atomPositions[10].moleculeId, 2);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 2);
  EXPECT_EQ(atomPositions[13].moleculeId, 3);
  EXPECT_EQ(atomPositions[14].moleculeId, 3);
  EXPECT_EQ(atomPositions[15].moleculeId, 3);

  EXPECT_EQ(atomPositions[0].componentId, 0);
  EXPECT_EQ(atomPositions[1].componentId, 0);
  EXPECT_EQ(atomPositions[2].componentId, 0);
  EXPECT_EQ(atomPositions[3].componentId, 0);
  EXPECT_EQ(atomPositions[4].componentId, 1);
  EXPECT_EQ(atomPositions[5].componentId, 1);
  EXPECT_EQ(atomPositions[6].componentId, 1);
  EXPECT_EQ(atomPositions[7].componentId, 1);
  EXPECT_EQ(atomPositions[8].componentId, 1);
  EXPECT_EQ(atomPositions[9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);
  EXPECT_EQ(atomPositions[13].componentId, 1);
  EXPECT_EQ(atomPositions[14].componentId, 1);
  EXPECT_EQ(atomPositions[15].componentId, 1);

  // delete CO2 1
  std::vector<Atom>::iterator iterator_CO2_2 = system.iteratorForMolecule(1, 1);
  system.deleteMolecule(1, 1, {iterator_CO2_2, 3});

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 4);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 4);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 0);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 3);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 3);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 0);

  // refresh span
  atomPositions = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomPositions[0].moleculeId, 0);
  EXPECT_EQ(atomPositions[1].moleculeId, 1);
  EXPECT_EQ(atomPositions[2].moleculeId, 2);
  EXPECT_EQ(atomPositions[3].moleculeId, 3);
  EXPECT_EQ(atomPositions[4].moleculeId, 0);
  EXPECT_EQ(atomPositions[5].moleculeId, 0);
  EXPECT_EQ(atomPositions[6].moleculeId, 0);
  EXPECT_EQ(atomPositions[7].moleculeId, 1);
  EXPECT_EQ(atomPositions[8].moleculeId, 1);
  EXPECT_EQ(atomPositions[9].moleculeId, 1);
  EXPECT_EQ(atomPositions[10].moleculeId, 2);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 2);

  EXPECT_EQ(atomPositions[0].componentId, 0);
  EXPECT_EQ(atomPositions[1].componentId, 0);
  EXPECT_EQ(atomPositions[2].componentId, 0);
  EXPECT_EQ(atomPositions[3].componentId, 0);
  EXPECT_EQ(atomPositions[4].componentId, 1);
  EXPECT_EQ(atomPositions[5].componentId, 1);
  EXPECT_EQ(atomPositions[6].componentId, 1);
  EXPECT_EQ(atomPositions[7].componentId, 1);
  EXPECT_EQ(atomPositions[8].componentId, 1);
  EXPECT_EQ(atomPositions[9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);

  // insert new Methane
  system.insertMolecule(0, Molecule(), {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 2, 0, 0)});

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 0);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 3);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 3);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 0);

  // refresh span
  atomPositions = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomPositions[0].moleculeId, 0);
  EXPECT_EQ(atomPositions[1].moleculeId, 1);
  EXPECT_EQ(atomPositions[2].moleculeId, 2);
  EXPECT_EQ(atomPositions[3].moleculeId, 3);
  EXPECT_EQ(atomPositions[4].moleculeId, 4);
  EXPECT_EQ(atomPositions[5].moleculeId, 0);
  EXPECT_EQ(atomPositions[6].moleculeId, 0);
  EXPECT_EQ(atomPositions[7].moleculeId, 0);
  EXPECT_EQ(atomPositions[8].moleculeId, 1);
  EXPECT_EQ(atomPositions[9].moleculeId, 1);
  EXPECT_EQ(atomPositions[10].moleculeId, 1);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 2);
  EXPECT_EQ(atomPositions[13].moleculeId, 2);

  EXPECT_EQ(atomPositions[0].componentId, 0);
  EXPECT_EQ(atomPositions[1].componentId, 0);
  EXPECT_EQ(atomPositions[2].componentId, 0);
  EXPECT_EQ(atomPositions[3].componentId, 0);
  EXPECT_EQ(atomPositions[4].componentId, 0);
  EXPECT_EQ(atomPositions[5].componentId, 1);
  EXPECT_EQ(atomPositions[6].componentId, 1);
  EXPECT_EQ(atomPositions[7].componentId, 1);
  EXPECT_EQ(atomPositions[8].componentId, 1);
  EXPECT_EQ(atomPositions[9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);
  EXPECT_EQ(atomPositions[13].componentId, 1);

  // delete Methane 2
  std::vector<Atom>::iterator iterator_methane3 = system.iteratorForMolecule(0, 2);
  system.deleteMolecule(0, 2, {iterator_methane3, 1});

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 4);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 4);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 0);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 3);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 3);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 0);

  // refresh span
  atomPositions = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomPositions[0].moleculeId, 0);
  EXPECT_EQ(atomPositions[1].moleculeId, 1);
  EXPECT_EQ(atomPositions[2].moleculeId, 2);
  EXPECT_EQ(atomPositions[3].moleculeId, 3);
  EXPECT_EQ(atomPositions[4].moleculeId, 0);
  EXPECT_EQ(atomPositions[5].moleculeId, 0);
  EXPECT_EQ(atomPositions[6].moleculeId, 0);
  EXPECT_EQ(atomPositions[7].moleculeId, 1);
  EXPECT_EQ(atomPositions[8].moleculeId, 1);
  EXPECT_EQ(atomPositions[9].moleculeId, 1);
  EXPECT_EQ(atomPositions[10].moleculeId, 2);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 2);

  EXPECT_EQ(atomPositions[0].componentId, 0);
  EXPECT_EQ(atomPositions[1].componentId, 0);
  EXPECT_EQ(atomPositions[2].componentId, 0);
  EXPECT_EQ(atomPositions[3].componentId, 0);
  EXPECT_EQ(atomPositions[4].componentId, 1);
  EXPECT_EQ(atomPositions[5].componentId, 1);
  EXPECT_EQ(atomPositions[6].componentId, 1);
  EXPECT_EQ(atomPositions[7].componentId, 1);
  EXPECT_EQ(atomPositions[8].componentId, 1);
  EXPECT_EQ(atomPositions[9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);

  // insert new CO2
  system.insertMolecule(
      1, Molecule(),
      {Atom(double3(0.0, 0.0, 1.149), -0.3256, 1.0, 0, 4, 1, 0), Atom(double3(0.0, 0.0, 0.0), 0.6512, 1.0, 0, 3, 1, 0),
       Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)});

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 4);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 4);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 0);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 4);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 4);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 0);

  // refresh span
  atomPositions = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomPositions[0].moleculeId, 0);
  EXPECT_EQ(atomPositions[1].moleculeId, 1);
  EXPECT_EQ(atomPositions[2].moleculeId, 2);
  EXPECT_EQ(atomPositions[3].moleculeId, 3);
  EXPECT_EQ(atomPositions[4].moleculeId, 0);
  EXPECT_EQ(atomPositions[5].moleculeId, 0);
  EXPECT_EQ(atomPositions[6].moleculeId, 0);
  EXPECT_EQ(atomPositions[7].moleculeId, 1);
  EXPECT_EQ(atomPositions[8].moleculeId, 1);
  EXPECT_EQ(atomPositions[9].moleculeId, 1);
  EXPECT_EQ(atomPositions[10].moleculeId, 2);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 2);
  EXPECT_EQ(atomPositions[13].moleculeId, 3);
  EXPECT_EQ(atomPositions[14].moleculeId, 3);
  EXPECT_EQ(atomPositions[15].moleculeId, 3);

  EXPECT_EQ(atomPositions[0].componentId, 0);
  EXPECT_EQ(atomPositions[1].componentId, 0);
  EXPECT_EQ(atomPositions[2].componentId, 0);
  EXPECT_EQ(atomPositions[3].componentId, 0);
  EXPECT_EQ(atomPositions[4].componentId, 1);
  EXPECT_EQ(atomPositions[5].componentId, 1);
  EXPECT_EQ(atomPositions[6].componentId, 1);
  EXPECT_EQ(atomPositions[7].componentId, 1);
  EXPECT_EQ(atomPositions[8].componentId, 1);
  EXPECT_EQ(atomPositions[9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);
  EXPECT_EQ(atomPositions[13].componentId, 1);
  EXPECT_EQ(atomPositions[14].componentId, 1);
  EXPECT_EQ(atomPositions[15].componentId, 1);

  // insert new Methane
  system.insertMolecule(0, Molecule(), {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 2, 0, 0)});

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 0);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 4);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 4);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 0);

  // refresh span
  atomPositions = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomPositions[0].moleculeId, 0);
  EXPECT_EQ(atomPositions[1].moleculeId, 1);
  EXPECT_EQ(atomPositions[2].moleculeId, 2);
  EXPECT_EQ(atomPositions[3].moleculeId, 3);
  EXPECT_EQ(atomPositions[4].moleculeId, 4);
  EXPECT_EQ(atomPositions[5].moleculeId, 0);
  EXPECT_EQ(atomPositions[6].moleculeId, 0);
  EXPECT_EQ(atomPositions[7].moleculeId, 0);
  EXPECT_EQ(atomPositions[8].moleculeId, 1);
  EXPECT_EQ(atomPositions[9].moleculeId, 1);
  EXPECT_EQ(atomPositions[10].moleculeId, 1);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 2);
  EXPECT_EQ(atomPositions[13].moleculeId, 2);
  EXPECT_EQ(atomPositions[14].moleculeId, 3);
  EXPECT_EQ(atomPositions[15].moleculeId, 3);
  EXPECT_EQ(atomPositions[16].moleculeId, 3);

  EXPECT_EQ(atomPositions[0].componentId, 0);
  EXPECT_EQ(atomPositions[1].componentId, 0);
  EXPECT_EQ(atomPositions[2].componentId, 0);
  EXPECT_EQ(atomPositions[3].componentId, 0);
  EXPECT_EQ(atomPositions[4].componentId, 0);
  EXPECT_EQ(atomPositions[5].componentId, 1);
  EXPECT_EQ(atomPositions[6].componentId, 1);
  EXPECT_EQ(atomPositions[7].componentId, 1);
  EXPECT_EQ(atomPositions[8].componentId, 1);
  EXPECT_EQ(atomPositions[9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);
  EXPECT_EQ(atomPositions[13].componentId, 1);
  EXPECT_EQ(atomPositions[14].componentId, 1);
  EXPECT_EQ(atomPositions[15].componentId, 1);
  EXPECT_EQ(atomPositions[16].componentId, 1);

  // delete CO2 0
  std::vector<Atom>::iterator iterator_CO2_3 = system.iteratorForMolecule(1, 1);
  system.deleteMolecule(1, 0, {iterator_CO2_3, 3});

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 0);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 3);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 3);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 0);

  // refresh span
  atomPositions = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomPositions[0].moleculeId, 0);
  EXPECT_EQ(atomPositions[1].moleculeId, 1);
  EXPECT_EQ(atomPositions[2].moleculeId, 2);
  EXPECT_EQ(atomPositions[3].moleculeId, 3);
  EXPECT_EQ(atomPositions[4].moleculeId, 4);
  EXPECT_EQ(atomPositions[5].moleculeId, 0);
  EXPECT_EQ(atomPositions[6].moleculeId, 0);
  EXPECT_EQ(atomPositions[7].moleculeId, 0);
  EXPECT_EQ(atomPositions[8].moleculeId, 1);
  EXPECT_EQ(atomPositions[9].moleculeId, 1);
  EXPECT_EQ(atomPositions[10].moleculeId, 1);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 2);
  EXPECT_EQ(atomPositions[13].moleculeId, 2);

  EXPECT_EQ(atomPositions[0].componentId, 0);
  EXPECT_EQ(atomPositions[1].componentId, 0);
  EXPECT_EQ(atomPositions[2].componentId, 0);
  EXPECT_EQ(atomPositions[3].componentId, 0);
  EXPECT_EQ(atomPositions[4].componentId, 0);
  EXPECT_EQ(atomPositions[5].componentId, 1);
  EXPECT_EQ(atomPositions[6].componentId, 1);
  EXPECT_EQ(atomPositions[7].componentId, 1);
  EXPECT_EQ(atomPositions[8].componentId, 1);
  EXPECT_EQ(atomPositions[9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);
  EXPECT_EQ(atomPositions[13].componentId, 1);

  // component 0: Methane
  EXPECT_EQ(atomPositions[0].componentId, 0);
  EXPECT_EQ(atomPositions[1].componentId, 0);
  EXPECT_EQ(atomPositions[2].componentId, 0);
  EXPECT_EQ(atomPositions[3].componentId, 0);
  EXPECT_EQ(atomPositions[4].componentId, 0);
  // component 1: CO2
  EXPECT_EQ(atomPositions[5].componentId, 1);
  EXPECT_EQ(atomPositions[6].componentId, 1);
  EXPECT_EQ(atomPositions[7].componentId, 1);
  EXPECT_EQ(atomPositions[8].componentId, 1);
  EXPECT_EQ(atomPositions[9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);
  EXPECT_EQ(atomPositions[13].componentId, 1);

  // component 0: Methane
  EXPECT_EQ(atomPositions[0].moleculeId, 0);
  EXPECT_EQ(atomPositions[1].moleculeId, 1);
  EXPECT_EQ(atomPositions[2].moleculeId, 2);
  EXPECT_EQ(atomPositions[3].moleculeId, 3);
  EXPECT_EQ(atomPositions[4].moleculeId, 4);
  // component 1: CO2
  EXPECT_EQ(atomPositions[5].moleculeId, 0);
  EXPECT_EQ(atomPositions[6].moleculeId, 0);
  EXPECT_EQ(atomPositions[7].moleculeId, 0);
  EXPECT_EQ(atomPositions[8].moleculeId, 1);
  EXPECT_EQ(atomPositions[9].moleculeId, 1);
  EXPECT_EQ(atomPositions[10].moleculeId, 1);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 2);
  EXPECT_EQ(atomPositions[13].moleculeId, 2);

  // component 0: Methane
  EXPECT_EQ(atomPositions[0].type, 2);
  EXPECT_EQ(atomPositions[1].type, 2);
  EXPECT_EQ(atomPositions[2].type, 2);
  EXPECT_EQ(atomPositions[3].type, 2);
  EXPECT_EQ(atomPositions[4].type, 2);
  // component 1: CO2
  EXPECT_EQ(atomPositions[5].type, 4);
  EXPECT_EQ(atomPositions[6].type, 3);
  EXPECT_EQ(atomPositions[7].type, 4);
  EXPECT_EQ(atomPositions[8].type, 4);
  EXPECT_EQ(atomPositions[9].type, 3);
  EXPECT_EQ(atomPositions[10].type, 4);
  EXPECT_EQ(atomPositions[11].type, 4);
  EXPECT_EQ(atomPositions[12].type, 3);
  EXPECT_EQ(atomPositions[13].type, 4);

  // component 0: Methane
  EXPECT_NEAR(atomPositions[0].charge, 0.0, 1e-6);
  EXPECT_NEAR(atomPositions[1].charge, 0.0, 1e-6);
  EXPECT_NEAR(atomPositions[2].charge, 0.0, 1e-6);
  EXPECT_NEAR(atomPositions[3].charge, 0.0, 1e-6);
  EXPECT_NEAR(atomPositions[4].charge, 0.0, 1e-6);
  // component 1: CO2
  EXPECT_NEAR(atomPositions[5].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomPositions[6].charge, 0.6512, 1e-6);
  EXPECT_NEAR(atomPositions[7].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomPositions[8].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomPositions[9].charge, 0.6512, 1e-6);
  EXPECT_NEAR(atomPositions[10].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomPositions[11].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomPositions[12].charge, 0.6512, 1e-6);
  EXPECT_NEAR(atomPositions[13].charge, -0.3256, 1e-6);
}
