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
  System system = System(0, forceField, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, 1.0, {}, {c}, {}, {20}, 5);

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 20uz);
  EXPECT_EQ(system.numberOfPseudoAtoms[0][2], 20uz);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();

  for (size_t i = 0; i != atomData.size(); ++i)
  {
    EXPECT_NEAR(atomData[i].charge, 0.0, 1e-6);
    EXPECT_NEAR(atomData[i].scalingVDW, 1.0, 1e-6);
    EXPECT_NEAR(atomData[i].scalingCoulomb, 1.0, 1e-6);
    EXPECT_EQ(atomData[i].moleculeId, i);
    EXPECT_EQ(atomData[i].type, 2);
    EXPECT_EQ(atomData[i].componentId, 0);
    EXPECT_EQ(atomData[i].groupId, 0);
  }
}

TEST(insertion_deletion, CO2_number_of_molecules_per_component)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Component c = TestFactories::makeCO2(forceField, 0, true);
  System system =
      System(0, forceField, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, 1.0, std::nullopt, {c}, {}, {3}, 5);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();

  for (size_t i = 0; i != atomData.size(); ++i)
  {
    EXPECT_NEAR(atomData[i].scalingVDW, 1.0, 1e-6);
    EXPECT_NEAR(atomData[i].scalingCoulomb, 1.0, 1e-6);
    EXPECT_EQ(atomData[i].componentId, 0);
    EXPECT_EQ(atomData[i].groupId, 0);
  }

  EXPECT_EQ(atomData[0].moleculeId, 0);
  EXPECT_EQ(atomData[1].moleculeId, 0);
  EXPECT_EQ(atomData[2].moleculeId, 0);
  EXPECT_EQ(atomData[3].moleculeId, 1);
  EXPECT_EQ(atomData[4].moleculeId, 1);
  EXPECT_EQ(atomData[5].moleculeId, 1);
  EXPECT_EQ(atomData[6].moleculeId, 2);
  EXPECT_EQ(atomData[7].moleculeId, 2);
  EXPECT_EQ(atomData[8].moleculeId, 2);

  EXPECT_EQ(atomData[0].type, 4);
  EXPECT_EQ(atomData[1].type, 3);
  EXPECT_EQ(atomData[2].type, 4);
  EXPECT_EQ(atomData[3].type, 4);
  EXPECT_EQ(atomData[4].type, 3);
  EXPECT_EQ(atomData[5].type, 4);
  EXPECT_EQ(atomData[6].type, 4);
  EXPECT_EQ(atomData[7].type, 3);
  EXPECT_EQ(atomData[8].type, 4);
}

TEST(insertion_deletion, CO2_Methane_number_of_molecules_per_component)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Framework f = TestFactories::makeMFI_Si(forceField, int3(2, 2, 2));
  Component methane = TestFactories::makeMethane(forceField, 0);
  Component co2 = TestFactories::makeCO2(forceField, 1, true);
  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {methane, co2}, {}, {5, 3}, 5);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();

  // component 0: Methane
  EXPECT_EQ(atomData[0].componentId, 0);
  EXPECT_EQ(atomData[1].componentId, 0);
  EXPECT_EQ(atomData[2].componentId, 0);
  EXPECT_EQ(atomData[3].componentId, 0);
  EXPECT_EQ(atomData[4].componentId, 0);
  // component 1: CO2
  EXPECT_EQ(atomData[5].componentId, 1);
  EXPECT_EQ(atomData[6].componentId, 1);
  EXPECT_EQ(atomData[7].componentId, 1);
  EXPECT_EQ(atomData[8].componentId, 1);
  EXPECT_EQ(atomData[9].componentId, 1);
  EXPECT_EQ(atomData[10].componentId, 1);
  EXPECT_EQ(atomData[11].componentId, 1);
  EXPECT_EQ(atomData[12].componentId, 1);
  EXPECT_EQ(atomData[13].componentId, 1);

  // component 0: Methane
  EXPECT_EQ(atomData[0].moleculeId, 0);
  EXPECT_EQ(atomData[1].moleculeId, 1);
  EXPECT_EQ(atomData[2].moleculeId, 2);
  EXPECT_EQ(atomData[3].moleculeId, 3);
  EXPECT_EQ(atomData[4].moleculeId, 4);
  // component 1: CO2
  EXPECT_EQ(atomData[5].moleculeId, 0);
  EXPECT_EQ(atomData[6].moleculeId, 0);
  EXPECT_EQ(atomData[7].moleculeId, 0);
  EXPECT_EQ(atomData[8].moleculeId, 1);
  EXPECT_EQ(atomData[9].moleculeId, 1);
  EXPECT_EQ(atomData[10].moleculeId, 1);
  EXPECT_EQ(atomData[11].moleculeId, 2);
  EXPECT_EQ(atomData[12].moleculeId, 2);
  EXPECT_EQ(atomData[13].moleculeId, 2);

  // component 0: Methane
  EXPECT_EQ(atomData[0].type, 2);
  EXPECT_EQ(atomData[1].type, 2);
  EXPECT_EQ(atomData[2].type, 2);
  EXPECT_EQ(atomData[3].type, 2);
  EXPECT_EQ(atomData[4].type, 2);
  // component 1: CO2
  EXPECT_EQ(atomData[5].type, 4);
  EXPECT_EQ(atomData[6].type, 3);
  EXPECT_EQ(atomData[7].type, 4);
  EXPECT_EQ(atomData[8].type, 4);
  EXPECT_EQ(atomData[9].type, 3);
  EXPECT_EQ(atomData[10].type, 4);
  EXPECT_EQ(atomData[11].type, 4);
  EXPECT_EQ(atomData[12].type, 3);
  EXPECT_EQ(atomData[13].type, 4);

  // component 0: Methane
  EXPECT_NEAR(atomData[0].charge, 0.0, 1e-6);
  EXPECT_NEAR(atomData[1].charge, 0.0, 1e-6);
  EXPECT_NEAR(atomData[2].charge, 0.0, 1e-6);
  EXPECT_NEAR(atomData[3].charge, 0.0, 1e-6);
  EXPECT_NEAR(atomData[4].charge, 0.0, 1e-6);
  // component 1: CO2
  EXPECT_NEAR(atomData[5].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomData[6].charge, 0.6512, 1e-6);
  EXPECT_NEAR(atomData[7].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomData[8].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomData[9].charge, 0.6512, 1e-6);
  EXPECT_NEAR(atomData[10].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomData[11].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomData[12].charge, 0.6512, 1e-6);
  EXPECT_NEAR(atomData[13].charge, -0.3256, 1e-6);
}

TEST(insertion_deletion, Dynamic_CO2_Methane_number_of_molecules_per_component)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Framework f = TestFactories::makeMFI_Si(forceField, int3(2, 2, 2));
  Component methane = TestFactories::makeMethane(forceField, 0);
  Component co2 = TestFactories::makeCO2(forceField, 1, true);
  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {methane, co2}, {}, {5, 3}, 5);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 3);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 3);

  // insert new CO2
  system.insertMolecule(1, Molecule(),
                        {Atom(double3(0.0, 0.0, 1.149), -0.3256, 1.0, 0, 4, 1, false, false),
                         Atom(double3(0.0, 0.0, 0.0), 0.6512, 1.0, 0, 3, 1, false, false),
                         Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, false, false)});

  // refresh span
  atomData = system.spanOfMoleculeAtoms();

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 0);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 4);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 4);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 0);

  EXPECT_EQ(atomData[0].moleculeId, 0);
  EXPECT_EQ(atomData[1].moleculeId, 1);
  EXPECT_EQ(atomData[2].moleculeId, 2);
  EXPECT_EQ(atomData[3].moleculeId, 3);
  EXPECT_EQ(atomData[4].moleculeId, 4);
  EXPECT_EQ(atomData[5].moleculeId, 0);
  EXPECT_EQ(atomData[6].moleculeId, 0);
  EXPECT_EQ(atomData[7].moleculeId, 0);
  EXPECT_EQ(atomData[8].moleculeId, 1);
  EXPECT_EQ(atomData[9].moleculeId, 1);
  EXPECT_EQ(atomData[10].moleculeId, 1);
  EXPECT_EQ(atomData[11].moleculeId, 2);
  EXPECT_EQ(atomData[12].moleculeId, 2);
  EXPECT_EQ(atomData[13].moleculeId, 2);
  EXPECT_EQ(atomData[14].moleculeId, 3);
  EXPECT_EQ(atomData[15].moleculeId, 3);
  EXPECT_EQ(atomData[16].moleculeId, 3);

  EXPECT_EQ(atomData[0].componentId, 0);
  EXPECT_EQ(atomData[1].componentId, 0);
  EXPECT_EQ(atomData[2].componentId, 0);
  EXPECT_EQ(atomData[3].componentId, 0);
  EXPECT_EQ(atomData[4].componentId, 0);
  EXPECT_EQ(atomData[5].componentId, 1);
  EXPECT_EQ(atomData[6].componentId, 1);
  EXPECT_EQ(atomData[7].componentId, 1);
  EXPECT_EQ(atomData[8].componentId, 1);
  EXPECT_EQ(atomData[9].componentId, 1);
  EXPECT_EQ(atomData[10].componentId, 1);
  EXPECT_EQ(atomData[11].componentId, 1);
  EXPECT_EQ(atomData[12].componentId, 1);
  EXPECT_EQ(atomData[13].componentId, 1);
  EXPECT_EQ(atomData[14].componentId, 1);
  EXPECT_EQ(atomData[15].componentId, 1);
  EXPECT_EQ(atomData[16].componentId, 1);

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
  atomData = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomData[0].moleculeId, 0);
  EXPECT_EQ(atomData[1].moleculeId, 1);
  EXPECT_EQ(atomData[2].moleculeId, 2);
  EXPECT_EQ(atomData[3].moleculeId, 3);
  EXPECT_EQ(atomData[4].moleculeId, 0);
  EXPECT_EQ(atomData[5].moleculeId, 0);
  EXPECT_EQ(atomData[6].moleculeId, 0);
  EXPECT_EQ(atomData[7].moleculeId, 1);
  EXPECT_EQ(atomData[8].moleculeId, 1);
  EXPECT_EQ(atomData[9].moleculeId, 1);
  EXPECT_EQ(atomData[10].moleculeId, 2);
  EXPECT_EQ(atomData[11].moleculeId, 2);
  EXPECT_EQ(atomData[12].moleculeId, 2);
  EXPECT_EQ(atomData[13].moleculeId, 3);
  EXPECT_EQ(atomData[14].moleculeId, 3);
  EXPECT_EQ(atomData[15].moleculeId, 3);

  EXPECT_EQ(atomData[0].componentId, 0);
  EXPECT_EQ(atomData[1].componentId, 0);
  EXPECT_EQ(atomData[2].componentId, 0);
  EXPECT_EQ(atomData[3].componentId, 0);
  EXPECT_EQ(atomData[4].componentId, 1);
  EXPECT_EQ(atomData[5].componentId, 1);
  EXPECT_EQ(atomData[6].componentId, 1);
  EXPECT_EQ(atomData[7].componentId, 1);
  EXPECT_EQ(atomData[8].componentId, 1);
  EXPECT_EQ(atomData[9].componentId, 1);
  EXPECT_EQ(atomData[10].componentId, 1);
  EXPECT_EQ(atomData[11].componentId, 1);
  EXPECT_EQ(atomData[12].componentId, 1);
  EXPECT_EQ(atomData[13].componentId, 1);
  EXPECT_EQ(atomData[14].componentId, 1);
  EXPECT_EQ(atomData[15].componentId, 1);

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
  atomData = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomData[0].moleculeId, 0);
  EXPECT_EQ(atomData[1].moleculeId, 1);
  EXPECT_EQ(atomData[2].moleculeId, 2);
  EXPECT_EQ(atomData[3].moleculeId, 0);
  EXPECT_EQ(atomData[4].moleculeId, 0);
  EXPECT_EQ(atomData[5].moleculeId, 0);
  EXPECT_EQ(atomData[6].moleculeId, 1);
  EXPECT_EQ(atomData[7].moleculeId, 1);
  EXPECT_EQ(atomData[8].moleculeId, 1);
  EXPECT_EQ(atomData[9].moleculeId, 2);
  EXPECT_EQ(atomData[10].moleculeId, 2);
  EXPECT_EQ(atomData[11].moleculeId, 2);
  EXPECT_EQ(atomData[12].moleculeId, 3);
  EXPECT_EQ(atomData[13].moleculeId, 3);
  EXPECT_EQ(atomData[14].moleculeId, 3);

  EXPECT_EQ(atomData[0].componentId, 0);
  EXPECT_EQ(atomData[1].componentId, 0);
  EXPECT_EQ(atomData[2].componentId, 0);
  EXPECT_EQ(atomData[3].componentId, 1);
  EXPECT_EQ(atomData[4].componentId, 1);
  EXPECT_EQ(atomData[5].componentId, 1);
  EXPECT_EQ(atomData[6].componentId, 1);
  EXPECT_EQ(atomData[7].componentId, 1);
  EXPECT_EQ(atomData[8].componentId, 1);
  EXPECT_EQ(atomData[9].componentId, 1);
  EXPECT_EQ(atomData[10].componentId, 1);
  EXPECT_EQ(atomData[11].componentId, 1);
  EXPECT_EQ(atomData[12].componentId, 1);
  EXPECT_EQ(atomData[13].componentId, 1);
  EXPECT_EQ(atomData[14].componentId, 1);

  // insert new CO2
  system.insertMolecule(1, Molecule(),
                        {Atom(double3(0.0, 0.0, 1.149), -0.3256, 1.0, 0, 4, 1, false, false),
                         Atom(double3(0.0, 0.0, 0.0), 0.6512, 1.0, 0, 3, 1, false, false),
                         Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, false, false)});

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 3);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 3);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 0);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 5);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 5);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 0);

  // refresh span
  atomData = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomData[0].moleculeId, 0);
  EXPECT_EQ(atomData[1].moleculeId, 1);
  EXPECT_EQ(atomData[2].moleculeId, 2);
  EXPECT_EQ(atomData[3].moleculeId, 0);
  EXPECT_EQ(atomData[4].moleculeId, 0);
  EXPECT_EQ(atomData[5].moleculeId, 0);
  EXPECT_EQ(atomData[6].moleculeId, 1);
  EXPECT_EQ(atomData[7].moleculeId, 1);
  EXPECT_EQ(atomData[8].moleculeId, 1);
  EXPECT_EQ(atomData[9].moleculeId, 2);
  EXPECT_EQ(atomData[10].moleculeId, 2);
  EXPECT_EQ(atomData[11].moleculeId, 2);
  EXPECT_EQ(atomData[12].moleculeId, 3);
  EXPECT_EQ(atomData[13].moleculeId, 3);
  EXPECT_EQ(atomData[14].moleculeId, 3);
  EXPECT_EQ(atomData[15].moleculeId, 4);
  EXPECT_EQ(atomData[16].moleculeId, 4);
  EXPECT_EQ(atomData[17].moleculeId, 4);

  EXPECT_EQ(atomData[0].componentId, 0);
  EXPECT_EQ(atomData[1].componentId, 0);
  EXPECT_EQ(atomData[2].componentId, 0);
  EXPECT_EQ(atomData[3].componentId, 1);
  EXPECT_EQ(atomData[4].componentId, 1);
  EXPECT_EQ(atomData[5].componentId, 1);
  EXPECT_EQ(atomData[6].componentId, 1);
  EXPECT_EQ(atomData[7].componentId, 1);
  EXPECT_EQ(atomData[8].componentId, 1);
  EXPECT_EQ(atomData[9].componentId, 1);
  EXPECT_EQ(atomData[10].componentId, 1);
  EXPECT_EQ(atomData[11].componentId, 1);
  EXPECT_EQ(atomData[12].componentId, 1);
  EXPECT_EQ(atomData[13].componentId, 1);
  EXPECT_EQ(atomData[14].componentId, 1);
  EXPECT_EQ(atomData[15].componentId, 1);
  EXPECT_EQ(atomData[16].componentId, 1);
  EXPECT_EQ(atomData[17].componentId, 1);

  // insert new Methane
  system.insertMolecule(0, Molecule(), {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 2, 0, false, false)});

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 4);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 4);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 0);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 5);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 5);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 0);

  // refresh span
  atomData = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomData[0].moleculeId, 0);
  EXPECT_EQ(atomData[1].moleculeId, 1);
  EXPECT_EQ(atomData[2].moleculeId, 2);
  EXPECT_EQ(atomData[3].moleculeId, 3);
  EXPECT_EQ(atomData[4].moleculeId, 0);
  EXPECT_EQ(atomData[5].moleculeId, 0);
  EXPECT_EQ(atomData[6].moleculeId, 0);
  EXPECT_EQ(atomData[7].moleculeId, 1);
  EXPECT_EQ(atomData[8].moleculeId, 1);
  EXPECT_EQ(atomData[9].moleculeId, 1);
  EXPECT_EQ(atomData[10].moleculeId, 2);
  EXPECT_EQ(atomData[11].moleculeId, 2);
  EXPECT_EQ(atomData[12].moleculeId, 2);
  EXPECT_EQ(atomData[13].moleculeId, 3);
  EXPECT_EQ(atomData[14].moleculeId, 3);
  EXPECT_EQ(atomData[15].moleculeId, 3);
  EXPECT_EQ(atomData[16].moleculeId, 4);
  EXPECT_EQ(atomData[17].moleculeId, 4);
  EXPECT_EQ(atomData[18].moleculeId, 4);

  EXPECT_EQ(atomData[0].componentId, 0);
  EXPECT_EQ(atomData[1].componentId, 0);
  EXPECT_EQ(atomData[2].componentId, 0);
  EXPECT_EQ(atomData[3].componentId, 0);
  EXPECT_EQ(atomData[4].componentId, 1);
  EXPECT_EQ(atomData[5].componentId, 1);
  EXPECT_EQ(atomData[6].componentId, 1);
  EXPECT_EQ(atomData[7].componentId, 1);
  EXPECT_EQ(atomData[8].componentId, 1);
  EXPECT_EQ(atomData[9].componentId, 1);
  EXPECT_EQ(atomData[10].componentId, 1);
  EXPECT_EQ(atomData[11].componentId, 1);
  EXPECT_EQ(atomData[12].componentId, 1);
  EXPECT_EQ(atomData[13].componentId, 1);
  EXPECT_EQ(atomData[14].componentId, 1);
  EXPECT_EQ(atomData[15].componentId, 1);
  EXPECT_EQ(atomData[16].componentId, 1);
  EXPECT_EQ(atomData[17].componentId, 1);
  EXPECT_EQ(atomData[18].componentId, 1);

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
  atomData = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomData[0].moleculeId, 0);
  EXPECT_EQ(atomData[1].moleculeId, 1);
  EXPECT_EQ(atomData[2].moleculeId, 2);
  EXPECT_EQ(atomData[3].moleculeId, 3);
  EXPECT_EQ(atomData[4].moleculeId, 0);
  EXPECT_EQ(atomData[5].moleculeId, 0);
  EXPECT_EQ(atomData[6].moleculeId, 0);
  EXPECT_EQ(atomData[7].moleculeId, 1);
  EXPECT_EQ(atomData[8].moleculeId, 1);
  EXPECT_EQ(atomData[9].moleculeId, 1);
  EXPECT_EQ(atomData[10].moleculeId, 2);
  EXPECT_EQ(atomData[11].moleculeId, 2);
  EXPECT_EQ(atomData[12].moleculeId, 2);
  EXPECT_EQ(atomData[13].moleculeId, 3);
  EXPECT_EQ(atomData[14].moleculeId, 3);
  EXPECT_EQ(atomData[15].moleculeId, 3);

  EXPECT_EQ(atomData[0].componentId, 0);
  EXPECT_EQ(atomData[1].componentId, 0);
  EXPECT_EQ(atomData[2].componentId, 0);
  EXPECT_EQ(atomData[3].componentId, 0);
  EXPECT_EQ(atomData[4].componentId, 1);
  EXPECT_EQ(atomData[5].componentId, 1);
  EXPECT_EQ(atomData[6].componentId, 1);
  EXPECT_EQ(atomData[7].componentId, 1);
  EXPECT_EQ(atomData[8].componentId, 1);
  EXPECT_EQ(atomData[9].componentId, 1);
  EXPECT_EQ(atomData[10].componentId, 1);
  EXPECT_EQ(atomData[11].componentId, 1);
  EXPECT_EQ(atomData[12].componentId, 1);
  EXPECT_EQ(atomData[13].componentId, 1);
  EXPECT_EQ(atomData[14].componentId, 1);
  EXPECT_EQ(atomData[15].componentId, 1);

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
  atomData = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomData[0].moleculeId, 0);
  EXPECT_EQ(atomData[1].moleculeId, 1);
  EXPECT_EQ(atomData[2].moleculeId, 2);
  EXPECT_EQ(atomData[3].moleculeId, 3);
  EXPECT_EQ(atomData[4].moleculeId, 0);
  EXPECT_EQ(atomData[5].moleculeId, 0);
  EXPECT_EQ(atomData[6].moleculeId, 0);
  EXPECT_EQ(atomData[7].moleculeId, 1);
  EXPECT_EQ(atomData[8].moleculeId, 1);
  EXPECT_EQ(atomData[9].moleculeId, 1);
  EXPECT_EQ(atomData[10].moleculeId, 2);
  EXPECT_EQ(atomData[11].moleculeId, 2);
  EXPECT_EQ(atomData[12].moleculeId, 2);

  EXPECT_EQ(atomData[0].componentId, 0);
  EXPECT_EQ(atomData[1].componentId, 0);
  EXPECT_EQ(atomData[2].componentId, 0);
  EXPECT_EQ(atomData[3].componentId, 0);
  EXPECT_EQ(atomData[4].componentId, 1);
  EXPECT_EQ(atomData[5].componentId, 1);
  EXPECT_EQ(atomData[6].componentId, 1);
  EXPECT_EQ(atomData[7].componentId, 1);
  EXPECT_EQ(atomData[8].componentId, 1);
  EXPECT_EQ(atomData[9].componentId, 1);
  EXPECT_EQ(atomData[10].componentId, 1);
  EXPECT_EQ(atomData[11].componentId, 1);
  EXPECT_EQ(atomData[12].componentId, 1);

  // insert new Methane
  system.insertMolecule(0, Molecule(), {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 2, 0, false, false)});

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 0);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 3);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 3);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 0);

  // refresh span
  atomData = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomData[0].moleculeId, 0);
  EXPECT_EQ(atomData[1].moleculeId, 1);
  EXPECT_EQ(atomData[2].moleculeId, 2);
  EXPECT_EQ(atomData[3].moleculeId, 3);
  EXPECT_EQ(atomData[4].moleculeId, 4);
  EXPECT_EQ(atomData[5].moleculeId, 0);
  EXPECT_EQ(atomData[6].moleculeId, 0);
  EXPECT_EQ(atomData[7].moleculeId, 0);
  EXPECT_EQ(atomData[8].moleculeId, 1);
  EXPECT_EQ(atomData[9].moleculeId, 1);
  EXPECT_EQ(atomData[10].moleculeId, 1);
  EXPECT_EQ(atomData[11].moleculeId, 2);
  EXPECT_EQ(atomData[12].moleculeId, 2);
  EXPECT_EQ(atomData[13].moleculeId, 2);

  EXPECT_EQ(atomData[0].componentId, 0);
  EXPECT_EQ(atomData[1].componentId, 0);
  EXPECT_EQ(atomData[2].componentId, 0);
  EXPECT_EQ(atomData[3].componentId, 0);
  EXPECT_EQ(atomData[4].componentId, 0);
  EXPECT_EQ(atomData[5].componentId, 1);
  EXPECT_EQ(atomData[6].componentId, 1);
  EXPECT_EQ(atomData[7].componentId, 1);
  EXPECT_EQ(atomData[8].componentId, 1);
  EXPECT_EQ(atomData[9].componentId, 1);
  EXPECT_EQ(atomData[10].componentId, 1);
  EXPECT_EQ(atomData[11].componentId, 1);
  EXPECT_EQ(atomData[12].componentId, 1);
  EXPECT_EQ(atomData[13].componentId, 1);

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
  atomData = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomData[0].moleculeId, 0);
  EXPECT_EQ(atomData[1].moleculeId, 1);
  EXPECT_EQ(atomData[2].moleculeId, 2);
  EXPECT_EQ(atomData[3].moleculeId, 3);
  EXPECT_EQ(atomData[4].moleculeId, 0);
  EXPECT_EQ(atomData[5].moleculeId, 0);
  EXPECT_EQ(atomData[6].moleculeId, 0);
  EXPECT_EQ(atomData[7].moleculeId, 1);
  EXPECT_EQ(atomData[8].moleculeId, 1);
  EXPECT_EQ(atomData[9].moleculeId, 1);
  EXPECT_EQ(atomData[10].moleculeId, 2);
  EXPECT_EQ(atomData[11].moleculeId, 2);
  EXPECT_EQ(atomData[12].moleculeId, 2);

  EXPECT_EQ(atomData[0].componentId, 0);
  EXPECT_EQ(atomData[1].componentId, 0);
  EXPECT_EQ(atomData[2].componentId, 0);
  EXPECT_EQ(atomData[3].componentId, 0);
  EXPECT_EQ(atomData[4].componentId, 1);
  EXPECT_EQ(atomData[5].componentId, 1);
  EXPECT_EQ(atomData[6].componentId, 1);
  EXPECT_EQ(atomData[7].componentId, 1);
  EXPECT_EQ(atomData[8].componentId, 1);
  EXPECT_EQ(atomData[9].componentId, 1);
  EXPECT_EQ(atomData[10].componentId, 1);
  EXPECT_EQ(atomData[11].componentId, 1);
  EXPECT_EQ(atomData[12].componentId, 1);

  // insert new CO2
  system.insertMolecule(1, Molecule(),
                        {Atom(double3(0.0, 0.0, 1.149), -0.3256, 1.0, 0, 4, 1, false, false),
                         Atom(double3(0.0, 0.0, 0.0), 0.6512, 1.0, 0, 3, 1, false, false),
                         Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, false, false)});

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 4);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 4);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 0);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 4);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 4);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 0);

  // refresh span
  atomData = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomData[0].moleculeId, 0);
  EXPECT_EQ(atomData[1].moleculeId, 1);
  EXPECT_EQ(atomData[2].moleculeId, 2);
  EXPECT_EQ(atomData[3].moleculeId, 3);
  EXPECT_EQ(atomData[4].moleculeId, 0);
  EXPECT_EQ(atomData[5].moleculeId, 0);
  EXPECT_EQ(atomData[6].moleculeId, 0);
  EXPECT_EQ(atomData[7].moleculeId, 1);
  EXPECT_EQ(atomData[8].moleculeId, 1);
  EXPECT_EQ(atomData[9].moleculeId, 1);
  EXPECT_EQ(atomData[10].moleculeId, 2);
  EXPECT_EQ(atomData[11].moleculeId, 2);
  EXPECT_EQ(atomData[12].moleculeId, 2);
  EXPECT_EQ(atomData[13].moleculeId, 3);
  EXPECT_EQ(atomData[14].moleculeId, 3);
  EXPECT_EQ(atomData[15].moleculeId, 3);

  EXPECT_EQ(atomData[0].componentId, 0);
  EXPECT_EQ(atomData[1].componentId, 0);
  EXPECT_EQ(atomData[2].componentId, 0);
  EXPECT_EQ(atomData[3].componentId, 0);
  EXPECT_EQ(atomData[4].componentId, 1);
  EXPECT_EQ(atomData[5].componentId, 1);
  EXPECT_EQ(atomData[6].componentId, 1);
  EXPECT_EQ(atomData[7].componentId, 1);
  EXPECT_EQ(atomData[8].componentId, 1);
  EXPECT_EQ(atomData[9].componentId, 1);
  EXPECT_EQ(atomData[10].componentId, 1);
  EXPECT_EQ(atomData[11].componentId, 1);
  EXPECT_EQ(atomData[12].componentId, 1);
  EXPECT_EQ(atomData[13].componentId, 1);
  EXPECT_EQ(atomData[14].componentId, 1);
  EXPECT_EQ(atomData[15].componentId, 1);

  // insert new Methane
  system.insertMolecule(0, Molecule(), {Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 2, 0, false, false)});

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 0);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 4);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 4);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 0);

  // refresh span
  atomData = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomData[0].moleculeId, 0);
  EXPECT_EQ(atomData[1].moleculeId, 1);
  EXPECT_EQ(atomData[2].moleculeId, 2);
  EXPECT_EQ(atomData[3].moleculeId, 3);
  EXPECT_EQ(atomData[4].moleculeId, 4);
  EXPECT_EQ(atomData[5].moleculeId, 0);
  EXPECT_EQ(atomData[6].moleculeId, 0);
  EXPECT_EQ(atomData[7].moleculeId, 0);
  EXPECT_EQ(atomData[8].moleculeId, 1);
  EXPECT_EQ(atomData[9].moleculeId, 1);
  EXPECT_EQ(atomData[10].moleculeId, 1);
  EXPECT_EQ(atomData[11].moleculeId, 2);
  EXPECT_EQ(atomData[12].moleculeId, 2);
  EXPECT_EQ(atomData[13].moleculeId, 2);
  EXPECT_EQ(atomData[14].moleculeId, 3);
  EXPECT_EQ(atomData[15].moleculeId, 3);
  EXPECT_EQ(atomData[16].moleculeId, 3);

  EXPECT_EQ(atomData[0].componentId, 0);
  EXPECT_EQ(atomData[1].componentId, 0);
  EXPECT_EQ(atomData[2].componentId, 0);
  EXPECT_EQ(atomData[3].componentId, 0);
  EXPECT_EQ(atomData[4].componentId, 0);
  EXPECT_EQ(atomData[5].componentId, 1);
  EXPECT_EQ(atomData[6].componentId, 1);
  EXPECT_EQ(atomData[7].componentId, 1);
  EXPECT_EQ(atomData[8].componentId, 1);
  EXPECT_EQ(atomData[9].componentId, 1);
  EXPECT_EQ(atomData[10].componentId, 1);
  EXPECT_EQ(atomData[11].componentId, 1);
  EXPECT_EQ(atomData[12].componentId, 1);
  EXPECT_EQ(atomData[13].componentId, 1);
  EXPECT_EQ(atomData[14].componentId, 1);
  EXPECT_EQ(atomData[15].componentId, 1);
  EXPECT_EQ(atomData[16].componentId, 1);

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
  atomData = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomData[0].moleculeId, 0);
  EXPECT_EQ(atomData[1].moleculeId, 1);
  EXPECT_EQ(atomData[2].moleculeId, 2);
  EXPECT_EQ(atomData[3].moleculeId, 3);
  EXPECT_EQ(atomData[4].moleculeId, 4);
  EXPECT_EQ(atomData[5].moleculeId, 0);
  EXPECT_EQ(atomData[6].moleculeId, 0);
  EXPECT_EQ(atomData[7].moleculeId, 0);
  EXPECT_EQ(atomData[8].moleculeId, 1);
  EXPECT_EQ(atomData[9].moleculeId, 1);
  EXPECT_EQ(atomData[10].moleculeId, 1);
  EXPECT_EQ(atomData[11].moleculeId, 2);
  EXPECT_EQ(atomData[12].moleculeId, 2);
  EXPECT_EQ(atomData[13].moleculeId, 2);

  EXPECT_EQ(atomData[0].componentId, 0);
  EXPECT_EQ(atomData[1].componentId, 0);
  EXPECT_EQ(atomData[2].componentId, 0);
  EXPECT_EQ(atomData[3].componentId, 0);
  EXPECT_EQ(atomData[4].componentId, 0);
  EXPECT_EQ(atomData[5].componentId, 1);
  EXPECT_EQ(atomData[6].componentId, 1);
  EXPECT_EQ(atomData[7].componentId, 1);
  EXPECT_EQ(atomData[8].componentId, 1);
  EXPECT_EQ(atomData[9].componentId, 1);
  EXPECT_EQ(atomData[10].componentId, 1);
  EXPECT_EQ(atomData[11].componentId, 1);
  EXPECT_EQ(atomData[12].componentId, 1);
  EXPECT_EQ(atomData[13].componentId, 1);

  // component 0: Methane
  EXPECT_EQ(atomData[0].componentId, 0);
  EXPECT_EQ(atomData[1].componentId, 0);
  EXPECT_EQ(atomData[2].componentId, 0);
  EXPECT_EQ(atomData[3].componentId, 0);
  EXPECT_EQ(atomData[4].componentId, 0);
  // component 1: CO2
  EXPECT_EQ(atomData[5].componentId, 1);
  EXPECT_EQ(atomData[6].componentId, 1);
  EXPECT_EQ(atomData[7].componentId, 1);
  EXPECT_EQ(atomData[8].componentId, 1);
  EXPECT_EQ(atomData[9].componentId, 1);
  EXPECT_EQ(atomData[10].componentId, 1);
  EXPECT_EQ(atomData[11].componentId, 1);
  EXPECT_EQ(atomData[12].componentId, 1);
  EXPECT_EQ(atomData[13].componentId, 1);

  // component 0: Methane
  EXPECT_EQ(atomData[0].moleculeId, 0);
  EXPECT_EQ(atomData[1].moleculeId, 1);
  EXPECT_EQ(atomData[2].moleculeId, 2);
  EXPECT_EQ(atomData[3].moleculeId, 3);
  EXPECT_EQ(atomData[4].moleculeId, 4);
  // component 1: CO2
  EXPECT_EQ(atomData[5].moleculeId, 0);
  EXPECT_EQ(atomData[6].moleculeId, 0);
  EXPECT_EQ(atomData[7].moleculeId, 0);
  EXPECT_EQ(atomData[8].moleculeId, 1);
  EXPECT_EQ(atomData[9].moleculeId, 1);
  EXPECT_EQ(atomData[10].moleculeId, 1);
  EXPECT_EQ(atomData[11].moleculeId, 2);
  EXPECT_EQ(atomData[12].moleculeId, 2);
  EXPECT_EQ(atomData[13].moleculeId, 2);

  // component 0: Methane
  EXPECT_EQ(atomData[0].type, 2);
  EXPECT_EQ(atomData[1].type, 2);
  EXPECT_EQ(atomData[2].type, 2);
  EXPECT_EQ(atomData[3].type, 2);
  EXPECT_EQ(atomData[4].type, 2);
  // component 1: CO2
  EXPECT_EQ(atomData[5].type, 4);
  EXPECT_EQ(atomData[6].type, 3);
  EXPECT_EQ(atomData[7].type, 4);
  EXPECT_EQ(atomData[8].type, 4);
  EXPECT_EQ(atomData[9].type, 3);
  EXPECT_EQ(atomData[10].type, 4);
  EXPECT_EQ(atomData[11].type, 4);
  EXPECT_EQ(atomData[12].type, 3);
  EXPECT_EQ(atomData[13].type, 4);

  // component 0: Methane
  EXPECT_NEAR(atomData[0].charge, 0.0, 1e-6);
  EXPECT_NEAR(atomData[1].charge, 0.0, 1e-6);
  EXPECT_NEAR(atomData[2].charge, 0.0, 1e-6);
  EXPECT_NEAR(atomData[3].charge, 0.0, 1e-6);
  EXPECT_NEAR(atomData[4].charge, 0.0, 1e-6);
  // component 1: CO2
  EXPECT_NEAR(atomData[5].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomData[6].charge, 0.6512, 1e-6);
  EXPECT_NEAR(atomData[7].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomData[8].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomData[9].charge, 0.6512, 1e-6);
  EXPECT_NEAR(atomData[10].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomData[11].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomData[12].charge, 0.6512, 1e-6);
  EXPECT_NEAR(atomData[13].charge, -0.3256, 1e-6);
}
