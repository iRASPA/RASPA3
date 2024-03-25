#include <gtest/gtest.h>

import <vector>;
import <tuple>;
import <algorithm>;
import <span>;
import <ranges>;

import double3;
import double3x3;

import atom;
import forcefield;
import framework;
import component;
import system;
import simulationbox;

TEST(insertion_deletion, methane_number_of_molecules_per_component)
{
  ForceField forceField = ForceField(
    { PseudoAtom("CH4", 16.04246, 0.0, 6, false) },
    { VDWParameters(158.5 / 1.2027242847, 3.72) },
    ForceField::MixingRule::Lorentz_Berthelot,
    12.0,
    true,
    false);
  Component c = Component(0,
    "methane",
    16.04246,
    190.564, 45599200, 0.01142,
    { Atom(double3(0.0, 0.0,  0.0),    0.0, 1.0, 0, 0, 0) }, 5, 21);

  System system = System(0, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, forceField, {}, { c }, { 20 }, 5);

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 20);
  EXPECT_EQ(system.numberOfPseudoAtoms[0][0], 20);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();

  for(size_t i = 0; i != atomPositions.size(); ++i)
  {
    EXPECT_NEAR(atomPositions[i].charge, 0.0, 1e-6);
    EXPECT_NEAR(atomPositions[i].scalingVDW, 1.0, 1e-6);
    EXPECT_NEAR(atomPositions[i].scalingCoulomb, 1.0, 1e-6);
    EXPECT_EQ(atomPositions[i].moleculeId, i);
    EXPECT_EQ(atomPositions[i].type, 0);
    EXPECT_EQ(atomPositions[i].componentId, 0);
    EXPECT_EQ(atomPositions[i].groupId, 0);
  }
}

TEST(insertion_deletion, CO2_number_of_molecules_per_component)
{
  ForceField forceField = ForceField(
  { PseudoAtom("Si",    28.0855,   2.05,  14, false),
    PseudoAtom("O",     15.999,   -1.025,  8, false),
    PseudoAtom("CH4",   16.04246,  0.0,    6, false),
    PseudoAtom("C_co2", 12.0,      0.6512, 6, false),
    PseudoAtom("O_co2", 15.9994,  -0.3256, 8, false),
  },
  { VDWParameters(22.0 / 1.2027242847, 2.30),
    VDWParameters(53.0 / 1.2027242847, 3.3),
    VDWParameters(158.5 / 1.2027242847, 3.72),
    VDWParameters(29.933 / 1.2027242847, 2.745),
    VDWParameters(85.671 / 1.2027242847, 3.017)
  },
  ForceField::MixingRule::Lorentz_Berthelot,
  12.0,
  true,
  false);

  Component c = Component(0,
    "CO2",
    43.9988,
    304.1282, 7377300.0, 0.22394,
    {
       Atom(double3(0.0, 0.0,  1.149), -0.3256, 1.0, 4, 1, 0),
       Atom(double3(0.0, 0.0,  0.0  ),  0.6512, 1.0, 3, 1, 0),
       Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 4, 1, 0)
    }, 5, 21);

  System system = System(0, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, forceField, { }, { c }, { 3 }, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();

  for(size_t i = 0; i != atomPositions.size(); ++i)
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
  ForceField forceField = ForceField(
  { PseudoAtom("Si",    28.0855,   2.05,  14, false),
    PseudoAtom("O",     15.999,   -1.025,  8, false),
    PseudoAtom("CH4",   16.04246,  0.0,    6, false),
    PseudoAtom("C_co2", 12.0,      0.6512, 6, false),
    PseudoAtom("O_co2", 15.9994,  -0.3256, 8, false),
  },
  { VDWParameters(22.0 / 1.2027242847, 2.30),
    VDWParameters(53.0 / 1.2027242847, 3.3),
    VDWParameters(158.5 / 1.2027242847, 3.72),
    VDWParameters(29.933 / 1.2027242847, 2.745),
    VDWParameters(85.671 / 1.2027242847, 3.017)
  },
  ForceField::MixingRule::Lorentz_Berthelot,
  12.0,
  true,
  false);

  Component methane = Component(0,
    "methane",
    16.04246,
    190.564, 45599200, 0.01142,
    { Atom(double3(0.0, 0.0,  0.0),    0.0, 1.0, 2, 0, 0) }, 5, 21);

  Component co2 = Component(1,
    "CO2",
    43.9988,
    304.1282, 7377300.0, 0.22394,
    {  // double3 position, double charge, double lambda, uint16_t type, uint8_t componentId, uint32_t moleculeId
       Atom(double3(0.0, 0.0,  1.149), -0.3256, 1.0, 4, 1, 0),
       Atom(double3(0.0, 0.0,  0.0  ),  0.6512, 1.0, 3, 1, 0),
       Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 4, 1, 0)
    }, 5, 21);

  System system = System(0, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, forceField, {}, { methane, co2 }, { 5, 3 }, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();

  // component 0: Methane
  EXPECT_EQ(atomPositions[ 0].componentId, 0);
  EXPECT_EQ(atomPositions[ 1].componentId, 0);
  EXPECT_EQ(atomPositions[ 2].componentId, 0);
  EXPECT_EQ(atomPositions[ 3].componentId, 0);
  EXPECT_EQ(atomPositions[ 4].componentId, 0);
  // component 1: CO2
  EXPECT_EQ(atomPositions[ 5].componentId, 1);
  EXPECT_EQ(atomPositions[ 6].componentId, 1);
  EXPECT_EQ(atomPositions[ 7].componentId, 1);
  EXPECT_EQ(atomPositions[ 8].componentId, 1);
  EXPECT_EQ(atomPositions[ 9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);
  EXPECT_EQ(atomPositions[13].componentId, 1);

  // component 0: Methane
  EXPECT_EQ(atomPositions[ 0].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 1].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 2].moleculeId, 2);
  EXPECT_EQ(atomPositions[ 3].moleculeId, 3);
  EXPECT_EQ(atomPositions[ 4].moleculeId, 4);
  // component 1: CO2
  EXPECT_EQ(atomPositions[ 5].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 6].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 7].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 8].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 9].moleculeId, 1);
  EXPECT_EQ(atomPositions[10].moleculeId, 1);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 2);
  EXPECT_EQ(atomPositions[13].moleculeId, 2);

  // component 0: Methane
  EXPECT_EQ(atomPositions[ 0].type, 2);
  EXPECT_EQ(atomPositions[ 1].type, 2);
  EXPECT_EQ(atomPositions[ 2].type, 2);
  EXPECT_EQ(atomPositions[ 3].type, 2);
  EXPECT_EQ(atomPositions[ 4].type, 2);
  // component 1: CO2
  EXPECT_EQ(atomPositions[ 5].type, 4);
  EXPECT_EQ(atomPositions[ 6].type, 3);
  EXPECT_EQ(atomPositions[ 7].type, 4);
  EXPECT_EQ(atomPositions[ 8].type, 4);
  EXPECT_EQ(atomPositions[ 9].type, 3);
  EXPECT_EQ(atomPositions[10].type, 4);
  EXPECT_EQ(atomPositions[11].type, 4);
  EXPECT_EQ(atomPositions[12].type, 3);
  EXPECT_EQ(atomPositions[13].type, 4);

  // component 0: Methane
  EXPECT_NEAR(atomPositions[ 0].charge, 0.0, 1e-6);
  EXPECT_NEAR(atomPositions[ 1].charge, 0.0, 1e-6);
  EXPECT_NEAR(atomPositions[ 2].charge, 0.0, 1e-6);
  EXPECT_NEAR(atomPositions[ 3].charge, 0.0, 1e-6);
  EXPECT_NEAR(atomPositions[ 4].charge, 0.0, 1e-6);
  // component 1: CO2
  EXPECT_NEAR(atomPositions[ 5].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomPositions[ 6].charge, 0.6512, 1e-6);
  EXPECT_NEAR(atomPositions[ 7].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomPositions[ 8].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomPositions[ 9].charge, 0.6512, 1e-6);
  EXPECT_NEAR(atomPositions[10].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomPositions[11].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomPositions[12].charge, 0.6512, 1e-6);
  EXPECT_NEAR(atomPositions[13].charge, -0.3256, 1e-6);
}

TEST(insertion_deletion, Dynamic_CO2_Methane_number_of_molecules_per_component)
{
  ForceField forceField = ForceField(
  { PseudoAtom("Si",    28.0855,   2.05,  14, false),
    PseudoAtom("O",     15.999,   -1.025,  8, false),
    PseudoAtom("CH4",   16.04246,  0.0,    6, false),
    PseudoAtom("C_co2", 12.0,      0.6512, 6, false),
    PseudoAtom("O_co2", 15.9994,  -0.3256, 8, false),
  },
  { VDWParameters(22.0 / 1.2027242847, 2.30),
    VDWParameters(53.0 / 1.2027242847, 3.3),
    VDWParameters(158.5 / 1.2027242847, 3.72),
    VDWParameters(29.933 / 1.2027242847, 2.745),
    VDWParameters(85.671 / 1.2027242847, 3.017)
  },
  ForceField::MixingRule::Lorentz_Berthelot,
  12.0,
  true,
  false);

  Component methane = Component(0,
    "methane",
    16.04246,
    190.564, 45599200, 0.01142,
    { Atom(double3(0.0, 0.0,  0.0),    0.0, 1.0, 2, 0, 0) }, 5, 21);

  Component co2 = Component(1,
    "CO2",
    43.9988,
    304.1282, 7377300.0, 0.22394,
    {  // double3 position, double charge, double lambda, uint16_t type, uint8_t componentId, uint32_t moleculeId
       Atom(double3(0.0, 0.0,  1.149), -0.3256, 1.0, 4, 1, 0),
       Atom(double3(0.0, 0.0,  0.0  ),  0.6512, 1.0, 3, 1, 0),
       Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 4, 1, 0)
    }, 5, 21);

  System system = System(0, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, forceField, {}, { methane, co2 }, { 5, 3 }, 5);

  //insert new CO2
  system.insertMolecule(1, {Atom(double3(0.0, 0.0,  1.149), -0.3256, 1.0, 4, 1, 0),
                            Atom(double3(0.0, 0.0,  0.0  ),  0.6512, 1.0, 3, 1, 0),
                            Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 4, 1, 0)});
  //delete Methane 2
  std::vector<Atom>::iterator iterator_methane1 = system.iteratorForMolecule(0, 2);
  system.deleteMolecule(0, 2, {iterator_methane1, 1});

  //delete Methane 2
  std::vector<Atom>::iterator iterator_methane2 = system.iteratorForMolecule(0, 2);
  system.deleteMolecule(0, 2, {iterator_methane2, 1});

  //insert new CO2
  system.insertMolecule(1, {Atom(double3(0.0, 0.0,  1.149), -0.3256, 1.0, 4, 1, 0),
                            Atom(double3(0.0, 0.0,  0.0  ),  0.6512, 1.0, 3, 1, 0),
                            Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 4, 1, 0)});


  //insert new Methane
  system.insertMolecule(0, {Atom(double3(0.0, 0.0,  0.0),    0.0, 1.0, 2, 0, 0) });

  //delete CO2 1
  std::vector<Atom>::iterator iterator_CO2 = system.iteratorForMolecule(1, 1);
  system.deleteMolecule(1, 1, {iterator_CO2, 3});

  //delete CO2 1
  std::vector<Atom>::iterator iterator_CO2_2 = system.iteratorForMolecule(1, 1);
  system.deleteMolecule(1, 1, {iterator_CO2_2, 3});

  //insert new Methane
  system.insertMolecule(0, {Atom(double3(0.0, 0.0,  0.0),    0.0, 1.0, 2, 0, 0) });



  //delete Methane 2
  std::vector<Atom>::iterator iterator_methane3 = system.iteratorForMolecule(0, 2);
  system.deleteMolecule(0, 2, {iterator_methane3, 1});

  //insert new CO2
  system.insertMolecule(1, {Atom(double3(0.0, 0.0,  1.149), -0.3256, 1.0, 4, 1, 0),
                            Atom(double3(0.0, 0.0,  0.0  ),  0.6512, 1.0, 3, 1, 0),
                            Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 4, 1, 0)});

  //insert new Methane
  system.insertMolecule(0, {Atom(double3(0.0, 0.0,  0.0),    0.0, 1.0, 2, 0, 0) });

  //delete CO2 0
  std::vector<Atom>::iterator iterator_CO2_3 = system.iteratorForMolecule(1, 1);
  system.deleteMolecule(1, 0, {iterator_CO2_3, 3});

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();

  // component 0: Methane
  EXPECT_EQ(atomPositions[ 0].componentId, 0);
  EXPECT_EQ(atomPositions[ 1].componentId, 0);
  EXPECT_EQ(atomPositions[ 2].componentId, 0);
  EXPECT_EQ(atomPositions[ 3].componentId, 0);
  EXPECT_EQ(atomPositions[ 4].componentId, 0);
  // component 1: CO2
  EXPECT_EQ(atomPositions[ 5].componentId, 1);
  EXPECT_EQ(atomPositions[ 6].componentId, 1);
  EXPECT_EQ(atomPositions[ 7].componentId, 1);
  EXPECT_EQ(atomPositions[ 8].componentId, 1);
  EXPECT_EQ(atomPositions[ 9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);
  EXPECT_EQ(atomPositions[13].componentId, 1);

  // component 0: Methane
  EXPECT_EQ(atomPositions[ 0].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 1].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 2].moleculeId, 2);
  EXPECT_EQ(atomPositions[ 3].moleculeId, 3);
  EXPECT_EQ(atomPositions[ 4].moleculeId, 4);
  // component 1: CO2
  EXPECT_EQ(atomPositions[ 5].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 6].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 7].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 8].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 9].moleculeId, 1);
  EXPECT_EQ(atomPositions[10].moleculeId, 1);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 2);
  EXPECT_EQ(atomPositions[13].moleculeId, 2);

  // component 0: Methane
  EXPECT_EQ(atomPositions[ 0].type, 2);
  EXPECT_EQ(atomPositions[ 1].type, 2);
  EXPECT_EQ(atomPositions[ 2].type, 2);
  EXPECT_EQ(atomPositions[ 3].type, 2);
  EXPECT_EQ(atomPositions[ 4].type, 2);
  // component 1: CO2
  EXPECT_EQ(atomPositions[ 5].type, 4);
  EXPECT_EQ(atomPositions[ 6].type, 3);
  EXPECT_EQ(atomPositions[ 7].type, 4);
  EXPECT_EQ(atomPositions[ 8].type, 4);
  EXPECT_EQ(atomPositions[ 9].type, 3);
  EXPECT_EQ(atomPositions[10].type, 4);
  EXPECT_EQ(atomPositions[11].type, 4);
  EXPECT_EQ(atomPositions[12].type, 3);
  EXPECT_EQ(atomPositions[13].type, 4);

  // component 0: Methane
  EXPECT_NEAR(atomPositions[ 0].charge, 0.0, 1e-6);
  EXPECT_NEAR(atomPositions[ 1].charge, 0.0, 1e-6);
  EXPECT_NEAR(atomPositions[ 2].charge, 0.0, 1e-6);
  EXPECT_NEAR(atomPositions[ 3].charge, 0.0, 1e-6);
  EXPECT_NEAR(atomPositions[ 4].charge, 0.0, 1e-6);
  // component 1: CO2
  EXPECT_NEAR(atomPositions[ 5].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomPositions[ 6].charge, 0.6512, 1e-6);
  EXPECT_NEAR(atomPositions[ 7].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomPositions[ 8].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomPositions[ 9].charge, 0.6512, 1e-6);
  EXPECT_NEAR(atomPositions[10].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomPositions[11].charge, -0.3256, 1e-6);
  EXPECT_NEAR(atomPositions[12].charge, 0.6512, 1e-6);
  EXPECT_NEAR(atomPositions[13].charge, -0.3256, 1e-6);
}
