#include <gtest/gtest.h>

#include <vector>
#include <tuple>
#include <algorithm>
#include <span>
#include <ranges>
#include <complex>

import int3;
import double3;
import double3x3;

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
  ForceField forceField = ForceField(
    { PseudoAtom("CH4", 16.04246, 0.0, 6, false) },
    { VDWParameters(158.5, 3.72) },
    ForceField::MixingRule::Lorentz_Berthelot,
    12.0,
    true,
    false);
  Component c = Component(0, forceField, "methane", 190.564, 45599200, 0.01142,
    { // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
      Atom(double3(0.0, 0.0,  0.0),    0.0, 1.0, 0, 0, 0, 0) 
    }, 5, 21);

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
  { VDWParameters(22.0, 2.30),
    VDWParameters(53.0, 3.3),
    VDWParameters(158.5, 3.72),
    VDWParameters(29.933, 2.745),
    VDWParameters(85.671, 3.017)
  },
  ForceField::MixingRule::Lorentz_Berthelot,
  12.0,
  true,
  false);

  Component c = Component(0, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
    { // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
      Atom(double3(0.0, 0.0,  1.149), -0.3256, 1.0, 0, 4, 1, 0),
      Atom(double3(0.0, 0.0,  0.0  ),  0.6512, 1.0, 0, 3, 1, 0),
      Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)
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
  { VDWParameters(22.0, 2.30),
    VDWParameters(53.0, 3.3),
    VDWParameters(158.5, 3.72),
    VDWParameters(29.933, 2.745),
    VDWParameters(85.671, 3.017)
  },
  ForceField::MixingRule::Lorentz_Berthelot,
  12.0,
  true,
  false);

   Framework f = Framework(0, forceField, "MFI_SI", SimulationBox(20.022, 19.899, 13.383),
    292,
    { // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
      Atom(double3(0.42238,  0.0565,  -0.33598), 2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.30716,  0.02772, -0.1893),  2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.27911,  0.06127,  0.0312),  2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.12215,  0.06298,  0.0267),  2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.07128,  0.02722, -0.18551), 2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.18641,  0.05896, -0.32818), 2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.42265, -0.1725,  -0.32718), 2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.30778, -0.13016, -0.18548), 2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.27554, -0.17279,  0.03109), 2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.12058, -0.1731,   0.02979), 2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.07044, -0.13037, -0.182),   2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.18706, -0.17327, -0.31933), 2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.3726,   0.0534,  -0.2442), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.3084,   0.0587,  -0.0789), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.2007,   0.0592,   0.0289), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.0969,   0.0611,  -0.0856), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.1149,   0.0541,  -0.2763), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.2435,   0.0553,  -0.246),  -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.3742,  -0.1561,  -0.2372), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.3085,  -0.1552,  -0.0728), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.198,   -0.1554,   0.0288), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.091,   -0.1614,  -0.0777), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.1169,  -0.1578,  -0.2694), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.2448,  -0.1594,  -0.2422), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.3047,  -0.051,   -0.1866), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.0768,  -0.0519,  -0.1769), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.4161,   0.1276,  -0.3896), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.4086,  -0.0017,  -0.4136), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.402,   -0.1314,  -0.4239), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.1886,   0.1298,  -0.3836), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.194,    0.0007,  -0.4082), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.1951,  -0.1291,  -0.419),  -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(-0.0037,  0.0502,  -0.208),  -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(-0.004,  -0.1528,  -0.2078), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.4192,  -0.25,    -0.354),  -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.1884,  -0.25,    -0.3538), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.2883,  -0.25,     0.0579), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.1085,  -0.25,     0.0611), -1.025, 1.0, 0, 1, 0, 0)
    },
    int3(2, 2, 2));


  Component methane = Component(0, forceField, "methane", 190.564, 45599200, 0.01142,
    { // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
      Atom(double3(0.0, 0.0,  0.0),    0.0, 1.0, 0, 2, 0, 0) 
    }, 5, 21);

  Component co2 = Component(1, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
    {  // double3 position, double charge, double lambda, uint16_t type, uint8_t componentId, uint32_t moleculeId
       Atom(double3(0.0, 0.0,  1.149), -0.3256, 1.0, 0, 4, 1, 0),
       Atom(double3(0.0, 0.0,  0.0  ),  0.6512, 1.0, 0, 3, 1, 0),
       Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)
    }, 5, 21);

  System system = System(0, std::nullopt, 300.0, 1e4, forceField, {f}, { methane, co2 }, { 5, 3 }, 5);

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
  { VDWParameters(22.0, 2.30),
    VDWParameters(53.0, 3.3),
    VDWParameters(158.5, 3.72),
    VDWParameters(29.933, 2.745),
    VDWParameters(85.671, 3.017)
  },
  ForceField::MixingRule::Lorentz_Berthelot,
  12.0,
  true,
  false);

   Framework f = Framework(0, forceField, "MFI_SI", SimulationBox(20.022, 19.899, 13.383),
    292,
    { // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
      Atom(double3(0.42238,  0.0565,  -0.33598), 2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.30716,  0.02772, -0.1893),  2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.27911,  0.06127,  0.0312),  2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.12215,  0.06298,  0.0267),  2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.07128,  0.02722, -0.18551), 2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.18641,  0.05896, -0.32818), 2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.42265, -0.1725,  -0.32718), 2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.30778, -0.13016, -0.18548), 2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.27554, -0.17279,  0.03109), 2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.12058, -0.1731,   0.02979), 2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.07044, -0.13037, -0.182),   2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.18706, -0.17327, -0.31933), 2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.3726,   0.0534,  -0.2442), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.3084,   0.0587,  -0.0789), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.2007,   0.0592,   0.0289), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.0969,   0.0611,  -0.0856), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.1149,   0.0541,  -0.2763), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.2435,   0.0553,  -0.246),  -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.3742,  -0.1561,  -0.2372), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.3085,  -0.1552,  -0.0728), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.198,   -0.1554,   0.0288), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.091,   -0.1614,  -0.0777), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.1169,  -0.1578,  -0.2694), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.2448,  -0.1594,  -0.2422), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.3047,  -0.051,   -0.1866), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.0768,  -0.0519,  -0.1769), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.4161,   0.1276,  -0.3896), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.4086,  -0.0017,  -0.4136), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.402,   -0.1314,  -0.4239), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.1886,   0.1298,  -0.3836), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.194,    0.0007,  -0.4082), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.1951,  -0.1291,  -0.419),  -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(-0.0037,  0.0502,  -0.208),  -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(-0.004,  -0.1528,  -0.2078), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.4192,  -0.25,    -0.354),  -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.1884,  -0.25,    -0.3538), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.2883,  -0.25,     0.0579), -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.1085,  -0.25,     0.0611), -1.025, 1.0, 0, 1, 0, 0)
    },
    int3(2, 2, 2));


  Component methane = Component(0, forceField, "methane", 190.564, 45599200, 0.01142,
    { 
      Atom(double3(0.0, 0.0,  0.0),    0.0, 1.0, 0, 2, 0, 0) 
    }, 5, 21);

  Component co2 = Component(1, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
    { // double3 position, double charge, double lambda, uint16_t type, uint8_t componentId, uint32_t moleculeId
      Atom(double3(0.0, 0.0,  1.149), -0.3256, 1.0, 0, 4, 1, 0),
      Atom(double3(0.0, 0.0,  0.0  ),  0.6512, 1.0, 0, 3, 1, 0),
      Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)
    }, 5, 21);

  System system = System(0, std::nullopt, 300.0, 1e4, forceField, {f}, { methane, co2 }, { 5, 3 }, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 3);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 3);

  // insert new CO2
  system.insertMolecule(1, Molecule(),
                           {Atom(double3(0.0, 0.0,  1.149), -0.3256, 1.0, 0, 4, 1, 0),
                            Atom(double3(0.0, 0.0,  0.0  ),  0.6512, 1.0, 0, 3, 1, 0),
                            Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)});

  // refresh span
  atomPositions = system.spanOfMoleculeAtoms();

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 0);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 4);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 4);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 0);

  EXPECT_EQ(atomPositions[ 0].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 1].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 2].moleculeId, 2);
  EXPECT_EQ(atomPositions[ 3].moleculeId, 3);
  EXPECT_EQ(atomPositions[ 4].moleculeId, 4);
  EXPECT_EQ(atomPositions[ 5].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 6].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 7].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 8].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 9].moleculeId, 1);
  EXPECT_EQ(atomPositions[10].moleculeId, 1);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 2);
  EXPECT_EQ(atomPositions[13].moleculeId, 2);
  EXPECT_EQ(atomPositions[14].moleculeId, 3);
  EXPECT_EQ(atomPositions[15].moleculeId, 3);
  EXPECT_EQ(atomPositions[16].moleculeId, 3);

  EXPECT_EQ(atomPositions[ 0].componentId, 0);
  EXPECT_EQ(atomPositions[ 1].componentId, 0);
  EXPECT_EQ(atomPositions[ 2].componentId, 0);
  EXPECT_EQ(atomPositions[ 3].componentId, 0);
  EXPECT_EQ(atomPositions[ 4].componentId, 0);
  EXPECT_EQ(atomPositions[ 5].componentId, 1);
  EXPECT_EQ(atomPositions[ 6].componentId, 1);
  EXPECT_EQ(atomPositions[ 7].componentId, 1);
  EXPECT_EQ(atomPositions[ 8].componentId, 1);
  EXPECT_EQ(atomPositions[ 9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);
  EXPECT_EQ(atomPositions[13].componentId, 1);
  EXPECT_EQ(atomPositions[14].componentId, 1);
  EXPECT_EQ(atomPositions[15].componentId, 1);
  EXPECT_EQ(atomPositions[16].componentId, 1);

  //delete Methane 2
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

  EXPECT_EQ(atomPositions[ 0].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 1].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 2].moleculeId, 2);
  EXPECT_EQ(atomPositions[ 3].moleculeId, 3);
  EXPECT_EQ(atomPositions[ 4].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 5].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 6].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 7].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 8].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 9].moleculeId, 1);
  EXPECT_EQ(atomPositions[10].moleculeId, 2);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 2);
  EXPECT_EQ(atomPositions[13].moleculeId, 3);
  EXPECT_EQ(atomPositions[14].moleculeId, 3);
  EXPECT_EQ(atomPositions[15].moleculeId, 3);

  EXPECT_EQ(atomPositions[ 0].componentId, 0);
  EXPECT_EQ(atomPositions[ 1].componentId, 0);
  EXPECT_EQ(atomPositions[ 2].componentId, 0);
  EXPECT_EQ(atomPositions[ 3].componentId, 0);
  EXPECT_EQ(atomPositions[ 4].componentId, 1);
  EXPECT_EQ(atomPositions[ 5].componentId, 1);
  EXPECT_EQ(atomPositions[ 6].componentId, 1);
  EXPECT_EQ(atomPositions[ 7].componentId, 1);
  EXPECT_EQ(atomPositions[ 8].componentId, 1);
  EXPECT_EQ(atomPositions[ 9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);
  EXPECT_EQ(atomPositions[13].componentId, 1);
  EXPECT_EQ(atomPositions[14].componentId, 1);
  EXPECT_EQ(atomPositions[15].componentId, 1);

  //delete Methane 2
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

  EXPECT_EQ(atomPositions[ 0].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 1].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 2].moleculeId, 2);
  EXPECT_EQ(atomPositions[ 3].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 4].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 5].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 6].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 7].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 8].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 9].moleculeId, 2);
  EXPECT_EQ(atomPositions[10].moleculeId, 2);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 3);
  EXPECT_EQ(atomPositions[13].moleculeId, 3);
  EXPECT_EQ(atomPositions[14].moleculeId, 3);

  EXPECT_EQ(atomPositions[ 0].componentId, 0);
  EXPECT_EQ(atomPositions[ 1].componentId, 0);
  EXPECT_EQ(atomPositions[ 2].componentId, 0);
  EXPECT_EQ(atomPositions[ 3].componentId, 1);
  EXPECT_EQ(atomPositions[ 4].componentId, 1);
  EXPECT_EQ(atomPositions[ 5].componentId, 1);
  EXPECT_EQ(atomPositions[ 6].componentId, 1);
  EXPECT_EQ(atomPositions[ 7].componentId, 1);
  EXPECT_EQ(atomPositions[ 8].componentId, 1);
  EXPECT_EQ(atomPositions[ 9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);
  EXPECT_EQ(atomPositions[13].componentId, 1);
  EXPECT_EQ(atomPositions[14].componentId, 1);


  //insert new CO2
  system.insertMolecule(1, Molecule(),
                           {Atom(double3(0.0, 0.0,  1.149), -0.3256, 1.0, 0, 4, 1, 0),
                            Atom(double3(0.0, 0.0,  0.0  ),  0.6512, 1.0, 0, 3, 1, 0),
                            Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)});

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 3);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 3);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 0);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 5);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 5);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 0);

  // refresh span
  atomPositions = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomPositions[ 0].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 1].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 2].moleculeId, 2);
  EXPECT_EQ(atomPositions[ 3].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 4].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 5].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 6].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 7].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 8].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 9].moleculeId, 2);
  EXPECT_EQ(atomPositions[10].moleculeId, 2);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 3);
  EXPECT_EQ(atomPositions[13].moleculeId, 3);
  EXPECT_EQ(atomPositions[14].moleculeId, 3);
  EXPECT_EQ(atomPositions[15].moleculeId, 4);
  EXPECT_EQ(atomPositions[16].moleculeId, 4);
  EXPECT_EQ(atomPositions[17].moleculeId, 4);

  EXPECT_EQ(atomPositions[ 0].componentId, 0);
  EXPECT_EQ(atomPositions[ 1].componentId, 0);
  EXPECT_EQ(atomPositions[ 2].componentId, 0);
  EXPECT_EQ(atomPositions[ 3].componentId, 1);
  EXPECT_EQ(atomPositions[ 4].componentId, 1);
  EXPECT_EQ(atomPositions[ 5].componentId, 1);
  EXPECT_EQ(atomPositions[ 6].componentId, 1);
  EXPECT_EQ(atomPositions[ 7].componentId, 1);
  EXPECT_EQ(atomPositions[ 8].componentId, 1);
  EXPECT_EQ(atomPositions[ 9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);
  EXPECT_EQ(atomPositions[13].componentId, 1);
  EXPECT_EQ(atomPositions[14].componentId, 1);
  EXPECT_EQ(atomPositions[15].componentId, 1);
  EXPECT_EQ(atomPositions[16].componentId, 1);
  EXPECT_EQ(atomPositions[17].componentId, 1);

  //insert new Methane
  system.insertMolecule(0, Molecule(), {Atom(double3(0.0, 0.0,  0.0),    0.0, 1.0, 0, 2, 0, 0) });

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 4);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 4);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 0);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 5);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 5);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 0);

  // refresh span
  atomPositions = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomPositions[ 0].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 1].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 2].moleculeId, 2);
  EXPECT_EQ(atomPositions[ 3].moleculeId, 3);
  EXPECT_EQ(atomPositions[ 4].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 5].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 6].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 7].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 8].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 9].moleculeId, 1);
  EXPECT_EQ(atomPositions[10].moleculeId, 2);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 2);
  EXPECT_EQ(atomPositions[13].moleculeId, 3);
  EXPECT_EQ(atomPositions[14].moleculeId, 3);
  EXPECT_EQ(atomPositions[15].moleculeId, 3);
  EXPECT_EQ(atomPositions[16].moleculeId, 4);
  EXPECT_EQ(atomPositions[17].moleculeId, 4);
  EXPECT_EQ(atomPositions[18].moleculeId, 4);

  EXPECT_EQ(atomPositions[ 0].componentId, 0);
  EXPECT_EQ(atomPositions[ 1].componentId, 0);
  EXPECT_EQ(atomPositions[ 2].componentId, 0);
  EXPECT_EQ(atomPositions[ 3].componentId, 0);
  EXPECT_EQ(atomPositions[ 4].componentId, 1);
  EXPECT_EQ(atomPositions[ 5].componentId, 1);
  EXPECT_EQ(atomPositions[ 6].componentId, 1);
  EXPECT_EQ(atomPositions[ 7].componentId, 1);
  EXPECT_EQ(atomPositions[ 8].componentId, 1);
  EXPECT_EQ(atomPositions[ 9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);
  EXPECT_EQ(atomPositions[13].componentId, 1);
  EXPECT_EQ(atomPositions[14].componentId, 1);
  EXPECT_EQ(atomPositions[15].componentId, 1);
  EXPECT_EQ(atomPositions[16].componentId, 1);
  EXPECT_EQ(atomPositions[17].componentId, 1);
  EXPECT_EQ(atomPositions[18].componentId, 1);

  //delete CO2 1
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

  EXPECT_EQ(atomPositions[ 0].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 1].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 2].moleculeId, 2);
  EXPECT_EQ(atomPositions[ 3].moleculeId, 3);
  EXPECT_EQ(atomPositions[ 4].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 5].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 6].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 7].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 8].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 9].moleculeId, 1);
  EXPECT_EQ(atomPositions[10].moleculeId, 2);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 2);
  EXPECT_EQ(atomPositions[13].moleculeId, 3);
  EXPECT_EQ(atomPositions[14].moleculeId, 3);
  EXPECT_EQ(atomPositions[15].moleculeId, 3);

  EXPECT_EQ(atomPositions[ 0].componentId, 0);
  EXPECT_EQ(atomPositions[ 1].componentId, 0);
  EXPECT_EQ(atomPositions[ 2].componentId, 0);
  EXPECT_EQ(atomPositions[ 3].componentId, 0);
  EXPECT_EQ(atomPositions[ 4].componentId, 1);
  EXPECT_EQ(atomPositions[ 5].componentId, 1);
  EXPECT_EQ(atomPositions[ 6].componentId, 1);
  EXPECT_EQ(atomPositions[ 7].componentId, 1);
  EXPECT_EQ(atomPositions[ 8].componentId, 1);
  EXPECT_EQ(atomPositions[ 9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);
  EXPECT_EQ(atomPositions[13].componentId, 1);
  EXPECT_EQ(atomPositions[14].componentId, 1);
  EXPECT_EQ(atomPositions[15].componentId, 1);


  //delete CO2 1
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

  EXPECT_EQ(atomPositions[ 0].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 1].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 2].moleculeId, 2);
  EXPECT_EQ(atomPositions[ 3].moleculeId, 3);
  EXPECT_EQ(atomPositions[ 4].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 5].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 6].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 7].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 8].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 9].moleculeId, 1);
  EXPECT_EQ(atomPositions[10].moleculeId, 2);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 2);

  EXPECT_EQ(atomPositions[ 0].componentId, 0);
  EXPECT_EQ(atomPositions[ 1].componentId, 0);
  EXPECT_EQ(atomPositions[ 2].componentId, 0);
  EXPECT_EQ(atomPositions[ 3].componentId, 0);
  EXPECT_EQ(atomPositions[ 4].componentId, 1);
  EXPECT_EQ(atomPositions[ 5].componentId, 1);
  EXPECT_EQ(atomPositions[ 6].componentId, 1);
  EXPECT_EQ(atomPositions[ 7].componentId, 1);
  EXPECT_EQ(atomPositions[ 8].componentId, 1);
  EXPECT_EQ(atomPositions[ 9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);

  //insert new Methane
  system.insertMolecule(0, Molecule(), {Atom(double3(0.0, 0.0,  0.0),    0.0, 1.0, 0, 2, 0, 0) });

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 0);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 3);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 3);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 0);

  // refresh span
  atomPositions = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomPositions[ 0].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 1].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 2].moleculeId, 2);
  EXPECT_EQ(atomPositions[ 3].moleculeId, 3);
  EXPECT_EQ(atomPositions[ 4].moleculeId, 4);
  EXPECT_EQ(atomPositions[ 5].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 6].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 7].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 8].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 9].moleculeId, 1);
  EXPECT_EQ(atomPositions[10].moleculeId, 1);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 2);
  EXPECT_EQ(atomPositions[13].moleculeId, 2);

  EXPECT_EQ(atomPositions[ 0].componentId, 0);
  EXPECT_EQ(atomPositions[ 1].componentId, 0);
  EXPECT_EQ(atomPositions[ 2].componentId, 0);
  EXPECT_EQ(atomPositions[ 3].componentId, 0);
  EXPECT_EQ(atomPositions[ 4].componentId, 0);
  EXPECT_EQ(atomPositions[ 5].componentId, 1);
  EXPECT_EQ(atomPositions[ 6].componentId, 1);
  EXPECT_EQ(atomPositions[ 7].componentId, 1);
  EXPECT_EQ(atomPositions[ 8].componentId, 1);
  EXPECT_EQ(atomPositions[ 9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);
  EXPECT_EQ(atomPositions[13].componentId, 1);


  //delete Methane 2
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

  EXPECT_EQ(atomPositions[ 0].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 1].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 2].moleculeId, 2);
  EXPECT_EQ(atomPositions[ 3].moleculeId, 3);
  EXPECT_EQ(atomPositions[ 4].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 5].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 6].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 7].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 8].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 9].moleculeId, 1);
  EXPECT_EQ(atomPositions[10].moleculeId, 2);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 2);

  EXPECT_EQ(atomPositions[ 0].componentId, 0);
  EXPECT_EQ(atomPositions[ 1].componentId, 0);
  EXPECT_EQ(atomPositions[ 2].componentId, 0);
  EXPECT_EQ(atomPositions[ 3].componentId, 0);
  EXPECT_EQ(atomPositions[ 4].componentId, 1);
  EXPECT_EQ(atomPositions[ 5].componentId, 1);
  EXPECT_EQ(atomPositions[ 6].componentId, 1);
  EXPECT_EQ(atomPositions[ 7].componentId, 1);
  EXPECT_EQ(atomPositions[ 8].componentId, 1);
  EXPECT_EQ(atomPositions[ 9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);

  //insert new CO2
  system.insertMolecule(1, Molecule(),
                           {Atom(double3(0.0, 0.0,  1.149), -0.3256, 1.0, 0, 4, 1, 0),
                            Atom(double3(0.0, 0.0,  0.0  ),  0.6512, 1.0, 0, 3, 1, 0),
                            Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)});

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 4);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 4);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 0);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 4);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 4);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 0);

  // refresh span
  atomPositions = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomPositions[ 0].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 1].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 2].moleculeId, 2);
  EXPECT_EQ(atomPositions[ 3].moleculeId, 3);
  EXPECT_EQ(atomPositions[ 4].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 5].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 6].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 7].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 8].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 9].moleculeId, 1);
  EXPECT_EQ(atomPositions[10].moleculeId, 2);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 2);
  EXPECT_EQ(atomPositions[13].moleculeId, 3);
  EXPECT_EQ(atomPositions[14].moleculeId, 3);
  EXPECT_EQ(atomPositions[15].moleculeId, 3);

  EXPECT_EQ(atomPositions[ 0].componentId, 0);
  EXPECT_EQ(atomPositions[ 1].componentId, 0);
  EXPECT_EQ(atomPositions[ 2].componentId, 0);
  EXPECT_EQ(atomPositions[ 3].componentId, 0);
  EXPECT_EQ(atomPositions[ 4].componentId, 1);
  EXPECT_EQ(atomPositions[ 5].componentId, 1);
  EXPECT_EQ(atomPositions[ 6].componentId, 1);
  EXPECT_EQ(atomPositions[ 7].componentId, 1);
  EXPECT_EQ(atomPositions[ 8].componentId, 1);
  EXPECT_EQ(atomPositions[ 9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);
  EXPECT_EQ(atomPositions[13].componentId, 1);
  EXPECT_EQ(atomPositions[14].componentId, 1);
  EXPECT_EQ(atomPositions[15].componentId, 1);

  //insert new Methane
  system.insertMolecule(0, Molecule(), {Atom(double3(0.0, 0.0,  0.0),    0.0, 1.0, 0, 2, 0, 0) });

  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 5);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[0], 0);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[1], 4);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[1], 4);
  EXPECT_EQ(system.numberOfFractionalMoleculesPerComponent[1], 0);

  // refresh span
  atomPositions = system.spanOfMoleculeAtoms();

  EXPECT_EQ(atomPositions[ 0].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 1].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 2].moleculeId, 2);
  EXPECT_EQ(atomPositions[ 3].moleculeId, 3);
  EXPECT_EQ(atomPositions[ 4].moleculeId, 4);
  EXPECT_EQ(atomPositions[ 5].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 6].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 7].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 8].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 9].moleculeId, 1);
  EXPECT_EQ(atomPositions[10].moleculeId, 1);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 2);
  EXPECT_EQ(atomPositions[13].moleculeId, 2);
  EXPECT_EQ(atomPositions[14].moleculeId, 3);
  EXPECT_EQ(atomPositions[15].moleculeId, 3);
  EXPECT_EQ(atomPositions[16].moleculeId, 3);

  EXPECT_EQ(atomPositions[ 0].componentId, 0);
  EXPECT_EQ(atomPositions[ 1].componentId, 0);
  EXPECT_EQ(atomPositions[ 2].componentId, 0);
  EXPECT_EQ(atomPositions[ 3].componentId, 0);
  EXPECT_EQ(atomPositions[ 4].componentId, 0);
  EXPECT_EQ(atomPositions[ 5].componentId, 1);
  EXPECT_EQ(atomPositions[ 6].componentId, 1);
  EXPECT_EQ(atomPositions[ 7].componentId, 1);
  EXPECT_EQ(atomPositions[ 8].componentId, 1);
  EXPECT_EQ(atomPositions[ 9].componentId, 1);
  EXPECT_EQ(atomPositions[10].componentId, 1);
  EXPECT_EQ(atomPositions[11].componentId, 1);
  EXPECT_EQ(atomPositions[12].componentId, 1);
  EXPECT_EQ(atomPositions[13].componentId, 1);
  EXPECT_EQ(atomPositions[14].componentId, 1);
  EXPECT_EQ(atomPositions[15].componentId, 1);
  EXPECT_EQ(atomPositions[16].componentId, 1);

  //delete CO2 0
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

  EXPECT_EQ(atomPositions[ 0].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 1].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 2].moleculeId, 2);
  EXPECT_EQ(atomPositions[ 3].moleculeId, 3);
  EXPECT_EQ(atomPositions[ 4].moleculeId, 4);
  EXPECT_EQ(atomPositions[ 5].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 6].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 7].moleculeId, 0);
  EXPECT_EQ(atomPositions[ 8].moleculeId, 1);
  EXPECT_EQ(atomPositions[ 9].moleculeId, 1);
  EXPECT_EQ(atomPositions[10].moleculeId, 1);
  EXPECT_EQ(atomPositions[11].moleculeId, 2);
  EXPECT_EQ(atomPositions[12].moleculeId, 2);
  EXPECT_EQ(atomPositions[13].moleculeId, 2);

  EXPECT_EQ(atomPositions[ 0].componentId, 0);
  EXPECT_EQ(atomPositions[ 1].componentId, 0);
  EXPECT_EQ(atomPositions[ 2].componentId, 0);
  EXPECT_EQ(atomPositions[ 3].componentId, 0);
  EXPECT_EQ(atomPositions[ 4].componentId, 0);
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
