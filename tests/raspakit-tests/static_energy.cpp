#include <gtest/gtest.h>

#include <algorithm>
#include <vector>
#include <cstddef>
#include <vector>
#include <span>
#include <complex>

import int3;
import double3;

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

TEST(StaticEnergy, Test_2_CO2_in_ITQ_29_1x1x1)
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
  Framework f = Framework(0, forceField, "ITQ-29", SimulationBox(11.8671, 11.8671, 11.8671),
    517,
    { // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
      Atom(double3(0.3683, 0.1847, 0),       2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.5,    0.2179, 0),      -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.2939, 0.2939, 0),      -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.3429, 0.1098, 0.1098), -1.025, 1.0, 0, 1, 0, 0)
    },
    int3(1, 1, 1));
  Component c = Component(1, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
    { // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
      Atom(double3(0.0, 0.0,  1.149), -0.3256, 1.0, 0, 4, 1, 0),
      Atom(double3(0.0, 0.0,  0.0  ),  0.6512, 1.0, 0, 3, 1, 0),
      Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)
    }, 5, 21);

  System system = System(0, std::nullopt, 300.0, 1e4, forceField, { f }, { c }, { 2 }, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  atomPositions[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  atomPositions[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  atomPositions[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  atomPositions[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  atomPositions[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  atomPositions[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);

  RunningEnergy energy;
  Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions, energy);
  Interactions::computeFrameworkMoleculeEnergy(system.forceField,system.simulationBox, frameworkAtoms, atomPositions, energy);

  EXPECT_NEAR(energy.frameworkMoleculeVDW * Units::EnergyToKelvin, -1545.62921755, 1e-6);
  EXPECT_NEAR(energy.frameworkMoleculeCharge * Units::EnergyToKelvin, -592.13188606, 1e-6);
  EXPECT_NEAR(energy.moleculeMoleculeVDW * Units::EnergyToKelvin, -242.42459776, 1e-6);
  EXPECT_NEAR(energy.moleculeMoleculeCharge * Units::EnergyToKelvin, 154.11883595, 1e-6);
}

TEST(StaticEnergy, Test_2_CO2_in_MFI_2x2x2_shifted)
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
  Component c = Component(1, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
    {  // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
       Atom(double3(0.0, 0.0,  1.149), -0.3256, 1.0, 0, 4, 1, 0),
       Atom(double3(0.0, 0.0,  0.0),    0.6512, 1.0, 0, 3, 1, 0),
       Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)
    }, 5, 21);

  System system = System(0, std::nullopt, 300.0, 1e4, forceField, { f } , { c }, { 2 }, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  atomPositions[0].position = double3(10.011, 4.97475 + 2.0,  1.149);
  atomPositions[1].position = double3(10.011, 4.97475 + 2.0,  0.0);
  atomPositions[2].position = double3(10.011, 4.97475 + 2.0, -1.149);
  atomPositions[3].position = double3(10.011, 4.97475 - 2.0, 1.149);
  atomPositions[4].position = double3(10.011, 4.97475 - 2.0, 0.0);
  atomPositions[5].position = double3(10.011, 4.97475 - 2.0, -1.149);
  
  RunningEnergy energy;
  Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions, energy);            
  Interactions::computeFrameworkMoleculeEnergy(system.forceField,system.simulationBox, frameworkAtoms, atomPositions, energy);

  EXPECT_NEAR(energy.frameworkMoleculeVDW * Units::EnergyToKelvin, -2525.36580663, 1e-6);
  EXPECT_NEAR(energy.frameworkMoleculeCharge * Units::EnergyToKelvin, 2167.45591472, 1e-6);
  EXPECT_NEAR(energy.moleculeMoleculeVDW * Units::EnergyToKelvin, -242.42459776, 1e-6);
  EXPECT_NEAR(energy.moleculeMoleculeCharge * Units::EnergyToKelvin, 154.11883595, 1e-6);
}

TEST(StaticEnergy, Test_2_CO2_in_MFI_2x2x2_truncated)
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
    false,
    true);
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
  Component c = Component(1, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
    {  // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
       Atom(double3(0.0, 0.0,  1.149), -0.3256, 1.0, 0, 4, 1, 0),
       Atom(double3(0.0, 0.0,  0.0),    0.6512, 1.0, 0, 3, 1, 0),
       Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)
    }, 5, 21);

  System system = System(0, std::nullopt, 300.0, 1e4, forceField, { f } ,{ c }, { 2 }, 5);

  std::span<Atom> moleculeAtomPositions = system.spanOfMoleculeAtoms();                                                  
  moleculeAtomPositions[0].position = double3(10.011, 4.97475 + 2.0, 1.149);
  moleculeAtomPositions[1].position = double3(10.011, 4.97475 + 2.0, 0.0);
  moleculeAtomPositions[2].position = double3(10.011, 4.97475 + 2.0, -1.149);
  moleculeAtomPositions[3].position = double3(10.011, 4.97475 - 2.0, 1.149);
  moleculeAtomPositions[4].position = double3(10.011, 4.97475 - 2.0, 0.0);
  moleculeAtomPositions[5].position = double3(10.011, 4.97475 - 2.0, -1.149);

  system.forceField.EwaldAlpha = 0.25;
  system.forceField.numberOfWaveVectors = int3(8, 8, 8);
  //double3 perpendicularWidths = system.simulationBox.perpendicularWidths();                                           
  //system.forceField.initializeEwaldParameters(perpendicularWidths);   

  RunningEnergy energy = system.computeTotalEnergies();

  EXPECT_NEAR(energy.frameworkMoleculeVDW * Units::EnergyToKelvin, -2657.36121975, 1e-6);
  EXPECT_NEAR(energy.frameworkMoleculeCharge * Units::EnergyToKelvin, 1971.00612979, 1e-6);
  EXPECT_NEAR(energy.moleculeMoleculeVDW * Units::EnergyToKelvin, -242.94298709, 1e-6);
  EXPECT_NEAR(energy.moleculeMoleculeCharge * Units::EnergyToKelvin, 162.41877650, 1e-6);
  //EXPECT_NEAR(energy.tail * Units::EnergyToKelvin, -127.81601515, 1e-6);
  EXPECT_NEAR(energy.tail * Units::EnergyToKelvin,-127.72803736223419, 1e-6);
  EXPECT_NEAR(energy.ewald * Units::EnergyToKelvin, -1197.23909965, 1e-6);
}
