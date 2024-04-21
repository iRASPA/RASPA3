#include <gtest/gtest.h>

#include <cstddef>
#include <algorithm>
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
import force_factor;
import running_energy;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;

TEST(Gradients, Test_CH4_in_Box_25x25x25)
{
  ForceField forceField = ForceField(
    { PseudoAtom("CH4", 16.04246, 0.0, 6, false) },
    { VDWParameters(158.5, 3.72) },
    ForceField::MixingRule::Lorentz_Berthelot,
    12.0,
    false,
    false);
  Component c = Component(0, forceField, "methane", 190.564, 45599200, 0.01142,
    { // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
      Atom(double3(0.0, 0.0,  0.0),    0.0, 1.0, 0, 0, 0, 0) 
    }, 5, 21);

  System system = System(0, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, forceField, { }, { c }, { 50 }, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  std::vector<Atom> atomPositions = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end());

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  [[maybe_unused]] ForceFactor factorFrameworkMolecular = 
    Interactions::computeFrameworkMoleculeGradient(system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms());
  [[maybe_unused]] ForceFactor factorInterMolecular = 
    Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox, system.spanOfMoleculeAtoms());

  double delta = 1e-5;
  double tolerance = 1e-6;
  double3 gradient;
  for (size_t i = 0; i < atomPositions.size(); ++i)
  {
    RunningEnergy x1, x2, y1, y2, z1, z2;

    // finite difference x
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x + 0.5 * delta;
    Interactions::computeFrameworkMoleculeEnergy(system.forceField,system.simulationBox, frameworkAtoms, atomPositions, x2);
    Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions, x2);

    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x - 0.5 * delta;
    Interactions::computeFrameworkMoleculeEnergy(system.forceField,system.simulationBox, frameworkAtoms, atomPositions, x1);
    Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions, x1);
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x;

    // finite difference y
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y + 0.5 * delta;
    Interactions::computeFrameworkMoleculeEnergy(system.forceField,system.simulationBox, frameworkAtoms, atomPositions, y2);
    Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions, y2);

    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y - 0.5 * delta;
    Interactions::computeFrameworkMoleculeEnergy(system.forceField,system.simulationBox, frameworkAtoms, atomPositions, y1);
    Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions, y1);
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y;

    // finite difference z
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z + 0.5 * delta;
    Interactions::computeFrameworkMoleculeEnergy(system.forceField,system.simulationBox, frameworkAtoms, atomPositions, z2);
    Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions, z2);

    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z - 0.5 * delta;
    Interactions::computeFrameworkMoleculeEnergy(system.forceField,system.simulationBox, frameworkAtoms, atomPositions, z1);
    Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions, z1);
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z;

    gradient.x = (x2.moleculeMoleculeVDW + x2.moleculeMoleculeCharge + x2.frameworkMoleculeVDW  + x2.frameworkMoleculeCharge 
                - (x1.moleculeMoleculeVDW + x1.moleculeMoleculeCharge + x1.frameworkMoleculeVDW + x1.frameworkMoleculeCharge)) / delta;
    gradient.y = (y2.moleculeMoleculeVDW + y2.moleculeMoleculeCharge + y2.frameworkMoleculeVDW + y2.frameworkMoleculeCharge 
                - (y1.moleculeMoleculeVDW + y1.moleculeMoleculeCharge + y1.frameworkMoleculeVDW + y1.frameworkMoleculeCharge)) / delta;
    gradient.z = (z2.moleculeMoleculeVDW + z2.moleculeMoleculeCharge + z2.frameworkMoleculeVDW + z2.frameworkMoleculeCharge 
                - (z1.moleculeMoleculeVDW + z1.moleculeMoleculeCharge + z1.frameworkMoleculeVDW + z1.frameworkMoleculeCharge)) / delta;

    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.x, gradient.x, tolerance);
    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.y, gradient.y, tolerance);
    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.z, gradient.z, tolerance);
  }
}


TEST(Gradients, Test_CO2_in_MFI_2x2x2)
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

  System system = System(0, std::nullopt, 300.0, 1e4, forceField, { f }, { c }, { 0, 10 }, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  std::vector<Atom> atomPositions = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end()); 

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  [[maybe_unused]] ForceFactor factor = 
    Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox, system.spanOfMoleculeAtoms());
  [[maybe_unused]] ForceFactor factor2 = 
    Interactions::computeFrameworkMoleculeGradient(system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms());

  double delta = 1e-5;
  double tolerance = 1e-4;
  double3 gradient;
  for (size_t i = 0; i < atomPositions.size(); ++i)
  {
    RunningEnergy x1, x2, y1, y2, z1, z2;

    // finite difference x
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x + 0.5 * delta;
    Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions, x2);
    Interactions::computeFrameworkMoleculeEnergy(system.forceField,system.simulationBox, frameworkAtoms, atomPositions, x2);

    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x - 0.5 * delta;
    Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions, x1);
    Interactions::computeFrameworkMoleculeEnergy(system.forceField,system.simulationBox, frameworkAtoms, atomPositions, x1);
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x;

    // finite difference y
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y + 0.5 * delta;
    Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions, y2);
    Interactions::computeFrameworkMoleculeEnergy(system.forceField,system.simulationBox, frameworkAtoms, atomPositions, y2);

    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y - 0.5 * delta;
    Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions, y1);
    Interactions::computeFrameworkMoleculeEnergy(system.forceField,system.simulationBox, frameworkAtoms, atomPositions, y1);
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y;

    // finite difference z
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z + 0.5 * delta;
    Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions, z2);
    Interactions::computeFrameworkMoleculeEnergy(system.forceField,system.simulationBox, frameworkAtoms, atomPositions, z2);

    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z - 0.5 * delta;
    Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions, z1);
    Interactions::computeFrameworkMoleculeEnergy(system.forceField,system.simulationBox, frameworkAtoms, atomPositions, z1);
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z;

    gradient.x = (x2.moleculeMoleculeVDW + x2.moleculeMoleculeCharge + x2.frameworkMoleculeVDW + x2.frameworkMoleculeCharge
      - (x1.moleculeMoleculeVDW + x1.moleculeMoleculeCharge + x1.frameworkMoleculeVDW + x1.frameworkMoleculeCharge)) / delta;
    gradient.y = (y2.moleculeMoleculeVDW + y2.moleculeMoleculeCharge + y2.frameworkMoleculeVDW + y2.frameworkMoleculeCharge
      - (y1.moleculeMoleculeVDW + y1.moleculeMoleculeCharge + y1.frameworkMoleculeVDW + y1.frameworkMoleculeCharge)) / delta;
    gradient.z = (z2.moleculeMoleculeVDW + z2.moleculeMoleculeCharge + z2.frameworkMoleculeVDW + z2.frameworkMoleculeCharge
      - (z1.moleculeMoleculeVDW + z1.moleculeMoleculeCharge + z1.frameworkMoleculeVDW + z1.frameworkMoleculeCharge)) / delta;

    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.x, gradient.x, tolerance);
    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.y, gradient.y, tolerance);
    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.z, gradient.z, tolerance);
  }
}

TEST(Gradients, Test_20_Na_Cl_in_Box_25x25x25)
{
  ForceField forceField = ForceField(
    { PseudoAtom("Si",   28.0855,   2.05,  14, false),
      PseudoAtom("O",    15.999,   -1.025,  8, false),
      PseudoAtom("CH4",  16.04246,  0.0,    6, false),
      PseudoAtom("Na+",  12.0,      0.0, 6, false),
      PseudoAtom("Cl-",  15.9994,   0.0, 8, false),
    },
    { VDWParameters(22.0, 2.30),
      VDWParameters(53.0, 3.3),
      VDWParameters(158.5, 3.72),
      VDWParameters(15.0966, 2.65755),
      VDWParameters(142.562, 3.51932)
    },
    ForceField::MixingRule::Lorentz_Berthelot,
    12.0,
    true,
    false);
  Component na = Component(0, forceField, "Na", 304.1282, 7377300.0, 0.22394,
    { // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
      Atom(double3(0.0, 0.0, 0.0), 1.0, 1.0, 0, 3, 0, 0),
    }, 5, 21);
  Component cl = Component(1, forceField, "Cl", 304.1282, 7377300.0, 0.22394,
    { // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
      Atom(double3(0.0, 0.0, 0.0), -1.0, 1.0, 0, 4, 1, 0),
    }, 5, 21);

  System system = System(0, SimulationBox(25.0, 25.0, 25.0), 300.0, 1e4, forceField, {}, { na, cl }, { 20, 20 }, 5);

  //std::fill(system.forceField.data.begin(), system.forceField.data.end(), VDWParameters(0.0, 1.0));

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::vector<Atom> atomPositions = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end());

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  for (size_t i = 0; i < 20; ++i)
  {
    atomPositions[i].charge = 1.0;
    system.atomPositions[i].charge = 1.0;
  }
  for (size_t i = 0; i < 20; ++i)
  {
    atomPositions[i + 20].charge = -1.0;
    system.atomPositions[i + 20].charge = -1.0;
  }

  [[maybe_unused]] ForceFactor factor = 
    Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox, system.spanOfMoleculeAtoms());

  double delta = 1e-5;
  double tolerance = 1e-4;
  double3 gradient;
  for (size_t i = 0; i < atomPositions.size(); ++i)
  {
    RunningEnergy x1, x2, y1, y2, z1, z2;

    // finite difference x
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x + 0.5 * delta;
    Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions, x2);

    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x - 0.5 * delta;
    Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions, x1);
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x;

    // finite difference y
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y + 0.5 * delta;
    Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions, y2);

    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y - 0.5 * delta;
    Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions, y1);
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y;

    // finite difference z
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z + 0.5 * delta;
    Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions, z2);

    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z - 0.5 * delta;
    Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, atomPositions, z1);
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z;

    gradient.x = (x2.moleculeMoleculeVDW + x2.moleculeMoleculeCharge - (x1.moleculeMoleculeVDW + x1.moleculeMoleculeCharge)) / delta;
    gradient.y = (y2.moleculeMoleculeVDW + y2.moleculeMoleculeCharge - (y1.moleculeMoleculeVDW + y1.moleculeMoleculeCharge)) / delta;
    gradient.z = (z2.moleculeMoleculeVDW + z2.moleculeMoleculeCharge - (z1.moleculeMoleculeVDW + z1.moleculeMoleculeCharge)) / delta;

    EXPECT_NEAR(system.atomPositions[i].gradient.x, gradient.x, tolerance);
    EXPECT_NEAR(system.atomPositions[i].gradient.y, gradient.y, tolerance);
    EXPECT_NEAR(system.atomPositions[i].gradient.z, gradient.z, tolerance);
  }
}
