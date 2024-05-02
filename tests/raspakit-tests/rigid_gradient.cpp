#include <gtest/gtest.h>

#include <algorithm>
#include <vector>
#include <cstddef>
#include <vector>
#include <span>
#include <complex>

import int3;
import double3;
import double3x3;
import simd_quatd;

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
import energy_status;

TEST(RigidGradient, Test_2_CO2_in_ITQ_29_2x2x2)
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
    11.8,
    true,
    false);
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);
  Framework f = Framework(0, forceField, "ITQ-29", SimulationBox(11.8671, 11.8671, 11.8671),
    517,
    { // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
      Atom(double3(0.3683, 0.1847, 0),       2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.5,    0.2179, 0),      -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.2939, 0.2939, 0),      -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.3429, 0.1098, 0.1098), -1.025, 1.0, 0, 1, 0, 0)
    },
    int3(2, 2, 2));
  Component c = Component(1, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
    { // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
      Atom(double3(0.0, 0.0,  1.149), -0.3256, 1.0, 0, 4, 1, 0),
      Atom(double3(0.0, 0.0,  0.0),    0.6512, 1.0, 0, 3, 1, 0),
      Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)
    }, 5, 21);

  System system = System(0, std::nullopt, 300.0, 1e4, forceField, { f }, { c }, { 2 }, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  atomPositions[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  atomPositions[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  atomPositions[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  atomPositions[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  atomPositions[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  atomPositions[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);

  system.moleculePositions[0].orientation = simd_quatd(0.0, 0.0, 0.0, 1.0);
  system.moleculePositions[1].orientation = simd_quatd(0.0, 0.0, 0.0, 1.0);

  for(Atom& atom: atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  system.computeTotalGradients();
  system.computeCenterOfMassAndQuaternionGradients();

  EXPECT_NEAR(atomPositions[0].gradient.x,     0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[0].gradient.y,   103.939706389550, 1e-4);
  EXPECT_NEAR(atomPositions[0].gradient.z,   -18.706660418853, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.x,     0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.y,  -574.210506496196, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.z,     0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.x,     0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.y,   103.939706389549, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.z,    18.706660418865, 1e-4);

  EXPECT_NEAR(atomPositions[3].gradient.x,     0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[3].gradient.y,  -103.939706389550, 1e-4);
  EXPECT_NEAR(atomPositions[3].gradient.z,   -18.706660418853, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.x,     0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.y,   574.210506496196, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.z,     0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.x,     0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.y,  -103.939706389549, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.z,    18.706660418865, 1e-4);

  EXPECT_NEAR(system.moleculePositions[0].gradient.x,     0.000000000000, 1e-4);
  EXPECT_NEAR(system.moleculePositions[0].gradient.y,  -366.33109372,     1e-4);
  EXPECT_NEAR(system.moleculePositions[0].gradient.z,     0.000000000000, 1e-4);

  EXPECT_NEAR(system.moleculePositions[1].gradient.x,     0.000000000000, 1e-4);
  EXPECT_NEAR(system.moleculePositions[1].gradient.y,   366.33109372,     1e-4);
  EXPECT_NEAR(system.moleculePositions[1].gradient.z,     0.000000000000, 1e-4);

  EXPECT_NEAR(system.moleculePositions[0].orientationGradient.ix,     0.000000000000, 1e-6);
  EXPECT_NEAR(system.moleculePositions[0].orientationGradient.iy,     0.000000000000, 1e-6);
  EXPECT_NEAR(system.moleculePositions[0].orientationGradient.iz,     0.000000000000, 1e-6);
  EXPECT_NEAR(system.moleculePositions[0].orientationGradient.r,      0.000000000000, 1e-6);

  EXPECT_NEAR(system.moleculePositions[1].orientationGradient.ix,     0.000000000000, 1e-6);
  EXPECT_NEAR(system.moleculePositions[1].orientationGradient.iy,     0.000000000000, 1e-6);
  EXPECT_NEAR(system.moleculePositions[1].orientationGradient.iz,     0.000000000000, 1e-6);
  EXPECT_NEAR(system.moleculePositions[1].orientationGradient.r,      0.000000000000, 1e-6);
}

TEST(RigidGradient, Test_2_CO2_in_ITQ_29_2x2x2_no_symmetry)
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
    11.8,
    true,
    false);
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);
  Framework f = Framework(0, forceField, "ITQ-29", SimulationBox(11.8671, 11.8671, 11.8671),
    517,
    { // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
      Atom(double3(0.3683, 0.1847, 0),       2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.5,    0.2179, 0),      -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.2939, 0.2939, 0),      -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.3429, 0.1098, 0.1098), -1.025, 1.0, 0, 1, 0, 0)
    },
    int3(2, 2, 2));
  Component c = Component(1, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
    { // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
      Atom(double3(0.0, 0.0,  1.149), -0.3256, 1.0, 0, 4, 1, 0),
      Atom(double3(0.0, 0.0,  0.0),    0.6512, 1.0, 0, 3, 1, 0),
      Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)
    }, 5, 21);

  System system = System(0, std::nullopt, 300.0, 1e4, forceField, { f }, { c }, { 2 }, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  atomPositions[0].position = double3(7.200346017629, 5.547320950450, 2.554032891860);
  atomPositions[1].position = double3(6.206010115191, 5.030510168430, 2.807811789148);
  atomPositions[2].position = double3(5.211674212752, 4.513699386409, 3.061590686436);
  atomPositions[3].position = double3(8.486922830523, 3.714749040941, 7.214116618614);
  atomPositions[4].position = double3(7.905757540326, 3.815476396960, 6.228062915605);
  atomPositions[5].position = double3(7.324592250129, 3.916203752979, 5.242009212595);

  system.moleculePositions[0].orientation = simd_quatd(0.44301474, 0.64356245, 0.55112425, 0.29296310);
  system.moleculePositions[1].orientation = simd_quatd(0.04909940, 0.26171955, 0.01333290, 0.96380203);

  for(Atom& atom: atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  system.computeTotalGradients();
  system.computeCenterOfMassAndQuaternionGradients();

  EXPECT_NEAR(atomPositions[0].gradient.x,  -401.918996903652, 1e-4);
  EXPECT_NEAR(atomPositions[0].gradient.y,  -250.086986676442, 1e-4);
  EXPECT_NEAR(atomPositions[0].gradient.z,    76.393683260493, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.x,  -112.917384315368, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.y,   135.881091839024, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.z,   191.752321499793, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.x,   649.466150035879, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.y,   183.761020135682, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.z,   542.179200465678, 1e-4);

  EXPECT_NEAR(atomPositions[3].gradient.x,   -30.785378731879, 1e-4);
  EXPECT_NEAR(atomPositions[3].gradient.y,    80.096820989863, 1e-4);
  EXPECT_NEAR(atomPositions[3].gradient.z,   113.734455039940, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.x,   -40.135560076449, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.y,    19.413856681000, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.z,  -302.257466181579, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.x,  -310.276701548065, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.y,   181.549370345199, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.z,   -33.308103959230, 1e-4);

  EXPECT_NEAR(system.moleculePositions[0].gradient.x,   134.62976882, 1e-4);
  EXPECT_NEAR(system.moleculePositions[0].gradient.y,    69.55512530, 1e-4);
  EXPECT_NEAR(system.moleculePositions[0].gradient.z,   810.32520523, 1e-4);

  EXPECT_NEAR(system.moleculePositions[1].gradient.x,  -381.19764036, 1e-4);
  EXPECT_NEAR(system.moleculePositions[1].gradient.y,   281.06004802, 1e-4);
  EXPECT_NEAR(system.moleculePositions[1].gradient.z,  -221.83111510, 1e-4);

  EXPECT_NEAR(system.moleculePositions[0].orientationGradient.ix,   454.91960800, 1e-4);
  EXPECT_NEAR(system.moleculePositions[0].orientationGradient.iy,   913.61679756, 1e-4);
  EXPECT_NEAR(system.moleculePositions[0].orientationGradient.iz, -1032.71684255, 1e-4);
  EXPECT_NEAR(system.moleculePositions[0].orientationGradient.r,   -752.14337829, 1e-4);

  EXPECT_NEAR(system.moleculePositions[1].orientationGradient.ix,   185.47928131, 1e-4);
  EXPECT_NEAR(system.moleculePositions[1].orientationGradient.iy,   361.21170189, 1e-4);
  EXPECT_NEAR(system.moleculePositions[1].orientationGradient.iz,   -33.44656884, 1e-4);
  EXPECT_NEAR(system.moleculePositions[1].orientationGradient.r,   -107.07297003, 1e-4);
}

