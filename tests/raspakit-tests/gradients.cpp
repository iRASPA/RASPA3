#include <gtest/gtest.h>

#include <cstddef>
#include <algorithm>
#include <vector>
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
import force_factor;
import running_energy;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;
import energy_status;

TEST(Gradients, Test_2_CO2_in_ITQ_29_2x2x2_inter)
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
      Atom(double3(0.0, 0.0,  1.149), 0.0, 1.0, 0, 4, 1, 0),
      Atom(double3(0.0, 0.0,  0.0  ), 0.0, 1.0, 0, 3, 1, 0),
      Atom(double3(0.0, 0.0, -1.149), 0.0, 1.0, 0, 4, 1, 0)
    }, 5, 21);

  System system = System(0, std::nullopt, 300.0, 1e4, forceField, { f }, { c }, { 2 }, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  atomPositions[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  atomPositions[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  atomPositions[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  atomPositions[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  atomPositions[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  atomPositions[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);

  std::pair<ForceFactor, ForceFactor> factorInterMolecular =
    Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox, system.spanOfMoleculeAtoms());

  EXPECT_NEAR(factorInterMolecular.first.energy  * Units::EnergyToKelvin,  -242.36960932, 1e-6);
  EXPECT_NEAR(factorInterMolecular.second.energy * Units::EnergyToKelvin,     0.00000000, 1e-6);

  EXPECT_NEAR(atomPositions[0].gradient.x,   0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[0].gradient.y,  90.951955774481, 1e-6);
  EXPECT_NEAR(atomPositions[0].gradient.z,  17.938271420558, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.x,   0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.y,  52.398564956433, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.z,   0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.x,   0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.y,  90.951955774481, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.z, -17.938271420558, 1e-6);

  EXPECT_NEAR(atomPositions[3].gradient.x,   0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[3].gradient.y, -90.951955774481, 1e-6);
  EXPECT_NEAR(atomPositions[3].gradient.z,  17.938271420558, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.x,   0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.y, -52.398564956433, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.z,   0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.x,   0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.y, -90.951955774481, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.z, -17.938271420558, 1e-6);

  for(Atom& atom: atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  std::pair<EnergyStatus, double3x3> pressureInfo = 
    Interactions::computeInterMolecularEnergyStrainDerivative(system.forceField, system.components, system.simulationBox,
                                                      system.spanOfMoleculeAtoms());

  EXPECT_NEAR(atomPositions[0].gradient.x,   0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[0].gradient.y,  90.951955774481, 1e-6);
  EXPECT_NEAR(atomPositions[0].gradient.z,  17.938271420558, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.x,   0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.y,  52.398564956433, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.z,   0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.x,   0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.y,  90.951955774481, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.z, -17.938271420558, 1e-6);

  EXPECT_NEAR(atomPositions[3].gradient.x,   0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[3].gradient.y, -90.951955774481, 1e-6);
  EXPECT_NEAR(atomPositions[3].gradient.z,  17.938271420558, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.x,   0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.y, -52.398564956433, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.z,   0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.x,   0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.y, -90.951955774481, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.z, -17.938271420558, 1e-6);
}

TEST(Gradients, Test_2_CO2_in_ITQ_29_2x2x2_framework_molecule)
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
      Atom(double3(0.0, 0.0,  1.149), 0.0, 1.0, 0, 4, 1, 0),
      Atom(double3(0.0, 0.0,  0.0  ), 0.0, 1.0, 0, 3, 1, 0),
      Atom(double3(0.0, 0.0, -1.149), 0.0, 1.0, 0, 4, 1, 0)
    }, 5, 21);

  System system = System(0, std::nullopt, 300.0, 1e4, forceField, { f }, { c }, { 2 }, 5);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  atomPositions[0].position = double3(5.93355, 7.93355, 7.08255);
  atomPositions[1].position = double3(5.93355, 7.93355, 5.93355);
  atomPositions[2].position = double3(5.93355, 7.93355, 4.78455);
  atomPositions[3].position = double3(5.93355, 3.93355, 7.08255);
  atomPositions[4].position = double3(5.93355, 3.93355, 5.93355);
  atomPositions[5].position = double3(5.93355, 3.93355, 4.78455);

  std::pair<ForceFactor, ForceFactor> factorFrameworkMolecular =
    Interactions::computeFrameworkMoleculeGradient(system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms());

  EXPECT_NEAR(factorFrameworkMolecular.first.energy  * Units::EnergyToKelvin, -1932.15586114, 1e-6);
  EXPECT_NEAR(factorFrameworkMolecular.second.energy * Units::EnergyToKelvin,     0.00000000, 1e-6);

  EXPECT_NEAR(atomPositions[0].gradient.x,    0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[0].gradient.y, -131.516544539514, 1e-6);
  EXPECT_NEAR(atomPositions[0].gradient.z,  -90.888952502176, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.x,    0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.y,  -48.847937726172, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.z,    0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.x,    0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.y, -131.516544539514, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.z,   90.888952502176, 1e-6);

  EXPECT_NEAR(atomPositions[3].gradient.x,    0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[3].gradient.y,  131.516544539514, 1e-6);
  EXPECT_NEAR(atomPositions[3].gradient.z,  -90.888952502176, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.x,    0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.y,   48.847937726172, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.z,    0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.x,    0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.y,  131.516544539514, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.z,   90.888952502176, 1e-6);

  for(Atom& atom: atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  std::pair<EnergyStatus, double3x3> pressureInfo =
    Interactions::computeFrameworkMoleculeEnergyStrainDerivative(system.forceField, system.frameworkComponents, system.components,
                                                       system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms());

  EXPECT_NEAR(atomPositions[0].gradient.x,    0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[0].gradient.y, -131.516544539514, 1e-6);
  EXPECT_NEAR(atomPositions[0].gradient.z,  -90.888952502176, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.x,    0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.y,  -48.847937726172, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.z,    0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.x,    0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.y, -131.516544539514, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.z,   90.888952502176, 1e-6);

  EXPECT_NEAR(atomPositions[3].gradient.x,    0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[3].gradient.y,  131.516544539514, 1e-6);
  EXPECT_NEAR(atomPositions[3].gradient.z,  -90.888952502176, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.x,    0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.y,   48.847937726172, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.z,    0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.x,    0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.y,  131.516544539514, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.z,   90.888952502176, 1e-6);
}

TEST(Gradients, Test_2_CO2_in_ITQ_29_2x2x2_NonEwald)
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

  std::pair<ForceFactor, ForceFactor> factorFrameworkMolecular =
    Interactions::computeFrameworkMoleculeGradient(system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms());
  std::pair<ForceFactor, ForceFactor> factorInterMolecular =
    Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox, system.spanOfMoleculeAtoms());

  EXPECT_NEAR(factorFrameworkMolecular.first.energy  * Units::EnergyToKelvin, -1932.15586114, 1e-6);
  EXPECT_NEAR(factorFrameworkMolecular.second.energy * Units::EnergyToKelvin,   554.41444763, 1e-6);
  EXPECT_NEAR(factorInterMolecular.first.energy  * Units::EnergyToKelvin,      -242.36960932, 1e-6);
  EXPECT_NEAR(factorInterMolecular.second.energy * Units::EnergyToKelvin,       162.41877650, 1e-6);

  EXPECT_NEAR(atomPositions[0].gradient.x,     0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[0].gradient.y,   466.705848885190, 1e-6);
  EXPECT_NEAR(atomPositions[0].gradient.z,   223.065047349758, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.x,     0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.y, -1259.011389396075, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.z,     0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.x,     0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.y,   466.705848885189, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.z,  -223.065047349755, 1e-6);

  EXPECT_NEAR(atomPositions[3].gradient.x,     0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[3].gradient.y,  -466.705848885190, 1e-6);
  EXPECT_NEAR(atomPositions[3].gradient.z,   223.065047349758, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.x,     0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.y,  1259.011389396075, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.z,     0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.x,     0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.y,  -466.705848885189, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.z,  -223.065047349755, 1e-6);

  for(Atom& atom: atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  std::pair<EnergyStatus, double3x3> pressureInfo1 =
    Interactions::computeFrameworkMoleculeEnergyStrainDerivative(system.forceField, system.frameworkComponents, system.components,
                                                       system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms());

  std::pair<EnergyStatus, double3x3> pressureInfo2 =
    Interactions::computeInterMolecularEnergyStrainDerivative(system.forceField, system.components, system.simulationBox,
                                                      system.spanOfMoleculeAtoms());

  EXPECT_NEAR(atomPositions[0].gradient.x,     0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[0].gradient.y,   466.705848885190, 1e-6);
  EXPECT_NEAR(atomPositions[0].gradient.z,   223.065047349758, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.x,     0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.y, -1259.011389396075, 1e-6);
  EXPECT_NEAR(atomPositions[1].gradient.z,     0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.x,     0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.y,   466.705848885189, 1e-6);
  EXPECT_NEAR(atomPositions[2].gradient.z,  -223.065047349755, 1e-6);

  EXPECT_NEAR(atomPositions[3].gradient.x,     0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[3].gradient.y,  -466.705848885190, 1e-6);
  EXPECT_NEAR(atomPositions[3].gradient.z,   223.065047349758, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.x,     0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.y,  1259.011389396075, 1e-6);
  EXPECT_NEAR(atomPositions[4].gradient.z,     0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.x,     0.000000000000, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.y,  -466.705848885189, 1e-6);
  EXPECT_NEAR(atomPositions[5].gradient.z,  -223.065047349755, 1e-6);
}

TEST(Gradients, Test_2_CO2_in_ITQ_29_2x2x2_Ewald)
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

  RunningEnergy energy, rigidenergy;
  Interactions::computeEwaldFourierRigidEnergy(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                        system.fixedFrameworkStoredEik,
                                        system.forceField, system.simulationBox,
                                        system.spanOfFrameworkAtoms(), rigidenergy);
  ForceFactor factorEwald =
    Interactions::computeEwaldFourierGradient(system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik,
                                              system.forceField, system.simulationBox, 
                                              system.components, system.numberOfMoleculesPerComponent,
                                              system.spanOfMoleculeAtoms());

  EXPECT_NEAR((factorEwald.energy - rigidenergy.ewald)  * Units::EnergyToKelvin, -759.67572774 + 38.02930863, 1e-4);


  EXPECT_NEAR(atomPositions[0].gradient.x,     0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[0].gradient.y,  -362.766142495638, 1e-4);
  EXPECT_NEAR(atomPositions[0].gradient.z,  -241.771707768611, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.x,     0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.y,   684.800882899876, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.z,     0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.x,     0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.y,  -362.766142495638, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.z,   241.771707768619, 1e-4);

  EXPECT_NEAR(atomPositions[3].gradient.x,     0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[3].gradient.y,   362.766142495638, 1e-4);
  EXPECT_NEAR(atomPositions[3].gradient.z,  -241.771707768611, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.x,     0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.y,  -684.800882899876, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.z,     0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.x,     0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.y,   362.766142495638, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.z,   241.771707768619, 1e-4);

  for(Atom& atom: atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  std::pair<EnergyStatus, double3x3> pressureInfo =
    Interactions::computeEwaldFourierEnergyStrainDerivative(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                            system.fixedFrameworkStoredEik, system.storedEik, system.forceField, system.simulationBox,
                                                            system.frameworkComponents, system.components, system.numberOfMoleculesPerComponent,
                                                            system.spanOfMoleculeAtoms());

  EXPECT_NEAR(atomPositions[0].gradient.x,     0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[0].gradient.y,  -362.766142495638, 1e-4);
  EXPECT_NEAR(atomPositions[0].gradient.z,  -241.771707768611, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.x,     0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.y,   684.800882899876, 1e-4);
  EXPECT_NEAR(atomPositions[1].gradient.z,     0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.x,     0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.y,  -362.766142495638, 1e-4);
  EXPECT_NEAR(atomPositions[2].gradient.z,   241.771707768619, 1e-4);

  EXPECT_NEAR(atomPositions[3].gradient.x,     0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[3].gradient.y,   362.766142495638, 1e-4);
  EXPECT_NEAR(atomPositions[3].gradient.z,  -241.771707768611, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.x,     0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.y,  -684.800882899876, 1e-4);
  EXPECT_NEAR(atomPositions[4].gradient.z,     0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.x,     0.000000000000, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.y,   362.766142495638, 1e-4);
  EXPECT_NEAR(atomPositions[5].gradient.z,   241.771707768619, 1e-4);
}

TEST(Gradients, Test_2_CO2_in_ITQ_29_2x2x2_Total)
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

  [[maybe_unused]] std::pair<ForceFactor, ForceFactor> factorFrameworkMolecular =
    Interactions::computeFrameworkMoleculeGradient(system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms());
  [[maybe_unused]] std::pair<ForceFactor, ForceFactor> factorInterMolecular =
    Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox, system.spanOfMoleculeAtoms());

  RunningEnergy energy, rigidenergy;
  Interactions::computeEwaldFourierRigidEnergy(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                        system.fixedFrameworkStoredEik,
                                        system.forceField, system.simulationBox,
                                        system.spanOfFrameworkAtoms(), rigidenergy);
  [[maybe_unused]] ForceFactor factorEwald =
    Interactions::computeEwaldFourierGradient(system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik,
                                              system.forceField, system.simulationBox, 
                                              system.components, system.numberOfMoleculesPerComponent,
                                              system.spanOfMoleculeAtoms());


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


  for(Atom& atom: atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }
  [[maybe_unused]] RunningEnergy gradientEnergy = system.computeTotalGradients();

  //EXPECT_NEAR(gradientEnergy.total()  * Units::EnergyToKelvin, -2179.338665434245, 1e-4);

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
}

/*
TEST(Gradients, Test_2_CO2_in_ITQ_29_1x1x1_numerical)
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
      //Atom(double3(0.0, 0.0,  1.149), 0.0, 1.0, 0, 4, 1, 0),
      //Atom(double3(0.0, 0.0,  0.0  ), 0.0, 1.0, 0, 3, 1, 0),
      //Atom(double3(0.0, 0.0, -1.149), 0.0, 1.0, 0, 4, 1, 0)
      Atom(double3(0.0, 0.0,  1.149), -0.3256, 1.0, 0, 4, 1, 0),
      Atom(double3(0.0, 0.0,  0.0  ),  0.6512, 1.0, 0, 3, 1, 0),
      Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)
    }, 5, 21);

  System system = System(0, std::nullopt, 300.0, 1e4, forceField, { f }, { c }, { 2 }, 5);

  //std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  //std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  //spanOfMoleculeAtoms[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  //spanOfMoleculeAtoms[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  //spanOfMoleculeAtoms[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  //spanOfMoleculeAtoms[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  //spanOfMoleculeAtoms[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  //spanOfMoleculeAtoms[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);

  //[[maybe_unused]] ForceFactor factorFrameworkMolecular = 
  //  Interactions::computeFrameworkMoleculeGradient(system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms());
  //[[maybe_unused]] ForceFactor factorInterMolecular = 
  //  Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox, system.spanOfMoleculeAtoms());

  //// make copy
  //std::vector<Atom> atomPositions = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end());

  //for (Atom& atom : atomPositions)
  //{
  //  atom.gradient = double3(0.0, 0.0, 0.0);
  //}

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  spanOfMoleculeAtoms[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);
  std::vector<Atom> atomPositions = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end());

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  //[[maybe_unused]] ForceFactor factorFrameworkMolecular = 
  //  Interactions::computeFrameworkMoleculeGradient(system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms());
  [[maybe_unused]] ForceFactor factorInterMolecular = 
    Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox, system.spanOfMoleculeAtoms());


  double delta = 1e-3;
  double tolerance = 1e-5;
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
*/

TEST(Gradients, Test_CO2_in_ITQ_29_1x1x1)
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
  Framework f = Framework(0, forceField, "ITQ-29", SimulationBox(11.8671, 11.8671, 11.8671),
    517,
    { // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
      Atom(double3(0.3683, 0.1847, 0),       2.05,  1.0, 0, 0, 0, 0),
      Atom(double3(0.5,    0.2179, 0),      -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.2939, 0.2939, 0),      -1.025, 1.0, 0, 1, 0, 0),
      Atom(double3(0.3429, 0.1098, 0.1098), -1.025, 1.0, 0, 1, 0, 0)
    },
    int3(2, 2, 2));
  Component c = Component(0, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
    {  // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId, uint8_t groupId
       Atom(double3(0.0, 0.0,  1.149), -0.3256, 1.0, 0, 4, 1, 0),
       Atom(double3(0.0, 0.0,  0.0),    0.6512, 1.0, 0, 3, 1, 0),
       Atom(double3(0.0, 0.0, -1.149), -0.3256, 1.0, 0, 4, 1, 0)
    }, 5, 21);

  System system = System(0, std::nullopt, 300.0, 1e4, forceField, { f }, { c }, { 2 }, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();

  spanOfMoleculeAtoms[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  std::vector<Atom> atomPositions = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end()); 

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  [[maybe_unused]] std::pair<ForceFactor, ForceFactor> factor = 
    Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox, system.spanOfMoleculeAtoms());
  [[maybe_unused]] std::pair<ForceFactor, ForceFactor> factor2 = 
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

  [[maybe_unused]] std::pair<ForceFactor, ForceFactor> factorFrameworkMolecular = 
    Interactions::computeFrameworkMoleculeGradient(system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms());
  [[maybe_unused]] std::pair<ForceFactor, ForceFactor> factorInterMolecular = 
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

  System system = System(0, std::nullopt, 300.0, 1e4, forceField, { f }, { c }, { 10 }, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();
  std::vector<Atom> atomPositions = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end()); 

  std::cout << "CHECK: " << spanOfMoleculeAtoms.size() << std::endl;

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  [[maybe_unused]] std::pair<ForceFactor, ForceFactor> factor = 
    Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox, system.spanOfMoleculeAtoms());
  [[maybe_unused]] std::pair<ForceFactor, ForceFactor> factor2 = 
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

  [[maybe_unused]] std::pair<ForceFactor, ForceFactor> factor = 
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
