#include <gtest/gtest.h>

import <vector>;
import <tuple>;
import <algorithm>;

import double3;

import forcefield;
import component;
import system;
import simulationbox;
import energy_factor;

TEST(dudlambda, Test_20_Na_Cl_in_Box_25x25x25)
{
  ForceField forceField = ForceField(
    { PseudoAtom("Si",   28.0855,   2.05,  14, false),
      PseudoAtom("O",    15.999,   -1.025,  8, false),
      PseudoAtom("CH4",  16.04246,  0.0,    6, false),
      PseudoAtom("Na+",  12.0,      0.0, 6, false),
      PseudoAtom("Cl-",  15.9994,   0.0, 8, false),
    },
    { VDWParameters(22.0 / 1.2027242847, 2.30),
      VDWParameters(53.0 / 1.2027242847, 3.3),
      VDWParameters(158.5 / 1.2027242847, 3.72),
      VDWParameters(15.0966 / 1.2027242847, 2.65755),
      VDWParameters(142.562 / 1.2027242847, 3.51932)
    },
    ForceField::MixingRule::Lorentz_Berthelot,
    12.0,
    true,
    false);
  Component na = Component(0,
    "Na",
    43.9988,
    SimulationBox(25.0, 25.0, 25.0),
    304.1282, 7377300.0, 0.22394,
    {
       Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 3, 0, 0),
    }, 5);
  Component cl = Component(1,
    "Cl",
    43.9988,
    SimulationBox(25.0, 25.0, 25.0),
    304.1282, 7377300.0, 0.22394,
    {
       Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 4, 1, 0),
    }, 5);

  System system = System(0, 300.0, 1e4, forceField, { na, cl }, { 20, 20 }, 5);

  std::fill(system.forceField.data.begin(), system.forceField.data.end(), VDWParameters(0.0, 1.0));

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

  
  
  
  
  

  system.atomPositions[8].scalingCoulomb = 0.45;
  system.atomPositions[2].scalingCoulomb = 0.5;
  system.atomPositions[12].scalingCoulomb = 0.4;
  system.atomPositions[4].scalingCoulomb = 0.6;
  system.atomPositions[14].scalingCoulomb = 0.7;
  EnergyFactor factor = system.computeInterMolecularGradient();

  double delta = 1e-6;
  double tolerance = 1e-4;

  system.atomPositions[8].scalingCoulomb = 0.45 + 0.5 * delta;
  system.atomPositions[2].scalingCoulomb = 0.5 + 0.5 * delta;
  system.atomPositions[12].scalingCoulomb = 0.4 + 0.5 * delta;
  system.atomPositions[4].scalingCoulomb = 0.6 + 0.5 * delta;
  system.atomPositions[14].scalingCoulomb = 0.7 + 0.5 * delta;
  EnergyFactor energyForward = system.computeInterMolecularGradient();

  system.atomPositions[8].scalingCoulomb = 0.45 - 0.5 * delta;
  system.atomPositions[2].scalingCoulomb = 0.5 - 0.5 * delta;
  system.atomPositions[12].scalingCoulomb = 0.4 - 0.5 * delta;
  system.atomPositions[4].scalingCoulomb = 0.6 - 0.5 * delta;
  system.atomPositions[14].scalingCoulomb = 0.7 - 0.5 * delta;
  EnergyFactor energyBackward = system.computeInterMolecularGradient();
  
  double dUdlambda = (energyForward.energy - energyBackward.energy) / delta;

  EXPECT_NEAR(factor.dUdlambda, dUdlambda, tolerance) << " ratio: " << factor.dUdlambda / dUdlambda << " " << dUdlambda / factor.dUdlambda;
}