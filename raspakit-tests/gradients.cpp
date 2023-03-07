#include <gtest/gtest.h>

import double3;

import forcefield;
import component;
import system;
import simulationbox;
import energy_factor;

TEST(Gradients, TestLennardJones)
{
  ForceField forceField = ForceField(
    { PseudoAtom("CH4", 16.04246, 0.0, 6, false) },
    { VDWParameters(158.5 / 1.2027242847, 3.72) },
    ForceField::MixingRule::Lorentz_Berthelot,
    12.0,
    false,
    false);
  Component c = Component(0,
    "methane",
    16.04246,
    SimulationBox(25.0, 25.0, 25.0),
    190.564, 45599200, 0.01142,
    { Atom(double3(0.0, 0.0,  0.0),    0.0, 1.0, 0, 0, 0) }, 5);

  System system = System(0, 300.0, 1e4, forceField, { c }, { 50 }, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::vector<Atom> atomPositions = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end());

  for (Atom& atom : atomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  EnergyFactor factor = system.computeInterMolecularGradient();

  double delta = 1e-5;
  double tolerance = 1e-6;
  double3 gradient;
  for (size_t i = 0; i < atomPositions.size(); ++i)
  {
    RunningEnergy x1, x2, y1, y2, z1, z2;

    // finite difference x
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x + 0.5 * delta;
    system.computeInterMolecularEnergy(system.simulationBox, std::span<const Atom>(atomPositions), x2);

    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x - 0.5 * delta;
    system.computeInterMolecularEnergy(system.simulationBox, std::span<const Atom>(atomPositions), x1);
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x;

    // finite difference y
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y + 0.5 * delta;
    system.computeInterMolecularEnergy(system.simulationBox, std::span<const Atom>(atomPositions), y2);

    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y - 0.5 * delta;
    system.computeInterMolecularEnergy(system.simulationBox, std::span<const Atom>(atomPositions), y1);
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y;

    // finite difference z
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z + 0.5 * delta;
    system.computeInterMolecularEnergy(system.simulationBox, std::span<const Atom>(atomPositions), z2);

    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z - 0.5 * delta;
    system.computeInterMolecularEnergy(system.simulationBox, std::span<const Atom>(atomPositions), z1);
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z;

    gradient.x = (x2.moleculeMoleculeVDW + x2.moleculeMoleculeCharge - (x1.moleculeMoleculeVDW + x1.moleculeMoleculeCharge)) / delta;
    gradient.y = (y2.moleculeMoleculeVDW + y2.moleculeMoleculeCharge - (y1.moleculeMoleculeVDW + y1.moleculeMoleculeCharge)) / delta;
    gradient.z = (z2.moleculeMoleculeVDW + z2.moleculeMoleculeCharge - (z1.moleculeMoleculeVDW + z1.moleculeMoleculeCharge)) / delta;

    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.x, gradient.x, tolerance);
    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.y, gradient.y, tolerance);
    EXPECT_NEAR(spanOfMoleculeAtoms[i].gradient.z, gradient.z, tolerance);
  }
}

