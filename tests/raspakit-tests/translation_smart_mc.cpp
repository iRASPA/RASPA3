#include <gtest/gtest.h>

import std;

import int3;
import double3;
import double3x3;
import units;
import atom;
import atom_dynamics;
import pseudo_atom;
import vdwparameters;
import forcefield;
import framework;
import component;
import system;
import simulationbox;
import running_energy;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;

// The force-biased translation move needs the total (center-of-mass) force on a single molecule, obtained cheaply
// from single-molecule kernels (framework + inter-molecular real space + reciprocal Ewald using the maintained
// structure factor). This test verifies that this fast-path force equals the exact force obtained by summing the
// full-system per-atom gradient over the atoms of that molecule.

static double3 sumMoleculeGradient(std::span<const AtomDynamics> dynamics, std::size_t begin, std::size_t count)
{
  double3 sum{};
  for (std::size_t i = begin; i != begin + count; ++i) sum += dynamics[i].gradient;
  return sum;
}

TEST(translation_smart_mc, single_molecule_force_matches_full_gradient)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

  Framework f = Framework::makeITQ29(forceField, int3(1, 1, 1));
  Component c = Component::makeCO2(forceField, 0, true);

  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {c}, {}, {2}, 5);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  atomData[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  atomData[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  atomData[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  atomData[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  atomData[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  atomData[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);

  // Populate the rigid-framework and total structure factors for the current configuration.
  system.precomputeTotalRigidEnergy();
  Interactions::computeEwaldFourierEnergy(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                          system.fixedFrameworkStoredEik, system.storedEik, system.forceField,
                                          system.simulationBox, system.components,
                                          system.numberOfMoleculesPerComponent, system.spanOfMoleculeAtoms());

  // Reference: full-system per-atom gradient.
  std::span<AtomDynamics> dynamics = system.spanOfMoleculeDynamics();
  for (AtomDynamics &d : dynamics) d.gradient = double3(0.0, 0.0, 0.0);

  Interactions::computeFrameworkMoleculeGradient(system.forceField, system.simulationBox,
                                                 system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(), dynamics,
                                                 system.interpolationGrids);
  Interactions::computeInterMolecularGradient(system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(),
                                              dynamics);
  Interactions::computeEwaldFourierGradient(system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik,
                                            system.fixedFrameworkStoredEik, system.forceField, system.simulationBox,
                                            system.components, system.numberOfMoleculesPerComponent,
                                            system.spanOfMoleculeAtoms(), dynamics);

  double3 referenceForceMol0 = -sumMoleculeGradient(dynamics, 0, 3);
  double3 referenceForceMol1 = -sumMoleculeGradient(dynamics, 3, 3);

  // Fast path used by the force-biased move for molecule 0.
  {
    std::span<const Atom> selected = std::span<const Atom>(system.spanOfMoleculeAtoms()).subspan(0, 3);
    std::vector<AtomDynamics> fastDynamics(3);
    Interactions::computeFrameworkMoleculeGradient(system.forceField, system.simulationBox,
                                                   system.spanOfFrameworkAtoms(), selected, fastDynamics,
                                                   system.interpolationGrids);
    Interactions::computeInterMolecularGradientMolecule(system.forceField, system.simulationBox,
                                                        system.spanOfMoleculeAtoms(), selected, fastDynamics);
    Interactions::computeEwaldFourierGradientSingleMolecule(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                            system.storedEik, system.forceField, system.simulationBox,
                                                            selected, fastDynamics);
    double3 fastForce{};
    for (const AtomDynamics &d : fastDynamics) fastForce += d.gradient;
    fastForce = -fastForce;

    EXPECT_NEAR(fastForce.x, referenceForceMol0.x, 1e-6);
    EXPECT_NEAR(fastForce.y, referenceForceMol0.y, 1e-6);
    EXPECT_NEAR(fastForce.z, referenceForceMol0.z, 1e-6);
  }

  // Fast path for molecule 1.
  {
    std::span<const Atom> selected = std::span<const Atom>(system.spanOfMoleculeAtoms()).subspan(3, 3);
    std::vector<AtomDynamics> fastDynamics(3);
    Interactions::computeFrameworkMoleculeGradient(system.forceField, system.simulationBox,
                                                   system.spanOfFrameworkAtoms(), selected, fastDynamics,
                                                   system.interpolationGrids);
    Interactions::computeInterMolecularGradientMolecule(system.forceField, system.simulationBox,
                                                        system.spanOfMoleculeAtoms(), selected, fastDynamics);
    Interactions::computeEwaldFourierGradientSingleMolecule(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                            system.storedEik, system.forceField, system.simulationBox,
                                                            selected, fastDynamics);
    double3 fastForce{};
    for (const AtomDynamics &d : fastDynamics) fastForce += d.gradient;
    fastForce = -fastForce;

    EXPECT_NEAR(fastForce.x, referenceForceMol1.x, 1e-6);
    EXPECT_NEAR(fastForce.y, referenceForceMol1.y, 1e-6);
    EXPECT_NEAR(fastForce.z, referenceForceMol1.z, 1e-6);
  }
}
