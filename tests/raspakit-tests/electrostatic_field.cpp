#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <cstddef>
#include <span>
#include <vector>

import int3;
import double3;
import double3x3;
import factory;
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
import gradient_factor;
import running_energy;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;
import energy_status;


// Test if E = F / q for framework-molecule interactions
TEST(electrostatic_field, Test_2_CO2_in_ITQ_29_2x2x2_NonEwald)
{ 
  double tolerance = 1e-6;

  ForceField forceField =  ForceField({{"Si", true, 28.0855, 2.05, 0.0, 14, false},
                     {"O", true, 15.999, -1.025, 0.0, 8, false},
                     {"C_co2", false, 12.0, 0.6512, 0.2, 6, false},
                     {"O_co2", false, 15.9994, -0.3256, 0.1, 8, false}},

                    {{0.0, 2.30},
                     {0.0, 3.30},
                     {0.0, 2.745},
                     {0.0, 3.017}},
                    ForceField::MixingRule::Lorentz_Berthelot, 11.8, 11.8, 11.8, true, false, true);
  
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);

  Framework f = TestFactories::makeITQ29(forceField, int3(2, 2, 2));

  Component CO2 = Component(0, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                   {Atom({0, 0,  1.149}, -0.3256, 1.0, 0, 4, 0, false, false),
                    Atom({0, 0,  0.000},  0.6512, 1.0, 0, 3, 0, false, false),
                    Atom({0, 0, -1.149}, -0.3256, 1.0, 0, 4, 0, false, false)},
                   5, 21);
  
  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {CO2}, {1}, 5);
  
  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms(); 
  atomPositions[0].position = double3(5.93355, 7.93355, 2.0 + 5.93355 + 1.149);
  atomPositions[1].position = double3(5.93355, 7.93355, 2.0 + 5.93355 + 0.0);
  atomPositions[2].position = double3(5.93355, 7.93355, 2.0 + 5.93355 - 1.149);
      
  //system.precomputeTotalRigidEnergy();
  //RunningEnergy factorEwald = Interactions::computeEwaldFourierGradient(
  //    system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.fixedFrameworkStoredEik,
  //    system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
  //    system.spanOfMoleculeAtoms());

  std::span<double3> moleculeElectricField = system.spanOfMoleculeElectricField();
  std::span<Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  std::fill(moleculeElectricField.begin(), moleculeElectricField.end(), double3(0.0, 0.0, 0.0));

  RunningEnergy energy1 = Interactions::computeFrameworkMoleculeElectricField(
        system.forceField, system.simulationBox, moleculeElectricField, frameworkAtoms, atomPositions);

  RunningEnergy energy2 = Interactions::computeFrameworkMoleculeGradient(
      system.forceField, system.simulationBox, frameworkAtoms, atomPositions, system.interpolationGrids);


  EXPECT_NEAR(energy1.frameworkMoleculeVDW, energy2.frameworkMoleculeVDW, tolerance);
  EXPECT_NEAR(energy1.frameworkMoleculeCharge, energy2.frameworkMoleculeCharge, tolerance);
  EXPECT_NEAR(moleculeElectricField[0].x, -atomPositions[0].gradient.x / atomPositions[0].charge, tolerance);
  EXPECT_NEAR(moleculeElectricField[0].y, -atomPositions[0].gradient.y / atomPositions[0].charge, tolerance);
  EXPECT_NEAR(moleculeElectricField[0].z, -atomPositions[0].gradient.z / atomPositions[0].charge, tolerance);
  EXPECT_NEAR(moleculeElectricField[1].x, -atomPositions[1].gradient.x / atomPositions[1].charge, tolerance);
  EXPECT_NEAR(moleculeElectricField[1].y, -atomPositions[1].gradient.y / atomPositions[1].charge, tolerance);
  EXPECT_NEAR(moleculeElectricField[1].z, -atomPositions[1].gradient.z / atomPositions[1].charge, tolerance);
  EXPECT_NEAR(moleculeElectricField[2].x, -atomPositions[2].gradient.x / atomPositions[2].charge, tolerance);
  EXPECT_NEAR(moleculeElectricField[2].y, -atomPositions[2].gradient.y / atomPositions[2].charge, tolerance);
  EXPECT_NEAR(moleculeElectricField[2].z, -atomPositions[2].gradient.z / atomPositions[2].charge, tolerance);
}

// Test if E = F / q for Ewald interactions
// True only if omitInterInteractions = true
// Otherwise, even a single molecule has interctions with its periodic images,
// leading to a small difference
TEST(electrostatic_field, Test_2_CO2_in_ITQ_29_2x2x2_Ewald)
{ 
  double tolerance = 1e-6;

  ForceField forceField =  ForceField({{"Si", true, 28.0855, 2.05, 0.0, 14, false},
                     {"O", true, 15.999, -1.025, 0.0, 8, false},
                     {"C_co2", false, 12.0, 0.6512, 0.2, 6, false},
                     {"O_co2", false, 15.9994, -0.3256, 0.1, 8, false}},

                    {{0.0, 2.30},
                     {0.0, 3.30},
                     {0.0, 2.745},
                     {0.0, 3.017}},
                    ForceField::MixingRule::Lorentz_Berthelot, 11.8, 11.8, 11.8, true, false, true);
  
  forceField.automaticEwald = false;
  forceField.EwaldAlpha = 0.25;
  forceField.numberOfWaveVectors = int3(8, 8, 8);
  forceField.omitInterInteractions = true;

  Framework f = TestFactories::makeITQ29(forceField, int3(2, 2, 2));

  Component CO2 = Component(0, forceField, "CO2", 304.1282, 7377300.0, 0.22394,
                   {Atom({0, 0,  1.149}, -0.3256, 1.0, 0, 4, 0, false, false),
                    Atom({0, 0,  0.000},  0.6512, 1.0, 0, 3, 0, false, false),
                    Atom({0, 0, -1.149}, -0.3256, 1.0, 0, 4, 0, false, false)},
                   5, 21);
  
  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {CO2}, {1}, 5);
  
  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms(); 
  atomPositions[0].position = double3(5.83355, 7.83355, 1.8 + 5.93355 + 1.149);
  atomPositions[1].position = double3(5.93355, 7.93355, 2.0 + 5.93355 + 0.0);
  atomPositions[2].position = double3(6.03355, 8.03355, 2.2 + 5.93355 - 1.149);
      

  std::span<double3> moleculeElectricField = system.spanOfMoleculeElectricField();

  std::fill(moleculeElectricField.begin(), moleculeElectricField.end(), double3(0.0, 0.0, 0.0));

  system.precomputeTotalRigidEnergy();
  RunningEnergy energy1 = Interactions::computeEwaldFourierElectricField(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.totalEik,
      system.forceField, system.simulationBox, moleculeElectricField, system.components, system.numberOfMoleculesPerComponent,
      system.spanOfMoleculeAtoms());

  RunningEnergy energy2 = Interactions::computeEwaldFourierGradient(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.fixedFrameworkStoredEik,
      system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
      system.spanOfMoleculeAtoms());

  EXPECT_NEAR(energy1.ewald_fourier, energy2.ewald_fourier, tolerance);
  EXPECT_NEAR(energy1.ewald_self, energy2.ewald_self, tolerance);
  EXPECT_NEAR(energy1.ewald_exclusion, energy2.ewald_exclusion, tolerance);
  EXPECT_NEAR(moleculeElectricField[0].x, -atomPositions[0].gradient.x / atomPositions[0].charge, tolerance);
  EXPECT_NEAR(moleculeElectricField[0].y, -atomPositions[0].gradient.y / atomPositions[0].charge, tolerance);
  EXPECT_NEAR(moleculeElectricField[0].z, -atomPositions[0].gradient.z / atomPositions[0].charge, tolerance);
  EXPECT_NEAR(moleculeElectricField[1].x, -atomPositions[1].gradient.x / atomPositions[1].charge, tolerance);
  EXPECT_NEAR(moleculeElectricField[1].y, -atomPositions[1].gradient.y / atomPositions[1].charge, tolerance);
  EXPECT_NEAR(moleculeElectricField[1].z, -atomPositions[1].gradient.z / atomPositions[1].charge, tolerance);
  EXPECT_NEAR(moleculeElectricField[2].x, -atomPositions[2].gradient.x / atomPositions[2].charge, tolerance);
  EXPECT_NEAR(moleculeElectricField[2].y, -atomPositions[2].gradient.y / atomPositions[2].charge, tolerance);
  EXPECT_NEAR(moleculeElectricField[2].z, -atomPositions[2].gradient.z / atomPositions[2].charge, tolerance);
}

/*
TEST(electrostatic_field, Test_2_CO2_in_ITQ_29_2x2x2)
{
  double tolerance = 1e-5;

  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Component c = TestFactories::makeCO2(forceField, 0, true);
  Framework f = TestFactories::makeITQ29(forceField, int3(2, 2, 2));

  forceField.computePolarization = true;
  forceField.omitInterPolarization = true;
  forceField.omitInterInteractions = true;
  forceField.omitEwaldFourier = false;

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {2}, 5);
  system.forceField.EwaldAlpha = 0.25;
  system.forceField.numberOfWaveVectors = int3(8, 8, 8);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtomPositions = system.spanOfFrameworkAtoms();
  std::span<double> moleculeElectricPotential = system.spanOfMoleculeElectricPotential();
  spanOfMoleculeAtoms[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[0].scalingCoulomb = 1.0;
  spanOfMoleculeAtoms[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[1].scalingCoulomb = 1.0;
  spanOfMoleculeAtoms[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[2].scalingCoulomb = 1.0;

  spanOfMoleculeAtoms[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[3].scalingCoulomb = 1.0;
  spanOfMoleculeAtoms[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[4].scalingCoulomb = 1.0;
  spanOfMoleculeAtoms[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[5].scalingCoulomb = 1.0;

  system.computeTotalElectricField();
  std::span<double3> electricField = system.spanOfMoleculeElectricField();

  // create copy
  std::vector<Atom> atomPositions = std::vector<Atom>(spanOfMoleculeAtoms.begin(), spanOfMoleculeAtoms.end());

  double delta = 1e-5;
  double3 gradient;
  for (size_t i = 0; i < atomPositions.size(); ++i)
  {
    double x1, x2, y1, y2, z1, z2;

    // finite difference x
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x + 0.5 * delta;
    std::fill(moleculeElectricPotential.begin(), moleculeElectricPotential.end(), 0.0);
    Interactions::computeFrameworkMoleculeElectricPotential(
        system.forceField, system.simulationBox, moleculeElectricPotential, frameworkAtomPositions, atomPositions);
    Interactions::computeEwaldFourierElectricPotential(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                       system.fixedFrameworkStoredEik, moleculeElectricPotential,
                                                       system.forceField, system.simulationBox, system.components,
                                                       system.numberOfMoleculesPerComponent, atomPositions);
    x2 = moleculeElectricPotential[i];

    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x - 0.5 * delta;
    std::fill(moleculeElectricPotential.begin(), moleculeElectricPotential.end(), 0.0);
    Interactions::computeFrameworkMoleculeElectricPotential(
        system.forceField, system.simulationBox, moleculeElectricPotential, frameworkAtomPositions, atomPositions);
    Interactions::computeEwaldFourierElectricPotential(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                       system.fixedFrameworkStoredEik, moleculeElectricPotential,
                                                       system.forceField, system.simulationBox, system.components,
                                                       system.numberOfMoleculesPerComponent, atomPositions);
    x1 = moleculeElectricPotential[i];
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x;

    // finite difference y
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y + 0.5 * delta;
    std::fill(moleculeElectricPotential.begin(), moleculeElectricPotential.end(), 0.0);
    Interactions::computeFrameworkMoleculeElectricPotential(
        system.forceField, system.simulationBox, moleculeElectricPotential, frameworkAtomPositions, atomPositions);
    Interactions::computeEwaldFourierElectricPotential(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                       system.fixedFrameworkStoredEik, moleculeElectricPotential,
                                                       system.forceField, system.simulationBox, system.components,
                                                       system.numberOfMoleculesPerComponent, atomPositions);
    y2 = moleculeElectricPotential[i];

    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y - 0.5 * delta;
    std::fill(moleculeElectricPotential.begin(), moleculeElectricPotential.end(), 0.0);
    Interactions::computeFrameworkMoleculeElectricPotential(
        system.forceField, system.simulationBox, moleculeElectricPotential, frameworkAtomPositions, atomPositions);
    Interactions::computeEwaldFourierElectricPotential(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                       system.fixedFrameworkStoredEik, moleculeElectricPotential,
                                                       system.forceField, system.simulationBox, system.components,
                                                       system.numberOfMoleculesPerComponent, atomPositions);
    y1 = moleculeElectricPotential[i];
    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y;

    // finite difference z
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z + 0.5 * delta;
    std::fill(moleculeElectricPotential.begin(), moleculeElectricPotential.end(), 0.0);
    Interactions::computeFrameworkMoleculeElectricPotential(
        system.forceField, system.simulationBox, moleculeElectricPotential, frameworkAtomPositions, atomPositions);
    Interactions::computeEwaldFourierElectricPotential(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                       system.fixedFrameworkStoredEik, moleculeElectricPotential,
                                                       system.forceField, system.simulationBox, system.components,
                                                       system.numberOfMoleculesPerComponent, atomPositions);
    z2 = moleculeElectricPotential[i];

    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z - 0.5 * delta;
    std::fill(moleculeElectricPotential.begin(), moleculeElectricPotential.end(), 0.0);
    Interactions::computeFrameworkMoleculeElectricPotential(
        system.forceField, system.simulationBox, moleculeElectricPotential, frameworkAtomPositions, atomPositions);
    Interactions::computeEwaldFourierElectricPotential(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                       system.fixedFrameworkStoredEik, moleculeElectricPotential,
                                                       system.forceField, system.simulationBox, system.components,
                                                       system.numberOfMoleculesPerComponent, atomPositions);
    z1 = moleculeElectricPotential[i];
    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z;

    gradient.x = -(x2 - x1) / delta;
    gradient.y = -(y2 - y1) / delta;
    gradient.z = -(z2 - z1) / delta;

    EXPECT_NEAR(electricField[i].x, gradient.x, tolerance);
    EXPECT_NEAR(electricField[i].y, gradient.y, tolerance);
    EXPECT_NEAR(electricField[i].z, gradient.z, tolerance);
  }
}
*/
