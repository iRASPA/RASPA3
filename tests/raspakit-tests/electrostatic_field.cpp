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

// Table 5, page 61, thesis D. Dubbeldam
TEST(electrostatic_field, Test_reference_system_1)
{
  ForceField forceField = ForceField(
      {
          PseudoAtom("t1", false, 1.0, 0.5, 0.0, 1, false),
          PseudoAtom("t2", false, 1.0, 1.5, 0.0, 1, false),
          PseudoAtom("t2", false, 1.0, -0.75, 0.0, 1, false),
          PseudoAtom("t4", false, 1.0, -1.25, 0.0, 1, false),
      },
      {VDWParameters(0.0, 1.0)}, ForceField::MixingRule::Lorentz_Berthelot, 50.0, 50.0, 50.0, true, false, true);

  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;
  Component c1 = Component(0, forceField, "t1", 0.0, 0.0, 0.0,
                           {
                               // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                               // uint8_t componentId, uint8_t groupId
                               Atom(double3(0.0, 0.0, 0.0), 0.5, 1.0, 0, 0, 0, 0),
                           },
                           5, 21);
  Component c2 = Component(1, forceField, "t2", 0.0, 0.0, 0.0,
                           {
                               // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                               // uint8_t componentId, uint8_t groupId
                               Atom(double3(0.0, 0.0, 0.0), 1.5, 1.0, 1, 1, 1, 0),
                           },
                           5, 21);
  Component c3 = Component(2, forceField, "t3", 0.0, 0.0, 0.0,
                           {
                               // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                               // uint8_t componentId, uint8_t groupId
                               Atom(double3(0.0, 0.0, 0.0), -0.75, 1.0, 2, 2, 2, 0),
                           },
                           5, 21);
  Component c4 = Component(3, forceField, "t4", 0.0, 0.0, 0.0,
                           {
                               // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                               // uint8_t componentId, uint8_t groupId
                               Atom(double3(0.0, 0.0, 0.0), -1.25, 1.0, 3, 3, 3, 0),
                           },
                           5, 21);

  System system =
      System(0, forceField, SimulationBox(100.0, 100.0, 100.0), 300.0, 1e4, 1.0, {}, {c1, c2, c3, c4}, {1, 1, 1, 1}, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  spanOfMoleculeAtoms[0].position = double3(-1.0, 0.0, 0.0);
  spanOfMoleculeAtoms[1].position = double3(1.0, 0.0, 0.0);
  spanOfMoleculeAtoms[2].position = double3(0.0, 1.0, 0.0);
  spanOfMoleculeAtoms[3].position = double3(0.0, -1.0, 0.0);

  system.computeTotalElectricField();
  std::span<double3> electricField = system.spanOfMoleculeElectricField();

  EXPECT_NEAR(electricField[0].x / Units::CoulombicConversionFactor, 0.332107, 1e-5);
  EXPECT_NEAR(electricField[0].y / Units::CoulombicConversionFactor, -0.17677, 1e-5);
  EXPECT_NEAR(electricField[0].z / Units::CoulombicConversionFactor, 0.0, 1e-5);

  EXPECT_NEAR(electricField[1].x / Units::CoulombicConversionFactor, -0.582107, 1e-5);
  EXPECT_NEAR(electricField[1].y / Units::CoulombicConversionFactor, -0.176777, 1e-5);
  EXPECT_NEAR(electricField[1].z / Units::CoulombicConversionFactor, 0.0, 1e-5);

  EXPECT_NEAR(electricField[2].x / Units::CoulombicConversionFactor, -0.353553, 1e-5);
  EXPECT_NEAR(electricField[2].y / Units::CoulombicConversionFactor, 0.394607, 1e-5);
  EXPECT_NEAR(electricField[2].z / Units::CoulombicConversionFactor, 0.0, 1e-5);

  EXPECT_NEAR(electricField[3].x / Units::CoulombicConversionFactor, -0.353553, 1e-5);
  EXPECT_NEAR(electricField[3].y / Units::CoulombicConversionFactor, -0.519607, 1e-5);
  EXPECT_NEAR(electricField[3].z / Units::CoulombicConversionFactor, 0.0, 1e-5);
}

// Table 5, page 61 thesis D. Dubbeldam
TEST(electrostatic_field, Test_reference_system_2)
{
  ForceField forceField = ForceField(
      {
          PseudoAtom("t1", false, 1.0, 0.5, 0.0, 1, false),
          PseudoAtom("t2", false, 1.0, -1.25, 0.0, 1, false),
          PseudoAtom("t2", false, 1.0, 1.5, 0.0, 1, false),
          PseudoAtom("t4", false, 1.0, -0.75, 0.0, 1, false),
      },
      {VDWParameters(0.0, 1.0)}, ForceField::MixingRule::Lorentz_Berthelot, 50.0, 50.0, 50.0, true, false, true);

  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;
  Component c1 = Component(
      0, forceField, "t1", 0.0, 0.0, 0.0,
      {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId,
       // uint8_t groupId
       Atom(double3(0.0, 0.0, 0.0), 0.5, 1.0, 0, 0, 0, 0), Atom(double3(0.0, 0.0, 0.0), -1.25, 1.0, 0, 1, 0, 0)},
      5, 21);
  Component c2 = Component(
      1, forceField, "t2", 0.0, 0.0, 0.0,
      {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId,
       // uint8_t groupId
       Atom(double3(0.0, 0.0, 0.0), 1.5, 1.0, 1, 2, 1, 0), Atom(double3(0.0, 0.0, 0.0), -0.75, 1.0, 1, 3, 1, 0)},
      5, 21);

  System system = System(0, forceField, SimulationBox(100.0, 100.0, 100.0), 300.0, 1e4, 1.0, {}, {c1, c2}, {1, 1}, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  spanOfMoleculeAtoms[0].position = double3(-1.0, 0.0, 0.0);
  spanOfMoleculeAtoms[1].position = double3(0.0, -1.0, 0.0);
  spanOfMoleculeAtoms[2].position = double3(1.0, 0.0, 0.0);
  spanOfMoleculeAtoms[3].position = double3(0.0, 1.0, 0.0);

  system.computeTotalElectricField();
  std::span<double3> electricField = system.spanOfMoleculeElectricField();

  EXPECT_NEAR(electricField[0].x / Units::CoulombicConversionFactor, -0.109835, 1e-5);
  EXPECT_NEAR(electricField[0].y / Units::CoulombicConversionFactor, 0.265165, 1e-5);
  EXPECT_NEAR(electricField[0].z / Units::CoulombicConversionFactor, 0.0, 1e-5);

  EXPECT_NEAR(electricField[1].x / Units::CoulombicConversionFactor, -0.53033, 1e-5);
  EXPECT_NEAR(electricField[1].y / Units::CoulombicConversionFactor, -0.34283, 1e-5);
  EXPECT_NEAR(electricField[1].z / Units::CoulombicConversionFactor, 0.0, 1e-5);

  EXPECT_NEAR(electricField[2].x / Units::CoulombicConversionFactor, -0.316942, 1e-5);
  EXPECT_NEAR(electricField[2].y / Units::CoulombicConversionFactor, -0.441941, 1e-5);
  EXPECT_NEAR(electricField[2].z / Units::CoulombicConversionFactor, 0.0, 1e-5);

  EXPECT_NEAR(electricField[3].x / Units::CoulombicConversionFactor, 0.176777, 1e-5);
  EXPECT_NEAR(electricField[3].y / Units::CoulombicConversionFactor, -0.135723, 1e-5);
  EXPECT_NEAR(electricField[3].z / Units::CoulombicConversionFactor, 0.0, 1e-5);
}

TEST(electrostatic_field, Test_2_CO2_in_ITQ_29_2x2x2)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Component c = TestFactories::makeCO2(forceField, 0, true);
  Framework f = TestFactories::makeITQ29(forceField, int3(2, 2, 2));

  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;
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
  double tolerance = 1e-4;
  double3 gradient;
  for (size_t i = 0; i < atomPositions.size(); ++i)
  {
    double x1, x2, y1, y2, z1, z2;

    // finite difference x
    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x + 0.5 * delta;
    std::fill(moleculeElectricPotential.begin(), moleculeElectricPotential.end(), 0.0);
    Interactions::computeFrameworkMoleculeElectricPotential(
        system.forceField, system.simulationBox, moleculeElectricPotential, frameworkAtomPositions, atomPositions);
    Interactions::computeInterMolecularElectricPotential(system.forceField, system.simulationBox,
                                                         moleculeElectricPotential, atomPositions);
    Interactions::computeEwaldFourierElectricPotential(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                       system.fixedFrameworkStoredEik, moleculeElectricPotential,
                                                       system.forceField, system.simulationBox, system.components,
                                                       system.numberOfMoleculesPerComponent, atomPositions);
    x2 = moleculeElectricPotential[i];

    atomPositions[i].position.x = spanOfMoleculeAtoms[i].position.x - 0.5 * delta;
    std::fill(moleculeElectricPotential.begin(), moleculeElectricPotential.end(), 0.0);
    Interactions::computeFrameworkMoleculeElectricPotential(
        system.forceField, system.simulationBox, moleculeElectricPotential, frameworkAtomPositions, atomPositions);
    Interactions::computeInterMolecularElectricPotential(system.forceField, system.simulationBox,
                                                         moleculeElectricPotential, atomPositions);
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
    Interactions::computeInterMolecularElectricPotential(system.forceField, system.simulationBox,
                                                         moleculeElectricPotential, atomPositions);
    Interactions::computeEwaldFourierElectricPotential(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                       system.fixedFrameworkStoredEik, moleculeElectricPotential,
                                                       system.forceField, system.simulationBox, system.components,
                                                       system.numberOfMoleculesPerComponent, atomPositions);
    y2 = moleculeElectricPotential[i];

    atomPositions[i].position.y = spanOfMoleculeAtoms[i].position.y - 0.5 * delta;
    std::fill(moleculeElectricPotential.begin(), moleculeElectricPotential.end(), 0.0);
    Interactions::computeFrameworkMoleculeElectricPotential(
        system.forceField, system.simulationBox, moleculeElectricPotential, frameworkAtomPositions, atomPositions);
    Interactions::computeInterMolecularElectricPotential(system.forceField, system.simulationBox,
                                                         moleculeElectricPotential, atomPositions);
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
    Interactions::computeInterMolecularElectricPotential(system.forceField, system.simulationBox,
                                                         moleculeElectricPotential, atomPositions);
    Interactions::computeEwaldFourierElectricPotential(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                       system.fixedFrameworkStoredEik, moleculeElectricPotential,
                                                       system.forceField, system.simulationBox, system.components,
                                                       system.numberOfMoleculesPerComponent, atomPositions);
    z2 = moleculeElectricPotential[i];

    atomPositions[i].position.z = spanOfMoleculeAtoms[i].position.z - 0.5 * delta;
    std::fill(moleculeElectricPotential.begin(), moleculeElectricPotential.end(), 0.0);
    Interactions::computeFrameworkMoleculeElectricPotential(
        system.forceField, system.simulationBox, moleculeElectricPotential, frameworkAtomPositions, atomPositions);
    Interactions::computeInterMolecularElectricPotential(system.forceField, system.simulationBox,
                                                         moleculeElectricPotential, atomPositions);
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

TEST(electrostatic_field, Test_CO2_in_ITQ_29_2x2x2_difference)
{
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);

  forceField.computePolarization = true;
  forceField.omitEwaldFourier = true;

  Framework f = TestFactories::makeITQ29(forceField, int3(2, 2, 2));
  Component c = TestFactories::makeCO2(forceField, 0, true);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {1}, 5);
  system.forceField.EwaldAlpha = 0.25;
  system.forceField.numberOfWaveVectors = int3(8, 8, 8);
  // system.forceField.omitEwaldFourier = true;

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  spanOfMoleculeAtoms[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[0].scalingCoulomb = 0.75;
  spanOfMoleculeAtoms[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[1].scalingCoulomb = 0.75;
  spanOfMoleculeAtoms[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[2].scalingCoulomb = 0.75;

  system.computeTotalElectricField();

  std::span<Atom> oldatoms = system.spanOfMolecule(0, 0);
  std::vector<Atom> newatoms = std::vector(oldatoms.begin(), oldatoms.end());

  newatoms[0].position = double3(5.93355 + 1.0, 7.93355, 5.93355 + 1.149);
  newatoms[0].scalingCoulomb = 0.65;
  newatoms[1].position = double3(5.93355 + 1.0, 7.93355, 5.93355 + 0.0);
  newatoms[1].scalingCoulomb = 0.65;
  newatoms[2].position = double3(5.93355 + 1.0, 7.93355, 5.93355 - 1.149);
  newatoms[2].scalingCoulomb = 0.65;

  std::span<double3> electricFieldMoleculeNew = system.spanElectricFieldNew(0, 0);

  std::span<const Atom> frameworkAtomPositions = system.spanOfFrameworkAtoms();
  [[maybe_unused]] std::optional<RunningEnergy> runningEnergy =
      Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                    frameworkAtomPositions, electricFieldMoleculeNew,
                                                                    newatoms, oldatoms);

  std::transform(system.electricField.begin(), system.electricField.end(), system.electricFieldNew.begin(),
                 system.electricFieldNew.begin(), std::plus<double3>());

  spanOfMoleculeAtoms[0].position = double3(5.93355 + 1.0, 7.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[0].scalingCoulomb = 0.65;
  spanOfMoleculeAtoms[1].position = double3(5.93355 + 1.0, 7.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[1].scalingCoulomb = 0.65;
  spanOfMoleculeAtoms[2].position = double3(5.93355 + 1.0, 7.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[2].scalingCoulomb = 0.65;

  system.computeTotalElectricField();

  std::span<double3> spanOfElectricField = system.spanOfMoleculeElectricField();
  std::span<double3> spanOfElectricFieldNew = system.spanOfMoleculeElectricFieldNew();
  double tolerance = 1e-5;
  for (size_t i = 0; i < spanOfElectricField.size(); ++i)
  {
    EXPECT_NEAR(spanOfElectricField[i].x, spanOfElectricFieldNew[i].x, tolerance);
    EXPECT_NEAR(spanOfElectricField[i].y, spanOfElectricFieldNew[i].y, tolerance);
    EXPECT_NEAR(spanOfElectricField[i].z, spanOfElectricFieldNew[i].z, tolerance);
  }
}

TEST(electrostatic_field, Test_2_CO2_in_ITQ_29_2x2x2_difference)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Component c = TestFactories::makeCO2(forceField, 0, true);
  Framework f = TestFactories::makeITQ29(forceField, int3(2, 2, 2));

  forceField.computePolarization = true;
  forceField.omitEwaldFourier = true;

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {2}, 5);
  system.forceField.EwaldAlpha = 0.25;
  system.forceField.numberOfWaveVectors = int3(8, 8, 8);
  // system.forceField.omitEwaldFourier = true;

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  spanOfMoleculeAtoms[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[0].scalingCoulomb = 0.9;
  spanOfMoleculeAtoms[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[1].scalingCoulomb = 0.9;
  spanOfMoleculeAtoms[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[2].scalingCoulomb = 0.9;
  spanOfMoleculeAtoms[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[3].scalingCoulomb = 0.45;
  spanOfMoleculeAtoms[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[4].scalingCoulomb = 0.45;
  spanOfMoleculeAtoms[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[5].scalingCoulomb = 0.45;

  system.computeTotalElectricField();

  std::span<Atom> oldatoms = system.spanOfMolecule(0, 0);
  std::vector<Atom> newatoms = std::vector(oldatoms.begin(), oldatoms.end());

  newatoms[0].position = double3(5.93355 + 1.0, 7.93355, 5.93355 + 1.149);
  newatoms[0].scalingCoulomb = 0.75;
  newatoms[1].position = double3(5.93355 + 1.0, 7.93355, 5.93355 + 0.0);
  newatoms[1].scalingCoulomb = 0.75;
  newatoms[2].position = double3(5.93355 + 1.0, 7.93355, 5.93355 - 1.149);
  newatoms[2].scalingCoulomb = 0.75;

  std::span<double3> electricFieldNew = system.spanOfMoleculeElectricFieldNew();
  std::span<double3> electricFieldMoleculeNew = system.spanElectricFieldNew(0, 0);

  std::span<const Atom> frameworkAtomPositions = system.spanOfFrameworkAtoms();
  [[maybe_unused]] std::optional<RunningEnergy> runningEnergyFramework =
      Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                    frameworkAtomPositions, electricFieldMoleculeNew,
                                                                    newatoms, oldatoms);
  [[maybe_unused]] std::optional<RunningEnergy> runningEnergyInter =
      Interactions::computeInterMolecularElectricFieldDifference(system.forceField, system.simulationBox,
                                                                 electricFieldNew, electricFieldMoleculeNew,
                                                                 spanOfMoleculeAtoms, newatoms, oldatoms);

  std::transform(system.electricField.begin(), system.electricField.end(), system.electricFieldNew.begin(),
                 system.electricFieldNew.begin(), std::plus<double3>());

  spanOfMoleculeAtoms[0].position = double3(5.93355 + 1.0, 7.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[0].scalingCoulomb = 0.75;
  spanOfMoleculeAtoms[1].position = double3(5.93355 + 1.0, 7.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[1].scalingCoulomb = 0.75;
  spanOfMoleculeAtoms[2].position = double3(5.93355 + 1.0, 7.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[2].scalingCoulomb = 0.75;
  spanOfMoleculeAtoms[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[3].scalingCoulomb = 0.45;
  spanOfMoleculeAtoms[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[4].scalingCoulomb = 0.45;
  spanOfMoleculeAtoms[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[5].scalingCoulomb = 0.45;

  system.computeTotalElectricField();

  std::span<double3> spanOfElectricField = system.spanOfMoleculeElectricField();
  std::span<double3> spanOfElectricFieldNew = system.spanOfMoleculeElectricFieldNew();
  double tolerance = 1e-5;
  for (size_t i = 0; i < spanOfElectricField.size(); ++i)
  {
    EXPECT_NEAR(spanOfElectricField[i].x, spanOfElectricFieldNew[i].x, tolerance);
    EXPECT_NEAR(spanOfElectricField[i].y, spanOfElectricFieldNew[i].y, tolerance);
    EXPECT_NEAR(spanOfElectricField[i].z, spanOfElectricFieldNew[i].z, tolerance);
  }
}

TEST(electrostatic_field, Test_2_CO2_in_ITQ_29_2x2x2_difference_non_Ewald)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Component c = TestFactories::makeCO2(forceField, 0, true);
  Framework f = TestFactories::makeITQ29(forceField, int3(2, 2, 2));

  forceField.computePolarization = true;
  forceField.omitEwaldFourier = true;

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {2}, 5);
  system.forceField.EwaldAlpha = 0.25;
  system.forceField.numberOfWaveVectors = int3(8, 8, 8);
  // system.forceField.omitEwaldFourier = true;

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  spanOfMoleculeAtoms[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[0].scalingCoulomb = 0.9;
  spanOfMoleculeAtoms[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[1].scalingCoulomb = 0.9;
  spanOfMoleculeAtoms[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[2].scalingCoulomb = 0.9;
  spanOfMoleculeAtoms[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[3].scalingCoulomb = 0.45;
  spanOfMoleculeAtoms[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[4].scalingCoulomb = 0.45;
  spanOfMoleculeAtoms[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[5].scalingCoulomb = 0.45;

  system.computeTotalElectricField();

  std::span<Atom> oldatoms = system.spanOfMolecule(0, 0);
  std::vector<Atom> newatoms = std::vector(oldatoms.begin(), oldatoms.end());

  newatoms[0].position = double3(5.93355 + 1.0, 7.93355, 5.93355 + 1.149);
  newatoms[0].scalingCoulomb = 0.75;
  newatoms[1].position = double3(5.93355 + 1.0, 7.93355, 5.93355 + 0.0);
  newatoms[1].scalingCoulomb = 0.75;
  newatoms[2].position = double3(5.93355 + 1.0, 7.93355, 5.93355 - 1.149);
  newatoms[2].scalingCoulomb = 0.75;

  std::span<double3> electricFieldNew = system.spanOfMoleculeElectricFieldNew();
  std::span<double3> electricFieldMoleculeNew = system.spanElectricFieldNew(0, 0);

  std::span<const Atom> frameworkAtomPositions = system.spanOfFrameworkAtoms();
  [[maybe_unused]] std::optional<RunningEnergy> runningEnergyFramework =
      Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                    frameworkAtomPositions, electricFieldMoleculeNew,
                                                                    newatoms, oldatoms);
  [[maybe_unused]] std::optional<RunningEnergy> runningEnergyInter =
      Interactions::computeInterMolecularElectricFieldDifference(system.forceField, system.simulationBox,
                                                                 electricFieldNew, electricFieldMoleculeNew,
                                                                 spanOfMoleculeAtoms, newatoms, oldatoms);

  // add the differences up to the old electr-field get the new electric-field
  std::transform(system.electricField.begin(), system.electricField.end(), system.electricFieldNew.begin(),
                 system.electricFieldNew.begin(), std::plus<double3>());

  spanOfMoleculeAtoms[0].position = double3(5.93355 + 1.0, 7.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[0].scalingCoulomb = 0.75;
  spanOfMoleculeAtoms[1].position = double3(5.93355 + 1.0, 7.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[1].scalingCoulomb = 0.75;
  spanOfMoleculeAtoms[2].position = double3(5.93355 + 1.0, 7.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[2].scalingCoulomb = 0.75;
  spanOfMoleculeAtoms[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[3].scalingCoulomb = 0.45;
  spanOfMoleculeAtoms[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[4].scalingCoulomb = 0.45;
  spanOfMoleculeAtoms[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[5].scalingCoulomb = 0.45;

  system.computeTotalElectricField();

  std::span<double3> spanOfElectricField = system.spanOfMoleculeElectricField();
  std::span<double3> spanOfElectricFieldNew = system.spanOfMoleculeElectricFieldNew();
  double tolerance = 1e-5;
  for (size_t i = 0; i < spanOfElectricField.size(); ++i)
  {
    EXPECT_NEAR(spanOfElectricField[i].x, spanOfElectricFieldNew[i].x, tolerance);
    EXPECT_NEAR(spanOfElectricField[i].y, spanOfElectricFieldNew[i].y, tolerance);
    EXPECT_NEAR(spanOfElectricField[i].z, spanOfElectricFieldNew[i].z, tolerance);
  }
}

TEST(electrostatic_field, Test_2_CO2_in_ITQ_29_2x2x2_difference_Ewald)
{
  ForceField forceField = TestFactories::makeDefaultFF(12.0, true, false, true);
  Component c = TestFactories::makeCO2(forceField, 0, true);
  Framework f = TestFactories::makeITQ29(forceField, int3(2, 2, 2));

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {2}, 5);
  system.forceField.EwaldAlpha = 0.25;
  system.forceField.numberOfWaveVectors = int3(8, 8, 8);
  system.forceField.omitInterPolarization = true;
  system.forceField.omitInterInteractions = true;

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  spanOfMoleculeAtoms[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[0].scalingCoulomb = 0.75;
  spanOfMoleculeAtoms[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[1].scalingCoulomb = 0.75;
  spanOfMoleculeAtoms[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[2].scalingCoulomb = 0.75;
  spanOfMoleculeAtoms[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[3].scalingCoulomb = 0.45;
  spanOfMoleculeAtoms[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[4].scalingCoulomb = 0.45;
  spanOfMoleculeAtoms[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[5].scalingCoulomb = 0.45;

  if (system.fixedFrameworkStoredEik.empty())
  {
    system.precomputeTotalRigidEnergy();
  }
  [[maybe_unused]] std::span<Atom> frameworkAtomPositions = system.spanOfFrameworkAtoms();
  std::span<Atom> moleculeAtomPositions = system.spanOfMoleculeAtoms();
  std::span<double3> moleculeElectricFieldBefore = system.spanOfMoleculeElectricField();

  std::fill(moleculeElectricFieldBefore.begin(), moleculeElectricFieldBefore.end(), double3(0.0, 0.0, 0.0));
  [[maybe_unused]] RunningEnergy energy_before = Interactions::computeEwaldFourierElectricField(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, moleculeElectricFieldBefore, system.components,
      system.numberOfMoleculesPerComponent, moleculeAtomPositions);

  // std::cout << "before: " << energy_before.ewald << std::endl;

  std::span<Atom> oldatoms = system.spanOfMolecule(0, 0);
  std::vector<Atom> newatoms = std::vector(oldatoms.begin(), oldatoms.end());

  newatoms[0].position = double3(5.93355 + 1.0, 7.93355, 5.93355 + 1.149);
  newatoms[0].scalingCoulomb = 0.75;
  newatoms[1].position = double3(5.93355 + 1.0, 7.93355, 5.93355 + 0.0);
  newatoms[1].scalingCoulomb = 0.75;
  newatoms[2].position = double3(5.93355 + 1.0, 7.93355, 5.93355 - 1.149);
  newatoms[2].scalingCoulomb = 0.75;

  std::vector<double3> moleculeElectricFieldDiff(moleculeElectricFieldBefore.size());
  std::fill(moleculeElectricFieldDiff.begin(), moleculeElectricFieldDiff.end(), double3(0.0, 0.0, 0.0));
  [[maybe_unused]] RunningEnergy energy_difference = Interactions::eletricFieldDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.totalEik, system.forceField, system.simulationBox, moleculeElectricFieldDiff, newatoms, oldatoms);

  spanOfMoleculeAtoms[0].position = double3(5.93355 + 1.0, 7.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[0].scalingCoulomb = 0.75;
  spanOfMoleculeAtoms[1].position = double3(5.93355 + 1.0, 7.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[1].scalingCoulomb = 0.75;
  spanOfMoleculeAtoms[2].position = double3(5.93355 + 1.0, 7.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[2].scalingCoulomb = 0.75;
  spanOfMoleculeAtoms[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[3].scalingCoulomb = 0.45;
  spanOfMoleculeAtoms[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[4].scalingCoulomb = 0.45;
  spanOfMoleculeAtoms[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[5].scalingCoulomb = 0.45;

  std::span<double3> moleculeElectricFieldAfter = system.spanOfMoleculeElectricFieldNew();
  std::fill(moleculeElectricFieldAfter.begin(), moleculeElectricFieldAfter.end(), double3(0.0, 0.0, 0.0));
  [[maybe_unused]] RunningEnergy energy_after = Interactions::computeEwaldFourierElectricField(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, moleculeElectricFieldAfter, system.components,
      system.numberOfMoleculesPerComponent, moleculeAtomPositions);
  // std::cout << "after: " << energy_after.ewald << std::endl;

  // std::cout << "Running energy difference: " << energy_difference.ewald << std::endl;
  // std::cout << "energy before - after: " << energy_after.ewald - energy_before.ewald << std::endl;

  // double tolerance = 1e-5;
  // for(size_t i = 0; i < moleculeElectricFieldAfter.size(); ++i)
  //{
  //   std::cout << "i: " << i << " " << moleculeElectricFieldBefore[i].x - moleculeElectricFieldAfter[i].x << " " <<
  //   moleculeElectricFieldDiff[i].x << std::endl; std::cout << "i: " << i << " " << moleculeElectricFieldBefore[i].y -
  //   moleculeElectricFieldAfter[i].y << " " << moleculeElectricFieldDiff[i].y << std::endl; std::cout << "i: " << i <<
  //   " " << moleculeElectricFieldBefore[i].z - moleculeElectricFieldAfter[i].z << " " <<
  //   moleculeElectricFieldDiff[i].z << std::endl;
  // }
}
