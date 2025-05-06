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

// Table 5, page 61 thesis D. Dubbeldam
TEST(electrostatic_potential, Test_reference_system_1)
{
  ForceField forceField = ForceField(
      {
          PseudoAtom("t1", false, 1.0, 0.5, 0.0, 1, false),
          PseudoAtom("t2", false, 1.0, 1.5, 0.0, 1, false),
          PseudoAtom("t2", false, 1.0, -0.75, 0.0, 1, false),
          PseudoAtom("t4", false, 1.0, -1.25, 0.0, 1, false),
      },
      {VDWParameters(0.0, 1.0)}, ForceField::MixingRule::Lorentz_Berthelot, 500.0, 500.0, 500.0, true, false, true);

  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
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

  System system = System(0, forceField, SimulationBox(1000.0, 1000.0, 1000.0), 300.0, 1e4, 1.0, {}, {c1, c2, c3, c4},
                         {1, 1, 1, 1}, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtomPositions = system.spanOfFrameworkAtoms();
  std::span<double> moleculeElectricPotential = system.spanOfMoleculeElectricPotential();
  spanOfMoleculeAtoms[0].position = double3(-1.0, 0.0, 0.0);
  spanOfMoleculeAtoms[1].position = double3(1.0, 0.0, 0.0);
  spanOfMoleculeAtoms[2].position = double3(0.0, 1.0, 0.0);
  spanOfMoleculeAtoms[3].position = double3(0.0, -1.0, 0.0);

  if (system.fixedFrameworkStoredEik.empty())
  {
    system.precomputeTotalRigidEnergy();
  }

  RunningEnergy energy =
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, spanOfMoleculeAtoms) +
      Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                   system.framework, frameworkAtomPositions, spanOfMoleculeAtoms) +
      Interactions::computeEwaldFourierEnergy(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                              system.fixedFrameworkStoredEik, system.storedEik, system.forceField,
                                              system.simulationBox, system.components,
                                              system.numberOfMoleculesPerComponent, spanOfMoleculeAtoms);

  system.computeTotalElectricPotential();

  // U = \sum_i q_i \phi_i
  double potentialEnergy{};
  for (size_t i = 0; i < 6; ++i)
  {
    potentialEnergy +=
        0.5 * spanOfMoleculeAtoms[i].scalingCoulomb * spanOfMoleculeAtoms[i].charge * moleculeElectricPotential[i];
  }

  EXPECT_NEAR(moleculeElectricPotential[0] / Units::CoulombicConversionFactor, -0.664214, 1e-5);
  EXPECT_NEAR(moleculeElectricPotential[1] / Units::CoulombicConversionFactor, -1.16421, 1e-5);
  EXPECT_NEAR(moleculeElectricPotential[2] / Units::CoulombicConversionFactor, 0.789214, 1e-5);
  EXPECT_NEAR(moleculeElectricPotential[3] / Units::CoulombicConversionFactor, 1.03921, 1e-5);

  EXPECT_NEAR(potentialEnergy * Units::EnergyToKelvin,
              (energy.frameworkMoleculeCharge + energy.moleculeMoleculeCharge + energy.ewald_fourier +
               energy.ewald_self + energy.ewald_exclusion) *
                  Units::EnergyToKelvin,
              1e-5);
  EXPECT_NEAR((energy.frameworkMoleculeCharge + energy.moleculeMoleculeCharge + energy.ewald_fourier +
               energy.ewald_self + energy.ewald_exclusion) /
                  Units::CoulombicConversionFactor,
              -1.98467, 1e-5);
}

// Table 5, page 61 thesis D. Dubbeldam
TEST(electrostatic_potential, Test_reference_system_2)
{
  ForceField forceField = ForceField(
      {
          PseudoAtom("t1", false, 1.0, 0.5, 0.0, 1, false),
          PseudoAtom("t2", false, 1.0, -1.25, 0.0, 1, false),
          PseudoAtom("t2", false, 1.0, 1.5, 0.0, 1, false),
          PseudoAtom("t4", false, 1.0, -0.75, 0.0, 1, false),
      },
      {VDWParameters(0.0, 1.0)}, ForceField::MixingRule::Lorentz_Berthelot, 500.0, 500.0, 500.0, true, false, true);

  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
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

  System system =
      System(0, forceField, SimulationBox(1000.0, 1000.0, 1000.0), 300.0, 1e4, 1.0, {}, {c1, c2}, {1, 1}, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtomPositions = system.spanOfFrameworkAtoms();
  std::span<double> moleculeElectricPotential = system.spanOfMoleculeElectricPotential();
  spanOfMoleculeAtoms[0].position = double3(-1.0, 0.0, 0.0);
  spanOfMoleculeAtoms[1].position = double3(0.0, -1.0, 0.0);
  spanOfMoleculeAtoms[2].position = double3(1.0, 0.0, 0.0);
  spanOfMoleculeAtoms[3].position = double3(0.0, 1.0, 0.0);

  if (system.fixedFrameworkStoredEik.empty())
  {
    system.precomputeTotalRigidEnergy();
  }

  RunningEnergy energy =
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, spanOfMoleculeAtoms) +
      Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                   system.framework, frameworkAtomPositions, spanOfMoleculeAtoms) +
      Interactions::computeEwaldFourierEnergy(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                              system.fixedFrameworkStoredEik, system.storedEik, system.forceField,
                                              system.simulationBox, system.components,
                                              system.numberOfMoleculesPerComponent, spanOfMoleculeAtoms);

  system.computeTotalElectricPotential();

  // U = \sum_i q_i \phi_i
  double potentialEnergy{};
  for (size_t i = 0; i < 6; ++i)
  {
    potentialEnergy +=
        0.5 * spanOfMoleculeAtoms[i].scalingCoulomb * spanOfMoleculeAtoms[i].charge * moleculeElectricPotential[i];
  }

  EXPECT_NEAR(moleculeElectricPotential[0] / Units::CoulombicConversionFactor, 0.21967, 1e-5);
  EXPECT_NEAR(moleculeElectricPotential[1] / Units::CoulombicConversionFactor, 0.68566, 1e-5);
  EXPECT_NEAR(moleculeElectricPotential[2] / Units::CoulombicConversionFactor, -0.633883, 1e-5);
  EXPECT_NEAR(moleculeElectricPotential[3] / Units::CoulombicConversionFactor, -0.271447, 1e-5);

  EXPECT_NEAR(potentialEnergy * Units::EnergyToKelvin,
              (energy.frameworkMoleculeCharge + energy.moleculeMoleculeCharge + energy.ewald_fourier +
               energy.ewald_self + energy.ewald_exclusion) *
                  Units::EnergyToKelvin,
              1e-4);

  // Check again, when net charges are implemented
  EXPECT_NEAR((energy.frameworkMoleculeCharge + energy.moleculeMoleculeCharge + energy.ewald_fourier +
               energy.ewald_self + energy.ewald_exclusion) /
                  Units::CoulombicConversionFactor,
              -0.745882, 1e-2);
}

TEST(electrostatic_potential, Test_2_CO2_in_ITQ_29_2x2x2)
{
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);

  forceField.computePolarization = true;
  forceField.omitInterPolarization = true;
  forceField.omitInterInteractions = true;

  Framework f = TestFactories::makeITQ29(forceField, int3(2, 2, 2));
  Component c = TestFactories::makeCO2(forceField, 0, true);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {2}, 5);
  system.forceField.EwaldAlpha = 0.25;
  system.forceField.numberOfWaveVectors = int3(8, 8, 8);
  // system.forceField.omitEwaldFourier = true;

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtomPositions = system.spanOfFrameworkAtoms();
  std::span<double> moleculeElectricPotential = system.spanOfMoleculeElectricPotential();
  spanOfMoleculeAtoms[0].position = double3(5.93355, 7.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[0].scalingCoulomb = 0.45;
  spanOfMoleculeAtoms[1].position = double3(5.93355, 7.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[1].scalingCoulomb = 0.45;
  spanOfMoleculeAtoms[2].position = double3(5.93355, 7.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[2].scalingCoulomb = 0.45;

  spanOfMoleculeAtoms[3].position = double3(5.93355, 3.93355, 5.93355 + 1.149);
  spanOfMoleculeAtoms[3].scalingCoulomb = 0.65;
  spanOfMoleculeAtoms[4].position = double3(5.93355, 3.93355, 5.93355 + 0.0);
  spanOfMoleculeAtoms[4].scalingCoulomb = 0.65;
  spanOfMoleculeAtoms[5].position = double3(5.93355, 3.93355, 5.93355 - 1.149);
  spanOfMoleculeAtoms[5].scalingCoulomb = 0.65;

  if (system.fixedFrameworkStoredEik.empty())
  {
    system.precomputeTotalRigidEnergy();
  }

  RunningEnergy energy =
      Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                   system.framework, frameworkAtomPositions, spanOfMoleculeAtoms) +
      Interactions::computeEwaldFourierEnergy(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                              system.fixedFrameworkStoredEik, system.storedEik, system.forceField,
                                              system.simulationBox, system.components,
                                              system.numberOfMoleculesPerComponent, spanOfMoleculeAtoms);

  system.computeTotalElectricPotential();

  // U = \sum_i q_i \phi_i
  double potentialEnergy{};
  for (size_t i = 0; i < 6; ++i)
  {
    potentialEnergy +=
        0.5 * spanOfMoleculeAtoms[i].scalingCoulomb * spanOfMoleculeAtoms[i].charge * moleculeElectricPotential[i];
  }

  EXPECT_NEAR(potentialEnergy * Units::EnergyToKelvin,
              (energy.frameworkMoleculeCharge + energy.moleculeMoleculeCharge + energy.ewald_fourier +
               energy.ewald_self + energy.ewald_exclusion) *
                  Units::EnergyToKelvin,
              1e-5);
}
