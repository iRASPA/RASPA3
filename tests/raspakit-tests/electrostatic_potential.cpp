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
TEST(electrostatic_potential, Test_reference_system_1_four_ions)
{
  ForceField forceField = ForceField(
      {
          // std::string name, bool framework, double mass, double charge, double polarizability, size_t atomicNumber,
          // bool printToPDB
          PseudoAtom("t1", false, 1.0, 0.5, 0.0, 1, false),
          PseudoAtom("t2", false, 1.0, 1.5, 0.0, 1, false),
          PseudoAtom("t2", false, 1.0, -0.75, 0.0, 1, false),
          PseudoAtom("t4", false, 1.0, -1.25, 0.0, 1, false),
      },
      {VDWParameters(0.0, 1.0), VDWParameters(0.0, 1.0), VDWParameters(0.0, 1.0), VDWParameters(0.0, 1.0)},
      ForceField::MixingRule::Lorentz_Berthelot, 500.0, 500.0, 500.0, true, false, true);

  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;

  Component c1 = Component(0, forceField, "t1", 0.0, 0.0, 0.0,
                           {
                               // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                               // uint8_t componentId, bool groupId, bool is_fractional
                               Atom(double3(0.0, 0.0, 0.0), 0.5, 1.0, 0, 0, 0, false, false),
                           },
                           {}, {}, 5, 21);
  Component c2 = Component(1, forceField, "t2", 0.0, 0.0, 0.0,
                           {
                               // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                               // uint8_t componentId, bool groupId, bool is_fractional
                               Atom(double3(0.0, 0.0, 0.0), 1.5, 1.0, 1, 1, 1, false, false),
                           },
                           {}, {}, 5, 21);
  Component c3 = Component(2, forceField, "t3", 0.0, 0.0, 0.0,
                           {
                               // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                               // uint8_t componentId, bool groupId, bool is_fractional
                               Atom(double3(0.0, 0.0, 0.0), -0.75, 1.0, 2, 2, 2, false, false),
                           },
                           {}, {}, 5, 21);
  Component c4 = Component(3, forceField, "t4", 0.0, 0.0, 0.0,
                           {
                               // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                               // uint8_t componentId, bool groupId, bool is_fractional
                               Atom(double3(0.0, 0.0, 0.0), -1.25, 1.0, 3, 3, 3, false, false),
                           },
                           {}, {}, 5, 21);

  System system = System(0, forceField, SimulationBox(1000.0, 1000.0, 1000.0), 300.0, 1e4, 1.0, {}, {c1, c2, c3, c4},
                         {}, {1, 1, 1, 1}, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtomPositions = system.spanOfFrameworkAtoms();
  std::span<double> moleculeElectrostaticPotential = system.spanOfMoleculeElectrostaticPotential();
  spanOfMoleculeAtoms[0].position = double3(-1.0, 0.0, 0.0);
  spanOfMoleculeAtoms[1].position = double3(1.0, 0.0, 0.0);
  spanOfMoleculeAtoms[2].position = double3(0.0, 1.0, 0.0);
  spanOfMoleculeAtoms[3].position = double3(0.0, -1.0, 0.0);

  // Initialize gradients to zero
  for (Atom& atom : spanOfMoleculeAtoms)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

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

  system.computeTotalElectrostaticPotential();

  EXPECT_NEAR(moleculeElectrostaticPotential[0] / Units::CoulombicConversionFactor, -0.664214, 1e-5);
  EXPECT_NEAR(moleculeElectrostaticPotential[1] / Units::CoulombicConversionFactor, -1.16421, 1e-5);
  EXPECT_NEAR(moleculeElectrostaticPotential[2] / Units::CoulombicConversionFactor, 0.789214, 1e-5);
  EXPECT_NEAR(moleculeElectrostaticPotential[3] / Units::CoulombicConversionFactor, 1.03921, 1e-5);

  EXPECT_NEAR((energy.frameworkMoleculeCharge + energy.moleculeMoleculeCharge + energy.ewald_fourier +
               energy.ewald_self + energy.ewald_exclusion) /
                  Units::CoulombicConversionFactor,
              -1.98467, 1e-5);

  // U = 0.5 \sum_i q_i \phi_i
  double potentialEnergy{};
  for (size_t i = 0; i < 4; ++i)
  {
    potentialEnergy += 0.5 * spanOfMoleculeAtoms[i].charge * moleculeElectrostaticPotential[i];
  }

  EXPECT_NEAR(potentialEnergy * Units::EnergyToKelvin,
              (energy.frameworkMoleculeCharge + energy.moleculeMoleculeCharge + energy.ewald_fourier +
               energy.ewald_self + energy.ewald_exclusion) *
                  Units::EnergyToKelvin,
              1e-5);

  std::span<double3> moleculeElectricField = system.spanOfMoleculeElectricField();
  std::fill(moleculeElectricField.begin(), moleculeElectricField.end(), double3(0.0, 0.0, 0.0));

  system.precomputeTotalRigidEnergy();

  [[maybe_unused]] RunningEnergy energy1 = Interactions::computeFrameworkMoleculeElectricField(
      system.forceField, system.simulationBox, moleculeElectricField, system.spanOfFrameworkAtoms(),
      system.spanOfMoleculeAtoms());

  [[maybe_unused]] RunningEnergy energy2 = Interactions::computeInterMolecularElectricField(
      system.forceField, system.simulationBox, moleculeElectricField, system.spanOfMoleculeAtoms());

  [[maybe_unused]] RunningEnergy energy3 = Interactions::computeEwaldFourierElectricField(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.totalEik,
      system.forceField, system.simulationBox, moleculeElectricField, system.components,
      system.numberOfMoleculesPerComponent, system.spanOfMoleculeAtoms());

  EXPECT_NEAR(moleculeElectricField[0].x / Units::CoulombicConversionFactor, 0.332107, 1e-5);
  EXPECT_NEAR(moleculeElectricField[0].y / Units::CoulombicConversionFactor, -0.176777, 1e-5);
  EXPECT_NEAR(moleculeElectricField[0].z / Units::CoulombicConversionFactor, 0.0, 1e-5);

  EXPECT_NEAR(moleculeElectricField[1].x / Units::CoulombicConversionFactor, -0.582107, 1e-5);
  EXPECT_NEAR(moleculeElectricField[1].y / Units::CoulombicConversionFactor, -0.176777, 1e-5);
  EXPECT_NEAR(moleculeElectricField[1].z / Units::CoulombicConversionFactor, 0.0, 1e-5);

  EXPECT_NEAR(moleculeElectricField[2].x / Units::CoulombicConversionFactor, -0.353553, 1e-5);
  EXPECT_NEAR(moleculeElectricField[2].y / Units::CoulombicConversionFactor, 0.394607, 1e-5);
  EXPECT_NEAR(moleculeElectricField[2].z / Units::CoulombicConversionFactor, 0.0, 1e-5);

  EXPECT_NEAR(moleculeElectricField[3].x / Units::CoulombicConversionFactor, -0.353553, 1e-5);
  EXPECT_NEAR(moleculeElectricField[3].y / Units::CoulombicConversionFactor, -0.519607, 1e-5);
  EXPECT_NEAR(moleculeElectricField[3].z / Units::CoulombicConversionFactor, 0.0, 1e-5);

  [[maybe_unused]] RunningEnergy energy4 = Interactions::computeFrameworkMoleculeGradient(
      system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
      system.interpolationGrids);

  [[maybe_unused]] RunningEnergy energy5 = Interactions::computeInterMolecularGradient(
      system.forceField, system.simulationBox, system.spanOfMoleculeAtoms());

  [[maybe_unused]] RunningEnergy energy6 = Interactions::computeEwaldFourierGradient(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.fixedFrameworkStoredEik,
      system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
      system.spanOfMoleculeAtoms());

  EXPECT_NEAR(spanOfMoleculeAtoms[0].gradient.x / Units::CoulombicConversionFactor, -0.166053, 1e-5);
  EXPECT_NEAR(spanOfMoleculeAtoms[0].gradient.y / Units::CoulombicConversionFactor, 0.0883883, 1e-5);
  EXPECT_NEAR(spanOfMoleculeAtoms[0].gradient.z / Units::CoulombicConversionFactor, 0.0, 1e-5);

  EXPECT_NEAR(spanOfMoleculeAtoms[1].gradient.x / Units::CoulombicConversionFactor, 0.87316, 1e-5);
  EXPECT_NEAR(spanOfMoleculeAtoms[1].gradient.y / Units::CoulombicConversionFactor, 0.265165, 1e-5);
  EXPECT_NEAR(spanOfMoleculeAtoms[1].gradient.z / Units::CoulombicConversionFactor, 0.0, 1e-5);

  EXPECT_NEAR(spanOfMoleculeAtoms[2].gradient.x / Units::CoulombicConversionFactor, -0.265165, 1e-5);
  EXPECT_NEAR(spanOfMoleculeAtoms[2].gradient.y / Units::CoulombicConversionFactor, 0.295955, 1e-5);
  EXPECT_NEAR(spanOfMoleculeAtoms[2].gradient.z / Units::CoulombicConversionFactor, 0.0, 1e-5);

  EXPECT_NEAR(spanOfMoleculeAtoms[3].gradient.x / Units::CoulombicConversionFactor, -0.441942, 1e-5);
  EXPECT_NEAR(spanOfMoleculeAtoms[3].gradient.y / Units::CoulombicConversionFactor, -0.649508, 1e-5);
  EXPECT_NEAR(spanOfMoleculeAtoms[3].gradient.z / Units::CoulombicConversionFactor, 0.0, 1e-5);
}

// Table 5, page 61 thesis D. Dubbeldam
TEST(electrostatic_potential, Test_reference_system_1_framework_molecule)
{
  ForceField forceField = ForceField(
      {
          // std::string name, bool framework, double mass, double charge, double polarizability, size_t atomicNumber,
          // bool printToPDB
          PseudoAtom("t1", false, 1.0, 0.5, 0.0, 1, false),
          PseudoAtom("t2", false, 1.0, 1.5, 0.0, 1, false),
          PseudoAtom("t2", false, 1.0, -0.75, 0.0, 1, false),
          PseudoAtom("t4", false, 1.0, -1.25, 0.0, 1, false),
      },
      {VDWParameters(0.0, 1.0), VDWParameters(0.0, 1.0), VDWParameters(0.0, 1.0), VDWParameters(0.0, 1.0)},
      ForceField::MixingRule::Lorentz_Berthelot, 500.0, 500.0, 500.0, true, false, true);

  forceField.computePolarization = true;
  forceField.omitInterPolarization = false;
  forceField.omitInterInteractions = false;

  Framework framework = Framework(0, forceField, "ions", SimulationBox(1000.0, 1000.0, 1000.0), 1,
                                  {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t
                                   // type, uint8_t componentId, bool groupId, bool is_fractional
                                   Atom({-1.0 / 1000.0, 0.0, 0.0}, 0.5, 1.0, 0, 0, 0, false, false),
                                   Atom({1.0 / 1000.0, 0.0, 0.0}, 1.5, 1.0, 0, 1, 0, false, false),
                                   Atom({0.0, 1.0 / 1000.0, 0.0}, -0.75, 1.0, 0, 2, 0, false, false)},
                                  {1, 1, 1});

  Component c4 = Component(0, forceField, "t4", 0.0, 0.0, 0.0,
                           {
                               // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                               // uint8_t componentId, bool groupId, bool is_fractional
                               Atom(double3(0.0, 0.0, 0.0), -1.25, 1.0, 0, 3, 0, false, false),
                           },
                           {}, {}, 5, 21);

  System system = System(0, forceField, SimulationBox(1000.0, 1000.0, 1000.0), 300.0, 1e4, 1.0, {framework}, {c4}, {},
                         {1, 1, 1, 1}, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::span<Atom> frameworkAtomPositions = system.spanOfFrameworkAtoms();
  frameworkAtomPositions[0].position = double3(-1.0, 0.0, 0.0);
  frameworkAtomPositions[1].position = double3(1.0, 0.0, 0.0);
  frameworkAtomPositions[2].position = double3(0.0, 1.0, 0.0);
  std::span<double> moleculeElectrostaticPotential = system.spanOfMoleculeElectrostaticPotential();
  spanOfMoleculeAtoms[0].position = double3(0.0, -1.0, 0.0);

  // Initialize gradients to zero
  for (Atom& atom : spanOfMoleculeAtoms)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  if (system.fixedFrameworkStoredEik.empty())
  {
    system.precomputeTotalRigidEnergy();
  }

  [[maybe_unused]] RunningEnergy energy =
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, spanOfMoleculeAtoms) +
      Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                   system.framework, frameworkAtomPositions, spanOfMoleculeAtoms) +
      Interactions::computeEwaldFourierEnergy(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                              system.fixedFrameworkStoredEik, system.storedEik, system.forceField,
                                              system.simulationBox, system.components,
                                              system.numberOfMoleculesPerComponent, spanOfMoleculeAtoms);

  system.computeTotalElectrostaticPotential();

  EXPECT_NEAR(moleculeElectrostaticPotential[0] / Units::CoulombicConversionFactor, 1.03921, 1e-5);

  std::span<double3> moleculeElectricField = system.spanOfMoleculeElectricField();
  std::fill(moleculeElectricField.begin(), moleculeElectricField.end(), double3(0.0, 0.0, 0.0));

  system.precomputeTotalRigidEnergy();

  [[maybe_unused]] RunningEnergy energy1 = Interactions::computeFrameworkMoleculeElectricField(
      system.forceField, system.simulationBox, moleculeElectricField, system.spanOfFrameworkAtoms(),
      system.spanOfMoleculeAtoms());

  [[maybe_unused]] RunningEnergy energy2 = Interactions::computeInterMolecularElectricField(
      system.forceField, system.simulationBox, moleculeElectricField, system.spanOfMoleculeAtoms());

  [[maybe_unused]] RunningEnergy energy3 = Interactions::computeEwaldFourierElectricField(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.totalEik,
      system.forceField, system.simulationBox, moleculeElectricField, system.components,
      system.numberOfMoleculesPerComponent, system.spanOfMoleculeAtoms());

  EXPECT_NEAR(moleculeElectricField[0].x / Units::CoulombicConversionFactor, -0.353553, 1e-5);
  EXPECT_NEAR(moleculeElectricField[0].y / Units::CoulombicConversionFactor, -0.519607, 1e-5);
  EXPECT_NEAR(moleculeElectricField[0].z / Units::CoulombicConversionFactor, 0.0, 1e-5);

  [[maybe_unused]] RunningEnergy energy4 = Interactions::computeFrameworkMoleculeGradient(
      system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
      system.interpolationGrids);

  [[maybe_unused]] RunningEnergy energy5 = Interactions::computeInterMolecularGradient(
      system.forceField, system.simulationBox, system.spanOfMoleculeAtoms());

  [[maybe_unused]] RunningEnergy energy6 = Interactions::computeEwaldFourierGradient(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.fixedFrameworkStoredEik,
      system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
      system.spanOfMoleculeAtoms());

  EXPECT_NEAR(spanOfMoleculeAtoms[0].gradient.x / Units::CoulombicConversionFactor, -0.441942, 1e-5);
  EXPECT_NEAR(spanOfMoleculeAtoms[0].gradient.y / Units::CoulombicConversionFactor, -0.649508, 1e-5);
  EXPECT_NEAR(spanOfMoleculeAtoms[0].gradient.z / Units::CoulombicConversionFactor, 0.0, 1e-5);
}

// Table 5, page 61 thesis D. Dubbeldam
TEST(electrostatic_potential, Test_reference_system_2_framework_molecule)
{
  ForceField forceField = ForceField(
      {
          // std::string name, bool framework, double mass, double charge, double polarizability, size_t atomicNumber,
          // bool printToPDB
          PseudoAtom("t1", false, 1.0, 0.5, 0.0, 1, false),
          PseudoAtom("t2", false, 1.0, -1.25, 0.0, 1, false),
          PseudoAtom("t2", false, 1.0, 1.5, 0.0, 1, false),
          PseudoAtom("t4", false, 1.0, -0.75, 0.0, 1, false),
      },
      {VDWParameters(0.0, 1.0), VDWParameters(0.0, 1.0), VDWParameters(0.0, 1.0), VDWParameters(0.0, 1.0)},
      ForceField::MixingRule::Lorentz_Berthelot, 500.0, 500.0, 500.0, true, false, true);

  forceField.computePolarization = false;
  forceField.omitInterPolarization = true;
  forceField.omitInterInteractions = true;

  Framework framework = Framework(0, forceField, "ions", SimulationBox(1000.0, 1000.0, 1000.0), 1,
                                  {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t
                                   // type, uint8_t componentId, bool groupId, bool is_fractional
                                   Atom({-1.0 / 1000.0, 0.0, 0.0}, 0.5, 1.0, 0, 0, 0, false, false),
                                   Atom({0.0, -1.0 / 1000.0, 0.0}, -1.25, 1.0, 0, 1, 0, false, false)},
                                  {1, 1, 1});

  Component c3 = Component(0, forceField, "t3", 0.0, 0.0, 0.0,
                           {
                               // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                               // uint8_t componentId, bool groupId, bool is_fractional
                               Atom(double3(0.0, 0.0, 0.0), 1.5, 1.0, 0, 2, 0, false, false),
                           },
                           {}, {}, 5, 21);

  Component c4 = Component(1, forceField, "t4", 0.0, 0.0, 0.0,
                           {
                               // double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                               // uint8_t componentId, bool groupId, bool is_fractional
                               Atom(double3(0.0, 0.0, 0.0), -0.75, 1.0, 0, 3, 0, false, false),
                           },
                           {}, {}, 5, 21);

  System system = System(0, forceField, SimulationBox(1000.0, 1000.0, 1000.0), 300.0, 1e4, 1.0, {framework}, {c3, c4},
                         {}, {1, 1, 1, 1}, 5);

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::span<Atom> frameworkAtomPositions = system.spanOfFrameworkAtoms();
  frameworkAtomPositions[0].position = double3(-1.0, 0.0, 0.0);
  frameworkAtomPositions[1].position = double3(0.0, -1.0, 0.0);
  std::span<double> moleculeElectrostaticPotential = system.spanOfMoleculeElectrostaticPotential();
  spanOfMoleculeAtoms[0].position = double3(1.0, 0.0, 0.0);
  spanOfMoleculeAtoms[1].position = double3(0.0, 1.0, 0.0);

  // Initialize gradients to zero
  for (Atom& atom : spanOfMoleculeAtoms)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  if (system.fixedFrameworkStoredEik.empty())
  {
    system.precomputeTotalRigidEnergy();
  }

  [[maybe_unused]] RunningEnergy energy =
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, spanOfMoleculeAtoms) +
      Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                   system.framework, frameworkAtomPositions, spanOfMoleculeAtoms) +
      Interactions::computeEwaldFourierEnergy(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                              system.fixedFrameworkStoredEik, system.storedEik, system.forceField,
                                              system.simulationBox, system.components,
                                              system.numberOfMoleculesPerComponent, spanOfMoleculeAtoms);

  EXPECT_NEAR((energy.frameworkMoleculeCharge + energy.moleculeMoleculeCharge + energy.ewald_fourier +
               energy.ewald_self + energy.ewald_exclusion) /
                  Units::CoulombicConversionFactor,
              -0.745882, 1e-3);

  system.computeTotalElectrostaticPotential();

  EXPECT_NEAR(moleculeElectrostaticPotential[0] / Units::CoulombicConversionFactor, -0.633883, 1e-2);
  EXPECT_NEAR(moleculeElectrostaticPotential[1] / Units::CoulombicConversionFactor, -0.271447, 1e-2);

  std::span<double3> moleculeElectricField = system.spanOfMoleculeElectricField();
  std::fill(moleculeElectricField.begin(), moleculeElectricField.end(), double3(0.0, 0.0, 0.0));

  system.precomputeTotalRigidEnergy();

  [[maybe_unused]] RunningEnergy energy1 = Interactions::computeFrameworkMoleculeElectricField(
      system.forceField, system.simulationBox, moleculeElectricField, system.spanOfFrameworkAtoms(),
      system.spanOfMoleculeAtoms());

  [[maybe_unused]] RunningEnergy energy2 = Interactions::computeInterMolecularElectricField(
      system.forceField, system.simulationBox, moleculeElectricField, system.spanOfMoleculeAtoms());

  [[maybe_unused]] RunningEnergy energy3 = Interactions::computeEwaldFourierElectricField(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.totalEik,
      system.forceField, system.simulationBox, moleculeElectricField, system.components,
      system.numberOfMoleculesPerComponent, system.spanOfMoleculeAtoms());

  EXPECT_NEAR(moleculeElectricField[0].x / Units::CoulombicConversionFactor, -0.316942, 1e-5);
  EXPECT_NEAR(moleculeElectricField[0].y / Units::CoulombicConversionFactor, -0.441942, 1e-5);
  EXPECT_NEAR(moleculeElectricField[0].z / Units::CoulombicConversionFactor, 0.0, 1e-5);

  EXPECT_NEAR(moleculeElectricField[1].x / Units::CoulombicConversionFactor, 0.176777, 1e-5);
  EXPECT_NEAR(moleculeElectricField[1].y / Units::CoulombicConversionFactor, -0.135723, 1e-5);
  EXPECT_NEAR(moleculeElectricField[1].z / Units::CoulombicConversionFactor, 0.0, 1e-5);

  [[maybe_unused]] RunningEnergy energy4 = Interactions::computeFrameworkMoleculeGradient(
      system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
      system.interpolationGrids);

  [[maybe_unused]] RunningEnergy energy5 = Interactions::computeInterMolecularGradient(
      system.forceField, system.simulationBox, system.spanOfMoleculeAtoms());

  [[maybe_unused]] RunningEnergy energy6 = Interactions::computeEwaldFourierGradient(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.fixedFrameworkStoredEik,
      system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
      system.spanOfMoleculeAtoms());

  EXPECT_NEAR(spanOfMoleculeAtoms[0].gradient.x / Units::CoulombicConversionFactor, 0.475413, 1e-5);
  EXPECT_NEAR(spanOfMoleculeAtoms[0].gradient.y / Units::CoulombicConversionFactor, 0.662913, 1e-5);
  EXPECT_NEAR(spanOfMoleculeAtoms[0].gradient.z / Units::CoulombicConversionFactor, 0.0, 1e-5);

  EXPECT_NEAR(spanOfMoleculeAtoms[1].gradient.x / Units::CoulombicConversionFactor, 0.132583, 1e-5);
  EXPECT_NEAR(spanOfMoleculeAtoms[1].gradient.y / Units::CoulombicConversionFactor, -0.101792, 1e-5);
  EXPECT_NEAR(spanOfMoleculeAtoms[1].gradient.z / Units::CoulombicConversionFactor, 0.0, 1e-5);
}

TEST(electrostatic_potential, Test_2_CO2_in_ITQ_29_2x2x2)
{
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);

  forceField.computePolarization = true;
  forceField.omitInterPolarization = true;
  forceField.omitInterInteractions = true;

  Framework f = TestFactories::makeITQ29(forceField, int3(2, 2, 2));
  Component c = TestFactories::makeCO2(forceField, 0, true);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {}, {1}, 5);
  system.forceField.EwaldAlpha = 0.25;
  system.forceField.numberOfWaveVectors = int3(8, 8, 8);
  system.forceField.omitEwaldFourier = true;

  std::span<Atom> spanOfMoleculeAtoms = system.spanOfMoleculeAtoms();
  std::span<const Atom> frameworkAtomPositions = system.spanOfFrameworkAtoms();
  std::span<double> moleculeElectrostaticPotential = system.spanOfMoleculeElectrostaticPotential();
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

  system.computeTotalElectrostaticPotential();

  // U = \sum_i q_i \phi_i
  // (Note: factor of 0.5 missing, since the framework electric potential is ignored)
  double potentialEnergy{};
  for (size_t i = 0; i < 3; ++i)
  {
    potentialEnergy +=
        spanOfMoleculeAtoms[i].scalingCoulomb * spanOfMoleculeAtoms[i].charge * moleculeElectrostaticPotential[i];
  }

  RunningEnergy energy =
      Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                   system.framework, frameworkAtomPositions, spanOfMoleculeAtoms) +
      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, spanOfMoleculeAtoms) +
      Interactions::computeEwaldFourierEnergy(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                              system.fixedFrameworkStoredEik, system.storedEik, system.forceField,
                                              system.simulationBox, system.components,
                                              system.numberOfMoleculesPerComponent, spanOfMoleculeAtoms);

  EXPECT_NEAR(potentialEnergy * Units::EnergyToKelvin,
              (energy.frameworkMoleculeCharge + energy.moleculeMoleculeCharge + energy.ewald_fourier +
               energy.ewald_self + energy.ewald_exclusion) *
                  Units::EnergyToKelvin,
              1e-5);
}
