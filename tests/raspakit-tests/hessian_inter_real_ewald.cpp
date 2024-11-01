#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <cstddef>
#include <span>
#include <vector>
#include <numbers>

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
import gradient_factor;
import hessian_factor;
import running_energy;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;
import energy_status;

TEST(hessian_inter_real_ewald, Test_gradient_cartesian_methane_in_CHA_triclinic_1x1x1)
{
  ForceField forceField = ForceField(
      {
          PseudoAtom("Si", true, 28.0855, 2.05, 0.0, 14, false),
          PseudoAtom("O", true, 15.999, -1.025, 0.0, 8, false),
          PseudoAtom("CH4", false, 16.04246, 0.0, 0.0, 6, false),
          PseudoAtom("C_co2", false, 12.0, 0.6512, 0.0, 6, false),
          PseudoAtom("O_co2", false, 15.9994, -0.3256, 0.0, 8, false),
      },
      {VDWParameters(22.0, 2.30), VDWParameters(53.0, 3.3), VDWParameters(158.5, 3.72), VDWParameters(29.933, 2.745),
       VDWParameters(85.671, 3.017)},
      ForceField::MixingRule::Lorentz_Berthelot, 11.8, 11.8, 11.8, true, false, true);

  Framework f = Framework(
      0, forceField, "CHA", SimulationBox(9.459, 9.459, 9.459, 
        94.07 * std::numbers::pi / 180.0, 94.07 * std::numbers::pi / 180.0, 94.07 * std::numbers::pi / 180.0), 1,
      {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId,
       // uint8_t groupId
        Atom(double3(0.10330000, 0.33310000, 0.87430000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.12570000, 0.89670000, 0.66690000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.66690000, 0.89670000, 0.12570000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.10330000, 0.87430000, 0.33310000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.87430000, 0.33310000, 0.10330000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.33310000, 0.87430000, 0.10330000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.12570000, 0.66690000, 0.89670000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.33310000, 0.10330000, 0.87430000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.66690000, 0.12570000, 0.89670000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.87430000, 0.10330000, 0.33310000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.89670000, 0.12570000, 0.66690000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.89670000, 0.66690000, 0.12570000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.00000000, 0.73350000, 0.26650000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.26650000, 0.73350000, 1.00000000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.26650000, 1.00000000, 0.73350000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.73350000, 1.00000000, 0.26650000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(1.00000000, 0.26650000, 0.73350000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.73350000, 0.26650000, 1.00000000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.15060000, 0.84940000, 0.50000000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.50000000, 0.84940000, 0.15060000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.15060000, 0.50000000, 0.84940000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.84940000, 0.15060000, 0.50000000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.84940000, 0.50000000, 0.15060000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.50000000, 0.15060000, 0.84940000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.25030000, 0.89300000, 0.25030000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.10700000, 0.74970000, 0.74970000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.25030000, 0.25030000, 0.89300000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.89300000, 0.25030000, 0.25030000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.74970000, 0.10700000, 0.74970000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.74970000, 0.74970000, 0.10700000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.02040000, 0.31930000, 0.02040000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.68070000, 0.97960000, 0.97960000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.02040000, 0.02040000, 0.31930000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.31930000, 0.02040000, 0.02040000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.97960000, 0.68070000, 0.97960000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.97960000, 0.97960000, 0.68070000), -1.025, 1.0, 0, 1, 0, 0)
      },
      int3(1, 1, 1));
  Component c = Component(0, forceField, "methane", 190.564, 45599200, 0.01142,
                          {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                           // uint8_t componentId, uint8_t groupId
                           Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 2, 0, 0)},
                          5, 21);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {0}, 5);

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  double delta = 1e-5;
  double tolerance = 1e-4;
  double3 numerical_gradient;
  double chargeB = 1.0;
  double3 posB = double3(5.0, 5.0, 5.0);

  {
    double3 s = system.simulationBox.inverseCell * posB;
    auto [_, reference_cartesian, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);
    double3 reference_fractional = system.simulationBox.cell.transpose() * reference_cartesian;

    // finite difference x
    posB = system.simulationBox.cell * double3(s.x + 0.5 * delta, s.y, s.z);
    auto [x2_energy, _, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x - 0.5 * delta, s.y, s.z);
    auto [x1_energy, _, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    // finite difference y
    posB = system.simulationBox.cell * double3(s.x, s.y + 0.5 * delta, s.z);
    auto [y2_energy, _, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y - 0.5 * delta, s.z);
    auto [y1_energy, _, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    // finite difference z
    posB = system.simulationBox.cell * double3(s.x, s.y, s.z + 0.5 * delta);
    auto [z2_energy, _, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z - 0.5 * delta);
    auto [z1_energy, _, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    numerical_gradient.x = (x2_energy - x1_energy) / delta;
    numerical_gradient.y = (y2_energy - y1_energy) / delta;
    numerical_gradient.z = (z2_energy - z1_energy) / delta;

    EXPECT_NEAR(reference_fractional.x, numerical_gradient.x, tolerance);
    EXPECT_NEAR(reference_fractional.y, numerical_gradient.y, tolerance);
    EXPECT_NEAR(reference_fractional.z, numerical_gradient.z, tolerance);
  }
}


TEST(hessian_inter_real_ewald, Test_gradient_fractional_methane_in_CHA_triclinic_1x1x1)
{
  ForceField forceField = ForceField(
      {
          PseudoAtom("Si", true, 28.0855, 2.05, 0.0, 14, false),
          PseudoAtom("O", true, 15.999, -1.025, 0.0, 8, false),
          PseudoAtom("CH4", false, 16.04246, 0.0, 0.0, 6, false),
          PseudoAtom("C_co2", false, 12.0, 0.6512, 0.0, 6, false),
          PseudoAtom("O_co2", false, 15.9994, -0.3256, 0.0, 8, false),
      },
      {VDWParameters(22.0, 2.30), VDWParameters(53.0, 3.3), VDWParameters(158.5, 3.72), VDWParameters(29.933, 2.745),
       VDWParameters(85.671, 3.017)},
      ForceField::MixingRule::Lorentz_Berthelot, 11.8, 11.8, 11.8, true, false, true);

  Framework f = Framework(
      0, forceField, "CHA", SimulationBox(9.459, 9.459, 9.459, 
        94.07 * std::numbers::pi / 180.0, 94.07 * std::numbers::pi / 180.0, 94.07 * std::numbers::pi / 180.0), 1,
      {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId,
       // uint8_t groupId
        Atom(double3(0.10330000, 0.33310000, 0.87430000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.12570000, 0.89670000, 0.66690000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.66690000, 0.89670000, 0.12570000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.10330000, 0.87430000, 0.33310000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.87430000, 0.33310000, 0.10330000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.33310000, 0.87430000, 0.10330000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.12570000, 0.66690000, 0.89670000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.33310000, 0.10330000, 0.87430000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.66690000, 0.12570000, 0.89670000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.87430000, 0.10330000, 0.33310000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.89670000, 0.12570000, 0.66690000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.89670000, 0.66690000, 0.12570000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.00000000, 0.73350000, 0.26650000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.26650000, 0.73350000, 1.00000000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.26650000, 1.00000000, 0.73350000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.73350000, 1.00000000, 0.26650000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(1.00000000, 0.26650000, 0.73350000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.73350000, 0.26650000, 1.00000000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.15060000, 0.84940000, 0.50000000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.50000000, 0.84940000, 0.15060000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.15060000, 0.50000000, 0.84940000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.84940000, 0.15060000, 0.50000000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.84940000, 0.50000000, 0.15060000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.50000000, 0.15060000, 0.84940000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.25030000, 0.89300000, 0.25030000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.10700000, 0.74970000, 0.74970000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.25030000, 0.25030000, 0.89300000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.89300000, 0.25030000, 0.25030000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.74970000, 0.10700000, 0.74970000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.74970000, 0.74970000, 0.10700000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.02040000, 0.31930000, 0.02040000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.68070000, 0.97960000, 0.97960000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.02040000, 0.02040000, 0.31930000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.31930000, 0.02040000, 0.02040000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.97960000, 0.68070000, 0.97960000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.97960000, 0.97960000, 0.68070000), -1.025, 1.0, 0, 1, 0, 0)
      },
      int3(1, 1, 1));
  Component c = Component(0, forceField, "methane", 190.564, 45599200, 0.01142,
                          {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                           // uint8_t componentId, uint8_t groupId
                           Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 2, 0, 0)},
                          5, 21);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {0}, 5);

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  double delta = 1e-5;
  double tolerance = 1e-4;
  double3 numerical_gradient;
  double chargeB = 1.0;
  double3 posB = double3(5.0, 5.0, 5.0);

  {
    double3 s = system.simulationBox.inverseCell * posB;
    auto [_, reference_cartesian, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);
    double3 reference_fractional = system.simulationBox.cell.transpose() * reference_cartesian;

    // finite difference x
    posB = system.simulationBox.cell * double3(s.x + 0.5 * delta, s.y, s.z);
    auto [x2_energy, _, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x - 0.5 * delta, s.y, s.z);
    auto [x1_energy, _, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    // finite difference y
    posB = system.simulationBox.cell * double3(s.x, s.y + 0.5 * delta, s.z);
    auto [y2_energy, _, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y - 0.5 * delta, s.z);
    auto [y1_energy, _, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    // finite difference z
    posB = system.simulationBox.cell * double3(s.x, s.y, s.z + 0.5 * delta);
    auto [z2_energy, _, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z - 0.5 * delta);
    auto [z1_energy, _, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    numerical_gradient.x = (x2_energy - x1_energy) / delta;
    numerical_gradient.y = (y2_energy - y1_energy) / delta;
    numerical_gradient.z = (z2_energy - z1_energy) / delta;

    EXPECT_NEAR(reference_fractional.x, numerical_gradient.x, tolerance);
    EXPECT_NEAR(reference_fractional.y, numerical_gradient.y, tolerance);
    EXPECT_NEAR(reference_fractional.z, numerical_gradient.z, tolerance);
  }
}


TEST(hessian_inter_real_ewald, Test_hessian_cartesian_methane_in_CHA_triclinic_1x1x1)
{
  ForceField forceField = ForceField(
      {
          PseudoAtom("Si", true, 28.0855, 2.05, 0.0, 14, false),
          PseudoAtom("O", true, 15.999, -1.025, 0.0, 8, false),
          PseudoAtom("CH4", false, 16.04246, 0.0, 0.0, 6, false),
          PseudoAtom("C_co2", false, 12.0, 0.6512, 0.0, 6, false),
          PseudoAtom("O_co2", false, 15.9994, -0.3256, 0.0, 8, false),
      },
      {VDWParameters(22.0, 2.30), VDWParameters(53.0, 3.3), VDWParameters(158.5, 3.72), VDWParameters(29.933, 2.745),
       VDWParameters(85.671, 3.017)},
      ForceField::MixingRule::Lorentz_Berthelot, 11.8, 11.8, 11.8, true, false, true);

  Framework f = Framework(
      0, forceField, "CHA", SimulationBox(9.459, 9.459, 9.459, 
        94.07 * std::numbers::pi / 180.0, 94.07 * std::numbers::pi / 180.0, 94.07 * std::numbers::pi / 180.0), 1,
      {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId,
       // uint8_t groupId
        Atom(double3(0.10330000, 0.33310000, 0.87430000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.12570000, 0.89670000, 0.66690000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.66690000, 0.89670000, 0.12570000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.10330000, 0.87430000, 0.33310000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.87430000, 0.33310000, 0.10330000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.33310000, 0.87430000, 0.10330000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.12570000, 0.66690000, 0.89670000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.33310000, 0.10330000, 0.87430000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.66690000, 0.12570000, 0.89670000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.87430000, 0.10330000, 0.33310000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.89670000, 0.12570000, 0.66690000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.89670000, 0.66690000, 0.12570000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.00000000, 0.73350000, 0.26650000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.26650000, 0.73350000, 1.00000000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.26650000, 1.00000000, 0.73350000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.73350000, 1.00000000, 0.26650000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(1.00000000, 0.26650000, 0.73350000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.73350000, 0.26650000, 1.00000000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.15060000, 0.84940000, 0.50000000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.50000000, 0.84940000, 0.15060000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.15060000, 0.50000000, 0.84940000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.84940000, 0.15060000, 0.50000000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.84940000, 0.50000000, 0.15060000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.50000000, 0.15060000, 0.84940000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.25030000, 0.89300000, 0.25030000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.10700000, 0.74970000, 0.74970000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.25030000, 0.25030000, 0.89300000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.89300000, 0.25030000, 0.25030000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.74970000, 0.10700000, 0.74970000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.74970000, 0.74970000, 0.10700000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.02040000, 0.31930000, 0.02040000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.68070000, 0.97960000, 0.97960000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.02040000, 0.02040000, 0.31930000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.31930000, 0.02040000, 0.02040000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.97960000, 0.68070000, 0.97960000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.97960000, 0.97960000, 0.68070000), -1.025, 1.0, 0, 1, 0, 0)
      },
      int3(1, 1, 1));
  Component c = Component(0, forceField, "methane", 190.564, 45599200, 0.01142,
                          {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                           // uint8_t componentId, uint8_t groupId
                           Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 2, 0, 0)},
                          5, 21);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {0}, 5);

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  double delta = 1e-5;
  double tolerance = 1e-4;
  double3x3 numerical_hessian;
  double chargeB = 1.0;
  double3 posB_reference = double3(5.0, 5.0, 5.0);

  {
    double3 posB;

    auto [_, _, reference_cartesian] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB_reference, chargeB, frameworkAtoms);

    // finite difference x
    posB = double3(posB_reference.x + 0.5 * delta, posB_reference.y, posB_reference.z);
    auto [_, x2_gradient, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    posB = double3(posB_reference.x - 0.5 * delta, posB_reference.y, posB_reference.z);
    auto [_, x1_gradient, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    // finite difference y
    posB = double3(posB_reference.x, posB_reference.y + 0.5 * delta, posB_reference.z);
    auto [_, y2_gradient, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    posB = double3(posB_reference.x, posB_reference.y - 0.5 * delta, posB_reference.z);
    auto [_, y1_gradient, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    // finite difference z
    posB = double3(posB_reference.x, posB_reference.y, posB_reference.z + 0.5 * delta);
    auto [_, z2_gradient, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    posB = double3(posB_reference.x, posB_reference.y, posB_reference.z - 0.5 * delta);
    auto [_, z1_gradient, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    numerical_hessian.ax = (x2_gradient.x - x1_gradient.x) / delta;
    numerical_hessian.ay = (x2_gradient.y - x1_gradient.y) / delta;
    numerical_hessian.az = (x2_gradient.z - x1_gradient.z) / delta;

    numerical_hessian.bx = (y2_gradient.x - y1_gradient.x) / delta;
    numerical_hessian.by = (y2_gradient.y - y1_gradient.y) / delta;
    numerical_hessian.bz = (y2_gradient.z - y1_gradient.z) / delta;

    numerical_hessian.cx = (z2_gradient.x - z1_gradient.x) / delta;
    numerical_hessian.cy = (z2_gradient.y - z1_gradient.y) / delta;
    numerical_hessian.cz = (z2_gradient.z - z1_gradient.z) / delta;

    EXPECT_NEAR(reference_cartesian.ax, numerical_hessian.ax, tolerance);
    EXPECT_NEAR(reference_cartesian.ay, numerical_hessian.ay, tolerance);
    EXPECT_NEAR(reference_cartesian.az, numerical_hessian.az, tolerance);

    EXPECT_NEAR(reference_cartesian.bx, numerical_hessian.bx, tolerance);
    EXPECT_NEAR(reference_cartesian.by, numerical_hessian.by, tolerance);
    EXPECT_NEAR(reference_cartesian.bz, numerical_hessian.bz, tolerance);

    EXPECT_NEAR(reference_cartesian.cx, numerical_hessian.cx, tolerance);
    EXPECT_NEAR(reference_cartesian.cy, numerical_hessian.cy, tolerance);
    EXPECT_NEAR(reference_cartesian.cz, numerical_hessian.cz, tolerance);
  }
}


TEST(hessian_inter_real_ewald, Test_hessian_fractional_methane_in_CHA_triclinic_1x1x1)
{
  ForceField forceField = ForceField(
      {
          PseudoAtom("Si", true, 28.0855, 2.05, 0.0, 14, false),
          PseudoAtom("O", true, 15.999, -1.025, 0.0, 8, false),
          PseudoAtom("CH4", false, 16.04246, 0.0, 0.0, 6, false),
          PseudoAtom("C_co2", false, 12.0, 0.6512, 0.0, 6, false),
          PseudoAtom("O_co2", false, 15.9994, -0.3256, 0.0, 8, false),
      },
      {VDWParameters(22.0, 2.30), VDWParameters(53.0, 3.3), VDWParameters(158.5, 3.72), VDWParameters(29.933, 2.745),
       VDWParameters(85.671, 3.017)},
      ForceField::MixingRule::Lorentz_Berthelot, 11.8, 11.8, 11.8, true, false, true);

  Framework f = Framework(
      0, forceField, "CHA", SimulationBox(9.459, 9.459, 9.459, 
        94.07 * std::numbers::pi / 180.0, 94.07 * std::numbers::pi / 180.0, 94.07 * std::numbers::pi / 180.0), 1,
      {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type, uint8_t componentId,
       // uint8_t groupId
        Atom(double3(0.10330000, 0.33310000, 0.87430000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.12570000, 0.89670000, 0.66690000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.66690000, 0.89670000, 0.12570000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.10330000, 0.87430000, 0.33310000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.87430000, 0.33310000, 0.10330000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.33310000, 0.87430000, 0.10330000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.12570000, 0.66690000, 0.89670000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.33310000, 0.10330000, 0.87430000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.66690000, 0.12570000, 0.89670000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.87430000, 0.10330000, 0.33310000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.89670000, 0.12570000, 0.66690000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.89670000, 0.66690000, 0.12570000), 2.05, 1.0, 0, 0, 0, 0),
        Atom(double3(0.00000000, 0.73350000, 0.26650000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.26650000, 0.73350000, 1.00000000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.26650000, 1.00000000, 0.73350000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.73350000, 1.00000000, 0.26650000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(1.00000000, 0.26650000, 0.73350000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.73350000, 0.26650000, 1.00000000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.15060000, 0.84940000, 0.50000000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.50000000, 0.84940000, 0.15060000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.15060000, 0.50000000, 0.84940000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.84940000, 0.15060000, 0.50000000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.84940000, 0.50000000, 0.15060000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.50000000, 0.15060000, 0.84940000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.25030000, 0.89300000, 0.25030000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.10700000, 0.74970000, 0.74970000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.25030000, 0.25030000, 0.89300000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.89300000, 0.25030000, 0.25030000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.74970000, 0.10700000, 0.74970000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.74970000, 0.74970000, 0.10700000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.02040000, 0.31930000, 0.02040000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.68070000, 0.97960000, 0.97960000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.02040000, 0.02040000, 0.31930000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.31930000, 0.02040000, 0.02040000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.97960000, 0.68070000, 0.97960000), -1.025, 1.0, 0, 1, 0, 0),
        Atom(double3(0.97960000, 0.97960000, 0.68070000), -1.025, 1.0, 0, 1, 0, 0)
      },
      int3(1, 1, 1));
  Component c = Component(0, forceField, "methane", 190.564, 45599200, 0.01142,
                          {// double3 position, double charge, double lambda, uint32_t moleculeId, uint16_t type,
                           // uint8_t componentId, uint8_t groupId
                           Atom(double3(0.0, 0.0, 0.0), 0.0, 1.0, 0, 2, 0, 0)},
                          5, 21);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {0}, 5);

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  double delta = 1e-5;
  double tolerance = 1e-4;
  double3x3 numerical_hessian;
  double chargeB = 1.0;
  double3 posB = double3(5.0, 5.0, 5.0);

  {
    double3 s = system.simulationBox.inverseCell * posB;
    auto [_, _, reference_cartesian] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);
    double3x3 reference_fractional = system.simulationBox.cell.transpose() * reference_cartesian * system.simulationBox.cell;

    double test = system.simulationBox.cell.ax * system.simulationBox.cell.ax * reference_cartesian.ax +
                  system.simulationBox.cell.ay * system.simulationBox.cell.ax * reference_cartesian.bx +
                  system.simulationBox.cell.az * system.simulationBox.cell.ax * reference_cartesian.cx +
                  system.simulationBox.cell.ax * system.simulationBox.cell.ay * reference_cartesian.ay +
                  system.simulationBox.cell.ay * system.simulationBox.cell.ay * reference_cartesian.by +
                  system.simulationBox.cell.az * system.simulationBox.cell.ay * reference_cartesian.cy +
                  system.simulationBox.cell.ax * system.simulationBox.cell.az * reference_cartesian.az +
                  system.simulationBox.cell.ay * system.simulationBox.cell.az * reference_cartesian.bz +
                  system.simulationBox.cell.az * system.simulationBox.cell.az * reference_cartesian.cz;


    // finite difference x
    posB = system.simulationBox.cell * double3(s.x + 0.5 * delta, s.y, s.z);
    auto [_, x2_gradient, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x - 0.5 * delta, s.y, s.z);
    auto [_, x1_gradient, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    // finite difference y
    posB = system.simulationBox.cell * double3(s.x, s.y + 0.5 * delta, s.z);
    auto [_, y2_gradient, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y - 0.5 * delta, s.z);
    auto [_, y1_gradient, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    // finite difference z
    posB = system.simulationBox.cell * double3(s.x, s.y, s.z + 0.5 * delta);
    auto [_, z2_gradient, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z - 0.5 * delta);
    auto [_, z1_gradient, _] = Interactions::calculateHessianAtPositionCoulomb(system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    double3 x2_gradient_corr = system.simulationBox.cell.transpose() * x2_gradient;
    double3 x1_gradient_corr = system.simulationBox.cell.transpose() * x1_gradient;

    double3 y2_gradient_corr = system.simulationBox.cell.transpose() * y2_gradient;
    double3 y1_gradient_corr = system.simulationBox.cell.transpose() * y1_gradient;

    double3 z2_gradient_corr = system.simulationBox.cell.transpose() * z2_gradient;
    double3 z1_gradient_corr = system.simulationBox.cell.transpose() * z1_gradient;

    numerical_hessian.ax = (x2_gradient_corr.x - x1_gradient_corr.x) / delta;
    numerical_hessian.ay = (x2_gradient_corr.y - x1_gradient_corr.y) / delta;
    numerical_hessian.az = (x2_gradient_corr.z - x1_gradient_corr.z) / delta;

    numerical_hessian.bx = (y2_gradient_corr.x - y1_gradient_corr.x) / delta;
    numerical_hessian.by = (y2_gradient_corr.y - y1_gradient_corr.y) / delta;
    numerical_hessian.bz = (y2_gradient_corr.z - y1_gradient_corr.z) / delta;

    numerical_hessian.cx = (z2_gradient_corr.x - z1_gradient_corr.x) / delta;
    numerical_hessian.cy = (z2_gradient_corr.y - z1_gradient_corr.y) / delta;
    numerical_hessian.cz = (z2_gradient_corr.z - z1_gradient_corr.z) / delta;

    EXPECT_NEAR(reference_fractional.ax, numerical_hessian.ax, tolerance);
    EXPECT_NEAR(reference_fractional.ay, numerical_hessian.ay, tolerance);
    EXPECT_NEAR(reference_fractional.az, numerical_hessian.az, tolerance);

    EXPECT_NEAR(reference_fractional.bx, numerical_hessian.bx, tolerance);
    EXPECT_NEAR(reference_fractional.by, numerical_hessian.by, tolerance);
    EXPECT_NEAR(reference_fractional.bz, numerical_hessian.bz, tolerance);

    EXPECT_NEAR(reference_fractional.cx, numerical_hessian.cx, tolerance);
    EXPECT_NEAR(reference_fractional.cy, numerical_hessian.cy, tolerance);
    EXPECT_NEAR(reference_fractional.cz, numerical_hessian.cz, tolerance);
  }
}

