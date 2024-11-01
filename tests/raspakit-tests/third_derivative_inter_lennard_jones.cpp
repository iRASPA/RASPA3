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
import third_derivative_factor;
import running_energy;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;
import energy_status;

TEST(third_derivative_inter_lennard_jones, Test_gradient_cartesian_methane_in_CHA_triclinic_1x1x1)
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
  size_t typeB = 2;
  double3 posB = double3(5.0, 5.0, 5.0);

  {
    double3 s = system.simulationBox.inverseCell * posB;
    auto [_, reference_cartesian, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);
    double3 reference_fractional = system.simulationBox.cell.transpose() * reference_cartesian;

    // finite difference x
    posB = system.simulationBox.cell * double3(s.x + 0.5 * delta, s.y, s.z);
    auto [x2_energy, _, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x - 0.5 * delta, s.y, s.z);
    auto [x1_energy, _, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    // finite difference y
    posB = system.simulationBox.cell * double3(s.x, s.y + 0.5 * delta, s.z);
    auto [y2_energy, _, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y - 0.5 * delta, s.z);
    auto [y1_energy, _, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    // finite difference z
    posB = system.simulationBox.cell * double3(s.x, s.y, s.z + 0.5 * delta);
    auto [z2_energy, _, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z - 0.5 * delta);
    auto [z1_energy, _, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    numerical_gradient.x = (x2_energy - x1_energy) / delta;
    numerical_gradient.y = (y2_energy - y1_energy) / delta;
    numerical_gradient.z = (z2_energy - z1_energy) / delta;

    EXPECT_NEAR(reference_fractional.x, numerical_gradient.x, tolerance);
    EXPECT_NEAR(reference_fractional.y, numerical_gradient.y, tolerance);
    EXPECT_NEAR(reference_fractional.z, numerical_gradient.z, tolerance);
  }
}


TEST(third_derivative_inter_lennard_jones, Test_gradient_fractional_methane_in_CHA_triclinic_1x1x1)
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
  size_t typeB = 2;
  double3 posB = double3(5.0, 5.0, 5.0);

  {
    double3 s = system.simulationBox.inverseCell * posB;
    auto [_, reference_cartesian, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);
    double3 reference_fractional = system.simulationBox.cell.transpose() * reference_cartesian;

    // finite difference x
    posB = system.simulationBox.cell * double3(s.x + 0.5 * delta, s.y, s.z);
    auto [x2_energy, _, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x - 0.5 * delta, s.y, s.z);
    auto [x1_energy, _, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    // finite difference y
    posB = system.simulationBox.cell * double3(s.x, s.y + 0.5 * delta, s.z);
    auto [y2_energy, _, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y - 0.5 * delta, s.z);
    auto [y1_energy, _, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    // finite difference z
    posB = system.simulationBox.cell * double3(s.x, s.y, s.z + 0.5 * delta);
    auto [z2_energy, _, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z - 0.5 * delta);
    auto [z1_energy, _, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    numerical_gradient.x = (x2_energy - x1_energy) / delta;
    numerical_gradient.y = (y2_energy - y1_energy) / delta;
    numerical_gradient.z = (z2_energy - z1_energy) / delta;

    EXPECT_NEAR(reference_fractional.x, numerical_gradient.x, tolerance);
    EXPECT_NEAR(reference_fractional.y, numerical_gradient.y, tolerance);
    EXPECT_NEAR(reference_fractional.z, numerical_gradient.z, tolerance);
  }
}


TEST(third_derivative_inter_lennard_jones, Test_hessian_cartesian_methane_in_CHA_triclinic_1x1x1)
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
  size_t typeB = 2;
  double3 posB_reference = double3(5.0, 5.0, 5.0);

  {
    double3 posB;

    auto [_, _, reference_cartesian, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB_reference, typeB, frameworkAtoms);

    // finite difference x
    posB = double3(posB_reference.x + 0.5 * delta, posB_reference.y, posB_reference.z);
    auto [_, x2_gradient, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    posB = double3(posB_reference.x - 0.5 * delta, posB_reference.y, posB_reference.z);
    auto [_, x1_gradient, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    // finite difference y
    posB = double3(posB_reference.x, posB_reference.y + 0.5 * delta, posB_reference.z);
    auto [_, y2_gradient, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    posB = double3(posB_reference.x, posB_reference.y - 0.5 * delta, posB_reference.z);
    auto [_, y1_gradient, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    // finite difference z
    posB = double3(posB_reference.x, posB_reference.y, posB_reference.z + 0.5 * delta);
    auto [_, z2_gradient, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    posB = double3(posB_reference.x, posB_reference.y, posB_reference.z - 0.5 * delta);
    auto [_, z1_gradient, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

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


TEST(third_derivative_inter_lennard_jones, Test_hessian_fractional_methane_in_CHA_triclinic_1x1x1)
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
  size_t typeB = 2;
  double3 posB = double3(5.0, 5.0, 5.0);

  {
    double3 s = system.simulationBox.inverseCell * posB;
    auto [_, _, reference_cartesian, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);
    double3x3 reference_fractional = system.simulationBox.cell.transpose() * reference_cartesian * system.simulationBox.cell;

    // finite difference x
    posB = system.simulationBox.cell * double3(s.x + 0.5 * delta, s.y, s.z);
    auto [_, x2_gradient, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x - 0.5 * delta, s.y, s.z);
    auto [_, x1_gradient, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    // finite difference y
    posB = system.simulationBox.cell * double3(s.x, s.y + 0.5 * delta, s.z);
    auto [_, y2_gradient, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y - 0.5 * delta, s.z);
    auto [_, y1_gradient, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    // finite difference z
    posB = system.simulationBox.cell * double3(s.x, s.y, s.z + 0.5 * delta);
    auto [_, z2_gradient, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z - 0.5 * delta);
    auto [_, z1_gradient, _, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

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


TEST(third_derivative_inter_lennard_jones, Test_third_derivative_cartesian_methane_in_CHA_triclinic_1x1x1)
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
  size_t typeB = 2;
  double3 posB_reference = double3(5.0, 5.0, 5.0);

  {
    double3 posB;

    auto [_, _, _, reference_cartesian] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB_reference, typeB, frameworkAtoms);

    posB = double3(posB_reference.x + 0.5 * delta, posB_reference.y, posB_reference.z);
    auto [_, _, x2_hessian, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    posB = double3(posB_reference.x - 0.5 * delta, posB_reference.y, posB_reference.z);
    auto [_, _, x1_hessian, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);


    double numerical_third_derivative_m111 = (x2_hessian.ax - x1_hessian.ax) / delta;
    double numerical_third_derivative_m112 = (x2_hessian.ay - x1_hessian.ay) / delta;
    double numerical_third_derivative_m113 = (x2_hessian.az - x1_hessian.az) / delta;
    double numerical_third_derivative_m121 = (x2_hessian.bx - x1_hessian.bx) / delta;
    double numerical_third_derivative_m122 = (x2_hessian.by - x1_hessian.by) / delta;
    double numerical_third_derivative_m123 = (x2_hessian.bz - x1_hessian.bz) / delta;
    double numerical_third_derivative_m131 = (x2_hessian.cx - x1_hessian.cx) / delta;
    double numerical_third_derivative_m132 = (x2_hessian.cy - x1_hessian.cy) / delta;
    double numerical_third_derivative_m133 = (x2_hessian.cz - x1_hessian.cz) / delta;

    EXPECT_NEAR(reference_cartesian.m111, numerical_third_derivative_m111, tolerance);
    EXPECT_NEAR(reference_cartesian.m112, numerical_third_derivative_m112, tolerance);
    EXPECT_NEAR(reference_cartesian.m113, numerical_third_derivative_m113, tolerance);
    EXPECT_NEAR(reference_cartesian.m121, numerical_third_derivative_m121, tolerance);
    EXPECT_NEAR(reference_cartesian.m122, numerical_third_derivative_m122, tolerance);
    EXPECT_NEAR(reference_cartesian.m123, numerical_third_derivative_m123, tolerance);
    EXPECT_NEAR(reference_cartesian.m131, numerical_third_derivative_m131, tolerance);
    EXPECT_NEAR(reference_cartesian.m132, numerical_third_derivative_m132, tolerance);
    EXPECT_NEAR(reference_cartesian.m133, numerical_third_derivative_m133, tolerance);


    posB = double3(posB_reference.x, posB_reference.y + 0.5 * delta, posB_reference.z);
    auto [_, _, y2_hessian, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    posB = double3(posB_reference.x, posB_reference.y - 0.5 * delta, posB_reference.z);
    auto [_, _, y1_hessian, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    double numerical_third_derivative_m211 = (y2_hessian.ax - y1_hessian.ax) / delta;
    double numerical_third_derivative_m212 = (y2_hessian.ay - y1_hessian.ay) / delta;
    double numerical_third_derivative_m213 = (y2_hessian.az - y1_hessian.az) / delta;
    double numerical_third_derivative_m221 = (y2_hessian.bx - y1_hessian.bx) / delta;
    double numerical_third_derivative_m222 = (y2_hessian.by - y1_hessian.by) / delta;
    double numerical_third_derivative_m223 = (y2_hessian.bz - y1_hessian.bz) / delta;
    double numerical_third_derivative_m231 = (y2_hessian.cx - y1_hessian.cx) / delta;
    double numerical_third_derivative_m232 = (y2_hessian.cy - y1_hessian.cy) / delta;
    double numerical_third_derivative_m233 = (y2_hessian.cz - y1_hessian.cz) / delta;

    EXPECT_NEAR(reference_cartesian.m211, numerical_third_derivative_m211, tolerance);
    EXPECT_NEAR(reference_cartesian.m212, numerical_third_derivative_m212, tolerance);
    EXPECT_NEAR(reference_cartesian.m213, numerical_third_derivative_m213, tolerance);
    EXPECT_NEAR(reference_cartesian.m221, numerical_third_derivative_m221, tolerance);
    EXPECT_NEAR(reference_cartesian.m222, numerical_third_derivative_m222, tolerance);
    EXPECT_NEAR(reference_cartesian.m223, numerical_third_derivative_m223, tolerance);
    EXPECT_NEAR(reference_cartesian.m231, numerical_third_derivative_m231, tolerance);
    EXPECT_NEAR(reference_cartesian.m232, numerical_third_derivative_m232, tolerance);
    EXPECT_NEAR(reference_cartesian.m233, numerical_third_derivative_m233, tolerance);


    posB = double3(posB_reference.x, posB_reference.y, posB_reference.z + 0.5 * delta);
    auto [_, _, z2_hessian, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    posB = double3(posB_reference.x, posB_reference.y, posB_reference.z - 0.5 * delta);
    auto [_, _, z1_hessian, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    double numerical_third_derivative_m311 = (z2_hessian.ax - z1_hessian.ax) / delta;
    double numerical_third_derivative_m312 = (z2_hessian.ay - z1_hessian.ay) / delta;
    double numerical_third_derivative_m313 = (z2_hessian.az - z1_hessian.az) / delta;
    double numerical_third_derivative_m321 = (z2_hessian.bx - z1_hessian.bx) / delta;
    double numerical_third_derivative_m322 = (z2_hessian.by - z1_hessian.by) / delta;
    double numerical_third_derivative_m323 = (z2_hessian.bz - z1_hessian.bz) / delta;
    double numerical_third_derivative_m331 = (z2_hessian.cx - z1_hessian.cx) / delta;
    double numerical_third_derivative_m332 = (z2_hessian.cy - z1_hessian.cy) / delta;
    double numerical_third_derivative_m333 = (z2_hessian.cz - z1_hessian.cz) / delta;

    EXPECT_NEAR(reference_cartesian.m311, numerical_third_derivative_m311, tolerance);
    EXPECT_NEAR(reference_cartesian.m312, numerical_third_derivative_m312, tolerance);
    EXPECT_NEAR(reference_cartesian.m313, numerical_third_derivative_m313, tolerance);
    EXPECT_NEAR(reference_cartesian.m321, numerical_third_derivative_m321, tolerance);
    EXPECT_NEAR(reference_cartesian.m322, numerical_third_derivative_m322, tolerance);
    EXPECT_NEAR(reference_cartesian.m323, numerical_third_derivative_m323, tolerance);
    EXPECT_NEAR(reference_cartesian.m331, numerical_third_derivative_m331, tolerance);
    EXPECT_NEAR(reference_cartesian.m332, numerical_third_derivative_m332, tolerance);
    EXPECT_NEAR(reference_cartesian.m333, numerical_third_derivative_m333, tolerance);
  }
}

TEST(third_derivative_inter_lennard_jones, Test_third_derivative_fractional_methane_in_CHA_triclinic_1x1x1)
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

  double delta = 1e-6;
  double tolerance = 1e-4;
  double3x3 numerical_hessian;
  size_t typeB = 2;
  double3 posB = double3(5.0, 5.0, 5.0);

  {
    double3 s = system.simulationBox.inverseCell * posB;
    auto [_, _, _, reference_cartesian] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    double reference_fractional = 
      system.simulationBox.cell.ax * system.simulationBox.cell.bx * system.simulationBox.cell.cx * reference_cartesian.m111 +
      system.simulationBox.cell.ax * system.simulationBox.cell.bx * system.simulationBox.cell.cy * reference_cartesian.m112 +
      system.simulationBox.cell.ax * system.simulationBox.cell.bx * system.simulationBox.cell.cz * reference_cartesian.m113 +
      system.simulationBox.cell.ax * system.simulationBox.cell.by * system.simulationBox.cell.cx * reference_cartesian.m121 +
      system.simulationBox.cell.ax * system.simulationBox.cell.by * system.simulationBox.cell.cy * reference_cartesian.m122 +
      system.simulationBox.cell.ax * system.simulationBox.cell.by * system.simulationBox.cell.cz * reference_cartesian.m132;

    /* cell.ay = cell.az = cell.bz = 0
      system.simulationBox.cell.ax * system.simulationBox.cell.bz * system.simulationBox.cell.cx * reference_cartesian.m131 +
      system.simulationBox.cell.ax * system.simulationBox.cell.bz * system.simulationBox.cell.cy * reference_cartesian.m132 +
      system.simulationBox.cell.ax * system.simulationBox.cell.bz * system.simulationBox.cell.cz * reference_cartesian.m133 +

      system.simulationBox.cell.ay * system.simulationBox.cell.bx * system.simulationBox.cell.cx * reference_cartesian.m211 +
      system.simulationBox.cell.ay * system.simulationBox.cell.bx * system.simulationBox.cell.cy * reference_cartesian.m212 +
      system.simulationBox.cell.ay * system.simulationBox.cell.bx * system.simulationBox.cell.cz * reference_cartesian.m213 +
      system.simulationBox.cell.ay * system.simulationBox.cell.by * system.simulationBox.cell.cx * reference_cartesian.m221 +
      system.simulationBox.cell.ay * system.simulationBox.cell.by * system.simulationBox.cell.cy * reference_cartesian.m222 +
      system.simulationBox.cell.ay * system.simulationBox.cell.by * system.simulationBox.cell.cz * reference_cartesian.m223 +
      system.simulationBox.cell.ay * system.simulationBox.cell.bz * system.simulationBox.cell.cx * reference_cartesian.m231 +
      system.simulationBox.cell.ay * system.simulationBox.cell.bz * system.simulationBox.cell.cy * reference_cartesian.m232 +
      system.simulationBox.cell.ay * system.simulationBox.cell.bz * system.simulationBox.cell.cz * reference_cartesian.m233 +

      system.simulationBox.cell.az * system.simulationBox.cell.bx * system.simulationBox.cell.cx * reference_cartesian.m311 +
      system.simulationBox.cell.az * system.simulationBox.cell.bx * system.simulationBox.cell.cy * reference_cartesian.m312 +
      system.simulationBox.cell.az * system.simulationBox.cell.bx * system.simulationBox.cell.cz * reference_cartesian.m313 +
      system.simulationBox.cell.az * system.simulationBox.cell.by * system.simulationBox.cell.cx * reference_cartesian.m321 +
      system.simulationBox.cell.az * system.simulationBox.cell.by * system.simulationBox.cell.cy * reference_cartesian.m322 +
      system.simulationBox.cell.az * system.simulationBox.cell.by * system.simulationBox.cell.cz * reference_cartesian.m323 +
      system.simulationBox.cell.az * system.simulationBox.cell.bz * system.simulationBox.cell.cx * reference_cartesian.m331 +
      system.simulationBox.cell.az * system.simulationBox.cell.bz * system.simulationBox.cell.cy * reference_cartesian.m332 +
      system.simulationBox.cell.az * system.simulationBox.cell.bz * system.simulationBox.cell.cz * reference_cartesian.m333;
      */

    posB = system.simulationBox.cell * double3(s.x + 0.5 * delta, s.y, s.z);
    auto [_, _, x2_hessian, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x - 0.5 * delta, s.y, s.z);
    auto [_, _, x1_hessian, _] = Interactions::calculateThirdDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    double3x3 x2_hessian_corr = system.simulationBox.cell.transpose() * x2_hessian * system.simulationBox.cell;
    double3x3 x1_hessian_corr = system.simulationBox.cell.transpose() * x1_hessian * system.simulationBox.cell;

    double numerical_third_derivative_m123 = (x2_hessian_corr.bz - x1_hessian_corr.bz) / delta;

    EXPECT_NEAR(numerical_third_derivative_m123, reference_fractional, tolerance);
  }
}
