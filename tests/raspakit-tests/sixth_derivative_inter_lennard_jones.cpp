#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <cstddef>
#include <span>
#include <vector>
#include <numbers>
#include <format>
#include <print>

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
import sixth_derivative_factor;
import running_energy;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;
import energy_status;

TEST(sixth_derivative_inter_lennard_jones, Test_gradient_cartesian_methane_in_CHA_triclinic_1x1x1)
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
    auto [_, reference_cartesian, _, _, _, _, _] = Interactions::calculateSixthDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);
    double3 reference_fractional = system.simulationBox.cell.transpose() * double3(reference_cartesian[0], reference_cartesian[1], reference_cartesian[2]);

    // finite difference x
    posB = system.simulationBox.cell * double3(s.x + 0.5 * delta, s.y, s.z);
    auto [x2_energy, _, _, _, _, _, _] = Interactions::calculateSixthDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x - 0.5 * delta, s.y, s.z);
    auto [x1_energy, _, _, _, _, _, _] = Interactions::calculateSixthDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    // finite difference y
    posB = system.simulationBox.cell * double3(s.x, s.y + 0.5 * delta, s.z);
    auto [y2_energy, _, _, _, _, _, _] = Interactions::calculateSixthDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y - 0.5 * delta, s.z);
    auto [y1_energy, _, _, _, _, _, _] = Interactions::calculateSixthDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    // finite difference z
    posB = system.simulationBox.cell * double3(s.x, s.y, s.z + 0.5 * delta);
    auto [z2_energy, _, _, _, _, _, _] = Interactions::calculateSixthDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z - 0.5 * delta);
    auto [z1_energy, _, _, _, _, _, _] = Interactions::calculateSixthDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

    numerical_gradient.x = (x2_energy - x1_energy) / delta;
    numerical_gradient.y = (y2_energy - y1_energy) / delta;
    numerical_gradient.z = (z2_energy - z1_energy) / delta;

    EXPECT_NEAR(reference_fractional.x, numerical_gradient.x, tolerance);
    EXPECT_NEAR(reference_fractional.y, numerical_gradient.y, tolerance);
    EXPECT_NEAR(reference_fractional.z, numerical_gradient.z, tolerance);
  }
}



TEST(sixth_derivative_inter_lennard_jones, Test_hessian_cartesian_methane_in_CHA_triclinic_1x1x1)
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

    auto [_, _, reference_cartesian, _, _, _, _] = Interactions::calculateSixthDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB_reference, typeB, frameworkAtoms);

    std::array<std::array<double,3>,3> numerical_second_derivative{};
    for(size_t i = 0; i != 3; ++i)
    {
      // finite difference x
      posB = posB_reference + double3(i==0 ? 0.5 * delta : 0.0, i==1 ? 0.5 * delta : 0.0, i==2 ? 0.5 * delta : 0.0);
      auto [_, gradient_plus, _, _, _, _, _] = Interactions::calculateSixthDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

      posB = posB_reference - double3(i==0 ? 0.5 * delta : 0.0, i==1 ? 0.5 * delta : 0.0, i==2 ? 0.5 * delta : 0.0);
      auto [_, gradient_minus, _, _, _, _, _] = Interactions::calculateSixthDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);


      for(size_t j = 0; j != 3; ++j)
      {
        numerical_second_derivative[i][j] = (gradient_plus[j] - gradient_minus[j]) / delta;
      }
    }

    for(size_t i = 0; i != 3; ++i)
    {
      for(size_t j = 0; j != 3; ++j)
      {
        EXPECT_NEAR(reference_cartesian[i][j], numerical_second_derivative[i][j], tolerance);
      }
    }

  }
}

TEST(sixth_derivative_inter_lennard_jones, Test_third_derivative_cartesian_methane_in_CHA_triclinic_1x1x1)
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

    auto [_, _, _, reference_cartesian, _, _, _] = Interactions::calculateSixthDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB_reference, typeB, frameworkAtoms);

    std::array<std::array<std::array<double,3>,3>,3> numerical_third_derivative{};
    for(size_t i = 0; i != 3; ++i)
    {
      // finite difference x
      posB = posB_reference + double3(i==0 ? 0.5 * delta : 0.0, i==1 ? 0.5 * delta : 0.0, i==2 ? 0.5 * delta : 0.0);
      auto [_, _, hessian_plus, _, _, _, _] = Interactions::calculateSixthDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

      posB = posB_reference - double3(i==0 ? 0.5 * delta : 0.0, i==1 ? 0.5 * delta : 0.0, i==2 ? 0.5 * delta : 0.0);
      auto [_, _, hessian_minus, _, _, _, _] = Interactions::calculateSixthDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);


      for(size_t j = 0; j != 3; ++j)
      {
        for(size_t k = 0; k != 3; ++k)
        {
          numerical_third_derivative[i][j][k] = (hessian_plus[j][k] - hessian_minus[j][k]) / delta;
        }
      }
    }

    for(size_t i = 0; i != 3; ++i)
    {
      for(size_t j = 0; j != 3; ++j)
      {
        for(size_t k = 0; k != 3; ++k)
        {
          EXPECT_NEAR(reference_cartesian[i][j][k], numerical_third_derivative[i][j][k], tolerance);
        }
      }
    }

  }
}

TEST(sixth_derivative_inter_lennard_jones, Test_fourth_derivative_cartesian_methane_in_CHA_triclinic_1x1x1)
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

    auto [_, _, _, _, reference_cartesian, _, _] = Interactions::calculateSixthDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB_reference, typeB, frameworkAtoms);

    std::array<std::array<std::array<std::array<double,3>,3>,3>,3> numerical_fourth_derivative{};
    for(size_t i = 0; i != 3; ++i)
    {
      // finite difference x
      posB = posB_reference + double3(i==0 ? 0.5 * delta : 0.0, i==1 ? 0.5 * delta : 0.0, i==2 ? 0.5 * delta : 0.0);
      auto [_, _, _, third_plus, _, _, _] = Interactions::calculateSixthDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

      posB = posB_reference - double3(i==0 ? 0.5 * delta : 0.0, i==1 ? 0.5 * delta : 0.0, i==2 ? 0.5 * delta : 0.0);
      auto [_, _, _, third_minus, _, _, _] = Interactions::calculateSixthDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);


      for(size_t j = 0; j != 3; ++j)
      {
        for(size_t k = 0; k != 3; ++k)
        {
          for(size_t l = 0; l != 3; ++l)
          {
            numerical_fourth_derivative[i][j][k][l] = (third_plus[j][k][l] - third_minus[j][k][l]) / delta;
          }
        }
      }
    }

    for(size_t i = 0; i != 3; ++i)
    {
      for(size_t j = 0; j != 3; ++j)
      {
        for(size_t k = 0; k != 3; ++k)
        {
          for(size_t l = 0; l != 3; ++l)
          {
            EXPECT_NEAR(reference_cartesian[i][j][k][l], numerical_fourth_derivative[i][j][k][l], tolerance) << std::format("{} {} {} {}\n",i,j,k,l);
          }
        }
      }
    }

  }
}

TEST(sixth_derivative_inter_lennard_jones, Test_fifth_derivative_cartesian_methane_in_CHA_triclinic_1x1x1)
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

    auto [_, _, _, _, _, reference_cartesian,  _] = Interactions::calculateSixthDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB_reference, typeB, frameworkAtoms);

    std::array<std::array<std::array<std::array<std::array<double,3>,3>,3>,3>,3> numerical_fifth_derivative{};
    for(size_t i = 0; i != 3; ++i)
    {
      // finite difference x
      posB = posB_reference + double3(i==0 ? 0.5 * delta : 0.0, i==1 ? 0.5 * delta : 0.0, i==2 ? 0.5 * delta : 0.0);
      auto [_, _, _, _, fourth_plus, _, _] = Interactions::calculateSixthDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

      posB = posB_reference - double3(i==0 ? 0.5 * delta : 0.0, i==1 ? 0.5 * delta : 0.0, i==2 ? 0.5 * delta : 0.0);
      auto [_, _, _, _, fourth_minus, _, _] = Interactions::calculateSixthDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);


      for(size_t j = 0; j != 3; ++j)
      {
        for(size_t k = 0; k != 3; ++k)
        {
          for(size_t l = 0; l != 3; ++l)
          {
            for(size_t m = 0; m != 3; ++m)
            {
              numerical_fifth_derivative[i][j][k][l][m] = (fourth_plus[j][k][l][m] - fourth_minus[j][k][l][m]) / delta;
            }
          }
        }
      }
    }

    for(size_t i = 0; i != 3; ++i)
    {
      for(size_t j = 0; j != 3; ++j)
      {
        for(size_t k = 0; k != 3; ++k)
        {
          for(size_t l = 0; l != 3; ++l)
          {
            for(size_t m = 0; m != 3; ++m)
            {
              EXPECT_NEAR(reference_cartesian[i][j][k][l][m], numerical_fifth_derivative[i][j][k][l][m], tolerance) << std::format("{} {} {} {} {}\n",i,j,k,l,m);
            }
          }
        }
      }
    }

  }
}

TEST(sixth_derivative_inter_lennard_jones, Test_sixth_derivative_cartesian_methane_in_CHA_triclinic_1x1x1)
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
      ForceField::MixingRule::Lorentz_Berthelot, 21.8, 21.8, 21.8, false, false, true);

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

  double delta = 1e-4;
  double tolerance = 1e-4;
  double3x3 numerical_hessian;
  size_t typeB = 2;
  double3 posB_reference = double3(5.1, 4.25, 5.4);

  {
    double3 posB;

    auto [_, _, _, _, _, _, reference_cartesian] = Interactions::calculateSixthDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB_reference, typeB, frameworkAtoms);

    std::array<std::array<std::array<std::array<std::array<std::array<double,3>,3>,3>,3>,3>,3> numerical_sixth_derivative{};
    for(size_t i = 0; i != 3; ++i)
    {
      // finite difference x
      posB = posB_reference + double3(i==0 ? 0.5 * delta : 0.0, i==1 ? 0.5 * delta : 0.0, i==2 ? 0.5 * delta : 0.0);
      auto [_, _, _, _, _, fifth_plus, _] = Interactions::calculateSixthDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);

      posB = posB_reference - double3(i==0 ? 0.5 * delta : 0.0, i==1 ? 0.5 * delta : 0.0, i==2 ? 0.5 * delta : 0.0);
      auto [_, _, _, _, _, fifth_minus, _] = Interactions::calculateSixthDerivativeAtPositionVDW(system.forceField, system.simulationBox, posB, typeB, frameworkAtoms);


      for(size_t j = 0; j != 3; ++j)
      {
        for(size_t k = 0; k != 3; ++k)
        {
          for(size_t l = 0; l != 3; ++l)
          {
            for(size_t m = 0; m != 3; ++m)
            {
              for(size_t n = 0; n != 3; ++n)
              {
                numerical_sixth_derivative[i][j][k][l][m][n] = (fifth_plus[j][k][l][m][n] - fifth_minus[j][k][l][m][n]) / delta;
              }
            }
          }
        }
      }
    }

    for(size_t i = 0; i != 3; ++i)
    {
      for(size_t j = 0; j != 3; ++j)
      {
        for(size_t k = 0; k != 3; ++k)
        {
          for(size_t l = 0; l != 3; ++l)
          {
            for(size_t m = 0; m != 3; ++m)
            {
              for(size_t n = 0; n != 3; ++n)
              {
                EXPECT_NEAR(reference_cartesian[i][j][k][l][m][n], numerical_sixth_derivative[i][j][k][l][m][n], tolerance) << std::format("{} {} {} {} {} {}\n",i,j,k,l,m,n);
              }
            }
          }
        }
      }
    }

  }
}
