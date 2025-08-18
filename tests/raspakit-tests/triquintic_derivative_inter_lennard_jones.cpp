#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <cstddef>
#include <format>
#include <numbers>
#include <print>
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
import hessian_factor;
import triquintic_derivative_factor;
import running_energy;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_framework_molecule_grid;
import interactions_ewald;
import energy_status;

// Cartesian derivatives
// =====================

TEST(sixth_derivative_inter_lennard_jones, Test_gradient_cartesian_methane_in_CHA_triclinic_1x1x1)
{
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);
  Framework f = TestFactories::makeCHA(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeMethane(forceField, 0);
  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {}, {0}, 5);

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  double delta = 1e-5;
  double tolerance = 1e-4;
  double3 numerical_gradient;
  size_t typeB = 2;
  double3 posB = double3(5.0, 5.0, 5.0);

  {
    double3 s = system.simulationBox.inverseCell * posB;
    auto [e, reference_cartesian, d2, d3, d4, d5, d6] = Interactions::calculateTriquinticDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);
    double3 reference_fractional = system.simulationBox.cell.transpose() *
                                   double3(reference_cartesian[0], reference_cartesian[1], reference_cartesian[2]);

    // finite difference x
    posB = system.simulationBox.cell * double3(s.x + 0.5 * delta, s.y, s.z);
    auto [x2_energy, d1x2, d2x2, d3x2, d4x2, d5x2, d6x2] = Interactions::calculateTriquinticDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x - 0.5 * delta, s.y, s.z);
    auto [x1_energy, d1x1, d2x1, d3x1, d4x1, d5x1, d6x1] = Interactions::calculateTriquinticDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    // finite difference y
    posB = system.simulationBox.cell * double3(s.x, s.y + 0.5 * delta, s.z);
    auto [y2_energy, d1y2, d2y2, d3y2, d4y2, d5y2, d6y2] = Interactions::calculateTriquinticDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y - 0.5 * delta, s.z);
    auto [y1_energy, d1y1, d2y1, d3y1, d4y1, d5y1, d6y1] = Interactions::calculateTriquinticDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    // finite difference z
    posB = system.simulationBox.cell * double3(s.x, s.y, s.z + 0.5 * delta);
    auto [z2_energy, d1z2, d2z2, d3z2, d4z2, d5z2, d6z2] = Interactions::calculateTriquinticDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z - 0.5 * delta);
    auto [z1_energy, d1z1, d2z1, d3z1, d4z1, d5z1, d6z1] = Interactions::calculateTriquinticDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

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
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);
  Framework f = TestFactories::makeCHA(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeMethane(forceField, 0);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {}, {0}, 5);

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  double delta = 1e-5;
  double tolerance = 1e-4;
  double3x3 numerical_hessian;
  size_t typeB = 2;
  double3 posB_reference = double3(5.0, 5.0, 5.0);

  {
    double3 posB;

    auto [e, d1, reference_cartesian, d3, d4, d5, d6] = Interactions::calculateTriquinticDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB_reference, typeB,
        frameworkAtoms);

    std::array<std::array<double, 3>, 3> numerical_second_derivative{};
    for (size_t i = 0; i != 3; ++i)
    {
      // finite difference x
      posB =
          posB_reference + double3(i == 0 ? 0.5 * delta : 0.0, i == 1 ? 0.5 * delta : 0.0, i == 2 ? 0.5 * delta : 0.0);
      auto [ep, gradient_plus, d2p, d3p, d4p, d5p, d6p] = Interactions::calculateTriquinticDerivativeAtPosition(
          ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
          frameworkAtoms);

      posB =
          posB_reference - double3(i == 0 ? 0.5 * delta : 0.0, i == 1 ? 0.5 * delta : 0.0, i == 2 ? 0.5 * delta : 0.0);
      auto [em, gradient_minus, d2m, d3m, d4m, d5m, d6m] = Interactions::calculateTriquinticDerivativeAtPosition(
          ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
          frameworkAtoms);

      for (size_t j = 0; j != 3; ++j)
      {
        numerical_second_derivative[i][j] = (gradient_plus[j] - gradient_minus[j]) / delta;
      }
    }

    for (size_t i = 0; i != 3; ++i)
    {
      for (size_t j = 0; j != 3; ++j)
      {
        EXPECT_NEAR(reference_cartesian[i][j], numerical_second_derivative[i][j], tolerance);
      }
    }
  }
}

TEST(sixth_derivative_inter_lennard_jones, Test_third_derivative_cartesian_methane_in_CHA_triclinic_1x1x1)
{
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);
  Framework f = TestFactories::makeCHA(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeMethane(forceField, 0);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {}, {0}, 5);

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  double delta = 1e-5;
  double tolerance = 1e-4;
  double3x3 numerical_hessian;
  size_t typeB = 2;
  double3 posB_reference = double3(5.0, 5.0, 5.0);

  {
    double3 posB;

    auto [e, d1, d2, reference_cartesian, d4, d5, d6] = Interactions::calculateTriquinticDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB_reference, typeB,
        frameworkAtoms);

    std::array<std::array<std::array<double, 3>, 3>, 3> numerical_third_derivative{};
    for (size_t i = 0; i != 3; ++i)
    {
      // finite difference x
      posB =
          posB_reference + double3(i == 0 ? 0.5 * delta : 0.0, i == 1 ? 0.5 * delta : 0.0, i == 2 ? 0.5 * delta : 0.0);
      auto [ep, d1p, hessian_plus, d3p, d4p, d5p, d6p] = Interactions::calculateTriquinticDerivativeAtPosition(
          ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
          frameworkAtoms);

      posB =
          posB_reference - double3(i == 0 ? 0.5 * delta : 0.0, i == 1 ? 0.5 * delta : 0.0, i == 2 ? 0.5 * delta : 0.0);
      auto [em, d1m, hessian_minus, d3m, d4m, d5m, d6m] = Interactions::calculateTriquinticDerivativeAtPosition(
          ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
          frameworkAtoms);

      for (size_t j = 0; j != 3; ++j)
      {
        for (size_t k = 0; k != 3; ++k)
        {
          numerical_third_derivative[i][j][k] = (hessian_plus[j][k] - hessian_minus[j][k]) / delta;
        }
      }
    }

    for (size_t i = 0; i != 3; ++i)
    {
      for (size_t j = 0; j != 3; ++j)
      {
        for (size_t k = 0; k != 3; ++k)
        {
          EXPECT_NEAR(reference_cartesian[i][j][k], numerical_third_derivative[i][j][k], tolerance);
        }
      }
    }
  }
}

TEST(sixth_derivative_inter_lennard_jones, Test_fourth_derivative_cartesian_methane_in_CHA_triclinic_1x1x1)
{
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);
  Framework f = TestFactories::makeCHA(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeMethane(forceField, 0);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {}, {0}, 5);

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  double delta = 1e-5;
  double tolerance = 1e-4;
  double3x3 numerical_hessian;
  size_t typeB = 2;
  double3 posB_reference = double3(5.0, 5.0, 5.0);

  {
    double3 posB;

    auto [e, d1, d2, d3, reference_cartesian, d5, d6] = Interactions::calculateTriquinticDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB_reference, typeB,
        frameworkAtoms);

    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3> numerical_fourth_derivative{};
    for (size_t i = 0; i != 3; ++i)
    {
      // finite difference x
      posB =
          posB_reference + double3(i == 0 ? 0.5 * delta : 0.0, i == 1 ? 0.5 * delta : 0.0, i == 2 ? 0.5 * delta : 0.0);
      auto [ep, d1p, d2p, third_plus, d4p, d5p, d6p] = Interactions::calculateTriquinticDerivativeAtPosition(
          ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
          frameworkAtoms);

      posB =
          posB_reference - double3(i == 0 ? 0.5 * delta : 0.0, i == 1 ? 0.5 * delta : 0.0, i == 2 ? 0.5 * delta : 0.0);
      auto [em, d1m, d2m, third_minus, d4m, d5m, d6m] = Interactions::calculateTriquinticDerivativeAtPosition(
          ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
          frameworkAtoms);

      for (size_t j = 0; j != 3; ++j)
      {
        for (size_t k = 0; k != 3; ++k)
        {
          for (size_t l = 0; l != 3; ++l)
          {
            numerical_fourth_derivative[i][j][k][l] = (third_plus[j][k][l] - third_minus[j][k][l]) / delta;
          }
        }
      }
    }

    for (size_t i = 0; i != 3; ++i)
    {
      for (size_t j = 0; j != 3; ++j)
      {
        for (size_t k = 0; k != 3; ++k)
        {
          for (size_t l = 0; l != 3; ++l)
          {
            EXPECT_NEAR(reference_cartesian[i][j][k][l], numerical_fourth_derivative[i][j][k][l], tolerance)
                << std::format("{} {} {} {}\n", i, j, k, l);
          }
        }
      }
    }
  }
}

TEST(sixth_derivative_inter_lennard_jones, Test_fifth_derivative_cartesian_methane_in_CHA_triclinic_1x1x1)
{
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);
  Framework f = TestFactories::makeCHA(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeMethane(forceField, 0);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {}, {0}, 5);

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  double delta = 1e-5;
  double tolerance = 1e-4;
  double3x3 numerical_hessian;
  size_t typeB = 2;
  double3 posB_reference = double3(5.0, 5.0, 5.0);

  {
    double3 posB;

    auto [e, d1, d2, d3, d4, reference_cartesian, d6] = Interactions::calculateTriquinticDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB_reference, typeB,
        frameworkAtoms);

    std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3> numerical_fifth_derivative{};
    for (size_t i = 0; i != 3; ++i)
    {
      // finite difference x
      posB =
          posB_reference + double3(i == 0 ? 0.5 * delta : 0.0, i == 1 ? 0.5 * delta : 0.0, i == 2 ? 0.5 * delta : 0.0);
      auto [ep, d1p, d2p, d3p, fourth_plus, d5p, d6p] = Interactions::calculateTriquinticDerivativeAtPosition(
          ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
          frameworkAtoms);

      posB =
          posB_reference - double3(i == 0 ? 0.5 * delta : 0.0, i == 1 ? 0.5 * delta : 0.0, i == 2 ? 0.5 * delta : 0.0);
      auto [em, d1m, d2m, d3m, fourth_minus, d5m, d6m] = Interactions::calculateTriquinticDerivativeAtPosition(
          ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
          frameworkAtoms);

      for (size_t j = 0; j != 3; ++j)
      {
        for (size_t k = 0; k != 3; ++k)
        {
          for (size_t l = 0; l != 3; ++l)
          {
            for (size_t m = 0; m != 3; ++m)
            {
              numerical_fifth_derivative[i][j][k][l][m] = (fourth_plus[j][k][l][m] - fourth_minus[j][k][l][m]) / delta;
            }
          }
        }
      }
    }

    for (size_t i = 0; i != 3; ++i)
    {
      for (size_t j = 0; j != 3; ++j)
      {
        for (size_t k = 0; k != 3; ++k)
        {
          for (size_t l = 0; l != 3; ++l)
          {
            for (size_t m = 0; m != 3; ++m)
            {
              EXPECT_NEAR(reference_cartesian[i][j][k][l][m], numerical_fifth_derivative[i][j][k][l][m], tolerance)
                  << std::format("{} {} {} {} {}\n", i, j, k, l, m);
            }
          }
        }
      }
    }
  }
}

TEST(sixth_derivative_inter_lennard_jones, Test_sixth_derivative_cartesian_methane_in_CHA_triclinic_1x1x1)
{
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);
  Framework f = TestFactories::makeCHA(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeMethane(forceField, 0);
  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {}, {0}, 5);

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  double delta = 1e-4;
  double tolerance = 1e-4;
  double3x3 numerical_hessian;
  size_t typeB = 2;
  double3 posB_reference = double3(5.1, 4.25, 5.4);

  {
    double3 posB;

    auto [e, d1, d2, d3, d4, d5, reference_cartesian] = Interactions::calculateTriquinticDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB_reference, typeB,
        frameworkAtoms);

    std::array<std::array<std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>, 3>, 3>
        numerical_sixth_derivative{};
    for (size_t i = 0; i != 3; ++i)
    {
      // finite difference x
      posB =
          posB_reference + double3(i == 0 ? 0.5 * delta : 0.0, i == 1 ? 0.5 * delta : 0.0, i == 2 ? 0.5 * delta : 0.0);
      auto [ep, d1p, d2p, d3p, d4p, fifth_plus, d6p] = Interactions::calculateTriquinticDerivativeAtPosition(
          ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
          frameworkAtoms);

      posB =
          posB_reference - double3(i == 0 ? 0.5 * delta : 0.0, i == 1 ? 0.5 * delta : 0.0, i == 2 ? 0.5 * delta : 0.0);
      auto [em, d1m, d2m, d3m, d4m, fifth_minus, d6m] = Interactions::calculateTriquinticDerivativeAtPosition(
          ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
          frameworkAtoms);

      for (size_t j = 0; j != 3; ++j)
      {
        for (size_t k = 0; k != 3; ++k)
        {
          for (size_t l = 0; l != 3; ++l)
          {
            for (size_t m = 0; m != 3; ++m)
            {
              for (size_t n = 0; n != 3; ++n)
              {
                numerical_sixth_derivative[i][j][k][l][m][n] =
                    (fifth_plus[j][k][l][m][n] - fifth_minus[j][k][l][m][n]) / delta;
              }
            }
          }
        }
      }
    }

    for (size_t i = 0; i != 3; ++i)
    {
      for (size_t j = 0; j != 3; ++j)
      {
        for (size_t k = 0; k != 3; ++k)
        {
          for (size_t l = 0; l != 3; ++l)
          {
            for (size_t m = 0; m != 3; ++m)
            {
              for (size_t n = 0; n != 3; ++n)
              {
                EXPECT_NEAR(reference_cartesian[i][j][k][l][m][n], numerical_sixth_derivative[i][j][k][l][m][n],
                            tolerance)
                    << std::format("{} {} {} {} {} {}\n", i, j, k, l, m, n);
              }
            }
          }
        }
      }
    }
  }
}

// Fractional derivatives for the triquintic interpolation algorithm
// =================================================================

TEST(sixth_derivative_inter_lennard_jones, Test_Txxyyzz_fractional_methane_in_CHA_triclinic_1x1x1)
{
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);
  Framework framework = TestFactories::makeCHA(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeMethane(forceField, 0);
  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {framework}, {c}, {}, {0}, 5);

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  size_t typeB = 2;
  double3 posB_reference = double3(5.1, 4.25, 5.4);
  double3 posB;

  // first derivatives
  {
    double delta = 1e-6;
    double tolerance = 1e-3;

    std::array<double, 27> analytical = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB_reference, typeB,
        framework.simulationBox, frameworkAtoms);

    double3 s = system.simulationBox.inverseCell * posB_reference;

    posB = system.simulationBox.cell * double3(s.x + 0.5 * delta, s.y, s.z);
    std::array<double, 27> x_right = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x - 0.5 * delta, s.y, s.z);
    std::array<double, 27> x_left = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y + 0.5 * delta, s.z);
    std::array<double, 27> y_right = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y - 0.5 * delta, s.z);
    std::array<double, 27> y_left = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z + 0.5 * delta);
    std::array<double, 27> z_right = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z - 0.5 * delta);
    std::array<double, 27> z_left = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    double numerical_gradient_x = (x_right[0] - x_left[0]) / delta;
    double numerical_gradient_y = (y_right[0] - y_left[0]) / delta;
    double numerical_gradient_z = (z_right[0] - z_left[0]) / delta;

    EXPECT_NEAR(analytical[1], numerical_gradient_x, tolerance);  // dU / dx
    EXPECT_NEAR(analytical[2], numerical_gradient_y, tolerance);  // dU / dy
    EXPECT_NEAR(analytical[3], numerical_gradient_z, tolerance);  // dU / dz
  }

  // second derivatives
  {
    double delta = 1e-6;
    double tolerance = 1e-3;

    std::array<double, 27> analytical = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB_reference, typeB,
        framework.simulationBox, frameworkAtoms);

    double3 s = system.simulationBox.inverseCell * posB_reference;

    posB = system.simulationBox.cell * double3(s.x + 0.5 * delta, s.y, s.z);
    std::array<double, 27> x_right = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x - 0.5 * delta, s.y, s.z);
    std::array<double, 27> x_left = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y + 0.5 * delta, s.z);
    std::array<double, 27> y_right = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y - 0.5 * delta, s.z);
    std::array<double, 27> y_left = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z + 0.5 * delta);
    std::array<double, 27> z_right = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z - 0.5 * delta);
    std::array<double, 27> z_left = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    double numerical_gradient_xx = (x_right[1] - x_left[1]) / delta;
    double numerical_gradient_xy = (y_right[1] - y_left[1]) / delta;
    double numerical_gradient_xz = (z_right[1] - z_left[1]) / delta;
    double numerical_gradient_yy = (y_right[2] - y_left[2]) / delta;
    double numerical_gradient_yz = (y_right[3] - y_left[3]) / delta;
    double numerical_gradient_zz = (z_right[3] - z_left[3]) / delta;

    EXPECT_NEAR(analytical[4], numerical_gradient_xx, tolerance);  // d^2U / dx^2
    EXPECT_NEAR(analytical[5], numerical_gradient_xy, tolerance);  // d^2U / dx dy
    EXPECT_NEAR(analytical[6], numerical_gradient_xz, tolerance);  // d^2U / dx dz
    EXPECT_NEAR(analytical[7], numerical_gradient_yy, tolerance);  // d^2U / dy^2
    EXPECT_NEAR(analytical[8], numerical_gradient_yz, tolerance);  // d^2U / dy dz
    EXPECT_NEAR(analytical[9], numerical_gradient_zz, tolerance);  // d^2U / dz
  }

  // third derivatives
  {
    double delta = 1e-6;
    double tolerance = 1e-3;

    std::array<double, 27> analytical = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB_reference, typeB,
        framework.simulationBox, frameworkAtoms);

    double3 s = system.simulationBox.inverseCell * posB_reference;

    posB = system.simulationBox.cell * double3(s.x, s.y + 0.5 * delta, s.z);
    std::array<double, 27> y_right = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y - 0.5 * delta, s.z);
    std::array<double, 27> y_left = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z + 0.5 * delta);
    std::array<double, 27> z_right = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z - 0.5 * delta);
    std::array<double, 27> z_left = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    double numerical_gradient_xxy = (y_right[4] - y_left[4]) / delta;
    double numerical_gradient_xxz = (z_right[4] - z_left[4]) / delta;
    double numerical_gradient_xyy = (y_right[5] - y_left[5]) / delta;
    double numerical_gradient_xyz = (z_right[5] - z_left[5]) / delta;
    double numerical_gradient_yyz = (z_right[7] - z_left[7]) / delta;
    double numerical_gradient_xzz = (z_right[6] - z_left[6]) / delta;
    double numerical_gradient_yzz = (z_right[8] - z_left[8]) / delta;

    EXPECT_NEAR(analytical[10], numerical_gradient_xxy, tolerance);  // d^3U / dx dx dy
    EXPECT_NEAR(analytical[11], numerical_gradient_xxz, tolerance);  // d^3U / dx dx dy
    EXPECT_NEAR(analytical[12], numerical_gradient_xyy, tolerance);  // d^3U / dx dy dy
    EXPECT_NEAR(analytical[13], numerical_gradient_xyz, tolerance);  // d^3U / dx dy dz
    EXPECT_NEAR(analytical[14], numerical_gradient_yyz, tolerance);  // d^3U / dy dy dz
    EXPECT_NEAR(analytical[15], numerical_gradient_xzz, tolerance);  // d^3U / dx dz dz
    EXPECT_NEAR(analytical[16], numerical_gradient_yzz, tolerance);  // d^3U / dy dz dz
  }

  // fourth derivatives
  {
    double delta = 1e-6;
    double tolerance = 1e-3;

    std::array<double, 27> analytical = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB_reference, typeB,
        framework.simulationBox, frameworkAtoms);

    double3 s = system.simulationBox.inverseCell * posB_reference;

    posB = system.simulationBox.cell * double3(s.x + 0.5 * delta, s.y, s.z);
    std::array<double, 27> x_right = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x - 0.5 * delta, s.y, s.z);
    std::array<double, 27> x_left = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y + 0.5 * delta, s.z);
    std::array<double, 27> y_right = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y - 0.5 * delta, s.z);
    std::array<double, 27> y_left = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z + 0.5 * delta);
    std::array<double, 27> z_right = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z - 0.5 * delta);
    std::array<double, 27> z_left = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    double numerical_gradient_xxyy = (y_right[10] - y_left[10]) / delta;
    double numerical_gradient_xxzz = (x_right[15] - x_left[15]) / delta;
    double numerical_gradient_yyzz = (z_right[14] - z_left[14]) / delta;
    double numerical_gradient_xxyz = (x_right[13] - x_left[13]) / delta;
    double numerical_gradient_xyyz = (y_right[13] - y_left[13]) / delta;
    double numerical_gradient_xyzz = (z_right[13] - z_left[13]) / delta;

    EXPECT_NEAR(analytical[17], numerical_gradient_xxyy, tolerance);  // d^4U / dx dx dy dy
    EXPECT_NEAR(analytical[18], numerical_gradient_xxzz, tolerance);  // d^4U / dx dx dz dz
    EXPECT_NEAR(analytical[19], numerical_gradient_yyzz, tolerance);  // d^4U / dy dy dz dz
    EXPECT_NEAR(analytical[20], numerical_gradient_xxyz, tolerance);  // d^4U / dx dx dy dz
    EXPECT_NEAR(analytical[21], numerical_gradient_xyyz, tolerance);  // d^4U / dx dy dy dz
    EXPECT_NEAR(analytical[22], numerical_gradient_xyzz, tolerance);  // d^4U / dx dy dz dz
  }

  // fifth derivatives
  {
    double delta = 1e-6;
    double tolerance = 3e-1;

    std::array<double, 27> analytical = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB_reference, typeB,
        framework.simulationBox, frameworkAtoms);

    double3 s = system.simulationBox.inverseCell * posB_reference;

    posB = system.simulationBox.cell * double3(s.x + 0.5 * delta, s.y, s.z);
    std::array<double, 27> x_right = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x - 0.5 * delta, s.y, s.z);
    std::array<double, 27> x_left = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y + 0.5 * delta, s.z);
    std::array<double, 27> y_right = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y - 0.5 * delta, s.z);
    std::array<double, 27> y_left = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z + 0.5 * delta);
    std::array<double, 27> z_right = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z - 0.5 * delta);
    std::array<double, 27> z_left = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    double numerical_gradient_xxyyz = (z_right[17] - z_left[17]) / delta;
    double numerical_gradient_xxyzz = (y_right[18] - y_left[18]) / delta;
    double numerical_gradient_xyyzz = (x_right[19] - x_left[19]) / delta;

    EXPECT_NEAR(analytical[23], numerical_gradient_xxyyz, tolerance);  // d^5U / dx dx dy dy dz
    EXPECT_NEAR(analytical[24], numerical_gradient_xxyzz, tolerance);  // d^5U / dx dx dy dz dz
    EXPECT_NEAR(analytical[25], numerical_gradient_xyyzz, tolerance);  // d^5U / dx dy dy dz dz
  }

  // sixth derivative
  {
    double delta = 1e-5;
    double tolerance = 1.0;

    std::array<double, 27> analytical = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB_reference, typeB,
        framework.simulationBox, frameworkAtoms);

    double3 s = system.simulationBox.inverseCell * posB_reference;

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z + 0.5 * delta);
    std::array<double, 27> z_right = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z - 0.5 * delta);
    std::array<double, 27> z_left = Interactions::calculateTriquinticFractionalAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        framework.simulationBox, frameworkAtoms);

    double numerical_gradient_xxyyzz = (z_right[23] - z_left[23]) / delta;

    EXPECT_NEAR(analytical[26], numerical_gradient_xxyyzz, tolerance);  // d^6U / dx dx dy dy dz dz
  }
}

// Convert Fractional derivatives back to Cartesian
// ================================================

TEST(sixth_derivative_inter_lennard_jones, Test_fractional_to_Cartesian_methane_in_CHA_triclinic_1x1x1)
{
  double tolerance = 1e-4;

  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);
  Framework framework = TestFactories::makeCHA(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeMethane(forceField, 0);
  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {framework}, {c}, {}, {0}, 5);

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  size_t typeB = 2;
  double3 posB_reference = double3(5.1, 4.25, 5.4);

  std::array<double, 27> analyticalCartesian = Interactions::calculateTriquinticCartesianAtPosition(
      ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB_reference, typeB,
      frameworkAtoms);

  std::array<double, 27> analyticalFractional = Interactions::calculateTriquinticFractionalAtPosition(
      ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB_reference, typeB,
      framework.simulationBox, frameworkAtoms);

  double3 first_derivative_fractional{analyticalFractional[1], analyticalFractional[2], analyticalFractional[3]};
  double3x3 second_derivative_fractional{
      double3{analyticalFractional[4], analyticalFractional[5], analyticalFractional[6]},
      double3{analyticalFractional[5], analyticalFractional[7], analyticalFractional[8]},
      double3{analyticalFractional[6], analyticalFractional[8], analyticalFractional[9]}};

  double3 first_derivative_Cartesian{};
  double3x3 second_derivative_Cartesian{};

  first_derivative_Cartesian = double3x3::transpose(system.simulationBox.inverseCell) * first_derivative_fractional;

  second_derivative_Cartesian = double3x3::transpose(system.simulationBox.inverseCell) * second_derivative_fractional *
                                system.simulationBox.inverseCell;

  EXPECT_NEAR(first_derivative_Cartesian[0], analyticalCartesian[1], tolerance);
  EXPECT_NEAR(first_derivative_Cartesian[1], analyticalCartesian[2], tolerance);
  EXPECT_NEAR(first_derivative_Cartesian[2], analyticalCartesian[3], tolerance);

  EXPECT_NEAR(second_derivative_Cartesian[0][0], analyticalCartesian[4], tolerance);
  EXPECT_NEAR(second_derivative_Cartesian[0][1], analyticalCartesian[5], tolerance);
  EXPECT_NEAR(second_derivative_Cartesian[0][2], analyticalCartesian[6], tolerance);
  EXPECT_NEAR(second_derivative_Cartesian[1][1], analyticalCartesian[7], tolerance);
  EXPECT_NEAR(second_derivative_Cartesian[1][2], analyticalCartesian[8], tolerance);
  EXPECT_NEAR(second_derivative_Cartesian[2][2], analyticalCartesian[9], tolerance);

  first_derivative_Cartesian = double3{};
  second_derivative_Cartesian = double3x3{};
  for (const size_t &p : std::array<size_t, 3>{0, 1, 2})
  {
    first_derivative_Cartesian[0] += system.simulationBox.inverseCell[0][p] * first_derivative_fractional[p];
    first_derivative_Cartesian[1] += system.simulationBox.inverseCell[1][p] * first_derivative_fractional[p];
    first_derivative_Cartesian[2] += system.simulationBox.inverseCell[2][p] * first_derivative_fractional[p];

    for (const size_t &q : std::array<size_t, 3>{0, 1, 2})
    {
      second_derivative_Cartesian[0][0] += system.simulationBox.inverseCell[0][p] *
                                           system.simulationBox.inverseCell[0][q] * second_derivative_fractional[p][q];
      second_derivative_Cartesian[0][1] += system.simulationBox.inverseCell[0][p] *
                                           system.simulationBox.inverseCell[1][q] * second_derivative_fractional[p][q];
      second_derivative_Cartesian[0][2] += system.simulationBox.inverseCell[0][p] *
                                           system.simulationBox.inverseCell[2][q] * second_derivative_fractional[p][q];
      second_derivative_Cartesian[1][1] += system.simulationBox.inverseCell[1][p] *
                                           system.simulationBox.inverseCell[1][q] * second_derivative_fractional[p][q];
      second_derivative_Cartesian[1][2] += system.simulationBox.inverseCell[1][p] *
                                           system.simulationBox.inverseCell[2][q] * second_derivative_fractional[p][q];
      second_derivative_Cartesian[2][2] += system.simulationBox.inverseCell[2][p] *
                                           system.simulationBox.inverseCell[2][q] * second_derivative_fractional[p][q];
    }
  }

  EXPECT_NEAR(first_derivative_Cartesian[0], analyticalCartesian[1], tolerance);
  EXPECT_NEAR(first_derivative_Cartesian[1], analyticalCartesian[2], tolerance);
  EXPECT_NEAR(first_derivative_Cartesian[2], analyticalCartesian[3], tolerance);

  EXPECT_NEAR(second_derivative_Cartesian[0][0], analyticalCartesian[4], tolerance);
  EXPECT_NEAR(second_derivative_Cartesian[0][1], analyticalCartesian[5], tolerance);
  EXPECT_NEAR(second_derivative_Cartesian[0][2], analyticalCartesian[6], tolerance);
  EXPECT_NEAR(second_derivative_Cartesian[1][1], analyticalCartesian[7], tolerance);
  EXPECT_NEAR(second_derivative_Cartesian[1][2], analyticalCartesian[8], tolerance);
  EXPECT_NEAR(second_derivative_Cartesian[2][2], analyticalCartesian[9], tolerance);
}
