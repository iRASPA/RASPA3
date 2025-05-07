#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <cstddef>
#include <numbers>
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
import tricubic_derivative_factor;
import running_energy;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_framework_molecule_grid;
import interactions_ewald;
import energy_status;

TEST(third_derivative_inter_lennard_jones, Test_gradient_cartesian_methane_in_CHA_triclinic_1x1x1)
{
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);
  Framework f = TestFactories::makeCHA(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeMethane(forceField, 0);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {0}, 5);

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  double delta = 1e-5;
  double tolerance = 1e-4;
  double3 numerical_gradient;
  size_t typeB = 2;
  double3 posB = double3(5.0, 5.0, 5.0);

  {
    double3 s = system.simulationBox.inverseCell * posB;
    auto [e, reference_cartesian, d1, d2] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);
    double3 reference_fractional = system.simulationBox.cell.transpose() * reference_cartesian;

    // finite difference x
    posB = system.simulationBox.cell * double3(s.x + 0.5 * delta, s.y, s.z);
    auto [x2_energy, d1x2, d2x2, d3x2] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x - 0.5 * delta, s.y, s.z);
    auto [x1_energy, d1x1, d2x1, d3x1] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    // finite difference y
    posB = system.simulationBox.cell * double3(s.x, s.y + 0.5 * delta, s.z);
    auto [y2_energy, d1y2, d2y2, d3y2] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y - 0.5 * delta, s.z);
    auto [y1_energy, d1y1, d2y1, d3y1] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    // finite difference z
    posB = system.simulationBox.cell * double3(s.x, s.y, s.z + 0.5 * delta);
    auto [z2_energy, d1z2, d2z2, d3z2] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z - 0.5 * delta);
    auto [z1_energy, d1z1, d2z1, d3z1] = Interactions::calculateTricubicDerivativeAtPosition(
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

TEST(third_derivative_inter_lennard_jones, Test_gradient_fractional_methane_in_CHA_triclinic_1x1x1)
{
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);
  Framework f = TestFactories::makeCHA(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeMethane(forceField, 0);
  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {0}, 5);

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  double delta = 1e-5;
  double tolerance = 1e-4;
  double3 numerical_gradient;
  size_t typeB = 2;
  double3 posB = double3(5.0, 5.0, 5.0);

  {
    double3 s = system.simulationBox.inverseCell * posB;
    auto [e, reference_cartesian, d1, d2] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);
    double3 reference_fractional = system.simulationBox.cell.transpose() * reference_cartesian;

    // finite difference x
    posB = system.simulationBox.cell * double3(s.x + 0.5 * delta, s.y, s.z);
    auto [x2_energy, d1x2, d2x2, d3x2] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x - 0.5 * delta, s.y, s.z);
    auto [x1_energy, d1x1, d2x1, d3x1] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    // finite difference y
    posB = system.simulationBox.cell * double3(s.x, s.y + 0.5 * delta, s.z);
    auto [y2_energy, d1y2, d2y2, d3y2] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y - 0.5 * delta, s.z);
    auto [y1_energy, d1y1, d2y1, d3y1] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    // finite difference z
    posB = system.simulationBox.cell * double3(s.x, s.y, s.z + 0.5 * delta);
    auto [z2_energy, d1z2, d2z2, d3z2] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z - 0.5 * delta);
    auto [z1_energy, d1z1, d2z1, d3z1] = Interactions::calculateTricubicDerivativeAtPosition(
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

TEST(third_derivative_inter_lennard_jones, Test_hessian_cartesian_methane_in_CHA_triclinic_1x1x1)
{
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);
  Framework f = TestFactories::makeCHA(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeMethane(forceField, 0);
  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {0}, 5);

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  double delta = 1e-5;
  double tolerance = 1e-4;
  double3x3 numerical_hessian;
  size_t typeB = 2;
  double3 posB_reference = double3(5.0, 5.0, 5.0);

  {
    double3 posB;

    auto [e, d1, reference_cartesian, d3] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB_reference, typeB,
        frameworkAtoms);

    // finite difference x
    posB = double3(posB_reference.x + 0.5 * delta, posB_reference.y, posB_reference.z);
    auto [ex2, x2_gradient, d2x2, d3x2] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    posB = double3(posB_reference.x - 0.5 * delta, posB_reference.y, posB_reference.z);
    auto [ex1, x1_gradient, d2x1, d3x1] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    // finite difference y
    posB = double3(posB_reference.x, posB_reference.y + 0.5 * delta, posB_reference.z);
    auto [ey2, y2_gradient, d2y2, d3y2] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    posB = double3(posB_reference.x, posB_reference.y - 0.5 * delta, posB_reference.z);
    auto [ey1, y1_gradient, d2y1, d3y1] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    // finite difference z
    posB = double3(posB_reference.x, posB_reference.y, posB_reference.z + 0.5 * delta);
    auto [ez2, z2_gradient, d2z2, d3z2] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    posB = double3(posB_reference.x, posB_reference.y, posB_reference.z - 0.5 * delta);
    auto [ez1, z1_gradient, d2z1, d3z1] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    numerical_hessian.ax = (x2_gradient[0] - x1_gradient[0]) / delta;
    numerical_hessian.ay = (x2_gradient[1] - x1_gradient[1]) / delta;
    numerical_hessian.az = (x2_gradient[2] - x1_gradient[2]) / delta;

    numerical_hessian.bx = (y2_gradient[0] - y1_gradient[0]) / delta;
    numerical_hessian.by = (y2_gradient[1] - y1_gradient[1]) / delta;
    numerical_hessian.bz = (y2_gradient[2] - y1_gradient[2]) / delta;

    numerical_hessian.cx = (z2_gradient[0] - z1_gradient[0]) / delta;
    numerical_hessian.cy = (z2_gradient[1] - z1_gradient[1]) / delta;
    numerical_hessian.cz = (z2_gradient[2] - z1_gradient[2]) / delta;

    EXPECT_NEAR(reference_cartesian[0][0], numerical_hessian.ax, tolerance);
    EXPECT_NEAR(reference_cartesian[0][1], numerical_hessian.ay, tolerance);
    EXPECT_NEAR(reference_cartesian[0][2], numerical_hessian.az, tolerance);

    EXPECT_NEAR(reference_cartesian[1][0], numerical_hessian.bx, tolerance);
    EXPECT_NEAR(reference_cartesian[1][1], numerical_hessian.by, tolerance);
    EXPECT_NEAR(reference_cartesian[1][2], numerical_hessian.bz, tolerance);

    EXPECT_NEAR(reference_cartesian[2][0], numerical_hessian.cx, tolerance);
    EXPECT_NEAR(reference_cartesian[2][1], numerical_hessian.cy, tolerance);
    EXPECT_NEAR(reference_cartesian[2][2], numerical_hessian.cz, tolerance);
  }
}

TEST(third_derivative_inter_lennard_jones, Test_hessian_fractional_methane_in_CHA_triclinic_1x1x1)
{
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);
  Framework f = TestFactories::makeCHA(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeMethane(forceField, 0);
  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {0}, 5);

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  double delta = 1e-5;
  double tolerance = 1e-4;
  double3x3 numerical_hessian;
  size_t typeB = 2;
  double3 posB = double3(5.0, 5.0, 5.0);

  {
    double3 s = system.simulationBox.inverseCell * posB;
    auto [e, d1, reference_cartesian, d3] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);
    double3x3 reference_fractional =
        system.simulationBox.cell.transpose() * reference_cartesian * system.simulationBox.cell;

    // finite difference x
    posB = system.simulationBox.cell * double3(s.x + 0.5 * delta, s.y, s.z);
    auto [ex2, x2_gradient, d1x2, d2x2] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x - 0.5 * delta, s.y, s.z);
    auto [ex1, x1_gradient, d1x1, d2x1] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    // finite difference y
    posB = system.simulationBox.cell * double3(s.x, s.y + 0.5 * delta, s.z);
    auto [ey2, y2_gradient, d1y2, d2y2] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y - 0.5 * delta, s.z);
    auto [ey1, y1_gradient, d1y1, d2y1] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    // finite difference z
    posB = system.simulationBox.cell * double3(s.x, s.y, s.z + 0.5 * delta);
    auto [ez2, z2_gradient, d1z2, d2z2] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z - 0.5 * delta);
    auto [ez1, z1_gradient, d1z1, d2z1] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

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
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);
  Framework f = TestFactories::makeCHA(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeMethane(forceField, 0);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {0}, 5);

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  double delta = 1e-5;
  double tolerance = 1e-4;
  size_t typeB = 2;
  double3 posB_reference = double3(5.0, 5.0, 5.0);

  {
    double3 posB;

    auto [e, d1, d2, reference_cartesian] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB_reference, typeB,
        frameworkAtoms);

    posB = double3(posB_reference.x + 0.5 * delta, posB_reference.y, posB_reference.z);
    auto [ex2, d1x2, x2_hessian, d3x2] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    posB = double3(posB_reference.x - 0.5 * delta, posB_reference.y, posB_reference.z);
    auto [ex1, d1x1, x1_hessian, d3x1] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    double numerical_third_derivative_m111 = (x2_hessian[0][0] - x1_hessian[0][0]) / delta;
    double numerical_third_derivative_m112 = (x2_hessian[0][1] - x1_hessian[0][1]) / delta;
    double numerical_third_derivative_m113 = (x2_hessian[0][2] - x1_hessian[0][2]) / delta;
    double numerical_third_derivative_m121 = (x2_hessian[1][0] - x1_hessian[1][0]) / delta;
    double numerical_third_derivative_m122 = (x2_hessian[1][1] - x1_hessian[1][1]) / delta;
    double numerical_third_derivative_m123 = (x2_hessian[1][2] - x1_hessian[1][2]) / delta;
    double numerical_third_derivative_m131 = (x2_hessian[2][0] - x1_hessian[2][0]) / delta;
    double numerical_third_derivative_m132 = (x2_hessian[2][1] - x1_hessian[2][1]) / delta;
    double numerical_third_derivative_m133 = (x2_hessian[2][2] - x1_hessian[2][2]) / delta;

    EXPECT_NEAR(reference_cartesian[0][0][0], numerical_third_derivative_m111, tolerance);
    EXPECT_NEAR(reference_cartesian[0][0][1], numerical_third_derivative_m112, tolerance);
    EXPECT_NEAR(reference_cartesian[0][0][2], numerical_third_derivative_m113, tolerance);
    EXPECT_NEAR(reference_cartesian[0][1][0], numerical_third_derivative_m121, tolerance);
    EXPECT_NEAR(reference_cartesian[0][1][1], numerical_third_derivative_m122, tolerance);
    EXPECT_NEAR(reference_cartesian[0][1][2], numerical_third_derivative_m123, tolerance);
    EXPECT_NEAR(reference_cartesian[0][2][0], numerical_third_derivative_m131, tolerance);
    EXPECT_NEAR(reference_cartesian[0][2][1], numerical_third_derivative_m132, tolerance);
    EXPECT_NEAR(reference_cartesian[0][2][2], numerical_third_derivative_m133, tolerance);

    posB = double3(posB_reference.x, posB_reference.y + 0.5 * delta, posB_reference.z);
    auto [ey2, d1y2, y2_hessian, d3y2] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    posB = double3(posB_reference.x, posB_reference.y - 0.5 * delta, posB_reference.z);
    auto [ey1, d1y1, y1_hessian, d3y1] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    double numerical_third_derivative_m211 = (y2_hessian[0][0] - y1_hessian[0][0]) / delta;
    double numerical_third_derivative_m212 = (y2_hessian[0][1] - y1_hessian[0][1]) / delta;
    double numerical_third_derivative_m213 = (y2_hessian[0][2] - y1_hessian[0][2]) / delta;
    double numerical_third_derivative_m221 = (y2_hessian[1][0] - y1_hessian[1][0]) / delta;
    double numerical_third_derivative_m222 = (y2_hessian[1][1] - y1_hessian[1][1]) / delta;
    double numerical_third_derivative_m223 = (y2_hessian[1][2] - y1_hessian[1][2]) / delta;
    double numerical_third_derivative_m231 = (y2_hessian[2][0] - y1_hessian[2][0]) / delta;
    double numerical_third_derivative_m232 = (y2_hessian[2][1] - y1_hessian[2][1]) / delta;
    double numerical_third_derivative_m233 = (y2_hessian[2][2] - y1_hessian[2][2]) / delta;

    EXPECT_NEAR(reference_cartesian[1][0][0], numerical_third_derivative_m211, tolerance);
    EXPECT_NEAR(reference_cartesian[1][0][1], numerical_third_derivative_m212, tolerance);
    EXPECT_NEAR(reference_cartesian[1][0][2], numerical_third_derivative_m213, tolerance);
    EXPECT_NEAR(reference_cartesian[1][1][0], numerical_third_derivative_m221, tolerance);
    EXPECT_NEAR(reference_cartesian[1][1][1], numerical_third_derivative_m222, tolerance);
    EXPECT_NEAR(reference_cartesian[1][1][2], numerical_third_derivative_m223, tolerance);
    EXPECT_NEAR(reference_cartesian[1][2][0], numerical_third_derivative_m231, tolerance);
    EXPECT_NEAR(reference_cartesian[1][2][1], numerical_third_derivative_m232, tolerance);
    EXPECT_NEAR(reference_cartesian[1][2][2], numerical_third_derivative_m233, tolerance);

    posB = double3(posB_reference.x, posB_reference.y, posB_reference.z + 0.5 * delta);
    auto [ez2, d1z2, z2_hessian, d3z2] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    posB = double3(posB_reference.x, posB_reference.y, posB_reference.z - 0.5 * delta);
    auto [ez1, d1z1, z1_hessian, d3z1] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    double numerical_third_derivative_m311 = (z2_hessian[0][0] - z1_hessian[0][0]) / delta;
    double numerical_third_derivative_m312 = (z2_hessian[0][1] - z1_hessian[0][1]) / delta;
    double numerical_third_derivative_m313 = (z2_hessian[0][2] - z1_hessian[0][2]) / delta;
    double numerical_third_derivative_m321 = (z2_hessian[1][0] - z1_hessian[1][0]) / delta;
    double numerical_third_derivative_m322 = (z2_hessian[1][1] - z1_hessian[1][1]) / delta;
    double numerical_third_derivative_m323 = (z2_hessian[1][2] - z1_hessian[1][2]) / delta;
    double numerical_third_derivative_m331 = (z2_hessian[2][0] - z1_hessian[2][0]) / delta;
    double numerical_third_derivative_m332 = (z2_hessian[2][1] - z1_hessian[2][1]) / delta;
    double numerical_third_derivative_m333 = (z2_hessian[2][2] - z1_hessian[2][2]) / delta;

    EXPECT_NEAR(reference_cartesian[2][0][0], numerical_third_derivative_m311, tolerance);
    EXPECT_NEAR(reference_cartesian[2][0][1], numerical_third_derivative_m312, tolerance);
    EXPECT_NEAR(reference_cartesian[2][0][2], numerical_third_derivative_m313, tolerance);
    EXPECT_NEAR(reference_cartesian[2][1][0], numerical_third_derivative_m321, tolerance);
    EXPECT_NEAR(reference_cartesian[2][1][1], numerical_third_derivative_m322, tolerance);
    EXPECT_NEAR(reference_cartesian[2][1][2], numerical_third_derivative_m323, tolerance);
    EXPECT_NEAR(reference_cartesian[2][2][0], numerical_third_derivative_m331, tolerance);
    EXPECT_NEAR(reference_cartesian[2][2][1], numerical_third_derivative_m332, tolerance);
    EXPECT_NEAR(reference_cartesian[2][2][2], numerical_third_derivative_m333, tolerance);
  }
}

TEST(third_derivative_inter_lennard_jones, Test_third_derivative_fractional_methane_in_CHA_triclinic_1x1x1)
{
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);
  Framework f = TestFactories::makeCHA(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeMethane(forceField, 0);
  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {0}, 5);

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  double delta = 1e-6;
  double tolerance = 1e-4;
  double3x3 numerical_hessian;
  size_t typeB = 2;
  double3 posB = double3(5.0, 5.0, 5.0);

  {
    double3 s = system.simulationBox.inverseCell * posB;
    auto [e, d1, d2, reference_cartesian] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    double reference_fractional = system.simulationBox.cell.ax * system.simulationBox.cell.bx *
                                      system.simulationBox.cell.cx * reference_cartesian[0][0][0] +
                                  system.simulationBox.cell.ax * system.simulationBox.cell.bx *
                                      system.simulationBox.cell.cy * reference_cartesian[0][0][1] +
                                  system.simulationBox.cell.ax * system.simulationBox.cell.bx *
                                      system.simulationBox.cell.cz * reference_cartesian[0][0][2] +
                                  system.simulationBox.cell.ax * system.simulationBox.cell.by *
                                      system.simulationBox.cell.cx * reference_cartesian[0][1][0] +
                                  system.simulationBox.cell.ax * system.simulationBox.cell.by *
                                      system.simulationBox.cell.cy * reference_cartesian[0][1][1] +
                                  system.simulationBox.cell.ax * system.simulationBox.cell.by *
                                      system.simulationBox.cell.cz * reference_cartesian[0][2][1];

    /* cell.ay = cell.az = cell.bz = 0
      system.simulationBox.cell.ax * system.simulationBox.cell.bz * system.simulationBox.cell.cx *
      reference_cartesian.m131 + system.simulationBox.cell.ax * system.simulationBox.cell.bz *
      system.simulationBox.cell.cy * reference_cartesian.m132 + system.simulationBox.cell.ax *
      system.simulationBox.cell.bz * system.simulationBox.cell.cz * reference_cartesian.m133 +

      system.simulationBox.cell.ay * system.simulationBox.cell.bx * system.simulationBox.cell.cx *
      reference_cartesian.m211 + system.simulationBox.cell.ay * system.simulationBox.cell.bx *
      system.simulationBox.cell.cy * reference_cartesian.m212 + system.simulationBox.cell.ay *
      system.simulationBox.cell.bx * system.simulationBox.cell.cz * reference_cartesian.m213 +
      system.simulationBox.cell.ay * system.simulationBox.cell.by * system.simulationBox.cell.cx *
      reference_cartesian.m221 + system.simulationBox.cell.ay * system.simulationBox.cell.by *
      system.simulationBox.cell.cy * reference_cartesian.m222 + system.simulationBox.cell.ay *
      system.simulationBox.cell.by * system.simulationBox.cell.cz * reference_cartesian.m223 +
      system.simulationBox.cell.ay * system.simulationBox.cell.bz * system.simulationBox.cell.cx *
      reference_cartesian.m231 + system.simulationBox.cell.ay * system.simulationBox.cell.bz *
      system.simulationBox.cell.cy * reference_cartesian.m232 + system.simulationBox.cell.ay *
      system.simulationBox.cell.bz * system.simulationBox.cell.cz * reference_cartesian.m233 +

      system.simulationBox.cell.az * system.simulationBox.cell.bx * system.simulationBox.cell.cx *
      reference_cartesian.m311 + system.simulationBox.cell.az * system.simulationBox.cell.bx *
      system.simulationBox.cell.cy * reference_cartesian.m312 + system.simulationBox.cell.az *
      system.simulationBox.cell.bx * system.simulationBox.cell.cz * reference_cartesian.m313 +
      system.simulationBox.cell.az * system.simulationBox.cell.by * system.simulationBox.cell.cx *
      reference_cartesian.m321 + system.simulationBox.cell.az * system.simulationBox.cell.by *
      system.simulationBox.cell.cy * reference_cartesian.m322 + system.simulationBox.cell.az *
      system.simulationBox.cell.by * system.simulationBox.cell.cz * reference_cartesian.m323 +
      system.simulationBox.cell.az * system.simulationBox.cell.bz * system.simulationBox.cell.cx *
      reference_cartesian.m331 + system.simulationBox.cell.az * system.simulationBox.cell.bz *
      system.simulationBox.cell.cy * reference_cartesian.m332 + system.simulationBox.cell.az *
      system.simulationBox.cell.bz * system.simulationBox.cell.cz * reference_cartesian.m333;
      */

    posB = system.simulationBox.cell * double3(s.x + 0.5 * delta, s.y, s.z);
    auto [ex2, d1x2, x2_hessian, d3x2] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x - 0.5 * delta, s.y, s.z);
    auto [ex1, d1x1, x1_hessian, d3x1] = Interactions::calculateTricubicDerivativeAtPosition(
        ForceField::InterpolationGridType::LennardJones, system.forceField, system.simulationBox, posB, typeB,
        frameworkAtoms);

    double3x3 x2_hessian_corr = system.simulationBox.cell.transpose() * x2_hessian * system.simulationBox.cell;
    double3x3 x1_hessian_corr = system.simulationBox.cell.transpose() * x1_hessian * system.simulationBox.cell;

    double numerical_third_derivative_m123 = (x2_hessian_corr.bz - x1_hessian_corr.bz) / delta;

    EXPECT_NEAR(numerical_third_derivative_m123, reference_fractional, tolerance);
  }
}
