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
import running_energy;
import interactions_intermolecular;
import interactions_framework_molecule;
import interactions_ewald;
import energy_status;

TEST(hessian_inter_real_ewald, Test_gradient_cartesian_methane_in_CHA_triclinic_1x1x1)
{
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);
  Framework f = TestFactories::makeCHA(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeMethane(forceField, 0);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {}, {0}, 5);

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  double delta = 1e-5;
  double tolerance = 1e-4;
  double3 numerical_gradient;
  double chargeB = 1.0;
  double3 posB = double3(5.0, 5.0, 5.0);

  {
    double3 s = system.simulationBox.inverseCell * posB;
    auto [e, reference_cartesian, d2] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);
    double3 reference_fractional = system.simulationBox.cell.transpose() * reference_cartesian;

    // finite difference x
    posB = system.simulationBox.cell * double3(s.x + 0.5 * delta, s.y, s.z);
    auto [x2_energy, d1x2, d2x2] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x - 0.5 * delta, s.y, s.z);
    auto [x1_energy, d1x1, d2x1] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    // finite difference y
    posB = system.simulationBox.cell * double3(s.x, s.y + 0.5 * delta, s.z);
    auto [y2_energy, d1y2, d2y2] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y - 0.5 * delta, s.z);
    auto [y1_energy, d1y1, d2y1] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    // finite difference z
    posB = system.simulationBox.cell * double3(s.x, s.y, s.z + 0.5 * delta);
    auto [z2_energy, d1z2, d2z2] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z - 0.5 * delta);
    auto [z1_energy, d1z1, d2z1] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

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
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);
  Framework f = TestFactories::makeCHA(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeMethane(forceField, 0);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {}, {0}, 5);

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  double delta = 1e-5;
  double tolerance = 1e-4;
  double3 numerical_gradient;
  double chargeB = 1.0;
  double3 posB = double3(5.0, 5.0, 5.0);

  {
    double3 s = system.simulationBox.inverseCell * posB;
    auto [e, reference_cartesian, d2] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);
    double3 reference_fractional = system.simulationBox.cell.transpose() * reference_cartesian;

    // finite difference x
    posB = system.simulationBox.cell * double3(s.x + 0.5 * delta, s.y, s.z);
    auto [x2_energy, d1x2, d2x2] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x - 0.5 * delta, s.y, s.z);
    auto [x1_energy, d1x1, d2x1] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    // finite difference y
    posB = system.simulationBox.cell * double3(s.x, s.y + 0.5 * delta, s.z);
    auto [y2_energy, d1y2, d2y2] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y - 0.5 * delta, s.z);
    auto [y1_energy, d1y1, d2y1] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    // finite difference z
    posB = system.simulationBox.cell * double3(s.x, s.y, s.z + 0.5 * delta);
    auto [z2_energy, d1z2, d2z2] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z - 0.5 * delta);
    auto [z1_energy, d1z1, d2z1] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

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
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);
  Framework f = TestFactories::makeCHA(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeMethane(forceField, 0);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {}, {0}, 5);

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  double delta = 1e-5;
  double tolerance = 1e-4;
  double3x3 numerical_hessian;
  double chargeB = 1.0;
  double3 posB_reference = double3(5.0, 5.0, 5.0);

  {
    double3 posB;

    auto [e, d1, reference_cartesian] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB_reference, chargeB, frameworkAtoms);

    // finite difference x
    posB = double3(posB_reference.x + 0.5 * delta, posB_reference.y, posB_reference.z);
    auto [ex2, x2_gradient, d2x2] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    posB = double3(posB_reference.x - 0.5 * delta, posB_reference.y, posB_reference.z);
    auto [ex1, x1_gradient, d2x1] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    // finite difference y
    posB = double3(posB_reference.x, posB_reference.y + 0.5 * delta, posB_reference.z);
    auto [ey2, y2_gradient, d2y2] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    posB = double3(posB_reference.x, posB_reference.y - 0.5 * delta, posB_reference.z);
    auto [ey1, y1_gradient, d2y1] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    // finite difference z
    posB = double3(posB_reference.x, posB_reference.y, posB_reference.z + 0.5 * delta);
    auto [ez2, z2_gradient, d2z2] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    posB = double3(posB_reference.x, posB_reference.y, posB_reference.z - 0.5 * delta);
    auto [ez1, z1_gradient, d2z1] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

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
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);
  Framework f = TestFactories::makeCHA(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeMethane(forceField, 0);

  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {}, {0}, 5);

  std::span<const Atom> frameworkAtoms = system.spanOfFrameworkAtoms();

  double delta = 1e-5;
  double tolerance = 1e-4;
  double3x3 numerical_hessian;
  double chargeB = 1.0;
  double3 posB = double3(5.0, 5.0, 5.0);

  {
    double3 s = system.simulationBox.inverseCell * posB;
    auto [e, d1, reference_cartesian] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);
    double3x3 reference_fractional =
        system.simulationBox.cell.transpose() * reference_cartesian * system.simulationBox.cell;

    [[maybe_unused]] double test =
        system.simulationBox.cell.ax * system.simulationBox.cell.ax * reference_cartesian.ax +
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
    auto [ex2, x2_gradient, d2x2] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x - 0.5 * delta, s.y, s.z);
    auto [ex1, x1_gradient, d2x1] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    // finite difference y
    posB = system.simulationBox.cell * double3(s.x, s.y + 0.5 * delta, s.z);
    auto [ey2, y2_gradient, d2y2] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y - 0.5 * delta, s.z);
    auto [ey1, y1_gradient, d2y1] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    // finite difference z
    posB = system.simulationBox.cell * double3(s.x, s.y, s.z + 0.5 * delta);
    auto [ez2, z2_gradient, d2z2] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

    posB = system.simulationBox.cell * double3(s.x, s.y, s.z - 0.5 * delta);
    auto [ez1, z1_gradient, d2z1] = Interactions::calculateHessianAtPositionCoulomb(
        system.forceField, system.simulationBox, posB, chargeB, frameworkAtoms);

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
