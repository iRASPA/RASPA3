#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <cstddef>
#include <span>
#include <vector>

import int3;
import double3;
import double3x3;
import simd_quatd;
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
import integrators_update;
import integrators;
import integrators_compute;
import molecule;
import interpolation_energy_grid;

void DOUBLE3_EXPECT_NEAR(double3 a, double3 b, double tol)
{
  EXPECT_NEAR(a.x, b.x, tol);
  EXPECT_NEAR(a.y, b.y, tol);
  EXPECT_NEAR(a.z, b.z, tol);
}

void QUATD_EXPECT_NEAR(simd_quatd a, simd_quatd b, double tol)
{
  EXPECT_NEAR(a.r, b.r, tol);
  EXPECT_NEAR(a.ix, b.ix, tol);
  EXPECT_NEAR(a.iy, b.iy, tol);
  EXPECT_NEAR(a.iz, b.iz, tol);
}

void print_system(System& system)
{
  std::cout << "Energy (tr, rot)\n";
  std::cout << Integrators::computeTranslationalKineticEnergy(system.moleculePositions) << std::endl;
  std::cout << Integrators::computeRotationalKineticEnergy(system.moleculePositions, system.components) << std::endl;

  std::cout << "Center of Mass\n";
  std::cout << system.moleculePositions[0].centerOfMassPosition.to_string() << std::endl;
  std::cout << system.moleculePositions[1].centerOfMassPosition.to_string() << std::endl;

  std::cout << "Velocity\n";
  std::cout << system.moleculePositions[0].velocity.to_string() << std::endl;
  std::cout << system.moleculePositions[1].velocity.to_string() << std::endl;

  std::cout << "Gradient\n";
  std::cout << system.moleculePositions[0].gradient.to_string() << std::endl;
  std::cout << system.moleculePositions[1].gradient.to_string() << std::endl;

  std::cout << "Orientation\n";
  std::cout << system.moleculePositions[0].orientation.to_string() << std::endl;
  std::cout << system.moleculePositions[1].orientation.to_string() << std::endl;

  std::cout << "OrientationMomentum\n";
  std::cout << system.moleculePositions[0].orientationMomentum.to_string() << std::endl;
  std::cout << system.moleculePositions[1].orientationMomentum.to_string() << std::endl;

  std::cout << "OrientationGradient\n";
  std::cout << system.moleculePositions[0].orientationGradient.to_string() << std::endl;
  std::cout << system.moleculePositions[1].orientationGradient.to_string() << std::endl;
}

TEST(integrators, Test_2_CO2_in_ITQ_29_2x2x2_inter)
{
  ForceField forceField = TestFactories::makeDefaultFF(11.8, true, false, true);
  Framework f = TestFactories::makeITQ29(forceField, int3(1, 1, 1));
  Component c = TestFactories::makeCO2(forceField, 0, false);
  System system = System(0, forceField, std::nullopt, 300.0, 1e4, 1.0, {f}, {c}, {2}, 5);

  std::span<Molecule> moleculePositions(system.moleculePositions);

  moleculePositions[0].centerOfMassPosition = double3(5.93355, 7.93355, 5.93355);
  moleculePositions[1].centerOfMassPosition = double3(5.93355, 3.93355, 5.93355);
  moleculePositions[0].orientation = simd_quatd(0.0, 0.0, 0.0, 1.0);
  moleculePositions[1].orientation = simd_quatd(0.0, 0.0, 0.0, 1.0);

  Integrators::createCartesianPositions(moleculePositions, system.spanOfMoleculeAtoms(), system.components);

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  DOUBLE3_EXPECT_NEAR(atomPositions[0].position, double3(5.93355, 7.93355, 7.08255), 1e-6);
  DOUBLE3_EXPECT_NEAR(atomPositions[1].position, double3(5.93355, 7.93355, 5.93355), 1e-6);
  DOUBLE3_EXPECT_NEAR(atomPositions[2].position, double3(5.93355, 7.93355, 4.78455), 1e-6);
  DOUBLE3_EXPECT_NEAR(atomPositions[3].position, double3(5.93355, 3.93355, 7.08255), 1e-6);
  DOUBLE3_EXPECT_NEAR(atomPositions[4].position, double3(5.93355, 3.93355, 5.93355), 1e-6);
  DOUBLE3_EXPECT_NEAR(atomPositions[5].position, double3(5.93355, 3.93355, 4.78455), 1e-6);

  moleculePositions[0].velocity = double3(0.0, -1.0, 0.0);
  moleculePositions[0].orientationMomentum = simd_quatd(10.0, 0.0, 0.0, 0.0);
  moleculePositions[1].velocity = double3(0.0, 1.0, 0.0);
  moleculePositions[1].orientationMomentum = simd_quatd(-10.0, 0.0, 0.0, 0.0);

  Integrators::velocityVerlet(system.moleculePositions, system.spanOfMoleculeAtoms(), system.components,
                              system.timeStep, system.thermostat, system.spanOfFrameworkAtoms(), system.forceField,
                              system.simulationBox, system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                              system.totalEik, system.fixedFrameworkStoredEik, system.interpolationGrids,
                              system.numberOfMoleculesPerComponent);

  DOUBLE3_EXPECT_NEAR(system.moleculePositions[0].centerOfMassPosition, double3(5.933550, 7.933050, 5.933550), 1e-6);
  DOUBLE3_EXPECT_NEAR(system.moleculePositions[1].centerOfMassPosition, double3(5.933550, 3.934050, 5.933550), 1e-6);
  DOUBLE3_EXPECT_NEAR(system.moleculePositions[0].velocity, double3(0.000218, -0.999519, 0.00004), 1e-6);
  DOUBLE3_EXPECT_NEAR(system.moleculePositions[1].velocity, double3(0.000218, 0.999519, 0.00004), 1e-6);
  QUATD_EXPECT_NEAR(system.moleculePositions[0].orientation, simd_quatd(0.00003, 0.0, 0.0, 1.0), 1e-6);
  QUATD_EXPECT_NEAR(system.moleculePositions[1].orientation, simd_quatd(-0.00003, 0.0, 0.0, 1.0), 1e-6);
  QUATD_EXPECT_NEAR(system.moleculePositions[0].orientationMomentum, simd_quatd(10.000011, 0.0, 0.0, -0.000296), 1e-6);
  QUATD_EXPECT_NEAR(system.moleculePositions[1].orientationMomentum, simd_quatd(-10.000011, 0.0, 0.0, -0.000296), 1e-6);

  Integrators::velocityVerlet(system.moleculePositions, system.spanOfMoleculeAtoms(), system.components,
                              system.timeStep, system.thermostat, system.spanOfFrameworkAtoms(), system.forceField,
                              system.simulationBox, system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                              system.totalEik, system.fixedFrameworkStoredEik, system.interpolationGrids,
                              system.numberOfMoleculesPerComponent);

  DOUBLE3_EXPECT_NEAR(system.moleculePositions[0].centerOfMassPosition, double3(5.933550, 7.932550, 5.933550), 1e-6);
  DOUBLE3_EXPECT_NEAR(system.moleculePositions[1].centerOfMassPosition, double3(5.933550, 3.934550, 5.933550), 1e-6);
  DOUBLE3_EXPECT_NEAR(system.moleculePositions[0].velocity, double3(0.000653, -0.998560, 0.000121), 1e-6);
  DOUBLE3_EXPECT_NEAR(system.moleculePositions[1].velocity, double3(0.000653, 0.998560, 0.000121), 1e-6);
  QUATD_EXPECT_NEAR(system.moleculePositions[0].orientation, simd_quatd(0.000059, 0.0, 0.0, 1.0), 1e-6);
  QUATD_EXPECT_NEAR(system.moleculePositions[1].orientation, simd_quatd(-0.000059, 0.0, 0.0, 1.0), 1e-6);
  QUATD_EXPECT_NEAR(system.moleculePositions[0].orientationMomentum, simd_quatd(10.000043, 0.0, 0.0, -0.000592), 1e-6);
  QUATD_EXPECT_NEAR(system.moleculePositions[1].orientationMomentum, simd_quatd(-10.000043, 0.0, 0.0, -0.000592), 1e-6);
}
