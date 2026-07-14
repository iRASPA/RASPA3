#include <gtest/gtest.h>

import std;

import int3;
import double3;
import double3x3;
import simd_quatd;
import units;
import atom;
import atom_dynamics;
import pseudo_atom;
import vdwparameters;
import forcefield;
import framework;
import component;
import connectivity_table;
import intra_molecular_potentials;
import bond_potential;
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
import thermostat;
import randomnumbers;

namespace
{
ForceField makeFlexibleFrameworkForceField()
{
  return ForceField({{"A", false, 2.0, 0.0, 0.0, 1, false}, {"B", false, 5.0, 0.0, 0.0, 2, false}},
                    {{0.0, 1.0}, {0.0, 1.0}}, ForceField::MixingRule::Lorentz_Berthelot, 8.0, 8.0, 8.0, false, false,
                    false);
}

Framework makeTwoAtomFramework(const ForceField& forceField, bool rigid)
{
  const SimulationBox box(20.0, 20.0, 20.0);
  std::vector<Atom> atoms{Atom({0.20, 0.25, 0.25}, 0.0, 1.0, 0, 0, 0, 0, true),
                          Atom({0.25, 0.25, 0.25}, 0.0, 1.0, 0, 1, 0, 0, true)};
  Framework framework(forceField, "two-atom-framework", box, 1, atoms, atoms, int3(1, 1, 1));
  framework.rigid = rigid;
  return framework;
}
}  // namespace

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
  std::cout << Integrators::computeTranslationalKineticEnergy(system.moleculeData, system.spanOfMoleculeAtoms(),
                                                              system.spanOfMoleculeDynamics(), system.components)
            << std::endl;
  std::cout << Integrators::computeRotationalKineticEnergy(system.moleculeData, system.components) << std::endl;

  std::cout << "Center of Mass\n";
  std::cout << system.moleculeData[0].centerOfMassPosition.to_string() << std::endl;
  std::cout << system.moleculeData[1].centerOfMassPosition.to_string() << std::endl;

  std::cout << "Velocity\n";
  std::cout << system.moleculeData[0].velocity.to_string() << std::endl;
  std::cout << system.moleculeData[1].velocity.to_string() << std::endl;

  std::cout << "Gradient\n";
  std::cout << system.moleculeData[0].gradient.to_string() << std::endl;
  std::cout << system.moleculeData[1].gradient.to_string() << std::endl;

  std::cout << "Orientation\n";
  std::cout << system.moleculeData[0].orientation.to_string() << std::endl;
  std::cout << system.moleculeData[1].orientation.to_string() << std::endl;

  std::cout << "OrientationMomentum\n";
  std::cout << system.moleculeData[0].orientationMomentum.to_string() << std::endl;
  std::cout << system.moleculeData[1].orientationMomentum.to_string() << std::endl;

  std::cout << "OrientationGradient\n";
  std::cout << system.moleculeData[0].orientationGradient.to_string() << std::endl;
  std::cout << system.moleculeData[1].orientationGradient.to_string() << std::endl;
}

TEST(integrators_flexible_framework, degrees_of_freedom_and_rigid_regression)
{
  ForceField forceField = makeFlexibleFrameworkForceField();
  System flexible(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1.0e4, 1.0,
                  makeTwoAtomFramework(forceField, false), {}, {}, {}, 5);
  System rigid(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1.0e4, 1.0,
               makeTwoAtomFramework(forceField, true), {}, {}, {}, 5);

  EXPECT_EQ(flexible.translationalDegreesOfFreedom, 6u);
  EXPECT_EQ(flexible.rotationalDegreesOfFreedom, 0u);
  EXPECT_EQ(rigid.translationalDegreesOfFreedom, 0u);
  EXPECT_EQ(rigid.rotationalDegreesOfFreedom, 0u);
  EXPECT_EQ(rigid.numberOfRigidFrameworkAtoms, 2u);
}

TEST(integrators_flexible_framework, mass_weighted_velocity_position_and_kinetic_energy)
{
  ForceField forceField = makeFlexibleFrameworkForceField();
  System system(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1.0e4, 1.0,
                makeTwoAtomFramework(forceField, false), {}, {}, {}, 5);
  std::span<Atom> atoms = system.spanOfFrameworkAtoms();
  std::span<AtomDynamics> dynamics = system.spanOfFrameworkDynamics();
  atoms[0].position = {4.0, 5.0, 6.0};
  atoms[1].position = {7.0, 8.0, 9.0};
  dynamics[0].velocity = {1.0, -2.0, 0.5};
  dynamics[1].velocity = {-0.5, 1.0, 2.0};
  dynamics[0].gradient = {4.0, -2.0, 6.0};
  dynamics[1].gradient = {-5.0, 10.0, 2.5};

  Integrators::updateVelocities({}, {}, {}, {}, 0.2, system.framework, atoms, dynamics, &system.forceField);
  DOUBLE3_EXPECT_NEAR(dynamics[0].velocity, double3(0.8, -1.9, 0.2), 1.0e-14);
  DOUBLE3_EXPECT_NEAR(dynamics[1].velocity, double3(-0.4, 0.8, 1.95), 1.0e-14);

  Integrators::updatePositions({}, {}, {}, {}, 0.25, system.framework, atoms, dynamics);
  DOUBLE3_EXPECT_NEAR(atoms[0].position, double3(4.2, 4.525, 6.05), 1.0e-14);
  DOUBLE3_EXPECT_NEAR(atoms[1].position, double3(6.9, 8.2, 9.4875), 1.0e-14);

  const double expected =
      0.5 * 2.0 * (0.8 * 0.8 + 1.9 * 1.9 + 0.2 * 0.2) + 0.5 * 5.0 * (0.4 * 0.4 + 0.8 * 0.8 + 1.95 * 1.95);
  EXPECT_NEAR(Integrators::computeTranslationalKineticEnergy({}, {}, {}, {}, system.framework, atoms, dynamics,
                                                             &system.forceField),
              expected, 1.0e-14);
}

TEST(integrators_flexible_framework, center_of_mass_drift_combines_framework_and_molecule_momentum)
{
  ForceField forceField = ForceField::makeZeoliteForceField(8.0, false, false, false);
  Framework framework = Framework::makeITQ29(forceField, int3(1, 1, 1));
  framework.rigid = false;
  Component methane = Component::makeMethane(forceField, 0);
  System system(forceField, std::nullopt, false, 300.0, 1.0e4, 1.0, framework, {methane}, {}, {1}, 5);

  system.moleculeData[0].velocity = {1.0, -0.5, 0.25};
  for (AtomDynamics& dynamics : system.spanOfFrameworkDynamics())
  {
    dynamics.velocity = {-0.2, 0.4, -0.1};
  }
  const double frameworkMass = *system.frameworkMass();
  const double moleculeMass = system.moleculeData[0].mass;
  const double3 initialMoleculeVelocity = system.moleculeData[0].velocity;
  const double3 initialFrameworkVelocity = system.spanOfFrameworkDynamics()[0].velocity;
  const double3 drift = (moleculeMass * initialMoleculeVelocity + frameworkMass * initialFrameworkVelocity) /
                        (moleculeMass + frameworkMass);

  Integrators::removeCenterOfMassVelocityDrift(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
      system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(), &system.forceField);

  DOUBLE3_EXPECT_NEAR(system.moleculeData[0].velocity, initialMoleculeVelocity - drift, 1.0e-14);
  DOUBLE3_EXPECT_NEAR(system.spanOfFrameworkDynamics()[0].velocity, initialFrameworkVelocity - drift, 1.0e-14);
  double3 finalMomentum = moleculeMass * system.moleculeData[0].velocity;
  for (std::size_t i = 0; i < system.spanOfFrameworkAtoms().size(); ++i)
  {
    const double mass = system.forceField.pseudoAtoms[system.spanOfFrameworkAtoms()[i].type].mass;
    finalMomentum += mass * system.spanOfFrameworkDynamics()[i].velocity;
  }
  DOUBLE3_EXPECT_NEAR(finalMomentum, double3{}, 1.0e-11);
}

TEST(integrators_flexible_framework, thermostat_scaling_and_target_degrees_of_freedom)
{
  ForceField forceField = makeFlexibleFrameworkForceField();
  System system(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1.0e4, 1.0,
                makeTwoAtomFramework(forceField, false), {}, {}, {}, 5);
  system.setThermostat(Thermostat(3, 1, 0.15));
  ASSERT_TRUE(system.thermostat.has_value());
  EXPECT_EQ(system.thermostat->translationalDegreesOfFreedom, 6u);
  EXPECT_EQ(system.thermostat->rotationalDegreesOfFreedom, 0u);

  std::span<AtomDynamics> dynamics = system.spanOfFrameworkDynamics();
  dynamics[0].velocity = {1.0, 2.0, 3.0};
  dynamics[1].velocity = {-2.0, 1.0, 0.5};
  Integrators::scaleVelocities({}, {}, {}, {}, {0.25, 7.0}, system.framework, dynamics);
  DOUBLE3_EXPECT_NEAR(dynamics[0].velocity, double3(0.25, 0.5, 0.75), 1.0e-14);
  DOUBLE3_EXPECT_NEAR(dynamics[1].velocity, double3(-0.5, 0.25, 0.125), 1.0e-14);

  const double kinetic = Integrators::computeTranslationalKineticEnergy(
      {}, {}, {}, {}, system.framework, system.spanOfFrameworkAtoms(), dynamics, &system.forceField);
  const double measuredTemperature =
      2.0 * kinetic / (Units::KB * static_cast<double>(system.translationalDegreesOfFreedom));
  EXPECT_NEAR(measuredTemperature, 2.0 * (0.5 * 2.0 * 0.875 + 0.5 * 5.0 * 0.328125) / (6.0 * Units::KB), 1.0e-12);
}

TEST(integrators_flexible_framework, nose_hoover_velocity_verlet_scales_framework_kinetic_energy)
{
  ForceField forceField = makeFlexibleFrameworkForceField();
  const auto makeSystem = [&]()
  {
    System system(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1.0e4, 1.0,
                  makeTwoAtomFramework(forceField, false), {}, {}, {}, 5);
    system.spanOfFrameworkAtoms()[0].position = {4.0, 5.0, 6.0};
    system.spanOfFrameworkAtoms()[1].position = {8.0, 9.0, 10.0};
    system.spanOfFrameworkDynamics()[0].velocity = {1.0, -2.0, 0.5};
    system.spanOfFrameworkDynamics()[1].velocity = {-0.4, 0.8, 1.2};
    return system;
  };
  System nve = makeSystem();
  System nvt = makeSystem();
  nvt.setThermostat(Thermostat(3, 1, 0.15));
  ASSERT_TRUE(nvt.thermostat.has_value());
  RandomNumber random(1729);
  nvt.thermostat->initialize(random);

  const RunningEnergy nveEnergy = Integrators::velocityVerlet(
      {}, {}, {}, {}, nve.timeStep, nve.thermostat, nve.spanOfFrameworkAtoms(), nve.forceField, nve.simulationBox,
      nve.eik_x, nve.eik_y, nve.eik_z, nve.eik_xy, nve.totalEik, nve.fixedFrameworkStoredEik, nve.interpolationGrids,
      {}, nve.framework, nve.spanOfFrameworkDynamics());
  const RunningEnergy nvtEnergy = Integrators::velocityVerlet(
      {}, {}, {}, {}, nvt.timeStep, nvt.thermostat, nvt.spanOfFrameworkAtoms(), nvt.forceField, nvt.simulationBox,
      nvt.eik_x, nvt.eik_y, nvt.eik_z, nvt.eik_xy, nvt.totalEik, nvt.fixedFrameworkStoredEik, nvt.interpolationGrids,
      {}, nvt.framework, nvt.spanOfFrameworkDynamics());

  DOUBLE3_EXPECT_NEAR(nve.spanOfFrameworkDynamics()[0].velocity, double3(1.0, -2.0, 0.5), 1.0e-14);
  DOUBLE3_EXPECT_NEAR(nve.spanOfFrameworkDynamics()[1].velocity, double3(-0.4, 0.8, 1.2), 1.0e-14);
  const double thermostatScale =
      nvt.spanOfFrameworkDynamics()[0].velocity.x / nve.spanOfFrameworkDynamics()[0].velocity.x;
  EXPECT_GT(std::abs(thermostatScale - 1.0), 1.0e-8);
  DOUBLE3_EXPECT_NEAR(nvt.spanOfFrameworkDynamics()[0].velocity,
                      thermostatScale * nve.spanOfFrameworkDynamics()[0].velocity, 1.0e-13);
  DOUBLE3_EXPECT_NEAR(nvt.spanOfFrameworkDynamics()[1].velocity,
                      thermostatScale * nve.spanOfFrameworkDynamics()[1].velocity, 1.0e-13);

  const double explicitFrameworkKinetic = Integrators::computeTranslationalKineticEnergy(
      {}, {}, {}, {}, nvt.framework, nvt.spanOfFrameworkAtoms(), nvt.spanOfFrameworkDynamics(), &nvt.forceField);
  EXPECT_GT(explicitFrameworkKinetic, 0.0);
  EXPECT_NEAR(nvtEnergy.translationalKineticEnergy, explicitFrameworkKinetic, 1.0e-13);
  EXPECT_NE(nvtEnergy.translationalKineticEnergy, nveEnergy.translationalKineticEnergy);
}

TEST(integrators_flexible_framework, harmonic_velocity_verlet_has_bounded_short_nve_drift)
{
  ForceField forceField = makeFlexibleFrameworkForceField();
  Framework framework = makeTwoAtomFramework(forceField, false);
  framework.intraMolecularPotentials.bonds = {BondPotential({0, 1}, BondType::Harmonic, {2000.0, 1.0})};
  framework.intraMolecularImageShifts.bonds = {{{int3{}, int3{}}}};
  System system(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1.0e4, 1.0, framework, {}, {}, {}, 5);
  system.spanOfFrameworkAtoms()[0].position = {5.0, 5.0, 5.0};
  system.spanOfFrameworkAtoms()[1].position = {6.1, 5.0, 5.0};
  system.spanOfFrameworkDynamics()[0].velocity = {0.1, 0.0, 0.0};
  system.spanOfFrameworkDynamics()[1].velocity = {-0.04, 0.0, 0.0};

  RunningEnergy initial = Integrators::updateGradients(
      {}, {}, {}, system.spanOfFrameworkAtoms(), system.forceField, system.simulationBox, {}, system.eik_x,
      system.eik_y, system.eik_z, system.eik_xy, system.totalEik, system.fixedFrameworkStoredEik,
      system.interpolationGrids, {}, system.framework, system.spanOfFrameworkDynamics());
  initial.translationalKineticEnergy =
      Integrators::computeTranslationalKineticEnergy({}, {}, {}, {}, system.framework, system.spanOfFrameworkAtoms(),
                                                     system.spanOfFrameworkDynamics(), &system.forceField);
  const double initialEnergy = initial.conservedEnergy();
  double maximumDrift = 0.0;

  for (std::size_t step = 0; step < 250; ++step)
  {
    const RunningEnergy current =
        Integrators::velocityVerlet({}, {}, {}, {}, 1.0e-4, system.thermostat, system.spanOfFrameworkAtoms(),
                                    system.forceField, system.simulationBox, system.eik_x, system.eik_y, system.eik_z,
                                    system.eik_xy, system.totalEik, system.fixedFrameworkStoredEik,
                                    system.interpolationGrids, {}, system.framework, system.spanOfFrameworkDynamics());
    maximumDrift = std::max(maximumDrift, std::abs(current.conservedEnergy() - initialEnergy));
  }

  EXPECT_GT(std::abs(system.spanOfFrameworkAtoms()[0].position.x - 5.0), 1.0e-5);
  EXPECT_LT(maximumDrift / std::max(1.0, std::abs(initialEnergy)), 2.0e-5);
}

TEST(integrators, Test_2_CO2_in_ITQ_29_2x2x2_inter)
{
  ForceField forceField = ForceField::makeZeoliteForceField(11.8, true, false, true);
  Framework f = Framework::makeITQ29(forceField, int3(1, 1, 1));
  Component c = Component::makeCO2(forceField, 0, false);
  System system = System(forceField, std::nullopt, false, 300.0, 1e4, 1.0, {f}, {c}, {}, {2}, 5);

  std::span<Molecule> moleculeData(system.moleculeData);

  moleculeData[0].centerOfMassPosition = double3(5.93355, 7.93355, 5.93355);
  moleculeData[1].centerOfMassPosition = double3(5.93355, 3.93355, 5.93355);
  moleculeData[0].orientation = simd_quatd(0.0, 0.0, 0.0, 1.0);
  moleculeData[1].orientation = simd_quatd(0.0, 0.0, 0.0, 1.0);

  Integrators::createCartesianPositions(moleculeData, system.spanOfMoleculeAtoms(), system.components);

  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  DOUBLE3_EXPECT_NEAR(atomData[0].position, double3(5.93355, 7.93355, 7.08255), 1e-6);
  DOUBLE3_EXPECT_NEAR(atomData[1].position, double3(5.93355, 7.93355, 5.93355), 1e-6);
  DOUBLE3_EXPECT_NEAR(atomData[2].position, double3(5.93355, 7.93355, 4.78455), 1e-6);
  DOUBLE3_EXPECT_NEAR(atomData[3].position, double3(5.93355, 3.93355, 7.08255), 1e-6);
  DOUBLE3_EXPECT_NEAR(atomData[4].position, double3(5.93355, 3.93355, 5.93355), 1e-6);
  DOUBLE3_EXPECT_NEAR(atomData[5].position, double3(5.93355, 3.93355, 4.78455), 1e-6);

  moleculeData[0].velocity = double3(0.0, -1.0, 0.0);
  moleculeData[0].orientationMomentum = simd_quatd(10.0, 0.0, 0.0, 0.0);
  moleculeData[1].velocity = double3(0.0, 1.0, 0.0);
  moleculeData[1].orientationMomentum = simd_quatd(-10.0, 0.0, 0.0, 0.0);

  Integrators::velocityVerlet(system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(),
                              system.components, system.timeStep, system.thermostat, system.spanOfFrameworkAtoms(),
                              system.forceField, system.simulationBox, system.eik_x, system.eik_y, system.eik_z,
                              system.eik_xy, system.totalEik, system.fixedFrameworkStoredEik, system.interpolationGrids,
                              system.numberOfMoleculesPerComponent);

  DOUBLE3_EXPECT_NEAR(system.moleculeData[0].centerOfMassPosition, double3(5.933550, 7.933050, 5.933550), 1e-6);
  DOUBLE3_EXPECT_NEAR(system.moleculeData[1].centerOfMassPosition, double3(5.933550, 3.934050, 5.933550), 1e-6);
  DOUBLE3_EXPECT_NEAR(system.moleculeData[0].velocity, double3(0.000218, -0.999519, 0.00004), 1e-6);
  DOUBLE3_EXPECT_NEAR(system.moleculeData[1].velocity, double3(0.000218, 0.999519, 0.00004), 1e-6);
  QUATD_EXPECT_NEAR(system.moleculeData[0].orientation, simd_quatd(0.00003, 0.0, 0.0, 1.0), 1e-6);
  QUATD_EXPECT_NEAR(system.moleculeData[1].orientation, simd_quatd(-0.00003, 0.0, 0.0, 1.0), 1e-6);
  QUATD_EXPECT_NEAR(system.moleculeData[0].orientationMomentum, simd_quatd(10.000011, 0.0, 0.0, -0.000296), 1e-6);
  QUATD_EXPECT_NEAR(system.moleculeData[1].orientationMomentum, simd_quatd(-10.000011, 0.0, 0.0, -0.000296), 1e-6);

  Integrators::velocityVerlet(system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(),
                              system.components, system.timeStep, system.thermostat, system.spanOfFrameworkAtoms(),
                              system.forceField, system.simulationBox, system.eik_x, system.eik_y, system.eik_z,
                              system.eik_xy, system.totalEik, system.fixedFrameworkStoredEik, system.interpolationGrids,
                              system.numberOfMoleculesPerComponent);

  DOUBLE3_EXPECT_NEAR(system.moleculeData[0].centerOfMassPosition, double3(5.933550, 7.932550, 5.933550), 1e-6);
  DOUBLE3_EXPECT_NEAR(system.moleculeData[1].centerOfMassPosition, double3(5.933550, 3.934550, 5.933550), 1e-6);
  DOUBLE3_EXPECT_NEAR(system.moleculeData[0].velocity, double3(0.000653, -0.998560, 0.000121), 1e-6);
  DOUBLE3_EXPECT_NEAR(system.moleculeData[1].velocity, double3(0.000653, 0.998560, 0.000121), 1e-6);
  QUATD_EXPECT_NEAR(system.moleculeData[0].orientation, simd_quatd(0.000059, 0.0, 0.0, 1.0), 1e-6);
  QUATD_EXPECT_NEAR(system.moleculeData[1].orientation, simd_quatd(-0.000059, 0.0, 0.0, 1.0), 1e-6);
  QUATD_EXPECT_NEAR(system.moleculeData[0].orientationMomentum, simd_quatd(10.000043, 0.0, 0.0, -0.000592), 1e-6);
  QUATD_EXPECT_NEAR(system.moleculeData[1].orientationMomentum, simd_quatd(-10.000043, 0.0, 0.0, -0.000592), 1e-6);
}
