#include <gtest/gtest.h>

import std;

import double3;
import units;
import atom;
import molecule;
import pseudo_atom;
import vdwparameters;
import forcefield;
import component;
import system;
import simulationbox;
import running_energy;
import thermostat;
import randomnumbers;
import integrators_compute;
import integrators_update;
import molecular_dynamics;
import mc_moves_parallel_tempering_swap;

namespace
{

ForceField makeLennardJonesForceField()
{
  return ForceField({{"X", false, 16.0, 0.0, 0.0, 6, false}}, {{120.0, 3.5}}, ForceField::MixingRule::Lorentz_Berthelot,
                    9.0, 9.0, 9.0, true, false, false);
}

// A canonical MD replica of 'numberOfParticles' single-site Lennard-Jones particles at 'temperature',
// with a Nose-Hoover chain and Maxwell-Boltzmann velocities, integrated for a few steps so that the
// running energies are those of the integrator (as in the replica-exchange driver).
System makeReplica(const ForceField &forceField, double temperature, std::size_t numberOfParticles,
                   std::size_t seed)
{
  const Component particle(forceField, "particle", 190.0, 4.6e6, 0.0,
                           {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21);
  System system(forceField, SimulationBox(22.0, 22.0, 22.0), false, temperature, 1.0e5, 1.0, {}, {particle}, {},
                {numberOfParticles}, 5);
  system.timeStep = 0.0005;

  RandomNumber random(seed);
  system.initializeGroupData();
  Integrators::initializeVelocities(random, system.moleculeData, system.spanOfMoleculeAtoms(),
                                    system.spanOfMoleculeDynamics(), system.components, system.temperature,
                                    system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(),
                                    &system.forceField, system.spanOfGroupData(), system.spanOfFrameworkGroupData());
  Integrators::removeCenterOfMassVelocityDrift(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
      system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(), &system.forceField,
      system.spanOfGroupData(), system.spanOfFrameworkGroupData());
  system.translationalCenterOfMassConstraint = 3;
  system.thermostat = Thermostat(temperature, system.timeStep, system.translationalDegreesOfFreedom,
                                 system.rotationalDegreesOfFreedom, 3, 1, 0.15);
  system.thermostat->translationalCenterOfMassConstraint = 3;
  system.thermostat->initialize(random);

  system.precomputeTotalGradients();
  for (std::size_t step = 0; step < 5; ++step)
  {
    system.runningEnergies = molecularDynamicsStep(system);
  }
  system.conservedEnergy = system.runningEnergies.conservedEnergy();
  system.referenceEnergy = system.conservedEnergy;
  return system;
}

double kineticEnergy(System &system)
{
  return Integrators::computeTranslationalKineticEnergy(
      system.moleculeData, system.spanOfMoleculeAtoms(), system.spanOfMoleculeDynamics(), system.components,
      system.framework, system.spanOfFrameworkAtoms(), system.spanOfFrameworkDynamics(), &system.forceField,
      system.spanOfGroupData(), system.spanOfFrameworkGroupData());
}

std::vector<double3> positions(const System &system)
{
  std::vector<double3> result;
  for (const Atom &atom : system.atomData) result.push_back(atom.position);
  return result;
}

std::vector<double3> velocities(const System &system)
{
  std::vector<double3> result;
  for (const Molecule &molecule : system.moleculeData) result.push_back(molecule.velocity);
  return result;
}

}  // namespace

// An accepted replica-exchange swap moves the configuration together with its momenta, rescales the
// momenta to the receiving replica's temperature (so the kinetic temperature of each replica is
// unchanged by the exchange), keeps the thermostats with their replicas, and leaves both replicas in
// a state the integrator can continue from (gradients consistent with the positions, running energies
// consistent with a recomputation, conserved-energy reference reset).
TEST(MD_PARALLEL_TEMPERING, accepted_swap_exchanges_configuration_and_rescales_momenta)
{
  const ForceField forceField = makeLennardJonesForceField();
  System systemA = makeReplica(forceField, 250.0, 24, 11);
  System systemB = makeReplica(forceField, 450.0, 24, 13);

  const double temperatureA = systemA.temperature;
  const double temperatureB = systemB.temperature;
  const std::vector<double> thermostatPositionsA = systemA.thermostat->thermostatPositionTranslation;
  const std::vector<double> thermostatPositionsB = systemB.thermostat->thermostatPositionTranslation;

  std::vector<double3> positionsA = positions(systemA);
  std::vector<double3> positionsB = positions(systemB);
  std::vector<double3> velocitiesA = velocities(systemA);
  std::vector<double3> velocitiesB = velocities(systemB);
  double kineticA = kineticEnergy(systemA);
  double kineticB = kineticEnergy(systemB);

  RandomNumber random(29);
  std::optional<std::pair<RunningEnergy, RunningEnergy>> accepted{};
  for (std::size_t attempt = 0; attempt < 512uz && !accepted.has_value(); ++attempt)
  {
    positionsA = positions(systemA);
    positionsB = positions(systemB);
    velocitiesA = velocities(systemA);
    velocitiesB = velocities(systemB);
    kineticA = kineticEnergy(systemA);
    kineticB = kineticEnergy(systemB);
    accepted = MC_Moves::ParallelTemperingSwapMolecularDynamics(random, systemA, systemB);
  }
  ASSERT_TRUE(accepted.has_value());

  // the thermodynamic state and the heat bath stay with the replica
  EXPECT_EQ(systemA.temperature, temperatureA);
  EXPECT_EQ(systemB.temperature, temperatureB);
  ASSERT_TRUE(systemA.thermostat.has_value());
  ASSERT_TRUE(systemB.thermostat.has_value());
  EXPECT_EQ(systemA.thermostat->temperature, temperatureA);
  EXPECT_EQ(systemB.thermostat->temperature, temperatureB);
  EXPECT_EQ(systemA.thermostat->thermostatPositionTranslation, thermostatPositionsA);
  EXPECT_EQ(systemB.thermostat->thermostatPositionTranslation, thermostatPositionsB);

  // the configuration migrated
  ASSERT_EQ(systemA.atomData.size(), positionsB.size());
  ASSERT_EQ(systemB.atomData.size(), positionsA.size());
  for (std::size_t i = 0; i < positionsB.size(); ++i)
  {
    EXPECT_EQ(systemA.atomData[i].position, positionsB[i]);
    EXPECT_EQ(systemB.atomData[i].position, positionsA[i]);
  }

  // the momenta migrated with it and were rescaled by sqrt(T_new / T_old)
  const double scaleA = std::sqrt(temperatureA / temperatureB);
  const double scaleB = std::sqrt(temperatureB / temperatureA);
  ASSERT_EQ(systemA.moleculeData.size(), velocitiesB.size());
  for (std::size_t i = 0; i < velocitiesB.size(); ++i)
  {
    EXPECT_NEAR(systemA.moleculeData[i].velocity.x, scaleA * velocitiesB[i].x, 1.0e-12);
    EXPECT_NEAR(systemA.moleculeData[i].velocity.y, scaleA * velocitiesB[i].y, 1.0e-12);
    EXPECT_NEAR(systemA.moleculeData[i].velocity.z, scaleA * velocitiesB[i].z, 1.0e-12);
    EXPECT_NEAR(systemB.moleculeData[i].velocity.x, scaleB * velocitiesA[i].x, 1.0e-12);
    EXPECT_NEAR(systemB.moleculeData[i].velocity.y, scaleB * velocitiesA[i].y, 1.0e-12);
    EXPECT_NEAR(systemB.moleculeData[i].velocity.z, scaleB * velocitiesA[i].z, 1.0e-12);
  }
  EXPECT_NEAR(kineticEnergy(systemA), kineticB * temperatureA / temperatureB, 1.0e-9 * kineticB);
  EXPECT_NEAR(kineticEnergy(systemB), kineticA * temperatureB / temperatureA, 1.0e-9 * kineticA);
  EXPECT_NEAR(systemA.runningEnergies.translationalKineticEnergy, kineticEnergy(systemA), 1.0e-9 * kineticB);
  EXPECT_NEAR(systemB.runningEnergies.translationalKineticEnergy, kineticEnergy(systemB), 1.0e-9 * kineticA);

  // the running energies are those of the swapped-in configuration
  for (System *system : {&systemA, &systemB})
  {
    const RunningEnergy recomputed = system->computeTotalEnergies();
    EXPECT_NEAR(system->runningEnergies.potentialEnergy(), recomputed.potentialEnergy(), 1.0e-6);
    EXPECT_NEAR(system->runningEnergies.NoseHooverEnergy, system->thermostat->getEnergy(), 1.0e-9);
    EXPECT_DOUBLE_EQ(system->referenceEnergy, system->runningEnergies.conservedEnergy());
    EXPECT_DOUBLE_EQ(system->conservedEnergy, system->referenceEnergy);
  }

  // the integrator can continue from the swapped state: the gradients belong to the new positions,
  // so the extended-system energy is conserved over the following steps
  for (System *system : {&systemA, &systemB})
  {
    const double reference = system->referenceEnergy;
    for (std::size_t step = 0; step < 50; ++step)
    {
      system->runningEnergies = molecularDynamicsStep(*system);
      EXPECT_LT(std::abs(system->runningEnergies.conservedEnergy() - reference) / std::abs(reference), 1.0e-4);
    }
    EXPECT_EQ(system->temperature, system == &systemA ? temperatureA : temperatureB);
  }
}

// A rejected exchange must leave both replicas untouched (positions, momenta and reference energy).
TEST(MD_PARALLEL_TEMPERING, rejected_swap_preserves_dynamical_state)
{
  const ForceField forceFieldA = makeLennardJonesForceField();
  // a different Hamiltonian makes the pair incompatible: the swap is rejected before any change
  const ForceField forceFieldB({{"X", false, 16.0, 0.0, 0.0, 6, false}}, {{60.0, 3.5}},
                               ForceField::MixingRule::Lorentz_Berthelot, 9.0, 9.0, 9.0, true, false, false);
  System systemA = makeReplica(forceFieldA, 250.0, 16, 17);
  System systemB = makeReplica(forceFieldB, 450.0, 16, 19);

  const std::vector<double3> positionsA = positions(systemA);
  const std::vector<double3> velocitiesA = velocities(systemA);
  const std::vector<double3> positionsB = positions(systemB);
  const std::vector<double3> velocitiesB = velocities(systemB);
  const double referenceA = systemA.referenceEnergy;
  const double referenceB = systemB.referenceEnergy;

  RandomNumber random(31);
  const std::size_t drawsBefore = random.count;
  EXPECT_FALSE(MC_Moves::ParallelTemperingSwapMolecularDynamics(random, systemA, systemB).has_value());
  EXPECT_EQ(random.count, drawsBefore);

  EXPECT_EQ(positions(systemA), positionsA);
  EXPECT_EQ(velocities(systemA), velocitiesA);
  EXPECT_EQ(positions(systemB), positionsB);
  EXPECT_EQ(velocities(systemB), velocitiesB);
  EXPECT_DOUBLE_EQ(systemA.referenceEnergy, referenceA);
  EXPECT_DOUBLE_EQ(systemB.referenceEnergy, referenceB);
}
