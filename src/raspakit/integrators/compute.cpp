module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <cstddef>
#include <optional>
#include <span>
#endif

module integrators_compute;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import molecule;
import atom;
import double3;
import component;
import simd_quatd;
import double3x3;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import running_energy;
import integrators_cputime;

double Integrators::computeTranslationalKineticEnergy(std::span<const Molecule> moleculePositions)
{
  std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();
  double energy{};
  for (const Molecule& molecule : moleculePositions)
  {
    // Accumulate kinetic energy: 0.5 * mass * velocity squared
    energy += 0.5 * molecule.mass * double3::dot(molecule.velocity, molecule.velocity);
  }
  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  // Update CPU time tracking for this function
  integratorsCPUTime.computeTranslationalKineticEnergy += end - begin;
  return energy;
}

double Integrators::computeRotationalKineticEnergy(std::span<const Molecule> moleculePositions,
                                                   const std::vector<Component> components)
{
  std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();
  double3 ang_vel;
  double energy{};

  for (const Molecule& molecule : moleculePositions)
  {
    // Get inertia and inverse inertia vectors for the molecule's component
    double3 inertiaVector = components[molecule.componentId].inertiaVector;
    double3 inverseInertiaVector = components[molecule.componentId].inverseInertiaVector;
    // Retrieve orientation momentum and orientation
    simd_quatd p = molecule.orientationMomentum;
    p.r = -p.r;
    simd_quatd q = molecule.orientation;
    simd_quatd pq = p * q;

    // Calculate angular velocity
    ang_vel = 0.5 * double3(pq.ix, pq.iy, pq.iz) * inverseInertiaVector;

    // Accumulate rotational kinetic energy: 0.5 * inertia * angular velocity squared
    energy += 0.5 * double3::dot(inertiaVector, ang_vel * ang_vel);
  }
  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  // Update CPU time tracking for this function
  integratorsCPUTime.computeRotationalKineticEnergy += end - begin;
  return energy;
}

double3 Integrators::computeCenterOfMass(std::span<const Molecule> moleculePositions)
{
  std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();
  double3 com{};
  double totalMass{};

  for (const Molecule& molecule : moleculePositions)
  {
    // Accumulate total mass and weighted positions
    totalMass += molecule.mass;
    com += molecule.mass * molecule.centerOfMassPosition;
  }
  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  // Update CPU time tracking for this function
  integratorsCPUTime.computeCenterOfMass += end - begin;

  // Return the center of mass position
  return com / totalMass;
}

// The velocity of the center of mass is the average velocity of all objects in the system weighted by their masses
double3 Integrators::computeCenterOfMassVelocity(std::span<const Molecule> moleculePositions)
{
  std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();
  double3 com_velocity{};
  double totalMass{};

  for (const Molecule& molecule : moleculePositions)
  {
    // Accumulate total mass and weighted velocities
    totalMass += molecule.mass;
    com_velocity += molecule.mass * molecule.velocity;
  }
  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  // Update CPU time tracking for this function
  integratorsCPUTime.computeCenterOfMassVelocity += end - begin;
  // Return the center of mass velocity
  return com_velocity / totalMass;
}

double3 Integrators::computeLinearMomentum(std::span<const Molecule> moleculePositions)
{
  std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();

  double3 com_momentum{};
  for (const Molecule& molecule : moleculePositions)
  {
    // Accumulate linear momentum: mass * velocity
    com_momentum += molecule.mass * molecule.velocity;
  }
  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  // Update CPU time tracking for this function
  integratorsCPUTime.computeLinearMomentum += end - begin;
  return com_momentum;
}
