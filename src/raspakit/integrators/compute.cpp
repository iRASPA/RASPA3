module;

#ifdef USE_LEGACY_HEADERS
#include <optional>
#include <span>
#endif

module integrators_compute;

#ifndef USE_LEGACY_HEADERS
import <span>;
import <optional>;
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

double Integrators::computeTranslationalKineticEnergy(std::span<const Molecule> moleculePositions)
{
  double energy{};
  for (const Molecule& molecule : moleculePositions)
  {
    energy += 0.5 * molecule.mass * double3::dot(molecule.velocity, molecule.velocity);
  }
  return energy;
}

double Integrators::computeRotationalKineticEnergy(std::span<const Molecule> moleculePositions,
                                                   const std::vector<Component> components)
{
  double3 ang_vel;
  double energy{};

  for (const Molecule& molecule : moleculePositions)
  {
    double3 inertiaVector = components[molecule.componentId].inertiaVector;
    double3 inverseInertiaVector = components[molecule.componentId].inverseInertiaVector;
    simd_quatd p = molecule.orientationMomentum;
    simd_quatd q = molecule.orientation;
    simd_quatd pq = p * q;

    ang_vel = 0.5 * double3(pq.ix, pq.iy, pq.iz) * inverseInertiaVector;
    energy += 0.5 * double3::dot(inertiaVector, ang_vel * ang_vel);
  }
  return energy;
}

double3 Integrators::computeCenterOfMass(std::span<const Molecule> moleculePositions)
{
  double3 com{};
  double totalMass{};

  for (const Molecule& molecule : moleculePositions)
  {
    totalMass += molecule.mass;
    com += molecule.mass * molecule.centerOfMassPosition;
  }
  return com / totalMass;
}

// The velocity of the center of mass is the average velocity of all objects in the system weighted by their masses
double3 Integrators::computeCenterOfMassVelocity(std::span<const Molecule> moleculePositions)
{
  double3 com_velocity{};
  double totalMass{};

  for (const Molecule& molecule : moleculePositions)
  {
    totalMass += molecule.mass;
    com_velocity += molecule.mass * molecule.velocity;
  }
  return com_velocity / totalMass;
}

double3 Integrators::computeLinearMomentum(std::span<const Molecule> moleculePositions)
{
  double3 com_momentum{};
  for (const Molecule& molecule : moleculePositions)
  {
    com_momentum += molecule.mass * molecule.velocity;
  }
  return com_momentum;
}
