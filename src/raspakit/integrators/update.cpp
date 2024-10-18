module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <complex>
#include <span>
#include <vector>
#endif

module integrators_update;

#ifndef USE_LEGACY_HEADERS
import <span>;
import <vector>;
import <complex>;
import <chrono>;
#endif

import molecule;
import double3;
import atom;
import component;
import simd_quatd;
import double3x3;
import rigid;
import forcefield;
import simulationbox;
import running_energy;
import forcefield;
import interactions_ewald;
import interactions_intermolecular;
import interactions_framework_molecule;
import integrators_cputime;
void Integrators::scaleVelocities(std::span<Molecule> moleculePositions, std::pair<double, double> scaling)
{
  std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();
  // Scale velocities and orientation momenta of each molecule
  for (Molecule& molecule : moleculePositions)
  {
    molecule.velocity = scaling.first * molecule.velocity;
    molecule.orientationMomentum = scaling.second * molecule.orientationMomentum;
  }
  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  integratorsCPUTime.scaleVelocities += end - begin;
}

void Integrators::updatePositions(std::span<Molecule> moleculePositions, double dt)
{
  std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();
  // Update the center of mass positions for each molecule
  for (Molecule& molecule : moleculePositions)
  {
    molecule.centerOfMassPosition = molecule.centerOfMassPosition + dt * molecule.velocity;
  }
  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  integratorsCPUTime.updatePositions += end - begin;
}

void Integrators::updateVelocities(std::span<Molecule> moleculePositions, double dt)
{
  std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();

  // Update velocities and orientation momenta based on gradients
  for (Molecule& molecule : moleculePositions)
  {
    molecule.velocity = molecule.velocity - 0.5 * dt * molecule.gradient * molecule.invMass;
    molecule.orientationMomentum = 0.5 * dt * molecule.orientationGradient;
  }
  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  integratorsCPUTime.updateVelocities += end - begin;
}

void Integrators::createCartesianPositions(std::span<const Molecule> moleculePositions,
                                           std::span<Atom> moleculeAtomPositions, std::vector<Component> components)
{
  std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();
  size_t index{};
  // Convert molecule positions and orientations to atom positions
  for (const Molecule& molecule : moleculePositions)
  {
    std::span<Atom> span = std::span(&moleculeAtomPositions[index], molecule.numberOfAtoms);
    simd_quatd q = molecule.orientation;
    double3x3 M = double3x3::buildRotationMatrixInverse(q);

    // Calculate positions of each atom in the molecule
    for (size_t i = 0; i != span.size(); i++)
    {
      span[i].position = molecule.centerOfMassPosition + M * components[molecule.componentId].atoms[i].position;
    }

    index += molecule.numberOfAtoms;
  }
  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  integratorsCPUTime.createCartesianPositions += end - begin;
}

void Integrators::noSquishFreeRotorOrderTwo(std::span<Molecule> moleculePositions,
                                            const std::vector<Component> components, double dt)
{
  std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();
  // Perform second-order NoSquish free rotor integration over 5 sub-steps
  for (size_t i = 0; i != 5; ++i)
  {
    noSquishRotate(moleculePositions, components, 3, 0.5 * dt / 5.0);
    noSquishRotate(moleculePositions, components, 2, 0.5 * dt / 5.0);
    noSquishRotate(moleculePositions, components, 1, dt / 5.0);
    noSquishRotate(moleculePositions, components, 2, 0.5 * dt / 5.0);
    noSquishRotate(moleculePositions, components, 3, 0.5 * dt / 5.0);
  }
  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  integratorsCPUTime.noSquishFreeRotorOrderTwo += end - begin;
}

void Integrators::noSquishRotate(std::span<Molecule> moleculePositions, const std::vector<Component> components,
                                 size_t k, double dt)
{
  std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();
  std::pair<simd_quatd, simd_quatd> pq;
  // Rotate each molecule using the NoSquish algorithm
  for (Molecule& molecule : moleculePositions)
  {
    double3 inverseInertiaVector = components[molecule.componentId].inverseInertiaVector;
    pq = std::make_pair(molecule.orientationMomentum, molecule.orientation);
    pq = Rigid::NoSquishRotate(k, dt, pq, inverseInertiaVector);
    molecule.orientationMomentum = pq.first;
    molecule.orientation = pq.second;
  }
  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  integratorsCPUTime.noSquishRotate += end - begin;
}

void Integrators::updateCenterOfMassAndQuaternionVelocities(std::span<Molecule> moleculePositions,
                                                            std::span<Atom> moleculeAtomPositions,
                                                            std::vector<Component> components)
{
  std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();
  size_t index{};
  // Update center of mass velocities and orientation momenta for each molecule
  for (Molecule& molecule : moleculePositions)
  {
    std::span<Atom> span = std::span(&moleculeAtomPositions[index], molecule.numberOfAtoms);

    double3 com_velocity;
    double3 angularMomentum;

    for (size_t i = 0; i != span.size(); i++)
    {
      double atomMass = components[molecule.componentId].definedAtoms[i].second;
      com_velocity += atomMass * span[i].velocity;

      double3 dr = span[i].position - molecule.centerOfMassPosition;
      angularMomentum += atomMass * double3::cross(dr, span[i].velocity);
    }
    molecule.velocity = com_velocity * molecule.invMass;

    double3 I = components[molecule.componentId].inertiaVector;
    double3 invI = components[molecule.componentId].inverseInertiaVector;
    simd_quatd q = molecule.orientation;
    double3x3 M = double3x3::buildRotationMatrix(q);
    double3 angularVelocity = (M * angularMomentum) * invI;

    // Update orientation momentum based on angular velocity and inertia
    molecule.orientationMomentum.ix =
        2.0 * (q.r * (I.x * angularVelocity.x) - q.iz * (I.y * angularVelocity.y) + q.iy * (I.z * angularVelocity.z));
    molecule.orientationMomentum.iy =
        2.0 * (q.iz * (I.x * angularVelocity.x) + q.r * (I.y * angularVelocity.y) - q.ix * (I.z * angularVelocity.z));
    molecule.orientationMomentum.iz =
        2.0 * (-q.iy * (I.x * angularVelocity.x) + q.ix * (I.y * angularVelocity.y) + q.r * (I.z * angularVelocity.z));
    molecule.orientationMomentum.r =
        2.0 * (-q.ix * (I.x * angularVelocity.x) - q.iy * (I.y * angularVelocity.y) - q.iz * (I.z * angularVelocity.z));

    index += molecule.numberOfAtoms;
  }
  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  integratorsCPUTime.updateCenterOfMassAndQuaternionVelocities += end - begin;
}

void Integrators::updateCenterOfMassAndQuaternionGradients(std::span<Molecule> moleculePositions,
                                                           std::span<Atom> moleculeAtomPositions,
                                                           std::vector<Component> components)
{
  std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();
  size_t index{};
  // Update gradients of center of mass and orientation for each molecule
  for (Molecule& molecule : moleculePositions)
  {
    std::span<Atom> span = std::span(&moleculeAtomPositions[index], molecule.numberOfAtoms);
    double3 com_gradient{};
    for (size_t i = 0; i != span.size(); i++)
    {
      com_gradient += span[i].gradient;
    }
    molecule.gradient = com_gradient;

    double3 torque{};
    simd_quatd orientation = molecule.orientation;
    double3x3 M = double3x3::buildRotationMatrix(orientation);
    for (size_t i = 0; i != span.size(); i++)
    {
      double atomMass = components[molecule.componentId].definedAtoms[i].second;
      double3 F = M * (span[i].gradient - com_gradient * atomMass * molecule.invMass);
      double3 dr = components[molecule.componentId].atoms[i].position;
      torque += double3::cross(F, dr);
    }

    // Compute orientation gradient based on torque
    molecule.orientationGradient.ix =
        -2.0 * (orientation.r * torque.x - orientation.iz * torque.y + orientation.iy * torque.z);
    molecule.orientationGradient.iy =
        -2.0 * (orientation.iz * torque.x + orientation.r * torque.y - orientation.ix * torque.z);
    molecule.orientationGradient.iz =
        -2.0 * (-orientation.iy * torque.x + orientation.ix * torque.y + orientation.r * torque.z);
    molecule.orientationGradient.r =
        -2.0 * (-orientation.ix * torque.x - orientation.iy * torque.y - orientation.iz * torque.z);

    index += molecule.numberOfAtoms;
  }
  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  integratorsCPUTime.updateCenterOfMassAndQuaternionGradients += end - begin;
}

RunningEnergy Integrators::updateGradients(
    std::span<Atom> moleculeAtomPositions, std::span<Atom> frameworkAtomPositions, const ForceField& forceField,
    const SimulationBox& simulationBox, const std::vector<Component> components,
    std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
    std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
    std::vector<std::pair<std::complex<double>, std::complex<double>>>& fixedFrameworkStoredEik,
    const std::vector<size_t> numberOfMoleculesPerComponent)
{
  std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();

  // Initialize gradients to zero
  for (Atom& atom : moleculeAtomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  // Compute gradients and energies due to interactions
  RunningEnergy frameworkMoleculeEnergy = Interactions::computeFrameworkMoleculeGradient(
      forceField, simulationBox, frameworkAtomPositions, moleculeAtomPositions);
  RunningEnergy intermolecularEnergy =
      Interactions::computeInterMolecularGradient(forceField, simulationBox, moleculeAtomPositions);
  RunningEnergy ewaldEnergy = Interactions::computeEwaldFourierGradient(
      eik_x, eik_y, eik_z, eik_xy, fixedFrameworkStoredEik, forceField, simulationBox, components,
      numberOfMoleculesPerComponent, moleculeAtomPositions);

  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  integratorsCPUTime.updateGradients += end - begin;
  return frameworkMoleculeEnergy + intermolecularEnergy + ewaldEnergy;
}
