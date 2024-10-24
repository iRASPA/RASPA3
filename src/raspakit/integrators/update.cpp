module;

#ifdef USE_LEGACY_HEADERS
#include <complex>
#include <span>
#include <vector>
#endif

module integrators_update;

#ifndef USE_LEGACY_HEADERS
import <span>;
import <vector>;
import <complex>;
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

void Integrators::scaleVelocities(std::span<Molecule> moleculePositions, std::pair<double, double> scaling)
{
  for (Molecule& molecule : moleculePositions)
  {
    molecule.velocity = scaling.first * molecule.velocity;
    molecule.orientationMomentum = scaling.second * molecule.orientationMomentum;
  }
}

void Integrators::updatePositions(std::span<Molecule> moleculePositions, double dt)
{
  for (Molecule& molecule : moleculePositions)
  {
    molecule.centerOfMassPosition = molecule.centerOfMassPosition + dt * molecule.velocity;
  }
}

void Integrators::updateVelocities(std::span<Molecule> moleculePositions, double dt)
{
  for (Molecule& molecule : moleculePositions)
  {
    molecule.velocity = molecule.velocity - 0.5 * dt * molecule.gradient * molecule.invMass;
    molecule.orientationMomentum -= 0.5 * dt * molecule.orientationGradient;
  }
}

void Integrators::createCartesianPositions(std::span<const Molecule> moleculePositions,
                                           std::span<Atom> moleculeAtomPositions, std::vector<Component> components)
{
  size_t index{};
  for (const Molecule& molecule : moleculePositions)
  {
    std::span<Atom> span = std::span(&moleculeAtomPositions[index], molecule.numberOfAtoms);
    simd_quatd q = molecule.orientation;
    double3x3 M = double3x3::buildRotationMatrixInverse(q);

    for (size_t i = 0; i != span.size(); i++)
    {
      span[i].position = molecule.centerOfMassPosition + M * components[molecule.componentId].atoms[i].position;
    }

    index += molecule.numberOfAtoms;
  }
}

void Integrators::noSquishFreeRotorOrderTwo(std::span<Molecule> moleculePositions,
                                            const std::vector<Component> components, double dt)
{
  for (size_t i = 0; i != 5; ++i)
  {
    noSquishRotate(moleculePositions, components, 3, 0.5 * dt / 5.0);
    noSquishRotate(moleculePositions, components, 2, 0.5 * dt / 5.0);
    noSquishRotate(moleculePositions, components, 1, dt / 5.0);
    noSquishRotate(moleculePositions, components, 2, 0.5 * dt / 5.0);
    noSquishRotate(moleculePositions, components, 3, 0.5 * dt / 5.0);
  }
}

void Integrators::noSquishRotate(std::span<Molecule> moleculePositions, const std::vector<Component> components,
                                 size_t k, double dt)
{
  std::pair<simd_quatd, simd_quatd> pq;
  for (Molecule& molecule : moleculePositions)
  {
    double3 inverseInertiaVector = components[molecule.componentId].inverseInertiaVector;
    pq = std::make_pair(molecule.orientationMomentum, molecule.orientation);
    pq = Rigid::NoSquishRotate(k, dt, pq, inverseInertiaVector);
    molecule.orientationMomentum = pq.first;
    molecule.orientation = pq.second;
  }
}

void Integrators::updateCenterOfMassAndQuaternionVelocities(std::span<Molecule> moleculePositions,
                                                            std::span<Atom> moleculeAtomPositions,
                                                            std::vector<Component> components)
{

  size_t index{};
  for (Molecule& molecule : moleculePositions)
  {
    std::span<Atom> span = std::span(&moleculeAtomPositions[index], molecule.numberOfAtoms);

    double3 com{};
    for (size_t i = 0; i != span.size(); i++)
    {
      double mass = components[molecule.componentId].definedAtoms[i].second;
      com += mass * span[i].position;
    }
    com = com * molecule.invMass;

    double3 com_velocity{};
    double3 angularMomentum{};
    for (size_t i = 0; i != span.size(); i++)
    {
      double mass = components[molecule.componentId].definedAtoms[i].second;
      com_velocity += mass * span[i].velocity;

      double3 dr = span[i].position - com;
      angularMomentum += mass * double3::cross(dr, span[i].velocity);
    }
    molecule.velocity = com_velocity * molecule.invMass;

    simd_quatd q = molecule.orientation;
    double3x3 M = double3x3::buildRotationMatrix(q);
    simd_quatd omega(0.0, M * angularMomentum);
    molecule.orientationMomentum = 2.0 * q * omega;

    index += molecule.numberOfAtoms;
  }
}

void Integrators::updateCenterOfMassAndQuaternionGradients(std::span<Molecule> moleculePositions,
                                                           std::span<Atom> moleculeAtomPositions,
                                                           std::vector<Component> components)
{
  size_t index{};
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
    simd_quatd q = molecule.orientation;
    double3x3 M = double3x3::buildRotationMatrix(q);
    for (size_t i = 0; i != span.size(); i++)
    {
      double atomMass = components[molecule.componentId].definedAtoms[i].second;
      double3 F = M * (span[i].gradient - com_gradient * atomMass * molecule.invMass);
      double3 dr = components[molecule.componentId].atoms[i].position;
      torque += double3::cross(F, dr);
    }

    molecule.orientationGradient = -2.0 * q * simd_quatd(0.0, torque);
    index += molecule.numberOfAtoms;
  }
}

RunningEnergy Integrators::updateGradients(
    std::span<Atom> moleculeAtomPositions, std::span<Atom> frameworkAtomPositions, const ForceField& forceField,
    const SimulationBox& simulationBox, const std::vector<Component> components,
    std::vector<std::complex<double>>& eik_x, std::vector<std::complex<double>>& eik_y,
    std::vector<std::complex<double>>& eik_z, std::vector<std::complex<double>>& eik_xy,
    std::vector<std::pair<std::complex<double>, std::complex<double>>>& fixedFrameworkStoredEik,
    const std::vector<size_t> numberOfMoleculesPerComponent)
{
  for (Atom& atom : moleculeAtomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  RunningEnergy frameworkMoleculeEnergy = Interactions::computeFrameworkMoleculeGradient(
      forceField, simulationBox, frameworkAtomPositions, moleculeAtomPositions);
  RunningEnergy intermolecularEnergy =
      Interactions::computeInterMolecularGradient(forceField, simulationBox, moleculeAtomPositions);
  RunningEnergy ewaldEnergy = Interactions::computeEwaldFourierGradient(
      eik_x, eik_y, eik_z, eik_xy, fixedFrameworkStoredEik, forceField, simulationBox, components,
      numberOfMoleculesPerComponent, moleculeAtomPositions);

  return frameworkMoleculeEnergy + intermolecularEnergy + ewaldEnergy;
}
