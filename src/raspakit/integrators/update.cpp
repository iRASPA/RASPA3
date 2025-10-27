module;

#ifdef USE_PRECOMPILED_HEADERS
#include "pch.h"
#endif

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <complex>
#include <cstddef>
#include <iostream>
#include <span>
#include <vector>
#endif

module integrators_update;

#ifdef USE_STD_IMPORT
import std;
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
import integrators_compute;
import randomnumbers;
import units;
import interpolation_energy_grid;

void Integrators::scaleVelocities(std::span<Molecule> moleculeData, std::pair<double, double> scaling)
{
  std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();
  // Scale velocities and orientation momenta of each molecule
  for (Molecule& molecule : moleculeData)
  {
    molecule.velocity = scaling.first * molecule.velocity;
    molecule.orientationMomentum = scaling.second * molecule.orientationMomentum;
  }
  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  integratorsCPUTime.scaleVelocities += end - begin;
}

void Integrators::removeCenterOfMassVelocityDrift(std::span<Molecule> moleculeData)
{
  std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();

  double3 totalVelocity = computeCenterOfMassVelocity(moleculeData);
  // Scale velocities and orientation momenta of each molecule
  for (Molecule& molecule : moleculeData)
  {
    molecule.velocity -= totalVelocity;
  }
  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  integratorsCPUTime.removeCenterOfMassVelocity += end - begin;
}

void Integrators::updatePositions(std::span<Molecule> moleculeData, double dt)
{
  std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();
  // Update the center of mass positions for each molecule
  for (Molecule& molecule : moleculeData)
  {
    molecule.centerOfMassPosition += dt * molecule.velocity;
  }
  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  integratorsCPUTime.updatePositions += end - begin;
}

void Integrators::updateVelocities(std::span<Molecule> moleculeData, double dt)
{
  std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();

  // Update velocities and orientation momenta based on gradients
  for (Molecule& molecule : moleculeData)
  {
    molecule.velocity -= dt * molecule.gradient * molecule.invMass;
    molecule.orientationMomentum -= dt * molecule.orientationGradient;
  }
  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  integratorsCPUTime.updateVelocities += end - begin;
}

void Integrators::initializeVelocities(RandomNumber& random, std::span<Molecule> moleculeData,
                                       const std::vector<Component> components, double temperature)
{
  for (Molecule& molecule : moleculeData)
  {
    // Draw reandom vector with variance kBT / m
    molecule.velocity = double3(random.Gaussian(), random.Gaussian(), random.Gaussian()) *
                        std::sqrt(Units::KB * temperature / molecule.mass);

    // Draw random quaternion with variance by kBT / I
    double3 I = components[molecule.componentId].inertiaVector;
    double3 invI = components[molecule.componentId].inverseInertiaVector;

    // factor sqrt(3/2) added to match correct mean
    double3 angularVelocity =
        double3(random.Gaussian(), random.Gaussian(), random.Gaussian()) * sqrt(invI * Units::KB * temperature);

    // rotate into orientation frame
    simd_quatd q = molecule.orientation;
    molecule.orientationMomentum.ix =
        2.0 * (q.r * (I.x * angularVelocity.x) - q.iz * (I.y * angularVelocity.y) + q.iy * (I.z * angularVelocity.z));
    molecule.orientationMomentum.iy =
        2.0 * (q.iz * (I.x * angularVelocity.x) + q.r * (I.y * angularVelocity.y) - q.ix * (I.z * angularVelocity.z));
    molecule.orientationMomentum.iz =
        2.0 * (-q.iy * (I.x * angularVelocity.x) + q.ix * (I.y * angularVelocity.y) + q.r * (I.z * angularVelocity.z));
    molecule.orientationMomentum.r =
        2.0 * (-q.ix * (I.x * angularVelocity.x) - q.iy * (I.y * angularVelocity.y) - q.iz * (I.z * angularVelocity.z));
  }
}

void Integrators::createCartesianPositions(std::span<const Molecule> moleculeData,
                                           std::span<Atom> moleculeAtomPositions, std::vector<Component> components)
{
  std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();
  std::size_t index{};
  // Convert molecule positions and orientations to atom positions
  for (const Molecule& molecule : moleculeData)
  {
    std::span<Atom> span = std::span(&moleculeAtomPositions[index], molecule.numberOfAtoms);
    simd_quatd q = molecule.orientation;

    // Calculate positions of each atom in the molecule
    for (std::size_t i = 0; i != span.size(); i++)
    {
      span[i].position = molecule.centerOfMassPosition + q * components[molecule.componentId].atoms[i].position;
    }
    index += molecule.numberOfAtoms;
  }
  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  integratorsCPUTime.createCartesianPositions += end - begin;
}

void Integrators::noSquishFreeRotorOrderTwo(std::span<Molecule> moleculeData, const std::vector<Component> components,
                                            double dt)
{
  std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();
  for (Molecule& molecule : moleculeData)
  {
    double3 inverseInertiaVector = components[molecule.componentId].inverseInertiaVector;
    std::pair<simd_quatd, simd_quatd> pq = std::make_pair(molecule.orientationMomentum, molecule.orientation);
    pq = Rigid::NoSquishFreeRotorOrderTwo(dt, pq, inverseInertiaVector);
    molecule.orientationMomentum = pq.first;
    molecule.orientation = pq.second;
  }
  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  integratorsCPUTime.noSquishFreeRotorOrderTwo += end - begin;
}

void Integrators::updateCenterOfMassAndQuaternionVelocities(std::span<Molecule> moleculeData,
                                                            std::span<Atom> moleculeAtomPositions,
                                                            std::vector<Component> components)
{
  std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();
  std::size_t index{};
  // Update center of mass velocities and orientation momenta for each molecule
  for (Molecule& molecule : moleculeData)
  {
    std::span<Atom> span = std::span(&moleculeAtomPositions[index], molecule.numberOfAtoms);

    double3 com{};
    for (std::size_t i = 0; i != span.size(); i++)
    {
      double mass = components[molecule.componentId].definedAtoms[i].second;
      com += mass * span[i].position;
    }
    com = com * molecule.invMass;

    double3 com_velocity{};
    double3 angularMomentum{};
    for (std::size_t i = 0; i != span.size(); i++)
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
  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  integratorsCPUTime.updateCenterOfMassAndQuaternionVelocities += end - begin;
}

void Integrators::updateCenterOfMassAndQuaternionGradients(std::span<Molecule> moleculeData,
                                                           std::span<Atom> moleculeAtomPositions,
                                                           std::vector<Component> components)
{
  std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();
  std::size_t index{};

  // Update gradients of center of mass and orientation for each molecule
  for (Molecule& molecule : moleculeData)
  {
    std::span<Atom> span = std::span(&moleculeAtomPositions[index], molecule.numberOfAtoms);
    double3 com_gradient{};
    for (std::size_t i = 0; i != span.size(); i++)
    {
      com_gradient += span[i].gradient;
    }
    molecule.gradient = com_gradient;

    double3 torque{};
    simd_quatd q = molecule.orientation;
    double3x3 M = double3x3::buildRotationMatrix(q);
    for (std::size_t i = 0; i != span.size(); i++)
    {
      double atomMass = components[molecule.componentId].definedAtoms[i].second;
      double3 F = M * (span[i].gradient - com_gradient * atomMass * molecule.invMass);
      double3 dr = components[molecule.componentId].atoms[i].position;
      torque += double3::cross(F, dr);
    }

    // Compute orientation gradient based on torque
    molecule.orientationGradient = -2.0 * q * simd_quatd(0.0, torque);

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
    std::vector<std::pair<std::complex<double>, std::complex<double>>>& totalEik,
    std::vector<std::pair<std::complex<double>, std::complex<double>>>& fixedFrameworkStoredEik,
    const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::vector<std::size_t> numberOfMoleculesPerComponent)
{
  std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();

  // Initialize gradients to zero
  for (Atom& atom : moleculeAtomPositions)
  {
    atom.gradient = double3(0.0, 0.0, 0.0);
  }

  // Compute gradients and energies due to interactions
  RunningEnergy frameworkMoleculeEnergy = Interactions::computeFrameworkMoleculeGradient(
      forceField, simulationBox, frameworkAtomPositions, moleculeAtomPositions, interpolationGrids);
  RunningEnergy intermolecularEnergy =
      Interactions::computeInterMolecularGradient(forceField, simulationBox, moleculeAtomPositions);
  RunningEnergy ewaldEnergy = Interactions::computeEwaldFourierGradient(
      eik_x, eik_y, eik_z, eik_xy, totalEik, fixedFrameworkStoredEik, forceField, simulationBox, components,
      numberOfMoleculesPerComponent, moleculeAtomPositions);

  RunningEnergy internal_energies{};
  for (const Component& component : components)
  {
    internal_energies += component.intraMolecularPotentials.computeInternalGradient(moleculeAtomPositions);
  }

  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  integratorsCPUTime.updateGradients += end - begin;
  return frameworkMoleculeEnergy + intermolecularEnergy + ewaldEnergy + internal_energies;
}
