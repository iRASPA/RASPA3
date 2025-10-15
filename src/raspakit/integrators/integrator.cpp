module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <complex>
#include <cstddef>
#include <iostream>
#include <optional>
#include <span>
#include <vector>
#endif

module integrators;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import molecule;
import atom;
import component;
import running_energy;
import thermostat;
import integrators_compute;
import integrators_update;
import integrators_cputime;
import interpolation_energy_grid;

RunningEnergy Integrators::velocityVerlet(
    std::span<Molecule> moleculeData, std::span<Atom> moleculeAtomPositions, const std::vector<Component> components,
    double dt, std::optional<Thermostat>& thermostat, std::span<Atom> frameworkAtomPositions,
    const ForceField& forceField, const SimulationBox& simulationBox, std::vector<std::complex<double>>& eik_x,
    std::vector<std::complex<double>>& eik_y, std::vector<std::complex<double>>& eik_z,
    std::vector<std::complex<double>>& eik_xy,
    std::vector<std::pair<std::complex<double>, std::complex<double>>>& totalEik,
    std::vector<std::pair<std::complex<double>, std::complex<double>>>& fixedFrameworkStoredEik,
    const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::vector<std::size_t> numberOfMoleculesPerComponent)
{
  // apply thermo for temperature control
  if (thermostat.has_value())
  {
    // Adjust velocities using Nose-Hoover thermostat
    double UKineticTranslation = Integrators::computeTranslationalKineticEnergy(moleculeData);
    double UKineticRotation = Integrators::computeRotationalKineticEnergy(moleculeData, components);
    std::pair<double, double> scaling = thermostat->NoseHooverNVT(UKineticTranslation, UKineticRotation);
    scaleVelocities(moleculeData, scaling);
  }

  // Start timing the integration step
  // NOTE: moved from first statement to here as workaround for parsing error in llvm 21.1.1
  std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();

  // evolve the positions a half timestep
  updateVelocities(moleculeData, 0.5 * dt);

  // evolve the positions a full timestep
  updatePositions(moleculeData, dt);

  // evolve the part of rigid bodies involving free rotation
  noSquishFreeRotorOrderTwo(moleculeData, components, dt);

  // create the Cartesian position from center of mass and orientation
  createCartesianPositions(moleculeData, moleculeAtomPositions, components);

  // compute the gradient on all the atoms
  RunningEnergy runningEnergies = updateGradients(
      moleculeAtomPositions, frameworkAtomPositions, forceField, simulationBox, components, eik_x, eik_y, eik_z, eik_xy,
      totalEik, fixedFrameworkStoredEik, interpolationGrids, numberOfMoleculesPerComponent);

  // compute the gradients on the center of mass and the orientation
  updateCenterOfMassAndQuaternionGradients(moleculeData, moleculeAtomPositions, components);

  // evolve the positions a half timestep
  updateVelocities(moleculeData, 0.5 * dt);

  // apply thermo for temperature control
  if (thermostat.has_value())
  {
    // Adjust velocities using Nose-Hoover thermostat
    double UKineticTranslation = computeTranslationalKineticEnergy(moleculeData);
    double UKineticRotation = computeRotationalKineticEnergy(moleculeData, components);
    std::pair<double, double> scaling = thermostat->NoseHooverNVT(UKineticTranslation, UKineticRotation);
    scaleVelocities(moleculeData, scaling);
  }

  // Update the running energies with current kinetic energies
  runningEnergies.translationalKineticEnergy = computeTranslationalKineticEnergy(moleculeData);
  runningEnergies.rotationalKineticEnergy = computeRotationalKineticEnergy(moleculeData, components);
  if (thermostat.has_value())
  {
    runningEnergies.NoseHooverEnergy = thermostat->getEnergy();
  }

  // Update the CPU time spent in the integrator
  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  integratorsCPUTime.velocityVerlet += end - begin;
  return runningEnergies;
}
