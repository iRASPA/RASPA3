module;

module integrators;

import std;

import molecule;
import atom;
import atom_dynamics;
import component;
import running_energy;
import thermostat;
import integrators_compute;
import integrators_update;
import integrators_cputime;
import interpolation_energy_grid;
import framework;

RunningEnergy Integrators::velocityVerlet(
    std::span<Molecule> moleculeData, std::span<Atom> moleculeAtomPositions, std::span<AtomDynamics> moleculeDynamics,
    const std::vector<Component> &components, double dt, std::optional<Thermostat>& thermostat,
    std::span<Atom> frameworkAtomPositions,
    const ForceField& forceField, const SimulationBox& simulationBox, std::vector<std::complex<double>>& eik_x,
    std::vector<std::complex<double>>& eik_y, std::vector<std::complex<double>>& eik_z,
    std::vector<std::complex<double>>& eik_xy,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& trialEik,
    std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>>& fixedFrameworkStoredEik,
    const std::vector<std::optional<InterpolationEnergyGrid>>& interpolationGrids,
    const std::vector<std::size_t> &numberOfMoleculesPerComponent, const std::optional<Framework>& framework,
    std::span<AtomDynamics> frameworkDynamics)
{
  // apply thermo for temperature control
  if (thermostat.has_value())
  {
    // Adjust velocities using Nose-Hoover thermostat
    double UKineticTranslation = computeTranslationalKineticEnergy(
        moleculeData, moleculeAtomPositions, moleculeDynamics, components, framework, frameworkAtomPositions,
        frameworkDynamics, &forceField);
    double UKineticRotation = computeRotationalKineticEnergy(moleculeData, components);
    std::pair<double, double> scaling = thermostat->NoseHooverNVT(UKineticTranslation, UKineticRotation);
    scaleVelocities(moleculeData, moleculeAtomPositions, moleculeDynamics, components, scaling, framework,
                    frameworkDynamics);
  }

  // Start timing the integration step
  // NOTE: moved from first statement to here as workaround for parsing error in llvm 21.1.1
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  // evolve the positions a half timestep
  updateVelocities(moleculeData, moleculeAtomPositions, moleculeDynamics, components, dt, framework,
                   frameworkAtomPositions, frameworkDynamics, &forceField);

  // evolve the positions a full timestep
  updatePositions(moleculeData, moleculeAtomPositions, moleculeDynamics, components, dt, framework,
                  frameworkAtomPositions, frameworkDynamics);

  // evolve the part of rigid bodies involving free rotation
  noSquishFreeRotorOrderTwo(moleculeData, components, dt);

  // create the Cartesian position from center of mass and orientation
  createCartesianPositions(moleculeData, moleculeAtomPositions, components);

  // compute the gradient on all the atoms
  RunningEnergy runningEnergies = updateGradients(
      moleculeData, moleculeAtomPositions, moleculeDynamics, frameworkAtomPositions, forceField, simulationBox,
      components, eik_x, eik_y, eik_z, eik_xy,
      trialEik, fixedFrameworkStoredEik, interpolationGrids, numberOfMoleculesPerComponent, framework,
      frameworkDynamics);

  // compute the gradients on the center of mass and the orientation
  updateCenterOfMassAndQuaternionGradients(moleculeData, moleculeAtomPositions, moleculeDynamics, components);

  // evolve the positions a half timestep
  updateVelocities(moleculeData, moleculeAtomPositions, moleculeDynamics, components, dt, framework,
                   frameworkAtomPositions, frameworkDynamics, &forceField);

  // apply thermo for temperature control
  if (thermostat.has_value())
  {
    // Adjust velocities using Nose-Hoover thermostat
    double UKineticTranslation = computeTranslationalKineticEnergy(
        moleculeData, moleculeAtomPositions, moleculeDynamics, components, framework, frameworkAtomPositions,
        frameworkDynamics, &forceField);
    double UKineticRotation = computeRotationalKineticEnergy(moleculeData, components);
    std::pair<double, double> scaling = thermostat->NoseHooverNVT(UKineticTranslation, UKineticRotation);
    scaleVelocities(moleculeData, moleculeAtomPositions, moleculeDynamics, components, scaling, framework,
                    frameworkDynamics);
  }

  // Update the running energies with current kinetic energies
  runningEnergies.translationalKineticEnergy = computeTranslationalKineticEnergy(
      moleculeData, moleculeAtomPositions, moleculeDynamics, components, framework, frameworkAtomPositions,
      frameworkDynamics, &forceField);
  runningEnergies.rotationalKineticEnergy = computeRotationalKineticEnergy(moleculeData, components);
  if (thermostat.has_value())
  {
    runningEnergies.NoseHooverEnergy = thermostat->getEnergy();
  }

  // Update the CPU time spent in the integrator
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  integratorsCPUTime.velocityVerlet += end - begin;
  return runningEnergies;
}
