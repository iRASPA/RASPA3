module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <iostream>
#include <optional>
#include <span>
#endif

module mc_moves_hybridmc;

#ifndef USE_LEGACY_HEADERS
import <chrono>;
import <span>;
import <optional>;
import <iostream>;
#endif

import double3;
import running_energy;
import randomnumbers;
import system;
import atom;
import molecule;
import integrators;
import integrators_update;
import integrators_compute;
import thermostat;
import units;

std::optional<RunningEnergy> MC_Moves::hybridMCMove(RandomNumber& random, System& system)
{
  // Initialize time points for performance measurement
  std::chrono::system_clock::time_point time_begin, time_end;

  system.mc_moves_statistics.hybridMC.counts += 1;
  system.mc_moves_statistics.hybridMC.totalCounts += 1;

  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  std::vector<Atom> moleculeAtomPositions(atomPositions.size());
  std::copy(atomPositions.begin(), atomPositions.end(), moleculeAtomPositions.begin());

  std::vector<Molecule> moleculePositions(system.moleculePositions);
  std::optional<Thermostat> thermostat(system.thermostat);
  double dt = system.timeStep;

  // Initialize velocities and remove center of mass velocity
  Integrators::initializeVelocities(random, moleculePositions, system.components, system.temperature);
  Integrators::removeCenterOfMassVelocity(moleculePositions);

  double uKinTrans = Integrators::computeTranslationalKineticEnergy(moleculePositions);
  double translationalTemperature = 2.0 * uKinTrans /
                                    (Units::KB * static_cast<double>(system.translationalDegreesOfFreedom -
                                                                     system.translationalCenterOfMassConstraint));

  double uKinRot = Integrators::computeRotationalKineticEnergy(moleculePositions, system.components);
  double rotationalTemperature = 2.0 * uKinRot / (Units::KB * static_cast<double>(system.rotationalDegreesOfFreedom));

  // Scale velocities based on temperature ratios
  std::pair<double, double> scaling(system.temperature / translationalTemperature,
                                    system.temperature / rotationalTemperature);
  Integrators::scaleVelocities(moleculePositions, scaling);

  // Recompute current energy
  system.runningEnergies =
      Integrators::updateGradients(system.spanOfMoleculeAtoms(), system.spanOfFrameworkAtoms(), system.forceField,
                                   system.simulationBox, system.components, system.eik_x, system.eik_y, system.eik_z,
                                   system.eik_xy, system.fixedFrameworkStoredEik, system.numberOfMoleculesPerComponent);

  RunningEnergy referenceEnergy = system.runningEnergies;
  referenceEnergy.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(moleculePositions);
  referenceEnergy.rotationalKineticEnergy =
      Integrators::computeRotationalKineticEnergy(moleculePositions, system.components);

  // Integrate for N steps
  time_begin = std::chrono::system_clock::now();
  RunningEnergy currentEnergy = referenceEnergy;
  for (size_t step = 0; step < system.numberOfHybridMCSteps; ++step)
  {
    currentEnergy = Integrators::velocityVerlet(
        moleculePositions, moleculeAtomPositions, system.components, dt, thermostat, system.spanOfFrameworkAtoms(),
        system.forceField, system.simulationBox, system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
        system.fixedFrameworkStoredEik, system.numberOfMoleculesPerComponent);
  }
  time_end = std::chrono::system_clock::now();

  system.mc_moves_cputime.hybridMCIntegration += (time_end - time_begin);
  system.mc_moves_statistics.hybridMC.constructed += 1;
  system.mc_moves_statistics.hybridMC.totalConstructed += 1;

  double drift = std::abs(currentEnergy.conservedEnergy() - referenceEnergy.conservedEnergy());

  // Accept or reject move based on energy difference
  if (random.uniform() < std::exp(-system.beta * drift))
  {
    std::cout << "accepted\n";
    system.mc_moves_statistics.hybridMC.accepted += 1;
    system.mc_moves_statistics.hybridMC.totalAccepted += 1;

    system.moleculePositions = moleculePositions;
    system.thermostat = thermostat;
    system.timeStep = dt;

    std::copy(moleculeAtomPositions.begin(), moleculeAtomPositions.end(), atomPositions.begin());
    system.spanOfMoleculeAtoms() = moleculeAtomPositions;

    Integrators::createCartesianPositions(system.moleculePositions, system.spanOfMoleculeAtoms(), system.components);
    return system.runningEnergies;
  }
  return std::nullopt;
}
