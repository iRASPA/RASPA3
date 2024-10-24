module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <optional>
#include <span>
#endif

module mc_moves_hybridmc;

#ifndef USE_LEGACY_HEADERS
import <chrono>;
import <span>;
import <optional>;
#endif

import running_energy;
import randomnumbers;
import system;
import atom;
import molecule;
import integrators;
import integrators_update;
import thermostat;


std::optional<RunningEnergy> MC_Moves::hybridMCMove(RandomNumber& random, System& system)
{
  std::chrono::system_clock::time_point time_begin, time_end;

  system.mc_moves_statistics.hybridMC.counts += 1;
  system.mc_moves_statistics.hybridMC.totalCounts += 1;

  // all copied data: moleculePositions, moleculeAtomPositions, thermostat, dt
  // all const data: components, forcefield, simulationbox, numberofmoleculespercomponents
  // all scratch data: eik_x, eik_y, eik_z, eik_xy, fixedFrameworkStoredEik
  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  std::vector<Atom> moleculeAtomPositions(atomPositions.size());
  std::copy(atomPositions.begin(), atomPositions.end(), moleculeAtomPositions.begin());

  std::vector<Molecule> moleculePositions(system.moleculePositions);
  std::optional<Thermostat> thermostat(system.thermostat);
  double dt = system.timeStep;

  // initialize the velocities according to Boltzmann distribution
  Integrators::initializeVelocities(random, moleculePositions, system.components, system.temperature);

  // integrate for N steps
  time_begin = std::chrono::system_clock::now();
  RunningEnergy currentEnergy;
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


  // accept or reject based on energy difference
  if (random.uniform() <
      std::exp(-system.beta * std::abs(currentEnergy.conservedEnergy() - system.runningEnergies.conservedEnergy())))
  {
    system.mc_moves_statistics.hybridMC.accepted += 1;
    system.mc_moves_statistics.hybridMC.totalAccepted += 1;

    system.moleculePositions = moleculePositions;
    system.thermostat = thermostat;
    system.timeStep = dt;

    std::copy(moleculeAtomPositions.begin(), moleculeAtomPositions.end(), atomPositions.begin());
    system.spanOfMoleculeAtoms() = moleculeAtomPositions;
    return system.runningEnergies;
  }
  return std::nullopt;
}