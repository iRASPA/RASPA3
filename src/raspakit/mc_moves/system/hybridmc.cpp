module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <cmath>
#include <cstddef>
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
import interactions_ewald;
import mc_moves_move_types;

std::optional<RunningEnergy> MC_Moves::hybridMCMove(RandomNumber& random, System& system)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  MoveTypes move = MoveTypes::HybridMC;

  system.mc_moves_statistics.addTrial(move);

  if (system.moleculePositions.size() <= 1)
  {
    return std::nullopt;
  }

  // all copied data: moleculePositions, moleculeAtomPositions, thermostat, dt
  // all const data: components, forcefield, simulationbox, numberofmoleculespercomponents, fixedFrameworkStoredEik
  // all scratch data: eik_x, eik_y, eik_z, eik_xy
  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  std::vector<Atom> moleculeAtomPositions(atomPositions.size());
  std::copy(atomPositions.begin(), atomPositions.end(), moleculeAtomPositions.begin());

  std::vector<Molecule> moleculePositions(system.moleculePositions);
  std::optional<Thermostat> thermostat(system.thermostat);

  // get Timestep from the max change
  double dt = system.mc_moves_statistics.getMaxChange(move);

  // initialize the velocities according to Boltzmann distribution
  // NOTE: it is important that the reference energy has the initial kinetic energies
  Integrators::initializeVelocities(random, moleculePositions, system.components, system.temperature);

  if (system.numberOfFrameworkAtoms == 0)
  {
    Integrators::removeCenterOfMassVelocityDrift(moleculePositions);
  }

  // before getting energy, recompute current energy
  system.precomputeTotalGradients();

  RunningEnergy referenceEnergy = system.runningEnergies;
  referenceEnergy.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(moleculePositions);
  referenceEnergy.rotationalKineticEnergy =
      Integrators::computeRotationalKineticEnergy(moleculePositions, system.components);
  RunningEnergy currentEnergy = referenceEnergy;

  // integrate for N steps
  time_begin = std::chrono::system_clock::now();
  for (size_t step = 0; step < system.numberOfHybridMCSteps; ++step)
  {
    currentEnergy = Integrators::velocityVerlet(
        moleculePositions, moleculeAtomPositions, system.components, dt, thermostat, system.spanOfFrameworkAtoms(),
        system.forceField, system.simulationBox, system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
        system.totalEik, system.fixedFrameworkStoredEik, system.numberOfMoleculesPerComponent);
  }
  time_end = std::chrono::system_clock::now();

  system.mc_moves_cputime[move]["Integration"] += (time_end - time_begin);
  system.mc_moves_statistics.addConstructed(move);

  double drift = std::abs(currentEnergy.conservedEnergy() - referenceEnergy.conservedEnergy());

  // accept or reject based on energy difference
  if (random.uniform() < std::exp(-system.beta * drift))
  {
    system.mc_moves_statistics.addAccepted(move);

    system.moleculePositions = moleculePositions;
    system.thermostat = thermostat;
    system.timeStep = dt;

    std::copy(moleculeAtomPositions.begin(), moleculeAtomPositions.end(), atomPositions.begin());
    system.spanOfMoleculeAtoms() = moleculeAtomPositions;

    Integrators::createCartesianPositions(system.moleculePositions, system.spanOfMoleculeAtoms(), system.components);
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
    return currentEnergy;
  }
  return std::nullopt;
}
