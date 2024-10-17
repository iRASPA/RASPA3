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

std::optional<RunningEnergy> MC_Moves::hybridMCMove(RandomNumber& random, System& system)
{
  std::chrono::system_clock::time_point time_begin, time_end;

  system.mc_moves_statistics.hybridMC.counts += 1;
  system.mc_moves_statistics.hybridMC.totalCounts += 1;

  // copy current data
  std::span<Atom> spanReferencePositions = system.spanOfMoleculeAtoms();
  std::vector<Atom> referenceAtomPositions(spanReferencePositions.size());
  std::copy(spanReferencePositions.begin(), spanReferencePositions.end(), referenceAtomPositions.begin());
  std::vector<Molecule> referenceMoleculePositions(system.moleculePositions);

  // eik_* fixedFrameworkStoredEik
  // frameworkMoleculeEnergy
  // intermolecularEnergy
  // ewaldEnergy
  //

  // initialize system for md move
  // system.createCartesianPositions();
  // system.initializeVelocities(random);
  // system.removeCenterOfMassVelocityDrift();
  // system.computeTotalGradients();
  // system.runningEnergies.translationalKineticEnergy = system.computeTranslationalKineticEnergy();
  // system.runningEnergies.rotationalKineticEnergy = Integrators::computeRotationalKineticEnergy(moleculePositions);

  // save current energy
  RunningEnergy referenceEnergy = system.runningEnergies;

  // integrate for n steps
  for (size_t step = 0; step < system.numberOfHybridMCSteps; ++step)
  {
    // system.integrate();
  }

  system.mc_moves_statistics.hybridMC.constructed += 1;
  system.mc_moves_statistics.hybridMC.totalConstructed += 1;

  // accept or reject based on energy difference
  if (random.uniform() <
      std::exp(-system.beta * (referenceEnergy.conservedEnergy() - system.runningEnergies.conservedEnergy())))
  {
    system.mc_moves_statistics.hybridMC.accepted += 1;
    system.mc_moves_statistics.hybridMC.totalAccepted += 1;
    return system.runningEnergies;
  }
  else
  {
    // reset everything to original state.
    std::copy(referenceAtomPositions.begin(), referenceAtomPositions.end(), spanReferencePositions.begin());
    std::copy(referenceMoleculePositions.begin(), referenceMoleculePositions.end(), system.moleculePositions.begin());
    system.runningEnergies = referenceEnergy;
  }

  return std::nullopt;
}