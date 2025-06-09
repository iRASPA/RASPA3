module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#endif

module mc_moves_volume;

#ifndef USE_LEGACY_HEADERS
import <complex>;
import <vector>;
import <array>;
import <tuple>;
import <optional>;
import <span>;
import <optional>;
import <tuple>;
import <algorithm>;
import <numeric>;
import <chrono>;
import <cmath>;
import <iostream>;
import <iomanip>;
#endif

import component;
import atom;
import molecule;
import int3;
import double3;
import double3x3;
import simd_quatd;
import simulationbox;
import cbmc;
import randomnumbers;
import system;
import energy_factor;
import energy_status;
import energy_status_inter;
import running_energy;
import property_lambda_probability_histogram;
import property_widom;
import averages;
import forcefield;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import mc_moves_move_types;

std::optional<RunningEnergy> MC_Moves::volumeMove(RandomNumber &random, System &system)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  MoveTypes move = MoveTypes::VolumeChange;

  // Update volume move counts
  system.mc_moves_statistics.addTrial(move);

  RunningEnergy oldTotalEnergy = system.runningEnergies;
  // Calculate the total number of molecules
  double numberOfMolecules = static_cast<double>(std::accumulate(system.numberOfIntegerMoleculesPerComponent.begin(),
                                                             system.numberOfIntegerMoleculesPerComponent.end(), 0));
  double oldVolume = system.simulationBox.volume;
  double maxVolumeChange = system.mc_moves_statistics.getMaxChange(move);

  // Propose a new volume change
  double newVolume = std::exp(std::log(oldVolume) + maxVolumeChange * (2.0 * random.uniform() - 1.0));
  // Compute scaling factor for box dimensions
  double scale = std::pow(newVolume / oldVolume, 1.0 / 3.0);

  SimulationBox newBox = system.simulationBox.scaled(scale);
  std::pair<std::vector<Molecule>, std::vector<Atom>> newPositions = system.scaledCenterOfMassPositions(scale);

  double cutOffFrameworkVDW_stored = system.forceField.cutOffFrameworkVDW;
  double cutOffMoleculeVDW_stored = system.forceField.cutOffMoleculeVDW;
  double cutOffCoulomb_stored = system.forceField.cutOffCoulomb;
  double ewald_alpha_stored = system.forceField.EwaldAlpha;
  int3 ewald_k_stored = system.forceField.numberOfWaveVectors;

  system.forceField.initializeAutomaticCutOff(newBox);
  system.forceField.initializeEwaldParameters(newBox);

  time_begin = std::chrono::system_clock::now();
  // Compute new intermolecular energy
  RunningEnergy newTotalInterEnergy =
      Interactions::computeInterMolecularEnergy(system.forceField, newBox, newPositions.second);
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime[move]["NonEwald"] += (time_end - time_begin);

  time_begin = std::chrono::system_clock::now();
  // Compute new tail corrections
  RunningEnergy newTotalTailEnergy =
      Interactions::computeInterMolecularTailEnergy(system.forceField, newBox, newPositions.second);
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);

  time_begin = std::chrono::system_clock::now();
  // Compute new Ewald Fourier energy
  RunningEnergy newTotalEwaldEnergy = Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.totalEik,
      system.forceField, newBox, system.components, system.numberOfMoleculesPerComponent, newPositions.second);
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);

  // Sum up all energy contributions
  RunningEnergy newTotalEnergy = newTotalInterEnergy + newTotalTailEnergy + newTotalEwaldEnergy;

  // Update constructed move counts
  system.mc_moves_statistics.addConstructed(move);

  // Apply acceptance/rejection rule
  if (random.uniform() < std::exp((numberOfMolecules + 1.0) * std::log(newVolume / oldVolume) -
                                  (system.pressure * (newVolume - oldVolume) +
                                   (newTotalEnergy.potentialEnergy() - oldTotalEnergy.potentialEnergy())) *
                                      system.beta))
  {
    // Move accepted: update system state
    system.mc_moves_statistics.addAccepted(move);

    system.simulationBox = newBox;
    std::copy(newPositions.first.begin(), newPositions.first.end(), system.moleculePositions.begin());
    std::copy(newPositions.second.begin(), newPositions.second.end(), system.atomPositions.begin());

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);

    return newTotalEnergy;
  }

  system.forceField.cutOffFrameworkVDW = cutOffFrameworkVDW_stored;
  system.forceField.cutOffMoleculeVDW = cutOffMoleculeVDW_stored;
  system.forceField.cutOffCoulomb = cutOffCoulomb_stored;
  system.forceField.EwaldAlpha = ewald_alpha_stored;
  system.forceField.numberOfWaveVectors = ewald_k_stored;

  return std::nullopt;
}
