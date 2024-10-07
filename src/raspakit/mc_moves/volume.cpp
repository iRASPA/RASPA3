module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
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
import move_statistics;
import mc_moves_probabilities_particles;
import mc_moves_probabilities_system;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;

std::optional<RunningEnergy> MC_Moves::volumeMove(RandomNumber &random, System &system)
{
  std::chrono::system_clock::time_point time_begin, time_end;

  system.mc_moves_statistics.volumeMove.counts += 1;
  system.mc_moves_statistics.volumeMove.totalCounts += 1;

  RunningEnergy oldTotalEnergy = system.runningEnergies;
  double numberOfMolecules = static_cast<double>(std::reduce(system.numberOfIntegerMoleculesPerComponent.begin(),
                                                             system.numberOfIntegerMoleculesPerComponent.end()));
  double oldVolume = system.simulationBox.volume;
  double maxVolumeChange = system.mc_moves_statistics.volumeMove.maxChange;
  double newVolume = std::exp(std::log(oldVolume) + maxVolumeChange * (2.0 * random.uniform() - 1.0));
  double scale = std::pow(newVolume / oldVolume, 1.0 / 3.0);

  SimulationBox newBox = system.simulationBox.scaled(scale);
  std::pair<std::vector<Molecule>, std::vector<Atom>> newPositions = system.scaledCenterOfMassPositions(scale);

  time_begin = std::chrono::system_clock::now();
  RunningEnergy newTotalInterEnergy =
      Interactions::computeInterMolecularEnergy(system.forceField, newBox, newPositions.second);
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime.volumeMoveNonEwald += (time_end - time_begin);

  time_begin = std::chrono::system_clock::now();
  RunningEnergy newTotalTailEnergy =
      Interactions::computeInterMolecularTailEnergy(system.forceField, newBox, newPositions.second);
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime.volumeMoveTail += (time_end - time_begin);

  time_begin = std::chrono::system_clock::now();
  RunningEnergy newTotalEwaldEnergy = Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.totalEik,
      system.forceField, newBox, system.components, system.numberOfMoleculesPerComponent, newPositions.second);
  time_end = std::chrono::system_clock::now();
  system.mc_moves_cputime.volumeMoveEwald += (time_end - time_begin);

  RunningEnergy newTotalEnergy = newTotalInterEnergy + newTotalTailEnergy + newTotalEwaldEnergy;

  system.mc_moves_statistics.volumeMove.constructed += 1;
  system.mc_moves_statistics.volumeMove.totalConstructed += 1;

  // apply acceptance/rejection rule
  if (random.uniform() < std::exp((numberOfMolecules + 1.0) * std::log(newVolume / oldVolume) -
                                  (system.pressure * (newVolume - oldVolume) +
                                   (newTotalEnergy.potentialEnergy() - oldTotalEnergy.potentialEnergy())) *
                                      system.beta))
  {
    system.mc_moves_statistics.volumeMove.accepted += 1;
    system.mc_moves_statistics.volumeMove.totalAccepted += 1;

    system.simulationBox = newBox;
    std::copy(newPositions.first.begin(), newPositions.first.end(), system.moleculePositions.begin());
    std::copy(newPositions.second.begin(), newPositions.second.end(), system.atomPositions.begin());

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);

    return newTotalEnergy;
  }

  return std::nullopt;
}
