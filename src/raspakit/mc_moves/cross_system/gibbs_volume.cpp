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

module mc_moves_gibbs_volume;

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
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import mc_moves_move_types;

std::optional<std::pair<RunningEnergy, RunningEnergy>> MC_Moves::GibbsVolumeMove(RandomNumber &random, System &systemA,
                                                                                 System &systemB)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  MoveTypes move = MoveTypes::GibbsVolume;

  systemA.mc_moves_statistics.addTrial(move);
  systemB.mc_moves_statistics.addTrial(move);

  // determine New box-volumes leaving the total volume constant
  double oldVolumeA = systemA.simulationBox->volume;
  double maxVolumeChangeA = systemA.mc_moves_statistics.getMaxChange(move);
  double oldVolumeB = systemB.simulationBox->volume;
  double totalVolume = oldVolumeA + oldVolumeB;
  double expdv = std::exp(std::log(oldVolumeA / oldVolumeB) + maxVolumeChangeA * (2.0 * random.uniform() - 1.0));
  double newVolumeA = expdv * totalVolume / (1.0 + expdv);
  double newVolumeB = totalVolume - newVolumeA;

  // Scale systemA to new volume and calculate new positions
  RunningEnergy oldTotalEnergyA = systemA.runningEnergies;
  double numberOfMoleculesA = static_cast<double>(std::reduce(systemA.numberOfIntegerMoleculesPerComponent.begin(),
                                                              systemA.numberOfIntegerMoleculesPerComponent.end()));
  double scaleA = std::pow(newVolumeA / oldVolumeA, 1.0 / 3.0);
  SimulationBox newBoxA = systemA.simulationBox->scaled(scaleA);
  std::pair<std::vector<Molecule>, std::vector<Atom>> newPositionsA = systemA.scaledCenterOfMassPositions(scaleA);

  // Compute new intermolecular energy for systemA
  time_begin = std::chrono::system_clock::now();
  RunningEnergy newTotalInterEnergyA =
      Interactions::computeInterMolecularEnergy(*systemA.forceField, newBoxA, newPositionsA.second);
  time_end = std::chrono::system_clock::now();
  systemA.mc_moves_cputime[move]["NonEwald"] += (time_end - time_begin);

  // Compute new tail corrections for systemA
  time_begin = std::chrono::system_clock::now();
  RunningEnergy newTotalTailEnergyA =
      Interactions::computeInterMolecularTailEnergy(*systemA.forceField, newBoxA, newPositionsA.second);
  time_end = std::chrono::system_clock::now();
  systemA.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);

  // Compute new Ewald Fourier energy for systemA
  time_begin = std::chrono::system_clock::now();
  RunningEnergy newTotalEwaldEnergyA = Interactions::computeEwaldFourierEnergy(
      systemA.eik_x, systemA.eik_y, systemA.eik_z, systemA.eik_xy, systemA.fixedFrameworkStoredEik, systemA.totalEik,
      *systemA.forceField, newBoxA, systemA.components, systemA.numberOfMoleculesPerComponent, newPositionsA.second);
  time_end = std::chrono::system_clock::now();
  systemA.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);

  // Update energy and statistics for systemA
  RunningEnergy newTotalEnergyA = newTotalInterEnergyA + newTotalTailEnergyA + newTotalEwaldEnergyA;

  systemA.mc_moves_statistics.addConstructed(move);

  // Scale systemB to new volume and calculate new positions
  RunningEnergy oldTotalEnergyB = systemB.runningEnergies;
  double numberOfMoleculesB = static_cast<double>(std::reduce(systemB.numberOfIntegerMoleculesPerComponent.begin(),
                                                              systemB.numberOfIntegerMoleculesPerComponent.end()));
  double scaleB = std::pow(newVolumeB / oldVolumeB, 1.0 / 3.0);
  SimulationBox newBoxB = systemB.simulationBox->scaled(scaleB);
  std::pair<std::vector<Molecule>, std::vector<Atom>> newPositionsB = systemB.scaledCenterOfMassPositions(scaleB);

  // Compute new intermolecular energy for systemB
  time_begin = std::chrono::system_clock::now();
  RunningEnergy newTotalInterEnergyB =
      Interactions::computeInterMolecularEnergy(*systemB.forceField, newBoxB, newPositionsB.second);
  time_end = std::chrono::system_clock::now();
  systemA.mc_moves_cputime[move]["NonEwald"] += (time_end - time_begin);

  // Compute new tail corrections for systemB
  time_begin = std::chrono::system_clock::now();
  RunningEnergy newTotalTailEnergyB =
      Interactions::computeInterMolecularTailEnergy(*systemB.forceField, newBoxB, newPositionsB.second);
  time_end = std::chrono::system_clock::now();
  systemA.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);

  // Compute new Ewald Fourier energy for systemB
  time_begin = std::chrono::system_clock::now();
  RunningEnergy newTotalEwaldEnergyB = Interactions::computeEwaldFourierEnergy(
      systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy, systemB.fixedFrameworkStoredEik, systemB.totalEik,
      *systemB.forceField, newBoxB, systemB.components, systemB.numberOfMoleculesPerComponent, newPositionsB.second);
  time_end = std::chrono::system_clock::now();
  systemA.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);

  // Update energy and statistics for systemB
  RunningEnergy newTotalEnergyB = newTotalInterEnergyB + newTotalTailEnergyB + newTotalEwaldEnergyB;

  systemB.mc_moves_statistics.addConstructed(move);

  // Calculate total energy change deltaU
  double deltaU = (newTotalEnergyA.potentialEnergy() - oldTotalEnergyA.potentialEnergy()) +
                  (newTotalEnergyB.potentialEnergy() - oldTotalEnergyB.potentialEnergy());

  // apply acceptance/rejection rule
  if (random.uniform() < std::exp(-systemA.beta * deltaU +
                                  (static_cast<double>(numberOfMoleculesA + 1.0) * std::log(newVolumeA / oldVolumeA)) +
                                  (static_cast<double>(numberOfMoleculesB + 1.0) * std::log(newVolumeB / oldVolumeB))))
  {
    // Accept the move: update systems A and B with new configurations and energies
    systemA.mc_moves_statistics.addAccepted(move);

    systemA.simulationBox = std::make_shared<SimulationBox>(std::move(newBoxA));
    std::copy(newPositionsA.first.begin(), newPositionsA.first.end(), systemA.moleculePositions.begin());
    std::copy(newPositionsA.second.begin(), newPositionsA.second.end(), systemA.atomPositions.begin());
    Interactions::acceptEwaldMove(*systemA.forceField, systemA.storedEik, systemA.totalEik);

    systemB.mc_moves_statistics.addAccepted(move);

    systemB.simulationBox = std::make_shared<SimulationBox>(std::move(newBoxB));
    std::copy(newPositionsB.first.begin(), newPositionsB.first.end(), systemB.moleculePositions.begin());
    std::copy(newPositionsB.second.begin(), newPositionsB.second.end(), systemB.atomPositions.begin());
    Interactions::acceptEwaldMove(*systemB.forceField, systemB.storedEik, systemB.totalEik);

    return std::make_pair(newTotalEnergyA, newTotalEnergyB);
  }

  return std::nullopt;
}
