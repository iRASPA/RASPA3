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

module mc_moves_parallel_tempering_swap;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import component;
import atom;
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

std::optional<std::pair<RunningEnergy, RunningEnergy>> MC_Moves::ParallelTemperingSwap(RandomNumber &random,
                                                                                       System &systemA, System &systemB)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  MoveTypes move = MoveTypes::ParallelTempering;

  // Update swap move counts for both systems
  systemA.mc_moves_statistics.addTrial(move);

  double acc = 0.0;

  if (systemA.forceField != systemB.forceField)
  {
    // Compute energy of system A using system B's force field
    time_begin = std::chrono::system_clock::now();
    RunningEnergy systemAHamiltonianB =
        Interactions::computeInterMolecularEnergy(systemB.forceField, systemA.simulationBox, systemA.atomData);
    RunningEnergy systemBHamiltonianA =
        Interactions::computeInterMolecularEnergy(systemA.forceField, systemB.simulationBox, systemB.atomData);
    time_end = std::chrono::system_clock::now();

    systemA.mc_moves_cputime[move]["Energy"] += (time_end - time_begin);

    // Calculate acceptance probability when force fields differ
    acc = std::exp(-systemA.beta * (systemBHamiltonianA.potentialEnergy() - systemA.runningEnergies.potentialEnergy()) -
                   systemB.beta * (systemAHamiltonianB.potentialEnergy() - systemB.runningEnergies.potentialEnergy()));
  }
  else
  {
    // Calculate acceptance probability when force fields are the same
    acc = std::exp((systemB.beta - systemA.beta) *
                   (systemB.runningEnergies.potentialEnergy() - systemA.runningEnergies.potentialEnergy()));
  }

  if (systemA.pressure != systemB.pressure)
  {
    /// Ref: "Hyper-parallel tempering Monte Carlo: Application to the Lennard-Jones fluid and the
    /// restricted primitive model",  G. Yan and J.J. de Pablo, JCP, 111(21): 9509-9516, 1999

    // Adjust acceptance probability for pressure differences
    time_begin = std::chrono::system_clock::now();
    acc *= std::pow(systemB.pressure / systemA.pressure,
                    systemB.loadings.totalNumberOfMolecules - systemA.loadings.totalNumberOfMolecules);
    time_end = std::chrono::system_clock::now();

    systemA.mc_moves_cputime[move]["Fugacity"] += (time_end - time_begin);
  }

  // Update constructed move counts for both systems
  systemA.mc_moves_statistics.addConstructed(move);

  // Apply acceptance/rejection rule
  if (random.uniform() < acc)
  {
    // Update accepted move counts for both systems
    systemA.mc_moves_statistics.addAccepted(move);

    // Swap configurations and properties between systems
    std::swap(systemA.atomData, systemB.atomData);
    std::swap(systemA.simulationBox, systemB.simulationBox);
    std::swap(systemA.numberOfMoleculesPerComponent, systemB.numberOfMoleculesPerComponent);
    std::swap(systemA.numberOfIntegerMoleculesPerComponent, systemB.numberOfIntegerMoleculesPerComponent);
    std::swap(systemA.numberOfPseudoAtoms, systemB.numberOfPseudoAtoms);
    std::swap(systemA.totalNumberOfPseudoAtoms, systemB.totalNumberOfPseudoAtoms);
    std::swap(systemA.runningEnergies, systemB.runningEnergies);
    std::swap(systemA.averageEnergies, systemB.averageEnergies);
    std::swap(systemA.mc_moves_probabilities, systemB.mc_moves_probabilities);
    std::swap(systemA.mc_moves_statistics, systemB.mc_moves_statistics);
    std::swap(systemA.mc_moves_cputime, systemB.mc_moves_cputime);

    return std::make_pair(systemA.runningEnergies, systemB.runningEnergies);
  }

  return std::nullopt;
}
