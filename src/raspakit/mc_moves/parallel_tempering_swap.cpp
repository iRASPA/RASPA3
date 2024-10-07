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

module mc_moves_parallel_tempering_swap;

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
import mc_moves_probabilities_system;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;

std::optional<std::pair<RunningEnergy, RunningEnergy>> MC_Moves::ParallelTemperingSwap(RandomNumber &random,
                                                                                       System &systemA, System &systemB)
{
  std::chrono::system_clock::time_point time_begin, time_end;

  systemA.mc_moves_statistics.ParallelTemperingSwap.counts += 1;
  systemA.mc_moves_statistics.ParallelTemperingSwap.totalCounts += 1;
  systemB.mc_moves_statistics.ParallelTemperingSwap.counts += 1;
  systemB.mc_moves_statistics.ParallelTemperingSwap.totalCounts += 1;

  double acc = 0.0;

  if (systemA.forceField != systemB.forceField)
  {
    time_begin = std::chrono::system_clock::now();
    RunningEnergy systemAHamiltonianB =
        Interactions::computeInterMolecularEnergy(systemB.forceField, systemA.simulationBox, systemA.atomPositions);
    time_end = std::chrono::system_clock::now();
    systemA.mc_moves_cputime.ParallelTemperingSwapEnergy += (time_end - time_begin);

    time_begin = std::chrono::system_clock::now();
    RunningEnergy systemBHamiltonianA =
        Interactions::computeInterMolecularEnergy(systemA.forceField, systemB.simulationBox, systemB.atomPositions);
    time_end = std::chrono::system_clock::now();
    systemB.mc_moves_cputime.ParallelTemperingSwapEnergy += (time_end - time_begin);

    acc = std::exp(-systemA.beta * (systemBHamiltonianA.potentialEnergy() - systemA.runningEnergies.potentialEnergy()) -
                   systemB.beta * (systemAHamiltonianB.potentialEnergy() - systemB.runningEnergies.potentialEnergy()));
  }
  else
  {
    acc = std::exp((systemB.beta - systemA.beta) *
                   (systemB.runningEnergies.potentialEnergy() - systemA.runningEnergies.potentialEnergy()));
  }

  if (systemA.pressure != systemB.pressure)
  {
    /// Ref: "Hyper-parallel tempering Monte Carlo: Appliation to the Lennard-Jones fluid and the
    /// restricted primitive model",  G. Yan and J.J. de Pablo, JCP, 111(21): 9509-9516, 1999
    time_begin = std::chrono::system_clock::now();
    acc *= std::pow(systemB.pressure / systemA.pressure,
                    systemB.loadings.totalNumberOfMolecules - systemA.loadings.totalNumberOfMolecules);
    time_end = std::chrono::system_clock::now();
    systemA.mc_moves_cputime.ParallelTemperingSwapFugacity += (time_end - time_begin);
    systemB.mc_moves_cputime.ParallelTemperingSwapFugacity += (time_end - time_begin);
  }

  systemA.mc_moves_statistics.ParallelTemperingSwap.constructed += 1;
  systemA.mc_moves_statistics.ParallelTemperingSwap.totalConstructed += 1;
  systemB.mc_moves_statistics.ParallelTemperingSwap.constructed += 1;
  systemB.mc_moves_statistics.ParallelTemperingSwap.totalConstructed += 1;

  // apply acceptance/rejection rule
  if (random.uniform() < acc)
  {
    systemA.mc_moves_statistics.ParallelTemperingSwap.accepted += 1;
    systemA.mc_moves_statistics.ParallelTemperingSwap.totalAccepted += 1;

    systemB.mc_moves_statistics.ParallelTemperingSwap.accepted += 1;
    systemB.mc_moves_statistics.ParallelTemperingSwap.totalAccepted += 1;

    std::swap(systemA.atomPositions, systemB.atomPositions);
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
    std::swap(systemA.mc_moves_count, systemB.mc_moves_count);

    return std::make_pair(systemA.runningEnergies, systemB.runningEnergies);
  }

  return std::nullopt;
}
