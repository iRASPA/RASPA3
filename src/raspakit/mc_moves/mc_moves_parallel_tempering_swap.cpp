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

  RunningEnergy systemAHamiltonianB;
  time_begin = std::chrono::system_clock::now();
  Interactions::computeInterMolecularEnergy(systemB.forceField, systemA.simulationBox, systemA.atomPositions,
                                            systemAHamiltonianB);
  time_end = std::chrono::system_clock::now();
  systemA.mc_moves_cputime.ParallelTemperingSwapEnergy += (time_end - time_begin);
  systemA.mc_moves_statistics.ParallelTemperingSwap.constructed += 1;
  systemA.mc_moves_statistics.ParallelTemperingSwap.totalConstructed += 1;

  RunningEnergy systemBHamiltonianA;
  time_begin = std::chrono::system_clock::now();
  Interactions::computeInterMolecularEnergy(systemA.forceField, systemB.simulationBox, systemB.atomPositions,
                                            systemBHamiltonianA);
  time_end = std::chrono::system_clock::now();
  systemB.mc_moves_cputime.ParallelTemperingSwapEnergy += (time_end - time_begin);
  systemB.mc_moves_statistics.ParallelTemperingSwap.constructed += 1;
  systemB.mc_moves_statistics.ParallelTemperingSwap.totalConstructed += 1;

  double acc = std::exp(-systemA.beta * (systemBHamiltonianA.total() - systemA.runningEnergies.total()) -
                        systemB.beta * (systemAHamiltonianB.total() - systemB.runningEnergies.total()));

  if (systemA.init_pressure != systemB.init_pressure)
  {
    /// Ref: "Hyper-parallel tempering Monte Carlo: Appliation to the Lennard-Jones fluid and the
    /// restricted primitive model",  G. Yan and J.J. de Pablo, JCP, 111(21): 9509-9516, 1999
    acc *= std::pow(systemB.pressure / systemA.pressure,
                    systemB.loadings.totalNumberOfMolecules - systemA.loadings.totalNumberOfMolecules);
  }

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

    return std::make_pair(systemAHamiltonianB, systemBHamiltonianA);
  }

  return std::nullopt;
}
