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
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#endif

module mc_moves_insertion_cbmc;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import component;
import molecule;
import atom;
import double3;
import double3x3;
import simd_quatd;
import simulationbox;
import cbmc;
import cbmc_chain_data;
import randomnumbers;
import system;
import energy_factor;
import energy_status;
import energy_status_inter;
import property_lambda_probability_histogram;
import property_widom;
import averages;
import running_energy;
import forcefield;
import transition_matrix;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import interactions_polarization;
import mc_moves_move_types;

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::insertionMoveCBMC(RandomNumber& random, System& system,
                                                                             std::size_t selectedComponent)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  MoveTypes move = MoveTypes::SwapCBMC;
  Component& component = system.components[selectedComponent];

  // Update move counts statistics for swap insertion move
  std::size_t selectedMolecule = system.numberOfMoleculesPerComponent[selectedComponent];
  component.mc_moves_statistics.addTrial(move, 0);

  // Extract cutoff distances and growth type for the selected component
  double cutOffFrameworkVDW = system.forceField.cutOffFrameworkVDW;
  double cutOffMoleculeVDW = system.forceField.cutOffMoleculeVDW;
  double cutOffCoulomb = system.forceField.cutOffCoulomb;
  Component::GrowType growType = component.growType;

  // Attempt to grow a new molecule using CBMC
  time_begin = std::chrono::system_clock::now();
  std::optional<ChainData> growData = CBMC::growMoleculeSwapInsertion(
      random, component, system.hasExternalField, system.components, system.forceField, system.simulationBox,
      system.interpolationGrids, system.framework, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
      system.beta, growType, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, selectedComponent, selectedMolecule,
      1.0, false, false, system.numberOfTrialDirections);
  time_end = std::chrono::system_clock::now();

  // Update CPU time statistics for the non-Ewald part of the move
  component.mc_moves_cputime[move]["NonEwald"] += (time_end - time_begin);
  system.mc_moves_cputime[move]["NonEwald"] += (time_end - time_begin);

  // If growth failed, reject the move
  if (!growData) return {std::nullopt, double3(0.0, 1.0, 0.0)};

  std::span<const Atom> newMolecule = std::span(growData->atom.begin(), growData->atom.end());

  // Check if the new molecule is inside blocked pockets
  if (system.insideBlockedPockets(system.components[selectedComponent], newMolecule))
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  // Update statistics for successfully constructed molecules
  system.components[selectedComponent].mc_moves_statistics.addConstructed(move, 0);

  // Compute energy difference due to Ewald Fourier components
  time_begin = std::chrono::system_clock::now();
  RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, newMolecule, {});
  time_end = std::chrono::system_clock::now();

  // Update CPU time statistics for the Ewald part of the move
  component.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);
  system.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);

  // Compute tail energy difference due to long-range corrections
  time_begin = std::chrono::system_clock::now();
  RunningEnergy tailEnergyDifference =
      Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                              system.spanOfMoleculeAtoms(), newMolecule, {}) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.spanOfFrameworkAtoms(), newMolecule, {});
  time_end = std::chrono::system_clock::now();

  // Update CPU time statistics for the tail corrections
  component.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);
  system.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);

  RunningEnergy polarizationDifference;
  if (system.forceField.computePolarization)
  {
    Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                  system.spanOfFrameworkAtoms(),
                                                                  growData->electricField, {}, growData->atom, {});

    Interactions::computeEwaldFourierElectricFieldDifference(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.totalEik, system.forceField, system.simulationBox, growData->electricField, {}, growData->atom, {});

    // Compute polarization energy difference
    polarizationDifference = Interactions::computePolarizationEnergyDifference(
        system.forceField, growData->electricField, {}, growData->atom, {});
  }

  // Calculate correction factor for Ewald energy difference
  double correctionFactorEwald =
      std::exp(-system.beta * (energyFourierDifference.potentialEnergy() + tailEnergyDifference.potentialEnergy() +
                               polarizationDifference.potentialEnergy()));

  // Compute the acceptance probability pre-factor
  double fugacity = component.fugacityCoefficient.value_or(1.0) * system.pressure;
  double idealGasRosenbluthWeight = component.idealGasRosenbluthWeight.value_or(1.0);
  double preFactor = correctionFactorEwald * system.beta * component.molFraction * fugacity *
                     system.simulationBox.volume /
                     double(1 + system.numberOfIntegerMoleculesPerComponent[selectedComponent]);

  // Calculate the acceptance probability Pacc
  double Pacc = preFactor * growData->RosenbluthWeight / idealGasRosenbluthWeight;

  std::size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];
  double biasTransitionMatrix = system.tmmc.biasFactor(oldN + 1, oldN);

  // Check if TMMC is enabled and macrostate limit is not exceeded
  if (system.tmmc.doTMMC)
  {
    std::size_t newN = oldN + 1;
    if (newN > system.tmmc.maxMacrostate)
    {
      return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
    }
  }

  // Apply acceptance/rejection criterion
  if (random.uniform() < biasTransitionMatrix * Pacc)
  {
    // Move accepted; update acceptance statistics
    component.mc_moves_statistics.addAccepted(move, 0);

    // Accept Ewald move and insert the new molecule into the system
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
    system.insertMoleculePolarization(selectedComponent, growData->molecule, growData->atom, growData->electricField);

    return {growData->energies + energyFourierDifference + tailEnergyDifference + polarizationDifference,
            double3(0.0, 1.0 - Pacc, Pacc)};
  };

  // Move rejected
  return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
}
