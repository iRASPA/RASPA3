module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#endif

module mc_moves_insertion_cbmc;

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
import <chrono>;
import <cmath>;
import <iostream>;
import <iomanip>;
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
import mc_moves_statistics;
import mc_moves_move_types;
import mc_moves_probabilities;

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::insertionMoveCBMC(RandomNumber& random, System& system,
                                                                             size_t selectedComponent)
{
  std::chrono::system_clock::time_point time_begin, time_end;

  // Update move counts statistics for swap insertion move
  size_t selectedMolecule = system.numberOfMoleculesPerComponent[selectedComponent];
  system.components[selectedComponent].mc_moves_statistics.addTrial(MoveTypes::SwapCBMC, 0);

  // Extract cutoff distances and growth type for the selected component
  double cutOffFrameworkVDW = system.forceField.cutOffFrameworkVDW;
  double cutOffMoleculeVDW = system.forceField.cutOffMoleculeVDW;
  double cutOffCoulomb = system.forceField.cutOffCoulomb;
  Component::GrowType growType = system.components[selectedComponent].growType;

  time_begin = std::chrono::system_clock::now();

  // Attempt to grow a new molecule using CBMC
  std::optional<ChainData> growData = CBMC::growMoleculeSwapInsertion(
      random, system.frameworkComponents, system.components[selectedComponent], system.hasExternalField,
      system.components, system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(),
      system.spanOfMoleculeAtoms(), system.beta, growType, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb,
      selectedComponent, selectedMolecule, 1.0, 0uz, system.numberOfTrialDirections);

  time_end = std::chrono::system_clock::now();

  // Update CPU time statistics for the non-Ewald part of the move
  system.components[selectedComponent].mc_moves_cputime.swapInsertionMoveCBMCNonEwald += (time_end - time_begin);
  system.mc_moves_cputime.swapInsertionMoveCBMCNonEwald += (time_end - time_begin);

  // If growth failed, reject the move
  if (!growData) return {std::nullopt, double3(0.0, 1.0, 0.0)};

  std::span<const Atom> newMolecule = std::span(growData->atom.begin(), growData->atom.end());

  // Check if the new molecule is inside blocked pockets
  if (system.insideBlockedPockets(system.components[selectedComponent], newMolecule))
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  // Update statistics for successfully constructed molecules
  system.components[selectedComponent].mc_moves_statistics.addConstructed(MoveTypes::SwapCBMC, 0);

  time_begin = std::chrono::system_clock::now();

  // Compute energy difference due to Ewald Fourier components
  RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, newMolecule, {});

  time_end = std::chrono::system_clock::now();

  // Update CPU time statistics for the Ewald part of the move
  system.components[selectedComponent].mc_moves_cputime.swapInsertionMoveCBMCEwald += (time_end - time_begin);
  system.mc_moves_cputime.swapInsertionMoveCBMCEwald += (time_end - time_begin);

  time_begin = std::chrono::system_clock::now();

  // Compute tail energy difference due to long-range corrections
  RunningEnergy tailEnergyDifference =
      Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                              system.spanOfMoleculeAtoms(), newMolecule, {}) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.spanOfFrameworkAtoms(), newMolecule, {});

  time_end = std::chrono::system_clock::now();

  // Update CPU time statistics for the tail corrections
  system.components[selectedComponent].mc_moves_cputime.swapInsertionMoveCBMCTail += (time_end - time_begin);
  system.mc_moves_cputime.swapInsertionMoveCBMCTail += (time_end - time_begin);

  // Calculate correction factor for Ewald energy difference
  double correctionFactorEwald =
      std::exp(-system.beta * (energyFourierDifference.potentialEnergy() + tailEnergyDifference.potentialEnergy()));

  // Compute the acceptance probability pre-factor
  double fugacity = system.components[selectedComponent].fugacityCoefficient.value_or(1.0) * system.pressure;
  double idealGasRosenbluthWeight = system.components[selectedComponent].idealGasRosenbluthWeight.value_or(1.0);
  double preFactor = correctionFactorEwald * system.beta * system.components[selectedComponent].molFraction * fugacity *
                     system.simulationBox.volume /
                     double(1 + system.numberOfIntegerMoleculesPerComponent[selectedComponent]);

  // Calculate the acceptance probability Pacc
  double Pacc = preFactor * growData->RosenbluthWeight / idealGasRosenbluthWeight;

  size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];
  double biasTransitionMatrix = system.tmmc.biasFactor(oldN + 1, oldN);

  // Check if TMMC is enabled and macrostate limit is not exceeded
  if (system.tmmc.doTMMC)
  {
    size_t newN = oldN + 1;
    if (newN > system.tmmc.maxMacrostate)
    {
      return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
    }
  }

  // Apply acceptance/rejection criterion
  if (random.uniform() < biasTransitionMatrix * Pacc)
  {
    // Move accepted; update acceptance statistics
    system.components[selectedComponent].mc_moves_statistics.addAccepted(MoveTypes::SwapCBMC, 0);

    // Accept Ewald move and insert the new molecule into the system
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
    system.insertMolecule(selectedComponent, growData->molecule, growData->atom);

    return {growData->energies + energyFourierDifference + tailEnergyDifference, double3(0.0, 1.0 - Pacc, Pacc)};
  };

  // Move rejected
  return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
}
