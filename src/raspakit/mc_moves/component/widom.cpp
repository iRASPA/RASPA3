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

module mc_moves_widom;

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
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import mc_moves_move_types;

std::pair<double, double> MC_Moves::WidomMove(RandomNumber& random, System& system, size_t selectedComponent)
{
  size_t selectedMolecule = system.numberOfMoleculesPerComponent[selectedComponent];
  MoveTypes move = MoveTypes::Widom;
  Component& component = system.components[selectedComponent];

  // Update move statistics for Widom insertion move.
  component.mc_moves_statistics.addTrial(move);

  double cutOffFrameworkVDW = system.forceField.cutOffFrameworkVDW;
  double cutOffMoleculeVDW = system.forceField.cutOffMoleculeVDW;
  double cutOffCoulomb = system.forceField.cutOffCoulomb;
  Component::GrowType growType = component.growType;

  // Record the start time for the CBMC molecule growth.
  std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
  // Attempt to grow a new molecule using Configurational Bias Monte Carlo (CBMC) insertion.
  std::optional<ChainData> growData = CBMC::growMoleculeSwapInsertion(
      random, system.frameworkComponents, component, system.hasExternalField, system.components, system.forceField,
      system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(), system.beta, growType,
      cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, selectedComponent, selectedMolecule, 1.0, 0uz,
      system.numberOfTrialDirections);
  // Record the end time for the CBMC molecule growth.
  std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
  // Update CPU time statistics for CBMC non-Ewald calculations.
  component.mc_moves_cputime[move]["NonEwald"] += (t2 - t1);
  system.mc_moves_cputime[move]["NonEwald"] += (t2 - t1);

  // If molecule growth failed, terminate the move.
  if (!growData) return {0.0, 0.0};

  [[maybe_unused]] std::span<const Atom> newMolecule = std::span(growData->atom.begin(), growData->atom.end());

  // Check if the new molecule is inside blocked pockets; if so, abort the move.
  if (system.insideBlockedPockets(component, newMolecule))
  {
    return {0.0, 0.0};
  }

  // Update statistics for successfully constructed molecules.
  component.mc_moves_statistics.addConstructed(move);

  // Record start time for Ewald Fourier energy difference calculation.
  std::chrono::system_clock::time_point u1 = std::chrono::system_clock::now();
  // Compute the energy difference in Ewald Fourier space due to the new molecule.
  RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, newMolecule, {});
  // Record end time for Ewald Fourier energy difference calculation.
  std::chrono::system_clock::time_point u2 = std::chrono::system_clock::now();
  // Update CPU time statistics for Ewald calculations.
  component.mc_moves_cputime[move]["Ewald"] += (u2 - u1);
  system.mc_moves_cputime[move]["Ewald"] += (u2 - u1);

  // Compute the tail corrections for the energy due to the new molecule.
  RunningEnergy tailEnergyDifference =
      Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                              system.spanOfMoleculeAtoms(), newMolecule, {}) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.spanOfFrameworkAtoms(), newMolecule, {});

  // Compute the correction factor from Ewald and tail energy differences.
  double correctionFactorEwald =
      // std::exp(-system.beta * (energyFourierDifference.potentialEnergy()));
      std::exp(-system.beta * (energyFourierDifference.potentialEnergy() + tailEnergyDifference.potentialEnergy()));

  double idealGasRosenbluthWeight = component.idealGasRosenbluthWeight.value_or(1.0);

  // Return the Widom insertion weight for the new molecule.
  return {correctionFactorEwald * growData->RosenbluthWeight / idealGasRosenbluthWeight, 0.0};
  // return {correctionFactorEwald * growData->RosenbluthWeight / idealGasRosenbluthWeight,
  // tailEnergyDifference.potentialEnergy()};
}
