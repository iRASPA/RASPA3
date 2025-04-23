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

module mc_moves_reinsertion;

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
import molecule;
import double3;
import double3x3;
import simd_quatd;
import simulationbox;
import cbmc;
import cbmc_chain_data;
import cbmc_interactions;
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

std::optional<RunningEnergy> MC_Moves::reinsertionMove(RandomNumber &random, System &system, size_t selectedComponent,
                                                       size_t selectedMolecule, Molecule &molecule,
                                                       std::span<Atom> molecule_atoms)
{
  // Variables to record timing for performance measurement.
  std::chrono::system_clock::time_point time_begin, time_end;
  MoveTypes move = MoveTypes::ReinsertionCBMC;
  Component &component = system.components[selectedComponent];

  // Increment move counts for reinsertion CBMC statistics.
  component.mc_moves_statistics.addTrial(move);

  // If no molecules of selected component are present, exit the move.
  if (system.numberOfMoleculesPerComponent[selectedComponent] == 0)
  {
    return std::nullopt;
  }

  // Determine cutoff distances based on whether dual cutoff is used.
  double cutOffFrameworkVDW =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffFrameworkVDW;
  double cutOffMoleculeVDW =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffMoleculeVDW;
  double cutOffCoulomb =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffCoulomb;

  time_begin = std::chrono::system_clock::now();
  // Attempt to grow the molecule using CBMC reinsertion.
  std::optional<ChainData> growData = CBMC::growMoleculeReinsertion(
      random, component, system.hasExternalField, system.components, system.forceField, system.simulationBox,
      system.interpolationGrids, system.framework, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
      system.beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, selectedComponent, selectedMolecule, molecule,
      molecule_atoms, system.numberOfTrialDirections);
  time_end = std::chrono::system_clock::now();
  // Record CPU time taken for the non-Ewald part of the move.
  component.mc_moves_cputime[move]["NonEwald"] += (time_end - time_begin);
  system.mc_moves_cputime[move]["NonEwald"] += (time_end - time_begin);

  // If growth was unsuccessful, exit the move.
  if (!growData) return std::nullopt;

  // Get the new molecule configuration.
  std::span<const Atom> newMolecule = std::span(growData->atom.begin(), growData->atom.end());

  // Check if the new molecule is inside blocked pockets; if so, exit the move.
  if (system.insideBlockedPockets(component, newMolecule))
  {
    return std::nullopt;
  }

  // Increment the constructed moves count.
  component.mc_moves_statistics.addConstructed(move);

  // Retrace the old molecule configuration using CBMC retracing.
  time_begin = std::chrono::system_clock::now();
  ChainData retraceData = CBMC::retraceMoleculeReinsertion(
      random, component, system.hasExternalField, system.components, system.forceField, system.simulationBox,
      system.interpolationGrids, system.framework, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
      system.beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, selectedComponent, selectedMolecule, molecule,
      molecule_atoms, growData->storedR, system.numberOfTrialDirections);
  time_end = std::chrono::system_clock::now();

  // Record CPU time taken for the retracing step.
  component.mc_moves_cputime[move]["NonEwald"] += (time_end - time_begin);
  system.mc_moves_cputime[move]["NonEwald"] += (time_end - time_begin);

  // Compute the energy difference in the Fourier space due to Ewald summation.
  time_begin = std::chrono::system_clock::now();
  RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, newMolecule, molecule_atoms);
  time_end = std::chrono::system_clock::now();
  // Record CPU time taken for the Ewald Fourier part of the move.
  component.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);
  system.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);

  double correctionFactorDualCutOff = 1.0;
  std::optional<RunningEnergy> energyNew;
  std::optional<RunningEnergy> energyOld;
  if (system.forceField.useDualCutOff)
  {
    // If dual cutoff is used, compute correction factor due to non-overlapping energies.
    energyNew = CBMC::computeExternalNonOverlappingEnergyDualCutOff(
        component, system.hasExternalField, system.forceField, system.simulationBox, system.interpolationGrids,
        system.framework, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
        system.forceField.cutOffFrameworkVDW, system.forceField.cutOffMoleculeVDW, system.forceField.cutOffCoulomb,
        growData->atom);
    energyOld = CBMC::computeExternalNonOverlappingEnergyDualCutOff(
        component, system.hasExternalField, system.forceField, system.simulationBox, system.interpolationGrids,
        system.framework, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
        system.forceField.cutOffFrameworkVDW, system.forceField.cutOffMoleculeVDW, system.forceField.cutOffCoulomb,
        retraceData.atom);
    correctionFactorDualCutOff =
        std::exp(-system.beta * (energyNew->potentialEnergy() - growData->energies.potentialEnergy() -
                                 (energyOld->potentialEnergy() - retraceData.energies.potentialEnergy())));
  }

  // Compute correction factor from the Fourier energy difference.
  double correctionFactorFourier = std::exp(-system.beta * energyFourierDifference.potentialEnergy());

  // Apply Metropolis acceptance criterion.
  if (random.uniform() <
      correctionFactorDualCutOff * correctionFactorFourier * growData->RosenbluthWeight / retraceData.RosenbluthWeight)
  {
    // Move is accepted; update statistics and state.
    component.mc_moves_statistics.addAccepted(move);

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
    std::copy(newMolecule.begin(), newMolecule.end(), molecule_atoms.begin());
    molecule = growData->molecule;

    if (system.forceField.useDualCutOff)
    {
      return (energyNew.value() - energyOld.value()) + energyFourierDifference;
    }

    return (growData->energies - retraceData.energies) + energyFourierDifference;
  };

  // Move is rejected.
  return std::nullopt;
}
