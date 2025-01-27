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
import move_statistics;
import mc_moves_probabilities_particles;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;

std::optional<RunningEnergy> MC_Moves::reinsertionMove(RandomNumber &random, System &system, size_t selectedComponent,
                                                       size_t selectedMolecule, Molecule &molecule,
                                                       std::span<Atom> molecule_atoms)
{
  // Variables to record timing for performance measurement.
  std::chrono::system_clock::time_point time_begin, time_end;

  // Increment move counts for reinsertion CBMC statistics.
  system.components[selectedComponent].mc_moves_statistics.reinsertionMove_CBMC.counts += 1;
  system.components[selectedComponent].mc_moves_statistics.reinsertionMove_CBMC.totalCounts += 1;

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
      random, system.frameworkComponents, system.components[selectedComponent], system.hasExternalField,
      system.components, system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(),
      system.spanOfMoleculeAtoms(), system.beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb,
      selectedComponent, selectedMolecule, molecule, molecule_atoms, system.numberOfTrialDirections);
  time_end = std::chrono::system_clock::now();
  // Record CPU time taken for the non-Ewald part of the move.
  system.components[selectedComponent].mc_moves_cputime.reinsertionMoveCBMCNonEwald += (time_end - time_begin);
  system.mc_moves_cputime.reinsertionMoveCBMCNonEwald += (time_end - time_begin);

  // If growth was unsuccessful, exit the move.
  if (!growData) return std::nullopt;

  // Get the new molecule configuration.
  std::span<const Atom> newMolecule = std::span(growData->atom.begin(), growData->atom.end());

  // Check if the new molecule is inside blocked pockets; if so, exit the move.
  if (system.insideBlockedPockets(system.components[selectedComponent], newMolecule))
  {
    return std::nullopt;
  }

  // Increment the constructed moves count.
  system.components[selectedComponent].mc_moves_statistics.reinsertionMove_CBMC.constructed += 1;
  system.components[selectedComponent].mc_moves_statistics.reinsertionMove_CBMC.totalConstructed += 1;

  time_begin = std::chrono::system_clock::now();
  // Retrace the old molecule configuration using CBMC retracing.
  ChainData retraceData = CBMC::retraceMoleculeReinsertion(
      random, system.frameworkComponents, system.components[selectedComponent], system.hasExternalField,
      system.components, system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(),
      system.spanOfMoleculeAtoms(), system.beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb,
      selectedComponent, selectedMolecule, molecule, molecule_atoms, growData->storedR, system.numberOfTrialDirections);
  time_end = std::chrono::system_clock::now();
  // Record CPU time taken for the retracing step.
  system.components[selectedComponent].mc_moves_cputime.reinsertionMoveCBMCNonEwald += (time_end - time_begin);
  system.mc_moves_cputime.reinsertionMoveCBMCNonEwald += (time_end - time_begin);

  time_begin = std::chrono::system_clock::now();
  // Compute the energy difference in the Fourier space due to Ewald summation.
  RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, newMolecule, molecule_atoms);
  time_end = std::chrono::system_clock::now();
  // Record CPU time taken for the Ewald Fourier part of the move.
  system.components[selectedComponent].mc_moves_cputime.reinsertionMoveCBMCEwald += (time_end - time_begin);
  system.mc_moves_cputime.reinsertionMoveCBMCEwald += (time_end - time_begin);

  double correctionFactorDualCutOff = 1.0;
  std::optional<RunningEnergy> energyNew;
  std::optional<RunningEnergy> energyOld;
  if (system.forceField.useDualCutOff)
  {
    // If dual cutoff is used, compute correction factor due to non-overlapping energies.
    energyNew = CBMC::computeExternalNonOverlappingEnergyDualCutOff(
        system.frameworkComponents, system.components[selectedComponent], system.hasExternalField, system.forceField,
        system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
        system.forceField.cutOffFrameworkVDW, system.forceField.cutOffMoleculeVDW, system.forceField.cutOffCoulomb,
        growData->atom);
    energyOld = CBMC::computeExternalNonOverlappingEnergyDualCutOff(
        system.frameworkComponents, system.components[selectedComponent], system.hasExternalField, system.forceField,
        system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
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
    system.components[selectedComponent].mc_moves_statistics.reinsertionMove_CBMC.accepted += 1;
    system.components[selectedComponent].mc_moves_statistics.reinsertionMove_CBMC.totalAccepted += 1;

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
