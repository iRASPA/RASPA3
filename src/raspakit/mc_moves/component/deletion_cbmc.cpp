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

module mc_moves_deletion_cbmc;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double3;
import double3x3;
import simd_quatd;
import component;
import atom;
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

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::deletionMoveCBMC(RandomNumber& random, System& system,
                                                                            std::size_t selectedComponent,
                                                                            std::size_t selectedMolecule)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  MoveTypes move = MoveTypes::SwapCBMC;
  Component& component = system.components[selectedComponent];

  // Increment the count of swap deletion moves for the selected component
  component.mc_moves_statistics.addTrial(move, 1);

  // Proceed only if there is at least one molecule of the selected component
  if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] > 0)
  {
    // Get a reference to the molecule being deleted
    std::span<Atom> molecule = system.spanOfMolecule(selectedComponent, selectedMolecule);
    std::copy(system.electricField.begin(), system.electricField.end(), system.electricFieldNew.begin());
    // std::span<double3> electricFieldMoleculeNew = system.spanElectricFieldNew(selectedComponent, selectedMolecule);

    // Retrieve cutoff distances from the force field
    double cutOffFrameworkVDW = system.forceField.cutOffFrameworkVDW;
    double cutOffMoleculeVDW = system.forceField.cutOffMoleculeVDW;
    double cutOffCoulomb = system.forceField.cutOffCoulomb;

    // Retrace the molecule for the swap deletion using CBMC algorithm
    time_begin = std::chrono::system_clock::now();
    ChainData retraceData = CBMC::retraceMoleculeSwapDeletion(
        random, system.components[selectedComponent], system.hasExternalField, system.components, system.forceField,
        system.simulationBox, system.interpolationGrids, system.framework, system.spanOfFrameworkAtoms(),
        system.spanOfMoleculeAtoms(), system.beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb,
        selectedComponent, selectedMolecule, molecule, 1.0, system.numberOfTrialDirections);
    time_end = std::chrono::system_clock::now();

    // Update the CPU time statistics for the non-Ewald part of the move
    component.mc_moves_cputime[move]["NonEwald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["NonEwald"] += (time_end - time_begin);

    // Compute the energy difference in Fourier space due to the deletion
    time_begin = std::chrono::system_clock::now();
    RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
        system.simulationBox, {}, molecule);
    time_end = std::chrono::system_clock::now();
    // Update the CPU time statistics for the Ewald part of the move
    component.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);

    // Compute the tail energy difference due to the deletion
    time_begin = std::chrono::system_clock::now();
    [[maybe_unused]] RunningEnergy tailEnergyDifference =
        Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                                system.spanOfMoleculeAtoms(), {}, molecule) +
        Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                   system.spanOfFrameworkAtoms(), {}, molecule);
    time_end = std::chrono::system_clock::now();
    // Update the CPU time statistics for the tail corrections
    component.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);

    // Update the constructed count for the move statistics
    component.mc_moves_statistics.addConstructed(move, 1);

    RunningEnergy polarizationDifference;
    if (system.forceField.computePolarization)
    {
      Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                    system.spanOfFrameworkAtoms(), {},
                                                                    retraceData.electricField, {}, retraceData.atom);

      Interactions::computeEwaldFourierElectricFieldDifference(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                               system.fixedFrameworkStoredEik, system.storedEik,
                                                               system.totalEik, system.forceField, system.simulationBox,
                                                               {}, retraceData.electricField, {}, retraceData.atom);

      // Compute polarization energy difference
      polarizationDifference = Interactions::computePolarizationEnergyDifference(
          system.forceField, {}, retraceData.electricField, {}, retraceData.atom);
    }

    // Calculate the correction factor for Ewald summation
    double correctionFactorEwald =
        std::exp(-system.beta * (energyFourierDifference.potentialEnergy() + tailEnergyDifference.potentialEnergy() +
                                 polarizationDifference.potentialEnergy()));

    // Compute acceptance probability factors
    double fugacity = component.fugacityCoefficient.value_or(1.0) * system.pressure;
    double idealGasRosenbluthWeight = component.idealGasRosenbluthWeight.value_or(1.0);
    double preFactor = correctionFactorEwald * double(system.numberOfIntegerMoleculesPerComponent[selectedComponent]) /
                       (system.beta * component.molFraction * fugacity * system.simulationBox.volume);
    double Pacc = preFactor * idealGasRosenbluthWeight / retraceData.RosenbluthWeight;
    std::size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];
    double biasTransitionMatrix = system.tmmc.biasFactor(oldN - 1, oldN);

    // Check if the new macrostate is within the allowed TMMC range
    if (system.tmmc.doTMMC)
    {
      std::size_t newN = oldN - 1;
      if (newN < system.tmmc.minMacrostate)
      {
        return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
      }
    }

    // Apply acceptance/rejection rule
    if (random.uniform() < biasTransitionMatrix * Pacc)
    {
      component.mc_moves_statistics.addAccepted(move, 1);

      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
      system.deleteMolecule(selectedComponent, selectedMolecule, molecule);

      return {retraceData.energies - energyFourierDifference - tailEnergyDifference - polarizationDifference,
              double3(Pacc, 1.0 - Pacc, 0.0)};
    };
    return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
  }

  // Return default values if no molecules are available for deletion
  return {std::nullopt, double3(0.0, 1.0, 0.0)};
}
