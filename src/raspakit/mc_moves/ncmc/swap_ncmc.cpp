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

module mc_moves_noneq_cbmc;

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
import mc_moves_move_types;

std::pair<std::optional<RunningEnergy>, double3> NonEqCBMC(RandomNumber& random, System& system,
                                                           size_t selectedComponent)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  MoveTypes move = MoveTypes::SwapNonEqCBMC;
  Component& component = system.components[selectedComponent];
  size_t oldN = system.numberOfMoleculesPerComponent[selectedComponent];

  // Update move counts statistics for swap insertion move
  component.mc_moves_statistics.addTrial(move, 0);

  // all copied data: moleculePositions, moleculeAtomPositions, thermostat, dt
  // all const data: components, forcefield, simulationbox, numberofmoleculespercomponents, fixedFrameworkStoredEik
  // all scratch data: eik_x, eik_y, eik_z, eik_xy
  std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
  std::vector<Atom> moleculeAtomPositions(atomPositions.size());
  std::copy(atomPositions.begin(), atomPositions.end(), moleculeAtomPositions.begin());

  std::vector<Molecule> moleculePositions(system.moleculePositions);
  std::optional<Thermostat> thermostat(system.thermostat);

  // get Timestep from the max change
  double dt = system.mc_moves_statistics.getMaxChange(move);

  // insertion / deletion without acceptance
  if (random.uniform() < 0.5)  // Insertion
  {
    // Attempt to grow a new molecule using CBMC
    time_begin = std::chrono::system_clock::now();
    std::optional<ChainData> growData = CBMC::growMoleculeSwapInsertion(
        random, system.frameworkComponents, component, system.hasExternalField, system.components, system.forceField,
        system.simulationBox, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(), system.beta, growType,
        cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, selectedComponent, selectedMolecule, 1.0, 0uz,
        system.numberOfTrialDirections);
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

    // Calculate correction factor for Ewald energy difference
    double correctionFactorEwald =
        std::exp(-system.beta * (energyFourierDifference.potentialEnergy() + tailEnergyDifference.potentialEnergy()));

    // Compute the acceptance probability pre-factor
    double fugacity = component.fugacityCoefficient.value_or(1.0) * system.pressure;
    double idealGasRosenbluthWeight = component.idealGasRosenbluthWeight.value_or(1.0);
    double preFactor = correctionFactorEwald * system.beta * component.molFraction * fugacity *
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
  }
  else
  {
    if (oldN == 0)
    {
      return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
    }

    // Get a reference to the molecule being deleted
    std::span<Atom> molecule = system.spanOfMolecule(selectedComponent, selectedMolecule);

    // Retrieve cutoff distances from the force field
    double cutOffFrameworkVDW = system.forceField.cutOffFrameworkVDW;
    double cutOffMoleculeVDW = system.forceField.cutOffMoleculeVDW;
    double cutOffCoulomb = system.forceField.cutOffCoulomb;

    // Retrace the molecule for the swap deletion using CBMC algorithm
    time_begin = std::chrono::system_clock::now();
    ChainData retraceData = CBMC::retraceMoleculeSwapDeletion(
        random, system.frameworkComponents, system.components[selectedComponent], system.hasExternalField,
        system.components, system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(),
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

    // Calculate the correction factor for Ewald summation
    double correctionFactorEwald =
        std::exp(-system.beta * (energyFourierDifference.potentialEnergy() + tailEnergyDifference.potentialEnergy()));

    // Compute acceptance probability factors
    double fugacity = component.fugacityCoefficient.value_or(1.0) * system.pressure;
    double idealGasRosenbluthWeight = component.idealGasRosenbluthWeight.value_or(1.0);
    double preFactor = correctionFactorEwald * double(system.numberOfIntegerMoleculesPerComponent[selectedComponent]) /
                       (system.beta * component.molFraction * fugacity * system.simulationBox.volume);
    double Pacc = preFactor * idealGasRosenbluthWeight / retraceData.RosenbluthWeight;
    size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];
    double biasTransitionMatrix = system.tmmc.biasFactor(oldN - 1, oldN);

    // Check if the new macrostate is within the allowed TMMC range
    if (system.tmmc.doTMMC)
    {
      size_t newN = oldN - 1;
      if (newN < system.tmmc.minMacrostate)
      {
        return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
      }
    }
  }

}