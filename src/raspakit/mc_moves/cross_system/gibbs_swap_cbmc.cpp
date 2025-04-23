module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#endif

module mc_moves_gibbs_swap_cbmc;

#ifndef USE_LEGACY_HEADERS
import <optional>;
import <span>;
import <chrono>;
import <vector>;
import <cmath>;
import <tuple>;
#endif

import randomnumbers;
import running_energy;
import system;
import atom;
import cbmc;
import cbmc_chain_data;
import energy_factor;
import energy_status;
import energy_status_inter;
import property_lambda_probability_histogram;
import property_widom;
import averages;
import running_energy;
import forcefield;
import component;
import simulationbox;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import mc_moves_move_types;

std::optional<std::pair<RunningEnergy, RunningEnergy>> MC_Moves::GibbsSwapMove_CBMC(RandomNumber& random,
                                                                                    System& systemA, System& systemB,
                                                                                    size_t selectedComponent)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  MoveTypes move = MoveTypes::GibbsSwapCBMC;
  Component& componentA = systemA.components[selectedComponent];
  Component& componentB = systemB.components[selectedComponent];

  // Return if there are no molecules of the selected component in system B
  if (systemB.numberOfIntegerMoleculesPerComponent[selectedComponent] == 0) return std::nullopt;

  // Index for the new molecule in system A
  size_t newMoleculeIndex = systemA.numberOfMoleculesPerComponent[selectedComponent];

  // Update move counts statistics for both systems
  componentA.mc_moves_statistics.addTrial(move);
  componentB.mc_moves_statistics.addTrial(move);

  // Retrieve cutoff distances and grow type from system A
  double cutOffFrameworkVDW = systemA.forceField.cutOffFrameworkVDW;
  double cutOffMoleculeVDW = systemA.forceField.cutOffMoleculeVDW;
  double cutOffCoulomb = systemA.forceField.cutOffCoulomb;
  Component::GrowType growType = componentA.growType;

  // std::vector<Atom> atoms = componentA.recenteredCopy(1.0,
  // systemA.numberOfMoleculesPerComponent[selectedComponent]);

  // Attempt to grow a new molecule in system A using CBMC insertion
  time_begin = std::chrono::system_clock::now();
  std::optional<ChainData> growData = CBMC::growMoleculeSwapInsertion(
      random, componentA, systemA.hasExternalField, systemA.components, systemA.forceField, systemA.simulationBox,
      systemA.interpolationGrids, systemA.framework, systemA.spanOfFrameworkAtoms(), systemA.spanOfMoleculeAtoms(),
      systemA.beta, growType, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, selectedComponent, newMoleculeIndex,
      1.0, 0uz, systemA.numberOfTrialDirections);
  time_end = std::chrono::system_clock::now();

  // Update CPU time statistics for CBMC insertion (non-Ewald part)
  componentA.mc_moves_cputime[move]["NonEwald"] += (time_end - time_begin);
  systemA.mc_moves_cputime[move]["NonEwald"] += (time_end - time_begin);

  if (!growData) return std::nullopt;  // Insertion failed, return

  // Get new molecule atoms
  std::span<const Atom> newMolecule = std::span(growData->atom.begin(), growData->atom.end());

  // Update statistics for successfully constructed molecules in system A
  componentA.mc_moves_statistics.addConstructed(move);

  // Compute Ewald Fourier energy difference for system A
  time_begin = std::chrono::system_clock::now();
  RunningEnergy energyFourierDifferenceA = Interactions::energyDifferenceEwaldFourier(
      systemA.eik_x, systemA.eik_y, systemA.eik_z, systemA.eik_xy, systemA.storedEik, systemA.totalEik,
      systemA.forceField, systemA.simulationBox, newMolecule, {});
  time_end = std::chrono::system_clock::now();

  // Update CPU time statistics for Ewald Fourier computation
  componentA.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);
  systemA.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);

  // Compute tail energy difference for system A
  time_begin = std::chrono::system_clock::now();
  [[maybe_unused]] RunningEnergy tailEnergyDifferenceA =
      Interactions::computeInterMolecularTailEnergyDifference(systemA.forceField, systemA.simulationBox,
                                                              systemA.spanOfMoleculeAtoms(), newMolecule, {}) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(systemA.forceField, systemA.simulationBox,
                                                                 systemA.spanOfFrameworkAtoms(), newMolecule, {});
  time_end = std::chrono::system_clock::now();

  // Update CPU time statistics for tail energy computation
  componentA.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);
  systemA.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);

  // Compute correction factor for Ewald energies in system A
  double correctionFactorEwaldA =
      std::exp(-systemA.beta * (energyFourierDifferenceA.potentialEnergy() + tailEnergyDifferenceA.potentialEnergy()));

  // Select a random molecule of the selected component from system B
  size_t selectedMolecule = systemB.randomMoleculeOfComponent(random, selectedComponent);
  std::span<Atom> molecule = systemB.spanOfMolecule(selectedComponent, selectedMolecule);

  // Retrace the selected molecule in system B for deletion using CBMC
  time_begin = std::chrono::system_clock::now();
  ChainData retraceData = CBMC::retraceMoleculeSwapDeletion(
      random, componentB, systemB.hasExternalField, systemB.components, systemB.forceField, systemB.simulationBox,
      systemB.interpolationGrids, systemB.framework, systemB.spanOfFrameworkAtoms(), systemB.spanOfMoleculeAtoms(),
      systemB.beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, selectedComponent, selectedMolecule, molecule,
      1.0, systemB.numberOfTrialDirections);
  time_end = std::chrono::system_clock::now();

  // Update CPU time statistics for CBMC deletion (non-Ewald part)
  componentA.mc_moves_cputime[move]["NonEwald"] += (time_end - time_begin);
  systemA.mc_moves_cputime[move]["NonEwald"] += (time_end - time_begin);

  // Compute Ewald Fourier energy difference for system B
  time_begin = std::chrono::system_clock::now();
  RunningEnergy energyFourierDifferenceB = Interactions::energyDifferenceEwaldFourier(
      systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy, systemB.storedEik, systemB.totalEik,
      systemB.forceField, systemB.simulationBox, {}, molecule);
  time_end = std::chrono::system_clock::now();

  // Update CPU time statistics for Ewald Fourier computation
  componentA.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);
  systemA.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);

  // Compute tail energy difference for system B
  time_begin = std::chrono::system_clock::now();
  [[maybe_unused]] RunningEnergy tailEnergyDifferenceB =
      Interactions::computeInterMolecularTailEnergyDifference(systemB.forceField, systemB.simulationBox,
                                                              systemB.spanOfMoleculeAtoms(), {}, molecule) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(systemB.forceField, systemB.simulationBox,
                                                                 systemB.spanOfFrameworkAtoms(), {}, molecule);
  time_end = std::chrono::system_clock::now();

  // Update CPU time statistics for tail energy computation
  componentA.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);
  systemA.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);

  // Update statistics for retraced molecules in system B
  componentB.mc_moves_statistics.addConstructed(move);

  // Compute correction factor for Ewald energies in system B
  double correctionFactorEwaldB =
      std::exp(systemB.beta * (energyFourierDifferenceB.potentialEnergy() + tailEnergyDifferenceB.potentialEnergy()));

  // Apply Metropolis acceptance criterion
  if (random.uniform() <
      (correctionFactorEwaldA * growData->RosenbluthWeight *
       static_cast<double>(systemB.numberOfIntegerMoleculesPerComponent[selectedComponent]) *
       systemA.simulationBox.volume) /
          (correctionFactorEwaldB * retraceData.RosenbluthWeight *
           (1.0 + static_cast<double>(systemA.numberOfIntegerMoleculesPerComponent[selectedComponent])) *
           systemB.simulationBox.volume))
  {
    // Update accepted move statistics for system A
    componentA.mc_moves_statistics.addAccepted(move);

    // Accept Ewald updates and insert the new molecule into system A
    Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.totalEik);
    systemA.insertMolecule(selectedComponent, growData->molecule, growData->atom);

    // Update accepted move statistics for system B
    componentB.mc_moves_statistics.addAccepted(move);

    // Accept Ewald updates and delete the selected molecule from system B
    Interactions::acceptEwaldMove(systemB.forceField, systemB.storedEik, systemB.totalEik);
    systemB.deleteMolecule(selectedComponent, selectedMolecule, molecule);

    // Return the energy differences for both systems
    return std::make_pair(growData->energies + energyFourierDifferenceA + tailEnergyDifferenceA,
                          -(retraceData.energies - energyFourierDifferenceB - tailEnergyDifferenceB));
  };

  return std::nullopt;
}
