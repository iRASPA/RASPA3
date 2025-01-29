module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <cmath>
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
import mc_moves_statistics;
import mc_moves_move_types;
import mc_moves_probabilities;

std::optional<std::pair<RunningEnergy, RunningEnergy>> MC_Moves::GibbsSwapMove_CBMC(RandomNumber& random,
                                                                                    System& systemA, System& systemB,
                                                                                    size_t selectedComponent)
{
  // Return if there are no molecules of the selected component in system B
  if (systemB.numberOfIntegerMoleculesPerComponent[selectedComponent] == 0) return std::nullopt;

  // Index for the new molecule in system A
  size_t newMoleculeIndex = systemA.numberOfMoleculesPerComponent[selectedComponent];

  // Update move counts statistics for both systems
  systemA.components[selectedComponent].mc_moves_statistics.addTrial(MoveTypes::GibbsSwapCBMC);
  systemB.components[selectedComponent].mc_moves_statistics.addTrial(MoveTypes::GibbsSwapCBMC);

  // Retrieve cutoff distances and grow type from system A
  double cutOffFrameworkVDW = systemA.forceField.cutOffFrameworkVDW;
  double cutOffMoleculeVDW = systemA.forceField.cutOffMoleculeVDW;
  double cutOffCoulomb = systemA.forceField.cutOffCoulomb;
  Component::GrowType growType = systemA.components[selectedComponent].growType;

  // std::vector<Atom> atoms = systemA.components[selectedComponent].recenteredCopy(1.0,
  //                                                         systemA.numberOfMoleculesPerComponent[selectedComponent]);
  std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();  // Start timing CBMC insertion

  // Attempt to grow a new molecule in system A using CBMC insertion
  std::optional<ChainData> growData = CBMC::growMoleculeSwapInsertion(
      random, systemA.frameworkComponents, systemA.components[selectedComponent], systemA.hasExternalField,
      systemA.components, systemA.forceField, systemA.simulationBox, systemA.spanOfFrameworkAtoms(),
      systemA.spanOfMoleculeAtoms(), systemA.beta, growType, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb,
      selectedComponent, newMoleculeIndex, 1.0, 0uz, systemA.numberOfTrialDirections);

  std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();  // End timing CBMC insertion

  // Update CPU time statistics for CBMC insertion (non-Ewald part)
  systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapMoveCBMCNonEwald += (t2 - t1);
  systemA.mc_moves_cputime.GibbsSwapMoveCBMCNonEwald += (t2 - t1);

  if (!growData) return std::nullopt;  // Insertion failed, return

  // Get new molecule atoms
  std::span<const Atom> newMolecule = std::span(growData->atom.begin(), growData->atom.end());

  // Update statistics for successfully constructed molecules in system A
  systemA.components[selectedComponent].mc_moves_statistics.addConstructed(MoveTypes::GibbsSwapCBMC);

  std::chrono::system_clock::time_point u1 =
      std::chrono::system_clock::now();  // Start timing Ewald Fourier computation

  // Compute Ewald Fourier energy difference for system A
  RunningEnergy energyFourierDifferenceA = Interactions::energyDifferenceEwaldFourier(
      systemA.eik_x, systemA.eik_y, systemA.eik_z, systemA.eik_xy, systemA.storedEik, systemA.totalEik,
      systemA.forceField, systemA.simulationBox, newMolecule, {});

  std::chrono::system_clock::time_point u2 = std::chrono::system_clock::now();  // End timing Ewald Fourier computation

  // Update CPU time statistics for Ewald Fourier computation
  systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapMoveCBMCEwald += (u2 - u1);
  systemA.mc_moves_cputime.GibbsSwapMoveCBMCEwald += (u2 - u1);

  std::chrono::system_clock::time_point v1 = std::chrono::system_clock::now();  // Start timing tail energy computation

  // Compute tail energy difference for system A
  [[maybe_unused]] RunningEnergy tailEnergyDifferenceA =
      Interactions::computeInterMolecularTailEnergyDifference(systemA.forceField, systemA.simulationBox,
                                                              systemA.spanOfMoleculeAtoms(), newMolecule, {}) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(systemA.forceField, systemA.simulationBox,
                                                                 systemA.spanOfFrameworkAtoms(), newMolecule, {});

  std::chrono::system_clock::time_point v2 = std::chrono::system_clock::now();  // End timing tail energy computation

  // Update CPU time statistics for tail energy computation
  systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapMoveCBMCTail += (v2 - v1);
  systemA.mc_moves_cputime.GibbsSwapMoveCBMCTail += (v2 - v1);

  // Compute correction factor for Ewald energies in system A
  double correctionFactorEwaldA =
      std::exp(-systemA.beta * (energyFourierDifferenceA.potentialEnergy() + tailEnergyDifferenceA.potentialEnergy()));

  // Select a random molecule of the selected component from system B
  size_t selectedMolecule = systemB.randomMoleculeOfComponent(random, selectedComponent);
  std::span<Atom> molecule = systemB.spanOfMolecule(selectedComponent, selectedMolecule);

  std::chrono::system_clock::time_point tB1 = std::chrono::system_clock::now();  // Start timing CBMC deletion

  // Retrace the selected molecule in system B for deletion using CBMC
  ChainData retraceData = CBMC::retraceMoleculeSwapDeletion(
      random, systemB.frameworkComponents, systemB.components[selectedComponent], systemB.hasExternalField,
      systemB.components, systemB.forceField, systemB.simulationBox, systemB.spanOfFrameworkAtoms(),
      systemB.spanOfMoleculeAtoms(), systemB.beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb,
      selectedComponent, selectedMolecule, molecule, 1.0, systemB.numberOfTrialDirections);

  std::chrono::system_clock::time_point tB2 = std::chrono::system_clock::now();  // End timing CBMC deletion

  // Update CPU time statistics for CBMC deletion (non-Ewald part)
  systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapMoveCBMCNonEwald += (tB2 - tB1);
  systemA.mc_moves_cputime.GibbsSwapMoveCBMCNonEwald += (tB2 - tB1);

  std::chrono::system_clock::time_point uB1 =
      std::chrono::system_clock::now();  // Start timing Ewald Fourier computation for system B

  // Compute Ewald Fourier energy difference for system B
  RunningEnergy energyFourierDifferenceB = Interactions::energyDifferenceEwaldFourier(
      systemB.eik_x, systemB.eik_y, systemB.eik_z, systemB.eik_xy, systemB.storedEik, systemB.totalEik,
      systemB.forceField, systemB.simulationBox, {}, molecule);

  std::chrono::system_clock::time_point uB2 =
      std::chrono::system_clock::now();  // End timing Ewald Fourier computation for system B

  // Update CPU time statistics for Ewald Fourier computation
  systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapMoveCBMCEwald += (uB2 - uB1);
  systemA.mc_moves_cputime.GibbsSwapMoveCBMCEwald += (uB2 - uB1);

  std::chrono::system_clock::time_point w1 =
      std::chrono::system_clock::now();  // Start timing tail energy computation for system B

  // Compute tail energy difference for system B
  [[maybe_unused]] RunningEnergy tailEnergyDifferenceB =
      Interactions::computeInterMolecularTailEnergyDifference(systemB.forceField, systemB.simulationBox,
                                                              systemB.spanOfMoleculeAtoms(), {}, molecule) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(systemB.forceField, systemB.simulationBox,
                                                                 systemB.spanOfFrameworkAtoms(), {}, molecule);

  std::chrono::system_clock::time_point w2 =
      std::chrono::system_clock::now();  // End timing tail energy computation for system B

  // Update CPU time statistics for tail energy computation
  systemA.components[selectedComponent].mc_moves_cputime.GibbsSwapMoveCBMCTail += (w2 - w1);
  systemA.mc_moves_cputime.GibbsSwapMoveCBMCTail += (w2 - w1);

  // Update statistics for retraced molecules in system B
  systemB.components[selectedComponent].mc_moves_statistics.addConstructed(MoveTypes::GibbsSwapCBMC);

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
    systemA.components[selectedComponent].mc_moves_statistics.addAccepted(MoveTypes::GibbsSwapCBMC);

    // Accept Ewald updates and insert the new molecule into system A
    Interactions::acceptEwaldMove(systemA.forceField, systemA.storedEik, systemA.totalEik);
    systemA.insertMolecule(selectedComponent, growData->molecule, growData->atom);

    // Update accepted move statistics for system B
    systemB.components[selectedComponent].mc_moves_statistics.addAccepted(MoveTypes::GibbsSwapCBMC);

    // Accept Ewald updates and delete the selected molecule from system B
    Interactions::acceptEwaldMove(systemB.forceField, systemB.storedEik, systemB.totalEik);
    systemB.deleteMolecule(selectedComponent, selectedMolecule, molecule);

    // Return the energy differences for both systems
    return std::make_pair(growData->energies + energyFourierDifferenceA + tailEnergyDifferenceA,
                          -(retraceData.energies - energyFourierDifferenceB - tailEnergyDifferenceB));
  };

  return std::nullopt;
}
