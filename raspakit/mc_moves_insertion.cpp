module;

module mc_moves_insertion;

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
import move_statistics;
import mc_moves_probabilities_particles;
import transition_matrix;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;

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


std::pair<std::optional<RunningEnergy>, double3> 
MC_Moves::insertionMove(RandomNumber &random, System& system, size_t selectedComponent)
{
  size_t selectedMolecule = system.numberOfMoleculesPerComponent[selectedComponent];
  system.components[selectedComponent].mc_moves_statistics.swapInsertionMove_CBMC.counts += 1;
  system.components[selectedComponent].mc_moves_statistics.swapInsertionMove_CBMC.totalCounts += 1;

  double cutOffVDW = system.forceField.cutOffVDW;
  double cutOffCoulomb = system.forceField.cutOffCoulomb;
  Component::GrowType growType = system.components[selectedComponent].growType;
  
  std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
  std::vector<Atom> atoms = 
    system.components[selectedComponent].recenteredCopy(1.0, system.numberOfMoleculesPerComponent[selectedComponent]);
  std::optional<ChainData> growData = 
    CBMC::growMoleculeSwapInsertion(random, system.components, system.forceField, system.simulationBox, 
                                    system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(), system.beta,
                                    growType, cutOffVDW, cutOffCoulomb, selectedComponent, selectedMolecule, 1.0, 
                                    atoms, system.numberOfTrialDirections);
  std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
  system.components[selectedComponent].mc_moves_cputime.swapInsertionMoveCBMCNonEwald += (t2 - t1);
  system.mc_moves_cputime.swapInsertionMoveCBMCNonEwald += (t2 - t1);
  if (!growData) return {std::nullopt, double3(0.0, 1.0, 0.0)};

  std::span<const Atom> newMolecule = std::span(growData->atom.begin(), growData->atom.end());


  
  system.components[selectedComponent].mc_moves_statistics.swapInsertionMove_CBMC.constructed += 1;
  system.components[selectedComponent].mc_moves_statistics.swapInsertionMove_CBMC.totalConstructed += 1;

  std::chrono::system_clock::time_point u1 = std::chrono::system_clock::now();
  RunningEnergy energyFourierDifference = 
    Interactions::energyDifferenceEwaldFourier(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                               system.storedEik, system.totalEik,
                                               system.forceField, system.simulationBox,
                                               newMolecule, {});
  std::chrono::system_clock::time_point u2 = std::chrono::system_clock::now();
  system.components[selectedComponent].mc_moves_cputime.swapInsertionMoveCBMCEwald += (u2 - u1);
  system.mc_moves_cputime.swapInsertionMoveCBMCEwald += (u2 - u1);

  std::chrono::system_clock::time_point v1 = std::chrono::system_clock::now();
  RunningEnergy tailEnergyDifference = 
    Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,                 
                                         system.spanOfMoleculeAtoms(), newMolecule, {}) +
    Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,            
                                         system.spanOfFrameworkAtoms(), newMolecule, {});
  std::chrono::system_clock::time_point v2 = std::chrono::system_clock::now();
  system.components[selectedComponent].mc_moves_cputime.swapInsertionMoveCBMCTail += (v2 - v1);
  system.mc_moves_cputime.swapInsertionMoveCBMCTail += (v2 - v1);

  double correctionFactorEwald = 
    std::exp(-system.beta * (energyFourierDifference.total() + tailEnergyDifference.total()));

  double idealGasRosenbluthWeight = system.components[selectedComponent].idealGasRosenbluthWeight.value_or(1.0);
  double preFactor = correctionFactorEwald * system.beta * system.components[selectedComponent].molFraction * 
                     system.pressure * system.simulationBox.volume /
                     double(1 + system.numberOfMoleculesPerComponent[selectedComponent]);
  double Pacc = preFactor * growData->RosenbluthWeight / idealGasRosenbluthWeight;
  size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];
  double biasTransitionMatrix = system.tmmc.biasFactor(oldN + 1, oldN);

  if(system.tmmc.doTMMC)
  {
    size_t newN = oldN + 1;
    if(newN > system.tmmc.maxMacrostate)
    {
      return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
    }
  }

  if (random.uniform() < biasTransitionMatrix * Pacc)
  {
    system.components[selectedComponent].mc_moves_statistics.swapInsertionMove_CBMC.accepted += 1;
    system.components[selectedComponent].mc_moves_statistics.swapInsertionMove_CBMC.totalAccepted += 1;

    //system.acceptEwaldMove();
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
    system.insertMolecule(selectedComponent, growData->atom);

    return {growData->energies + energyFourierDifference + tailEnergyDifference, double3(0.0, 1.0 - Pacc, Pacc)};
  };
  
  return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
}

