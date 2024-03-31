module;

module mc_moves_reinsertion;

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

import component;
import atom;
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


std::optional<RunningEnergy> 
MC_Moves::reinsertionMove(RandomNumber &random, System& system, size_t selectedComponent, size_t selectedMolecule, 
                          std::span<Atom> molecule)
{
  std::chrono::system_clock::time_point time_begin, time_end;

  system.components[selectedComponent].mc_moves_statistics.reinsertionMove_CBMC.counts += 1;
  system.components[selectedComponent].mc_moves_statistics.reinsertionMove_CBMC.totalCounts += 1;

  if (system.numberOfMoleculesPerComponent[selectedComponent] == 0) 
  {
    return std::nullopt;
  }

  double cutOffVDW = system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffVDW;
  double cutOffCoulomb = system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffCoulomb;

  time_begin = std::chrono::system_clock::now();
  std::optional<ChainData> growData = 
    CBMC::growMoleculeReinsertion(random, system.hasExternalField, system.components, system.forceField, system.simulationBox, 
                                  system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(), system.beta,
                                  cutOffVDW, cutOffCoulomb, selectedComponent, selectedMolecule, molecule, 
                                  system.numberOfTrialDirections);
  time_end = std::chrono::system_clock::now();
  system.components[selectedComponent].mc_moves_cputime.reinsertionMoveCBMCNonEwald += (time_end - time_begin);
  system.mc_moves_cputime.reinsertionMoveCBMCNonEwald += (time_end - time_begin);

  if (!growData) return std::nullopt;

  std::span<const Atom> newMolecule = std::span(growData->atom.begin(), growData->atom.end());

  system.components[selectedComponent].mc_moves_statistics.reinsertionMove_CBMC.constructed += 1;
  system.components[selectedComponent].mc_moves_statistics.reinsertionMove_CBMC.totalConstructed += 1;

  time_begin = std::chrono::system_clock::now();
  ChainData retraceData = 
    CBMC::retraceMoleculeReinsertion(random, system.hasExternalField, system.components, system.forceField, system.simulationBox, 
                                     system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(), system.beta,
                                     cutOffVDW, cutOffCoulomb, selectedComponent, selectedMolecule, molecule, 
                                     growData->storedR, system.numberOfTrialDirections);
  time_end = std::chrono::system_clock::now();
  system.components[selectedComponent].mc_moves_cputime.reinsertionMoveCBMCNonEwald += (time_end - time_begin);
  system.mc_moves_cputime.reinsertionMoveCBMCNonEwald += (time_end - time_begin);

  time_begin = std::chrono::system_clock::now();
  RunningEnergy energyFourierDifference = 
    Interactions::energyDifferenceEwaldFourier(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                               system.storedEik, system.totalEik,
                                               system.forceField, system.simulationBox,
                                               newMolecule, molecule);
  time_end = std::chrono::system_clock::now();
  system.components[selectedComponent].mc_moves_cputime.reinsertionMoveCBMCEwald += (time_end - time_begin);
  system.mc_moves_cputime.reinsertionMoveCBMCEwald += (time_end - time_begin);

  double correctionFactorDualCutOff = 1.0;
  std::optional<RunningEnergy> energyNew;
  std::optional<RunningEnergy> energyOld;
  if(system.forceField.useDualCutOff)
  {
    energyNew = CBMC::computeExternalNonOverlappingEnergyDualCutOff(system.hasExternalField, system.forceField, system.simulationBox, 
                               system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(), 
                               system.forceField.cutOffVDW, system.forceField.cutOffCoulomb, growData->atom);
    energyOld = CBMC::computeExternalNonOverlappingEnergyDualCutOff(system.hasExternalField, system.forceField, system.simulationBox, 
                               system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
                               system.forceField.cutOffVDW, system.forceField.cutOffCoulomb, retraceData.atom);
    correctionFactorDualCutOff = std::exp(-system.beta * (energyNew->total() - growData->energies.total() 
                                                       - (energyOld->total() - retraceData.energies.total())) );
  }

  double correctionFactorFourier = std::exp(-system.beta * energyFourierDifference.total());

  // apply acceptance/rejection rule
  if (random.uniform() < correctionFactorDualCutOff * correctionFactorFourier * 
                         growData->RosenbluthWeight / retraceData.RosenbluthWeight)
  {
    system.components[selectedComponent].mc_moves_statistics.reinsertionMove_CBMC.accepted += 1;
    system.components[selectedComponent].mc_moves_statistics.reinsertionMove_CBMC.totalAccepted += 1;

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
    std::copy(newMolecule.begin(), newMolecule.end(), molecule.begin());

    if(system.forceField.useDualCutOff)
    {
      return (energyNew.value() - energyOld.value()) + energyFourierDifference;
    }

    return (growData->energies - retraceData.energies) + energyFourierDifference;
  };

  return std::nullopt;
}
