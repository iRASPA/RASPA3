module;

module mc_moves_widom;

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
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;


std::optional<double> 
MC_Moves::WidomMove(RandomNumber &random, System& system, size_t selectedComponent)
{
  size_t selectedMolecule = system.numberOfMoleculesPerComponent[selectedComponent];
  system.components[selectedComponent].mc_moves_statistics.WidomMove_CBMC.counts += 1;

  double cutOffVDW = system.forceField.cutOffVDW;
  double cutOffCoulomb = system.forceField.cutOffCoulomb;
  Component::GrowType growType = system.components[selectedComponent].growType;
  
  std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
  std::vector<Atom> atoms = 
    system.components[selectedComponent].recenteredCopy(1.0, system.numberOfMoleculesPerComponent[selectedComponent]);
  std::optional<ChainData> growData = 
    CBMC::growMoleculeSwapInsertion(random, system.hasExternalField, system.components, system.forceField, system.simulationBox, 
                                    system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(), system.beta, 
                                    growType, cutOffVDW, cutOffCoulomb, selectedComponent, selectedMolecule, 1.0, 
                                    atoms, system.numberOfTrialDirections);
  std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
  system.components[selectedComponent].mc_moves_cputime.WidomMoveCBMCNonEwald += (t2 - t1);
  system.mc_moves_cputime.WidomMoveCBMCNonEwald += (t2 - t1);

  if (!growData) return std::nullopt;

  [[maybe_unused]] std::span<const Atom> newMolecule = std::span(growData->atom.begin(), growData->atom.end());
  

  system.components[selectedComponent].mc_moves_statistics.WidomMove_CBMC.constructed += 1;

  std::chrono::system_clock::time_point u1 = std::chrono::system_clock::now();
  RunningEnergy energyFourierDifference = 
    Interactions::energyDifferenceEwaldFourier(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                               system.storedEik, system.totalEik,
                                               system.forceField, system.simulationBox,
                                               newMolecule, {});
  std::chrono::system_clock::time_point u2 = std::chrono::system_clock::now();
  system.components[selectedComponent].mc_moves_cputime.WidomMoveCBMCEwald += (u2 - u1);
  system.mc_moves_cputime.WidomMoveCBMCEwald += (u2 - u1);

  RunningEnergy tailEnergyDifference =
    Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                         system.spanOfMoleculeAtoms(),newMolecule, {}) +
    Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,            
                                         system.spanOfFrameworkAtoms(), newMolecule, {});

  double correctionFactorEwald = std::exp(-system.beta * (energyFourierDifference.total() + 
                                                          tailEnergyDifference.total()));

  double idealGasRosenbluthWeight = system.components[selectedComponent].idealGasRosenbluthWeight.value_or(1.0);

  return  correctionFactorEwald * growData->RosenbluthWeight / idealGasRosenbluthWeight;
}

