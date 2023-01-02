module;

module mc_moves;

import component;
import atom;
import double3;
import double3x3;
import simd_quatd;
import simulationbox;
import cbmc;
import randomnumbers;
import system;
import energy_factor;
import energy_status;
import energy_status_inter;
import lambda;
import property_widom;
import averages;
import running_energy;
import forcefield;

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


std::optional<RunningEnergy> MC_Moves::insertionMove(System& system, size_t selectedComponent)
{
  size_t selectedMolecule = system.numberOfMoleculesPerComponent[selectedComponent];
  system.components[selectedComponent].statistics_SwapInsertionMove_CBMC.counts += 1;

  double cutOffVDW = system.forceField.cutOffVDW;
  double cutOffCoulomb = system.forceField.cutOffCoulomb;
  
  std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
  std::optional<ChainData> growData = system.growMoleculeSwapInsertion(cutOffVDW, cutOffCoulomb, selectedComponent, selectedMolecule, 1.0);
  if (!growData) return std::nullopt;

  std::span<const Atom> newMolecule = std::span(growData->atom.begin(), growData->atom.end());

  std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
  system.components[selectedComponent].cpuTime_SwapInsertionMove_CBMC_NonEwald += (t2 - t1);

  
  system.components[selectedComponent].statistics_SwapInsertionMove_CBMC.constructed += 1;

  std::chrono::system_clock::time_point u1 = std::chrono::system_clock::now();
  RunningEnergy energyFourierDifference = system.energyDifferenceEwaldFourier(system.storedEik, newMolecule, {});
  std::chrono::system_clock::time_point u2 = std::chrono::system_clock::now();
  system.components[selectedComponent].cpuTime_SwapInsertionMove_CBMC_Ewald += (u2 - u1);

  //RunningEnergy tailEnergyDifference = system.computeTailCorrectionVDWAddEnergy(selectedComponent) - 
  //                                     system.computeTailCorrectionVDWOldEnergy();
  RunningEnergy tailEnergyDifference;
  double correctionFactorEwald = std::exp(-system.simulationBox.Beta * (energyFourierDifference.total() + tailEnergyDifference.total()));


  double idealGasRosenbluthWeight = system.components[selectedComponent].idealGasRosenbluthWeight.value_or(1.0);
  double preFactor = correctionFactorEwald * system.simulationBox.Beta * system.components[selectedComponent].molFraction * 
                     system.simulationBox.pressure * system.simulationBox.volume /
                     double(1 + system.numberOfMoleculesPerComponent[selectedComponent]);

  if (RandomNumber::Uniform() < preFactor * growData->RosenbluthWeight / idealGasRosenbluthWeight)
  {
      system.components[selectedComponent].statistics_SwapInsertionMove_CBMC.accepted += 1;

      system.acceptEwaldMove();
      system.insertMolecule(selectedComponent, growData->atom);

      // Debug
      //assert(system.checkMoleculeIds());

      return growData->energies + energyFourierDifference + tailEnergyDifference;
  };
  
  return std::nullopt;
}

