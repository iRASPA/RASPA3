module;

module mc_moves_deletion;

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
import move_statistics;
import mc_moves_probabilities_particles;
import transition_matrix;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;


std::pair<std::optional<RunningEnergy>, double3> 
MC_Moves::deletionMove(RandomNumber &random, System& system, size_t selectedComponent, size_t selectedMolecule)
{
  system.components[selectedComponent].mc_moves_statistics.swapDeletionMove_CBMC.counts += 1;
  system.components[selectedComponent].mc_moves_statistics.swapDeletionMove_CBMC.totalCounts += 1;
  
  if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] > 0)
  {
    system.components[selectedComponent].mc_moves_statistics.swapDeletionMove_CBMC.constructed += 1;
    system.components[selectedComponent].mc_moves_statistics.swapDeletionMove_CBMC.totalConstructed += 1;

    std::span<Atom> molecule = system.spanOfMolecule(selectedComponent, selectedMolecule);

    double cutOffVDW = system.forceField.cutOffVDW;
    double cutOffCoulomb = system.forceField.cutOffCoulomb;

    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
    ChainData retraceData = 
      CBMC::retraceMoleculeSwapDeletion(random, system.components, system.forceField, system.simulationBox, 
                                        system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(), system.beta, 
                                        cutOffVDW, cutOffCoulomb, selectedComponent, selectedMolecule, molecule, 
                                        1.0, 0.0, system.numberOfTrialDirections);
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
    system.components[selectedComponent].mc_moves_cputime.swapDeletionMoveCBMCNonEwald += (t2 - t1);
    system.mc_moves_cputime.swapDeletionMoveCBMCNonEwald += (t2 - t1);

    std::chrono::system_clock::time_point u1 = std::chrono::system_clock::now();
    RunningEnergy energyFourierDifference = //system.energyDifferenceEwaldFourier(system.storedEik, {}, molecule);
      Interactions::energyDifferenceEwaldFourier(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                               system.storedEik, system.totalEik,
                                               system.forceField, system.simulationBox,
                                               {}, molecule);
    std::chrono::system_clock::time_point u2 = std::chrono::system_clock::now();
    system.components[selectedComponent].mc_moves_cputime.swapDeletionMoveCBMCEwald += (u2 - u1);
    system.mc_moves_cputime.swapDeletionMoveCBMCEwald += (u2 - u1);


    std::chrono::system_clock::time_point v1 = std::chrono::system_clock::now();
    [[maybe_unused]] RunningEnergy tailEnergyDifference = 
      Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,                           
                                         system.spanOfMoleculeAtoms(), {}, molecule) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,            
                                         system.spanOfFrameworkAtoms(), {}, molecule);
    std::chrono::system_clock::time_point v2 = std::chrono::system_clock::now();
    system.components[selectedComponent].mc_moves_cputime.swapDeletionMoveCBMCTail += (v2 - v1);
    system.mc_moves_cputime.swapDeletionMoveCBMCTail += (v2 - v1);

    double correctionFactorEwald = std::exp(-system.beta * (energyFourierDifference.total() + 
                                                            tailEnergyDifference.total()));

    double idealGasRosenbluthWeight = system.components[selectedComponent].idealGasRosenbluthWeight.value_or(1.0);
    double preFactor = correctionFactorEwald * double(system.numberOfMoleculesPerComponent[selectedComponent]) /
                       (system.beta * system.components[selectedComponent].molFraction * 
                        system.pressure * system.simulationBox.volume);
    double Pacc = preFactor * idealGasRosenbluthWeight / retraceData.RosenbluthWeight;
    size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];
    double biasTransitionMatrix = system.tmmc.biasFactor(oldN - 1, oldN);

    if(system.tmmc.doTMMC)
    {
      size_t newN = oldN - 1;
      if(newN < system.tmmc.minMacrostate)
      {
        return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
      }
    }

    if (random.uniform() < biasTransitionMatrix * Pacc)
    {
      system.components[selectedComponent].mc_moves_statistics.swapDeletionMove_CBMC.accepted += 1;
      system.components[selectedComponent].mc_moves_statistics.swapDeletionMove_CBMC.totalAccepted += 1;

      //system.acceptEwaldMove();
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
      system.deleteMolecule(selectedComponent, selectedMolecule, molecule);

      return {retraceData.energies - energyFourierDifference - tailEnergyDifference, double3(Pacc, 1.0 - Pacc, 0.0)};
    };
    return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
  }

  return {std::nullopt, double3(0.0, 1.0, 0.0)};
}
