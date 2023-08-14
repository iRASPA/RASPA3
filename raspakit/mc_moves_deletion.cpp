module;

module mc_moves;

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

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::deletionMove(System& system, size_t selectedComponent, size_t selectedMolecule)
{
  system.components[selectedComponent].mc_moves_probabilities.statistics_SwapDeletionMove_CBMC.counts += 1;
  
  if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] > 0)
  {
    system.components[selectedComponent].mc_moves_probabilities.statistics_SwapDeletionMove_CBMC.constructed += 1;

    std::span<Atom> molecule = system.spanOfMolecule(selectedComponent, selectedMolecule);

    double cutOffVDW = system.forceField.cutOffVDW;
    double cutOffCoulomb = system.forceField.cutOffCoulomb;

    std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
    ChainData retraceData = system.retraceMoleculeSwapDeletion(cutOffVDW, cutOffCoulomb, selectedComponent, selectedMolecule, molecule, 1.0, 0.0);
    std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
    system.components[selectedComponent].mc_moves_cputime.swapDeletionMoveCBMCNonEwald += (t2 - t1);
    system.mc_moves_cputime.swapDeletionMoveCBMCNonEwald += (t2 - t1);

    std::chrono::system_clock::time_point u1 = std::chrono::system_clock::now();
    RunningEnergy energyFourierDifference = system.energyDifferenceEwaldFourier(system.storedEik, {}, molecule);
    std::chrono::system_clock::time_point u2 = std::chrono::system_clock::now();
    system.components[selectedComponent].mc_moves_cputime.swapDeletionMoveCBMCEwald += (u2 - u1);
    system.mc_moves_cputime.swapDeletionMoveCBMCEwald += (u2 - u1);


    std::chrono::system_clock::time_point v1 = std::chrono::system_clock::now();
    [[maybe_unused]] RunningEnergy tailEnergyDifference = system.computeInterMolecularTailEnergyDifference({}, molecule) +
                                                          system.computeFrameworkMoleculeTailEnergyDifference({}, molecule);
    std::chrono::system_clock::time_point v2 = std::chrono::system_clock::now();
    system.components[selectedComponent].mc_moves_cputime.swapDeletionMoveCBMCTail += (v2 - v1);
    system.mc_moves_cputime.swapDeletionMoveCBMCTail += (v2 - v1);

    double correctionFactorEwald = std::exp(-system.beta * (energyFourierDifference.total() + tailEnergyDifference.total()));
    //double correctionFactorEwald = std::exp(-system.beta * (energyFourierDifference.total()));

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

    if (RandomNumber::Uniform() < biasTransitionMatrix * Pacc)
    {
      system.components[selectedComponent].mc_moves_probabilities.statistics_SwapDeletionMove_CBMC.accepted += 1;

      system.acceptEwaldMove();
      system.deleteMolecule(selectedComponent, selectedMolecule, molecule);

      return {retraceData.energies - energyFourierDifference - tailEnergyDifference, double3(Pacc, 1.0 - Pacc, 0.0)};
    };
    return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
  }

  return {std::nullopt, double3(0.0, 1.0, 0.0)};
}
