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



std::optional<RunningEnergy> MC_Moves::deletionMove(System& system, size_t selectedComponent, size_t selectedMolecule)
{
    system.components[selectedComponent].statistics_SwapDeletionMove_CBMC.counts += 1;
    
    if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] > 0)
    {
        system.components[selectedComponent].statistics_SwapDeletionMove_CBMC.constructed += 1;

        std::span<Atom> molecule = system.spanOfMolecule(selectedComponent, selectedMolecule);

        std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
        ChainData retraceData = system.retraceMoleculeSwapDeletion(selectedComponent, selectedMolecule, molecule, 1.0, 0.0);
        std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
        system.components[selectedComponent].cpuTime_SwapDeletionMove_CBMC_NonEwald += (t2 - t1);

        std::chrono::system_clock::time_point u1 = std::chrono::system_clock::now();
        RunningEnergy energyFourierDifference = system.energyDifferenceEwaldFourier(system.storedEik, {}, molecule);
        std::chrono::system_clock::time_point u2 = std::chrono::system_clock::now();
        system.components[selectedComponent].cpuTime_SwapDeletionMove_CBMC_Ewald += (u2 - u1);

        //EnergyStatus tailEnergyDifference = system.computeTailCorrectionVDWRemoveEnergy(selectedComponent) - 
        //                                    system.computeTailCorrectionVDWOldEnergy();
        RunningEnergy tailEnergyDifference;
        double correctionFactorEwald = std::exp(-system.simulationBox.Beta * (energyFourierDifference.total() + tailEnergyDifference.total()));

        double idealGasRosenbluthWeight = system.components[selectedComponent].idealGasRosenbluthWeight.value_or(1.0);
        double preFactor = correctionFactorEwald * double(system.numberOfMoleculesPerComponent[selectedComponent]) /
                           (system.simulationBox.Beta * system.components[selectedComponent].molFraction * 
                            system.simulationBox.pressure * system.simulationBox.volume);
        if (RandomNumber::Uniform() < preFactor * idealGasRosenbluthWeight / retraceData.RosenbluthWeight)
        {
            system.components[selectedComponent].statistics_SwapDeletionMove_CBMC.accepted += 1;

            system.acceptEwaldMove();
            system.deleteMolecule(selectedComponent, selectedMolecule, molecule);

            // Debug
            //assert(system.checkMoleculeIds());

            return retraceData.energies - energyFourierDifference - tailEnergyDifference;
        };
    }

    return std::nullopt;
}
