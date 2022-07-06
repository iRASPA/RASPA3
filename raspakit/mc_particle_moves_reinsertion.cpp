module;

module mc_particle_moves;

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


std::optional<EnergyStatus> MC_Particle_Moves::reinsertionMove(System& system, size_t selectedComponent, size_t selectedMolecule, std::span<Atom> molecule)
{
	system.components[selectedComponent].statistics_ReinsertionMove_CBMC.counts += 1;

	if (system.numberOfMoleculesPerComponent[selectedComponent] > 0)
	{
        std::chrono::system_clock::time_point t1 = std::chrono::system_clock::now();
		std::optional<ChainData> growData = system.growMoleculeReinsertion(selectedComponent, selectedMolecule, molecule);
        std::chrono::system_clock::time_point t2 = std::chrono::system_clock::now();
        system.components[selectedComponent].cpuTime_ReinsertionGrowMove_CBMC_NonEwald += (t2 - t1);

		if (!growData) return std::nullopt;

        std::span<const Atom> newMolecule = std::span(growData->atom.begin(), growData->atom.end());

		system.components[selectedComponent].statistics_ReinsertionMove_CBMC.constructed += 1;

        std::chrono::system_clock::time_point u1 = std::chrono::system_clock::now();
		ChainData retraceData = system.retraceMoleculeReinsertion(selectedComponent, selectedMolecule, molecule, growData->storedR);
        std::chrono::system_clock::time_point u2 = std::chrono::system_clock::now();
        system.components[selectedComponent].cpuTime_ReinsertionRetraceMove_CBMC_NonEwald += (u2 - u1);

        std::chrono::system_clock::time_point v1 = std::chrono::system_clock::now();
        EnergyStatus energyFourierDifference = system.energyDifferenceEwaldFourier(system.storedEik, newMolecule, molecule);
        std::chrono::system_clock::time_point v2 = std::chrono::system_clock::now();
        system.components[selectedComponent].cpuTime_ReinsertionMove_CBMC_Ewald += (v2 - v1);

        double correctionFactorFourier = std::exp(-system.simulationBox.Beta * energyFourierDifference.totalEnergy.energy);

		if (RandomNumber::Uniform() < correctionFactorFourier * growData->RosenbluthWeight / retraceData.RosenbluthWeight)
		{
			system.components[selectedComponent].statistics_ReinsertionMove_CBMC.accepted += 1;

            system.acceptEwaldMove();
			std::copy(newMolecule.begin(), newMolecule.end(), molecule.begin());

			return (growData->energies - retraceData.energies) + energyFourierDifference;
		};
	}

	return std::nullopt;
}
