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
		std::optional<ChainData> growData = system.growMoleculeReinsertion(selectedComponent, selectedMolecule, molecule);

		if (!growData) return std::nullopt;

        std::span<const Atom> newMolecule = std::span(growData->atom.begin(), growData->atom.end());

		system.components[selectedComponent].statistics_ReinsertionMove_CBMC.constructed += 1;

		ChainData retraceData = system.retraceMoleculeReinsertion(selectedComponent, selectedMolecule, molecule, growData->storedR);
        EnergyStatus energyFourierDifference = system.energyDifferenceEwaldFourier(system.storedEik, newMolecule, molecule);
        double correctionFactorFourier = std::exp(-system.simulationBox.Beta * energyFourierDifference.totalEnergy);

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
