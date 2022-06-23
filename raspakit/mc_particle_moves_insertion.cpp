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


std::optional<EnergyStatus> MC_Particle_Moves::insertionMove(System& system, size_t selectedComponent)
{
	size_t selectedMolecule = system.numberOfMoleculesPerComponent[selectedComponent];
	system.components[selectedComponent].statistics_SwapInsertionMove_CBMC.counts += 1;
	
	std::optional<ChainData> growData = system.growMoleculeSwapInsertion(selectedComponent, selectedMolecule, 1.0);
    std::span<const Atom> newMolecule = std::span(growData->atom.begin(), growData->atom.end());
	if (!growData) return std::nullopt;

	system.components[selectedComponent].statistics_SwapInsertionMove_CBMC.constructed += 1;

    EnergyStatus energyFourierDifference = system.energyDifferenceEwaldFourier(system.storedEik, newMolecule, {});
    EnergyStatus tailEnergyDifference = system.computeTailCorrectionVDWAddEnergy(selectedComponent) - 
                                        system.computeTailCorrectionVDWOldEnergy();
    double correctionFactorEwald = std::exp(-system.simulationBox.Beta * (energyFourierDifference.totalEnergy.energy + tailEnergyDifference.totalEnergy.energy));


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

