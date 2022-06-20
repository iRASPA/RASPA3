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



std::optional<double> MC_Particle_Moves::WidomMove(System& system, size_t selectedComponent)
{
	size_t selectedMolecule = system.numberOfMoleculesPerComponent[selectedComponent];
	system.components[selectedComponent].statistics_WidomMove_CBMC.counts += 1;
	
	std::optional<ChainData> growData = system.growMoleculeSwapInsertion(selectedComponent, selectedMolecule, 1.0);
    std::span<const Atom> newMolecule = std::span(growData->atom.begin(), growData->atom.end());
	if (!growData) return std::nullopt;

	system.components[selectedComponent].statistics_WidomMove_CBMC.constructed += 1;

    EnergyStatus energyFourierDifference = system.energyDifferenceEwaldFourier(system.storedEik, newMolecule, {});
    EnergyStatus tailEnergyDifference = system.computeTailCorrectionVDWAddEnergy(selectedComponent) - 
                                        system.computeTailCorrectionVDWOldEnergy();
    double correctionFactorEwald = std::exp(-system.simulationBox.Beta * (energyFourierDifference.totalEnergy + tailEnergyDifference.totalEnergy));


	double idealGasRosenbluthWeight = system.components[selectedComponent].idealGasRosenbluthWeight.value_or(1.0);

	return  correctionFactorEwald * growData->RosenbluthWeight / idealGasRosenbluthWeight;
}

