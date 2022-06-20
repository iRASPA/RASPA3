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


std::optional<EnergyStatus> MC_Particle_Moves::translationMove(System & system, size_t selectedComponent, std::span<Atom> molecule)
{
	double3 displacement{};
	double3 maxDisplacement = system.components[selectedComponent].statistics_TranslationMove.maxChange;
	size_t selectedDirection = size_t(3.0 * RandomNumber::Uniform());
	displacement[selectedDirection] = maxDisplacement[selectedDirection] * 2.0 * (RandomNumber::Uniform() - 0.5);
	system.components[selectedComponent].statistics_TranslationMove.counts[selectedDirection] += 1;

	std::vector<Atom> trialPositions(molecule.size());
	std::transform(molecule.begin(), molecule.end(), trialPositions.begin(),
		[&](Atom a) { a.position += displacement; return a; });
    std::span<Atom> newMolecule{trialPositions.begin(), trialPositions.end()};

    std::optional<EnergyStatus> frameworkMolecule = system.computeFrameworkMoleculeEnergyDifference(newMolecule, molecule);
	if (!frameworkMolecule.has_value()) return std::nullopt;

    std::optional<EnergyStatus> interMolecule = system.computeInterMolecularEnergyDifference(newMolecule, molecule);
	if (!interMolecule.has_value()) return std::nullopt;

    EnergyStatus ewaldFourierEnergy = system.energyDifferenceEwaldFourier(system.storedEik, newMolecule, molecule);
    EnergyStatus energyDifference = frameworkMolecule.value() + interMolecule.value() + ewaldFourierEnergy;

	system.components[selectedComponent].statistics_TranslationMove.constructed[selectedDirection] += 1;

    if (RandomNumber::Uniform() < std::exp(-system.simulationBox.Beta * energyDifference.totalEnergy))
	{
		system.components[selectedComponent].statistics_TranslationMove.accepted[selectedDirection] += 1;

        system.acceptEwaldMove();
		std::copy(trialPositions.cbegin(), trialPositions.cend(), molecule.begin());

		return energyDifference;
	};
	return std::nullopt;
}
