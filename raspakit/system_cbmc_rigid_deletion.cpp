module;

module system;

import randomnumbers;
import component;
import atom;
import double3;
import double3x3;
import simulationbox;
import energy_status;
import cbmc;
import cbmc_growing_status;
import forcefield;

import <vector>;
import <tuple>;
import <optional>;
import <span>;

import <iostream>;
import <algorithm>;
import <numeric>;


[[nodiscard]] ChainData System::retraceMoleculeSwapDeletion(size_t selectedComponent, [[maybe_unused]] size_t selectedMolecule, std::span<Atom> molecule, double scaling, double storedR) const noexcept
{
	size_t startingBead = components[selectedComponent].startingBead;

	const FirstBeadData firstBeadData = retraceMultipleFirstBeadSwapDeletion(molecule[startingBead], scaling, storedR);

	const ChainData rigidRotationData = retraceChain(startingBead, scaling, molecule);

	return ChainData(std::vector<Atom>(molecule.begin(), molecule.end()), firstBeadData.energies + rigidRotationData.energies, firstBeadData.RosenbluthWeight * rigidRotationData.RosenbluthWeight, 0.0);
}

[[nodiscard]] FirstBeadData System::retraceMultipleFirstBeadSwapDeletion(const Atom& atom, double scaling, [[maybe_unused]] double storedR) const noexcept
{
	std::vector<Atom> trialPositions(numberOfTrialDirections, atom);
    for(Atom &trialPosition: trialPositions) {trialPosition.setScaling(scaling);}
    std::for_each(trialPositions.begin() + 1, trialPositions.end(),
            [this](Atom& a) {a.position = simulationBox.randomPosition();});

	const std::vector<std::pair<Atom, EnergyStatus>> externalEnergies = computeExternalNonOverlappingEnergies(trialPositions);

	std::vector<double> logBoltmannFactors{};
	std::transform(std::begin(externalEnergies), std::end(externalEnergies), std::back_inserter(logBoltmannFactors),
		[this](const std::pair<Atom, EnergyStatus>& v) {return -simulationBox.Beta * v.second.totalEnergy; });

	double RosenbluthWeight = std::reduce(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
		[](const double& acc, const double& logBoltmannFactor) {return acc + std::exp(logBoltmannFactor); });

	return FirstBeadData(atom, externalEnergies[0].second, RosenbluthWeight / double(numberOfTrialDirections), 0.0);

}


[[nodiscard]] ChainData System::retraceChain(size_t startingBead, double scaling, std::span<Atom> molecule) const noexcept
{
	std::vector<Atom> trialPosition = std::vector<Atom>(molecule.begin(), molecule.end());
	std::for_each(trialPosition.begin(), trialPosition.end(), [&](Atom& a) {a.setScaling(scaling); });
	std::vector<std::vector<Atom>> trialPositions = { trialPosition };

	for (size_t i = 1; i < numberOfTrialDirections; ++i)
	{
		trialPositions.push_back(rotateRandomlyAround(trialPosition, startingBead));
	};

	const std::vector<std::pair<std::vector<Atom>, EnergyStatus>> externalEnergies = computeExternalNonOverlappingEnergies(trialPositions, std::make_signed_t<std::size_t>(startingBead));

	std::vector<double> logBoltmannFactors{};
	std::transform(std::begin(externalEnergies), std::end(externalEnergies),
		std::back_inserter(logBoltmannFactors), [this](const std::pair<std::vector<Atom>, EnergyStatus>& v) {return -simulationBox.Beta * v.second.totalEnergy; });

	double RosenbluthWeight = std::reduce(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
		[](const double& acc, const double& logBoltmannFactor) {return acc + std::exp(logBoltmannFactor); });

	return ChainData(trialPositions[0], externalEnergies[0].second, RosenbluthWeight / double(numberOfTrialDirections), 0.0);
}
