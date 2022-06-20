/*

module;

module cbmc_linear_chain;

import randomnumbers;
import component;
import atom;
import double3;
import double3x3;
import simulationbox;
import energy_status;

import <vector>;
import <tuple>;
import <optional>;
import <span>;

import <iostream>;
import <algorithm>;
import <numeric>;

[[nodiscard]] std::optional<ChainData> CBMCLinearChain::growMoleculeInsertion(int selectedComponent, int selectedMolecule, std::span<Atom> molecule, double scaling) const noexcept
{
	int start = system.components[selectedComponent].startingBead;

	std::optional<FirstBeadData> firstBeadData = growMultipleFirstBead(selectedComponent, selectedMolecule, molecule, scaling);

	if (!firstBeadData) return std::nullopt;

	return ChainData({ firstBeadData->atom}, firstBeadData->energies, firstBeadData->logRosenbluthWeight, firstBeadData->storedR);
}

[[nodiscard]] ChainData CBMCLinearChain::retraceMoleculeInsertion(int selectedComponent, int selectedMolecule, std::span<Atom> molecule, double storedR) const noexcept
{
	Atom atom{};
	FirstBeadData firstBeadData = retraceMultipleFirstBead(selectedComponent, selectedMolecule, *molecule.begin(), storedR);

	return ChainData({ firstBeadData.atom }, firstBeadData.energies, firstBeadData.logRosenbluthWeight, firstBeadData.storedR);
}

[[nodiscard]] std::optional<FirstBeadData> CBMCLinearChain::growMultipleFirstBead(int selectedComponent, int selectedMolecule, std::span<Atom> molecule, double scaling) const noexcept
{
	std::vector<Atom> trialPositions(numberOfTrialDirections);
	std::for_each(trialPositions.begin(), trialPositions.end(),
		[this, selectedComponent, selectedMolecule, scaling](Atom& a)
		{a.position = system.simulationBox.randomPosition();
		 a.type = system.components[selectedComponent].atomTypes[0];
		 a.componentId = selectedComponent;
		 a.moleculeId = selectedMolecule;
		 a.scaling = scaling; });

	const std::vector<EnergyStatus> externalEnergies = system.computeExternalEnergies(trialPositions);

	if (externalEnergies.empty()) return std::nullopt;

	std::vector<double> logBoltmannFactors{};
	std::transform(std::begin(externalEnergies), std::end(externalEnergies),
		std::back_inserter(logBoltmannFactors), [this](const EnergyStatus& v) {return -system.simulationBox.Beta * v.totalEnergy; });

	int selected = selectTrialPosition(logBoltmannFactors);

	double RosenbluthWeight = std::reduce(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
		[](const double& acc, const double &logBoltmannFactor) {return acc+std::exp(logBoltmannFactor); });

	if(RosenbluthWeight < minimumRosenbluthFactor) return std::nullopt;

	// r=w(n)-exp(-beta U[h_n]) Eq.16 from Esselink et al.
	double storedR = RosenbluthWeight - std::exp(logBoltmannFactors[selected]);

	return FirstBeadData(trialPositions[selected], externalEnergies[selected], RosenbluthWeight / double(numberOfTrialDirections), storedR);
}

[[nodiscard]] FirstBeadData CBMCLinearChain::retraceMultipleFirstBead(int selectedComponent, int selectedMolecule, const Atom& atom, double storedR) const noexcept
{
	std::vector<Atom> trialPositions(numberOfTrialDirections);

	std::for_each(trialPositions.begin(), trialPositions.end(),
		[this, selectedComponent, selectedMolecule, atom](Atom& a)
		{a.position = system.simulationBox.randomPosition();
		 a.type = system.components[selectedComponent].atomTypes[0];
		 a.componentId = selectedComponent;
		 a.moleculeId = selectedMolecule;
		 a.scaling = atom.scaling; });

		 trialPositions[0] = atom;

	const std::vector<EnergyStatus> externalEnergies = system.computeExternalEnergies(trialPositions);

	std::vector<double> logBoltmannFactors{};
	std::transform(std::begin(externalEnergies), std::end(externalEnergies), std::back_inserter(logBoltmannFactors),
				   [this](const EnergyStatus& v) {return -system.simulationBox.Beta * v.totalEnergy; });

	double RosenbluthWeight = std::reduce(logBoltmannFactors.begin(), logBoltmannFactors.end(), 0.0,
		[](const double& acc, const double& logBoltmannFactor) {return acc + std::exp(logBoltmannFactor); });

	// w(o)=exp(-beta u(o))+r  Eq. 18 from Esselink et al.
	return FirstBeadData(atom, externalEnergies[0], (RosenbluthWeight + storedR) / double(numberOfTrialDirections), 0.0);

}


// LogBoltzmannFactors are (-Beta U)
int CBMCLinearChain::selectTrialPosition(std::vector <double> LogBoltzmannFactors) const noexcept
{
	std::vector<double> ShiftedBoltzmannFactors(LogBoltzmannFactors.size());

	// Energies are always bounded from below [-U_max, infinity>
	// Find the lowest energy value, i.e. the largest value of (-Beta U)
	double largest_value = *std::max_element(LogBoltzmannFactors.begin(), LogBoltzmannFactors.end());

	// Standard trick: shift the Boltzmann factors down to avoid numerical problems
	// The largest value of 'ShiftedBoltzmannFactors' will be 1 (which corresponds to the lowest energy).
	double SumShiftedBoltzmannFactors = 0.0;
	for (size_t i = 0; i < LogBoltzmannFactors.size(); ++i)
	{
		ShiftedBoltzmannFactors[i] = exp(LogBoltzmannFactors[i] - largest_value);
		SumShiftedBoltzmannFactors += ShiftedBoltzmannFactors[i];
	}

	// select the Boltzmann factor
	int selected = 0;
	double cumw = ShiftedBoltzmannFactors[0];
	double ws = RandomNumber::Uniform() * SumShiftedBoltzmannFactors;
	while (cumw < ws)
		cumw += ShiftedBoltzmannFactors[++selected];

	return selected;
}

std::vector<Atom> CBMCLinearChain::filterNonOverlappingTrialPositions(std::vector <std::pair<Atom, EnergyStatus>> trialPositions) const noexcept
{
	std::vector <std::pair<Atom, EnergyStatus>> nonOverlappingAtomEnergies{};
	std::copy_if(std::begin(trialPositions), std::end(trialPositions),
		std::back_inserter(nonOverlappingAtomEnergies), [this](const std::pair<Atom, EnergyStatus>& v) {return v.second.totalEnergy < overlapCriteria; });

	std::vector <Atom> nonOverlappingAtoms{};
	std::transform(std::begin(nonOverlappingAtomEnergies), std::end(nonOverlappingAtomEnergies),
		std::back_inserter(nonOverlappingAtoms), [this](const std::pair<Atom, EnergyStatus>& v) {return v.first; });

	return nonOverlappingAtoms;
}

*/