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



// LogBoltzmannFactors are (-Beta U)
size_t System::selectTrialPosition(std::vector <double> LogBoltzmannFactors) const noexcept
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
	size_t selected = 0;
	double cumw = ShiftedBoltzmannFactors[0];
	double ws = RandomNumber::Uniform() * SumShiftedBoltzmannFactors;
	while (cumw < ws)
		cumw += ShiftedBoltzmannFactors[++selected];

	return selected;
}

std::vector<Atom> System::filterNonOverlappingTrialPositions(std::vector <std::pair<Atom, EnergyStatus>> trialPositions) const noexcept
{
	std::vector <std::pair<Atom, EnergyStatus>> nonOverlappingAtomEnergies{};
	std::copy_if(std::begin(trialPositions), std::end(trialPositions),
		std::back_inserter(nonOverlappingAtomEnergies), [this](const std::pair<Atom, EnergyStatus>& v) {return v.second.totalEnergy < forceField.overlapCriteria; });

	std::vector <Atom> nonOverlappingAtoms{};
	std::transform(std::begin(nonOverlappingAtomEnergies), std::end(nonOverlappingAtomEnergies),
		std::back_inserter(nonOverlappingAtoms), [](const std::pair<Atom, EnergyStatus>& v) {return v.first; });

	return nonOverlappingAtoms;
}

