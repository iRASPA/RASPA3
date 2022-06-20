export module cbmc_branched_chain;

/*
import component;
import atom;
import double3;
import double3x3;
import simulationbox;
import energy_status;
import system;

import <vector>;
import <tuple>;
import <optional>;
import <span>;

struct FirstBeadData
{
	Atom atom;
	EnergyStatus energies;
	double logRosenbluthWeight;
	double storedR;

	FirstBeadData(size_t size) : energies(size)
	{

	}

	FirstBeadData(Atom atom, EnergyStatus energies, double logRosenbluthWeight, double storedR) :
		atom(atom), energies(energies), logRosenbluthWeight(logRosenbluthWeight), storedR(storedR)
	{
	}
};

export struct ChainData
{
	std::vector<Atom> atom;
	EnergyStatus energies;
	double logRosenbluthWeight;
	double storedR;

	ChainData(std::vector<Atom> atom, EnergyStatus energies, double logRosenbluthWeight, double storedR) :
		atom(atom), energies(energies), logRosenbluthWeight(logRosenbluthWeight), storedR(storedR)
	{
	}
};

export struct CBMCLinearChain
{
	CBMCLinearChain(const System& system) :
		system(system),
		numberOfTrialDirections(16)
	{

	}

	[[nodiscard]] std::optional<ChainData> growMoleculeInsertion(int selectedComponent, int selectedMolecule, std::span<Atom> atoms, double scaling) const noexcept;
	[[nodiscard]] ChainData retraceMoleculeInsertion(int selectedComponent, int selectedMolecule, std::span<Atom> atoms, double storedR) const noexcept;

	[[nodiscard]] std::optional<FirstBeadData> growMultipleFirstBead(int selectedComponent, int selectedMolecule, std::span<Atom> atoms, double scaling) const noexcept;
	[[nodiscard]] FirstBeadData retraceMultipleFirstBead(int selectedComponent, int selectedMolecule, const Atom& atom, double storedR) const noexcept;

	int selectTrialPosition(std::vector <double> BoltzmannFactors) const noexcept;
	std::vector<Atom> filterNonOverlappingTrialPositions(std::vector <std::pair<Atom, EnergyStatus>> trialPositions) const noexcept;

	const System& system;
	int numberOfTrialDirections;
	double overlapCriteria = 1e6;
	double minimumRosenbluthFactor = 1e-150;
};
*/