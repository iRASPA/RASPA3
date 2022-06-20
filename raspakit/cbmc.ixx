export module cbmc;

import atom;
import energy_status;
import double3x3;
import double3;
import energy_status_intra;
import energy_status_inter;

import <vector>;

export struct FirstBeadData
{
	Atom atom;
	EnergyStatus energies;
	double RosenbluthWeight;
	double storedR;

	FirstBeadData(size_t size) noexcept : energies(size)
	{

	}

	FirstBeadData(Atom atom, EnergyStatus energies, double RosenbluthWeight, double storedR) noexcept :
		atom(atom), energies(energies), RosenbluthWeight(RosenbluthWeight), storedR(storedR)
	{
	}

	FirstBeadData() = delete;
	FirstBeadData(const FirstBeadData& a) noexcept = default;
	FirstBeadData& operator=(const FirstBeadData& a) noexcept = default;
	FirstBeadData(FirstBeadData&& a) noexcept = default;
	FirstBeadData& operator=(FirstBeadData&& a) noexcept = default;
	~FirstBeadData() noexcept = default;
};

export struct ChainData
{
	std::vector<Atom> atom;
	EnergyStatus energies;
	double RosenbluthWeight;
	double storedR;

	ChainData(std::vector<Atom> atom, EnergyStatus energies, double RosenbluthWeight, double storedR) noexcept :
		atom(atom), energies(energies), RosenbluthWeight(RosenbluthWeight), storedR(storedR)
	{
	}
	ChainData() = delete;
	ChainData(const ChainData& a) noexcept = default;
	ChainData& operator=(const ChainData& a) noexcept = default;
	ChainData(ChainData&& a) noexcept = default;
	ChainData& operator=(ChainData&& a) noexcept = default;
	~ChainData() noexcept = default;
};

export inline std::vector<Atom> rotateRandomlyAround(std::vector<Atom> atoms, size_t startingBead)
{
	double3x3 randomRotationMatrix = double3x3::randomRotationMatrix();
	std::vector<Atom> randomlyRotatedAtoms{};
	for (size_t i = 0; i < atoms.size(); ++i)
	{
		Atom b = atoms[i];
		b.position = atoms[startingBead].position + randomRotationMatrix * (b.position - atoms[startingBead].position);
		randomlyRotatedAtoms.push_back(b);
	}
	return randomlyRotatedAtoms;
}
