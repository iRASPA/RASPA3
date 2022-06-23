export module mc_particle_moves;

import component;
import atom;
import double3;
import double3x3;
import simulationbox;
import randomnumbers;
import system;
import energy_status;

import <vector>;
import <tuple>;
import <optional>;
import <span>;


export struct MC_Particle_Moves
{
	enum class ParticleMove : size_t
	{
		Translation = 0,
		Rotation = 1,
		CBMC_Reinsertion = 2,
		CBMC_Swap = 3,
		CBMC_Gibbs_Swap = 4,
		Volume = 5,
		IdentityChange = 6,
		CFCMC_Swap = 7,
		CFCMC_CBMC_Swap = 8,
        WidomMove = 9,
		NumberOfMoves = 10
	};


	MC_Particle_Moves() {}

	void performRandomMove(System& selectedSystem, size_t selectedComponent);
	void performRandomMoveProduction(System& selectedSystem, size_t selectedComponent, size_t currentBlock);

	std::optional<EnergyStatus> translationMove(System& system, size_t selectedComponent, std::span<Atom> molecule);
	std::optional<EnergyStatus> rotationMove(System& system, size_t selectedComponent, std::span<Atom> molecule);
	std::optional<EnergyStatus> reinsertionMove(System& system, size_t selectedComponent, size_t selectedMolecule, std::span<Atom> atoms);
	std::optional<EnergyStatus> insertionMove(System& system, size_t selectedComponent);
	std::optional<EnergyStatus> deletionMove(System& system, size_t selectedComponent, size_t selectedMolecule);
	std::optional<EnergyStatus> swapMove_CFCMC_CBMC(System& system, size_t selectedComponent, size_t selectedMolecule, bool insertionDisabled=false, bool deletionDisabled=false);
	std::optional<double> WidomMove(System& system, size_t selectedComponent);
    std::optional<EnergyStatus> WidomMove_CFCMC_CBMC(System& system, size_t selectedComponent);

	double energyOverlapCriteria = 1e6;
};
