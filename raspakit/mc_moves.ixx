export module mc_moves;

import component;
import atom;
import double3;
import double3x3;
import simulationbox;
import randomnumbers;
import system;
import energy_status;
import running_energy;

import <vector>;
import <tuple>;
import <optional>;
import <span>;


export struct MC_Moves
{
    MC_Moves() {}

    void performRandomMove(System& selectedSystem, size_t selectedComponent);
    void performRandomMoveProduction(System& selectedSystem, size_t selectedComponent, size_t currentBlock);

    std::optional<RunningEnergy> translationMove(System& system, size_t selectedComponent, std::span<Atom> molecule) const;
    std::optional<RunningEnergy> randomTranslationMove(System& system, size_t selectedComponent, std::span<Atom> molecule);
    std::optional<RunningEnergy> rotationMove(System& system, size_t selectedComponent, std::span<Atom> molecule);
    std::optional<RunningEnergy> randomRotationMove(System& system, size_t selectedComponent, std::span<Atom> molecule);
    std::optional<RunningEnergy> volumeMove([[maybe_unused]] System &system) const;
    std::optional<RunningEnergy> reinsertionMove(System& system, size_t selectedComponent, size_t selectedMolecule, std::span<Atom> atoms);
    std::optional<RunningEnergy> insertionMove(System& system, size_t selectedComponent);
    std::optional<RunningEnergy> deletionMove(System& system, size_t selectedComponent, size_t selectedMolecule);
    std::optional<RunningEnergy> swapMove_CFCMC_CBMC(System& system, size_t selectedComponent, size_t selectedMolecule, bool insertionDisabled=false, bool deletionDisabled=false);
    std::optional<double> WidomMove(System& system, size_t selectedComponent);

    double energyOverlapCriteria = 1e6;
    bool useDualCutOff{ true };
};
