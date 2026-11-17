module;

export module mc_moves_gibbs_swap_serial_common;

import std;

import running_energy;
import atom;
import system;

/**
 * \brief Shared bookkeeping of the serial Gibbs CFCMC swap moves (GibbsSwapCFCMC and GibbsSwapCBCFCMC).
 *
 * In serial Gibbs CFCMC every box holds one fractional molecule per swappable component, but for each component
 * only the box with 'System::containsTheFractionalMolecule(component)' has it active; in the other box it is a
 * switched-off placeholder ('lambda = 0', no dU/dlambda group). Every component has its own flag: a move of one
 * component transfers the activity of its own fractional molecule and leaves those of the other components where
 * they are. The position of a placeholder is never used: the fractional molecule that becomes active in a box is
 * either drawn anew there (the shuffle) or an integer molecule of that box (the interchange). These helpers compute
 * and apply the change of the activity of the selected component's fractional molecule.
 */
export namespace MC_Moves::GibbsSwapSerialCommon
{
/// The fractional molecule(s) of a box before a move: ((componentId, fractional index), atoms).
using FractionalSnapshot = std::vector<std::pair<std::pair<std::size_t, std::size_t>, std::vector<Atom>>>;

/// The snapshot holds the fractional molecule of 'selectedComponent' (the only one a move changes).
[[nodiscard]] FractionalSnapshot saveFractionalMolecules(const System& system, std::size_t selectedComponent);
void restoreFractionalMolecules(System& system, const FractionalSnapshot& snapshot);

/// Brick-CFCMC-style aggregated tail-correction difference of a sub-step; advances the threaded effective counts.
[[nodiscard]] RunningEnergy computeTailEnergyDifference(
    System& system, std::vector<double>& tailEffectiveCounts,
    std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups>& tailGroupCounts, std::span<const Atom> newAtoms,
    std::span<const Atom> oldAtoms);

/// External-field difference of a sub-step (zero without an external field); 'std::nullopt' when a new position
/// lies outside a confining geometry. The field scales with the fractional scaling of the atoms, so every sub-step
/// that changes positions or scalings carries it, like the framework and intermolecular terms.
[[nodiscard]] std::optional<RunningEnergy> computeExternalFieldEnergyDifference(const System& system,
                                                                                std::span<const Atom> newAtoms,
                                                                                std::span<const Atom> oldAtoms) noexcept;

/// Energy difference of the fractional molecule in 'beforeSnapshot' taking the scaling of the target activity state
/// ('activeLambda' when the box becomes the active one of the component, switched off otherwise), relative to
/// 'beforeSnapshot'. The first Ewald sub-step of the move in 'system': 'storedEik' -> 'trialEik'. The atoms in
/// memory are read as they are (the caller may already have replaced the molecule's atoms).
/// 'runningNetCharge' is the Coulomb-scaled net charge threaded across the sub-steps like the tail counts: the
/// net-charge correction expands around it, and it is advanced by the change of this sub-step.
[[nodiscard]] std::optional<RunningEnergy> computeFlagSwapFractionalEnergyDifference(
    System& system, std::vector<double>& tailEffectiveCounts,
    std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups>& tailGroupCounts, double& runningNetCharge,
    const FractionalSnapshot& beforeSnapshot, bool systemBecomesActive, double activeLambda) noexcept;

/// Applies the scaling of the target activity state to the fractional molecule of 'selectedComponent' in memory, and
/// refreshes its records in the system's cell list (which caches position, scaling, and group id) so that a
/// subsequent CBMC grow or retrace sees the edited state. 'restoreFractionalMolecules' refreshes the records the
/// same way.
void applyFlagSwapScaling(System& system, std::size_t selectedComponent, bool systemBecomesActive,
                          double activeLambda);

/// The fractional molecule of 'selectedComponent' now belongs to 'systemBecomingActive'; the flags of the other
/// components are untouched.
void transferFractionalMoleculeOwnership(System& systemBecomingInactive, System& systemBecomingActive,
                                         std::size_t selectedComponent);

/// Moves the switched-off fractional molecule (the placeholder) of 'selectedComponent' to the periodic image with
/// its centre inside the box of 'system'. A placeholder enters a box with the Cartesian coordinates it had in the
/// other box, and every volume move of its new box scales them: carried back and forth between two boxes whose
/// volumes fluctuate, the coordinates performed a random walk on a logarithmic scale and grew without bound (until
/// the stencil of the cell list overflowed). The position of a placeholder is no coordinate of the state (it has no
/// interactions and is never switched on in place), and the image is a lattice translation, so no energy changes.
void parkInactiveFractionalMolecule(System& system, std::size_t selectedComponent);
}  // namespace MC_Moves::GibbsSwapSerialCommon
