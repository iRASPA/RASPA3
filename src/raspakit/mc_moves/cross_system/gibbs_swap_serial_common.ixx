module;

export module mc_moves_gibbs_swap_serial_common;

import std;

import running_energy;
import atom;
import system;

/**
 * \brief Shared flag-swap bookkeeping of the serial Gibbs CFCMC swap moves (GibbsSwapCFCMC and GibbsSwapCBCFCMC).
 *
 * In serial Gibbs CFCMC every box holds one fractional molecule per swappable component, but only the box with
 * 'containsTheFractionalMolecule' has them active; in the other box they are switched off ('lambda = 0', no
 * dU/dlambda group). The flag is a property of the box, so when a move of one component transfers the activity to
 * the other box, the fractional molecules of all other swappable components follow: they are switched off in the box
 * that becomes inactive and on (at their own current lambda) in the box that becomes active. These helpers compute
 * and apply that combined change.
 */
export namespace MC_Moves::GibbsSwapSerialCommon
{
/// The fractional molecules of a box before a move: ((componentId, fractional index), atoms) per swappable component.
using FractionalSnapshot = std::vector<std::pair<std::pair<std::size_t, std::size_t>, std::vector<Atom>>>;

[[nodiscard]] FractionalSnapshot saveFractionalMolecules(const System& system);
void restoreFractionalMolecules(System& system, const FractionalSnapshot& snapshot);

/// Brick-CFCMC-style aggregated tail-correction difference of a sub-step; advances the threaded effective counts.
[[nodiscard]] RunningEnergy computeTailEnergyDifference(
    System& system, std::vector<double>& tailEffectiveCounts,
    std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups>& tailGroupCounts, std::span<const Atom> newAtoms,
    std::span<const Atom> oldAtoms);

/// Energy difference of all fractional molecules of 'system' taking the scaling of the target activity state (the
/// selected component at 'selectedActiveLambda', the others at their lambda in 'lambdaSourceSystem'), relative to
/// 'beforeSnapshot'. The first Ewald sub-step of the move in 'system': 'storedEik' -> 'trialEik'. The atoms in
/// memory are read as they are (the caller may already have replaced the selected molecule's atoms).
[[nodiscard]] std::optional<RunningEnergy> computeFlagSwapFractionalEnergyDifference(
    System& system, std::vector<double>& tailEffectiveCounts,
    std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups>& tailGroupCounts,
    const FractionalSnapshot& beforeSnapshot, std::size_t selectedComponent, bool systemBecomesActive,
    double selectedActiveLambda, const System& lambdaSourceSystem) noexcept;

/// The same for the fractional molecules of the components other than 'selectedComponent' only. Called after the
/// selected component's own sub-steps in 'system', so the Ewald difference chains on 'trialEik'.
[[nodiscard]] std::optional<RunningEnergy> computeFlagSwapOtherComponentsEnergyDifference(
    System& system, std::vector<double>& tailEffectiveCounts,
    std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups>& tailGroupCounts,
    const FractionalSnapshot& beforeSnapshot, std::size_t selectedComponent, bool systemBecomesActive,
    const System& lambdaSourceSystem) noexcept;

/// Applies the scaling of the target activity state to all fractional molecules of 'system' in memory, and refreshes
/// their records in the system's cell list (which caches position, scaling, and group id) so that a subsequent CBMC
/// grow or retrace sees the edited state. 'restoreFractionalMolecules' refreshes the records the same way.
void applyFlagSwapScaling(System& system, std::size_t selectedComponent, bool systemBecomesActive,
                          double selectedActiveLambda, const System& lambdaSourceSystem);

/// The lambda bins of the other components travel with the activity flag.
void syncOtherComponentLambdaBins(System& systemBecomingInactive, System& systemBecomingActive,
                                  std::size_t selectedComponent);
}  // namespace MC_Moves::GibbsSwapSerialCommon
