module;

export module dctst_moves;

import std;

import system;
import randomnumbers;
import running_energy;
import mc_moves_move_types;

/**
 * \brief The Monte Carlo moves of the molecules around the tagged one in the dcTST drivers.
 *
 * Both drivers sample a fixed-loading (NVT) ensemble with a constraint on one molecule: the other molecules
 * move with the single-molecule moves of their component, dispatched here. Every move the drivers do not
 * perform is refused up front ('validate', called by the input reader and by the constructors); a move that
 * is not in the list was performed as a rigid translation, which froze the conformations of flexible
 * molecules and booked its trials under the wrong move.
 */
export namespace DcTSTMoves
{
/// Whether the drivers perform 'move' on a molecule other than the tagged one.
[[nodiscard]] bool supported(Move::Types move);

/// Performs 'move' on the molecule; throws std::logic_error for a move that is not 'supported'.
std::optional<RunningEnergy> perform(RandomNumber& random, System& system, std::size_t selectedComponent,
                                     std::size_t selectedMolecule, Move::Types move);

/// Refuses every move of every component that the drivers do not perform (particle-number, identity-change,
/// fractional-molecule, volume, cross-system and multi-molecule moves among them), a polyatomic tagged molecule
/// that would only be translated ('taggedReinsertionProbability' zero: its orientation and conformation never
/// change), and,
/// with the occupancy constraint ('occupancyConstraint'), the moves whose undo the constraint does not cover.
void validate(const System& system, std::size_t taggedComponent, std::string_view driverName,
              double taggedReinsertionProbability, bool occupancyConstraint);
}  // namespace DcTSTMoves
