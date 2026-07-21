module;

export module mc_moves_group_swap;

import std;

import double3;
import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * \brief Conventional group insertion without distance bias.
 *
 * Inserts a group of molecules in a single move: one molecule of the selected component (the
 * central molecule) plus the satellite molecules listed in Component::groupComponentIds. The
 * central molecule is grown at a random position in the box; every satellite is grown one after
 * the other with its starting bead placed uniformly inside the sphere of radius
 * Component::maximumGroupDistance around the central molecule's starting bead (the distance is
 * drawn as r = R_max * cbrt(u), uniform in the sphere volume). Intended for inserting a group of
 * charged ions that is overall neutral (e.g. Ca2+ with two Cl- satellites).
 */
std::pair<std::optional<RunningEnergy>, double3> groupInsertionMove(RandomNumber& random, System& system,
                                                                    std::size_t selectedComponent);

/**
 * \brief Conventional group deletion; the exact reverse of groupInsertionMove.
 *
 * Deletes a group in a single move: the given integer molecule of the selected component plus, for
 * every satellite slot, a partner selected uniformly among the integer molecules of the slot's
 * component whose starting bead lies within R_max of the central molecule (sequentially, without
 * replacement). The per-slot candidate counts enter the acceptance rule so that detailed balance
 * with the insertion move is obeyed.
 */
std::pair<std::optional<RunningEnergy>, double3> groupDeletionMove(RandomNumber& random, System& system,
                                                                   std::size_t selectedComponent,
                                                                   std::size_t selectedMolecule);

/**
 * \brief Distance-biased group CBMC insertion (generalization of the Orkoulas-Panagiotopoulos
 *        ion-pair insertion to groups of molecules).
 *
 * Same construction as groupInsertionMove, but every molecule is grown with CBMC and the satellite
 * distances are drawn uniformly in r in [0, R_max], which introduces an explicit distance-bias
 * factor 3 r^2 / R_max^2 per satellite in the acceptance rule.
 */
std::pair<std::optional<RunningEnergy>, double3> groupInsertionMoveCBMC(RandomNumber& random, System& system,
                                                                        std::size_t selectedComponent);

/**
 * \brief Distance-biased group CBMC deletion; the exact reverse of groupInsertionMoveCBMC.
 */
std::pair<std::optional<RunningEnergy>, double3> groupDeletionMoveCBMC(RandomNumber& random, System& system,
                                                                       std::size_t selectedComponent,
                                                                       std::size_t selectedMolecule);
}  // namespace MC_Moves
