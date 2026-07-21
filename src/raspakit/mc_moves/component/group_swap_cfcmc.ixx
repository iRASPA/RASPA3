module;

export module mc_moves_group_swap_cfcmc;

import std;

import double3;
import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * \brief CFCMC group swap: insertion/deletion of a whole group of molecules through a single
 *        shared lambda coordinate.
 *
 * The group consists of one molecule of the selected (central) component plus one molecule per
 * entry of Component::groupComponentIds (satellites; a component may occur multiple times, e.g.
 * two Cl of CaCl2). Every group member holds its own fractional-molecule slot (allocated
 * separately from the GC-CFCMC and pair-CFCMC slots, so the move can be combined with those), and
 * all members are coupled to the lambdaGroupSwap histogram of the central component.
 *
 * Three sub-moves, exactly as in the ion-pair CFCMC move generalized to K = 1 + #satellites
 * molecules:
 *  - lambda change: all K fractional molecules are rescaled together to the new lambda.
 *  - insertion (lambda crosses 1): the K fractional molecules become integer and a new fractional
 *    group is inserted at independent random positions with the residual lambda.
 *  - deletion (lambda crosses 0): the K fractional molecules are switched off and K randomly
 *    selected integer molecules (drawn without replacement per component) become the new
 *    fractional group.
 *
 * The acceptance rules use one factor beta*f_c*V/(N_c + k) per member (k counts members of the
 * same component), which reduces to the pair-CFCMC rule for one satellite and keeps detailed
 * balance for repeated components.
 */
std::pair<std::optional<RunningEnergy>, double3> groupSwapMove_CFCMC(RandomNumber& random, System& system,
                                                                     std::size_t selectedComponent);

/**
 * \brief CB/CFCMC group swap: as groupSwapMove_CFCMC but the new fractional group is grown with
 *        CBMC at insertion and the removed group is retraced with CBMC at deletion (Rosenbluth
 *        weights divided by the ideal-gas weights enter the acceptance rule).
 */
std::pair<std::optional<RunningEnergy>, double3> groupSwapMove_CFCMC_CBMC(RandomNumber& random, System& system,
                                                                          std::size_t selectedComponent);
}  // namespace MC_Moves
