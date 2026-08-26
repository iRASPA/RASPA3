module;

export module mc_moves_identity_switch;

import std;

import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * \brief Performs a canonical CBMC identity-switch move (the MCCCS-MN "swatch" move).
 *
 * Selects one molecule of \p selectedComponent and one molecule of a partner component (drawn
 * from \c Component::identitySwitches) in the same box and exchanges their identities: a molecule
 * of the partner type is grown with CBMC at the first molecule's starting-bead position and vice
 * versa, while both old molecules are retraced.
 *
 * In contrast to \c MC_Moves::identityChangeMove (the semi-grand "swotch" move, which converts a
 * single molecule and is driven by the imposed fugacity ratio), this move exchanges a *pair* and
 * therefore conserves the number of molecules of every component.  Fugacities, ideal-gas
 * Rosenbluth weights and the number-of-molecules factors cancel identically, so the acceptance
 * rule contains only the Rosenbluth ratio of the two grows over the two retraces, corrected for
 * the terms the CBMC weights do not contain (Fourier-space Ewald, polarization, and the direct
 * interaction between the two exchanged molecules).  Tail corrections cancel exactly because the
 * multiset of pseudo-atom types in the box is unchanged.
 *
 * Because composition is conserved, the move imposes no thermodynamic constraint of its own: it
 * leaves every equilibrium average untouched and serves purely to accelerate sampling in dense
 * multicomponent systems, where insertion and deletion of whole molecules is inefficient.
 *
 * Requires \c Component::identitySwitches to be populated with the partner component indices.
 *
 * \param random Random number generator instance.
 * \param system The simulation system.
 * \param selectedComponent Index of the component selected for the move attempt.
 *
 * \return Energy difference if the move is accepted; \c std::nullopt otherwise.
 */
std::optional<RunningEnergy> identitySwitchMove(RandomNumber &random, System &system, std::size_t selectedComponent);
}  // namespace MC_Moves
