module;

export module mc_moves_tethered_proton_hop;

import std;

import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * \brief Tethered proton-hop move: relocate a proton among its discrete candidate sites.
 *
 * The selected component represents Broensted protons modelled as single-site cations. Every proton
 * (one integer molecule of the component) is tethered to a group of discrete candidate positions,
 * stored in Component::tetheredProtonHopSiteGroups (e.g. the four oxygens around one framework Al).
 * The move selects the candidate group belonging to the chosen proton, identifies the site the
 * proton currently occupies (the nearest candidate) and proposes a jump to one of the other sites of
 * the same group, drawn uniformly. The proposal is symmetric (P(i->j) = P(j->i) = 1/(n-1)), so the
 * move is accepted with the plain Metropolis rule on the total energy difference; no bias factor
 * enters.
 *
 * \param random Random-number generator.
 * \param system The system holding the framework, molecules and force field.
 * \param selectedComponent Index of the proton component.
 * \param selectedMolecule Index of the proton (integer molecule) to relocate; also selects the
 *        candidate group Component::tetheredProtonHopSiteGroups[selectedMolecule].
 * \return The energy difference on acceptance, std::nullopt otherwise.
 */
std::optional<RunningEnergy> tetheredProtonHopMove(RandomNumber& random, System& system,
                                                   std::size_t selectedComponent, std::size_t selectedMolecule);
}  // namespace MC_Moves
