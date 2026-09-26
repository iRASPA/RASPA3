module;

export module mc_moves_double_rebridging;

import std;

import randomnumbers;
import running_energy;
import system;

export namespace MC_Moves
{
/**
 * \brief Performs an intramolecular double rebridging (IDR) move on a single chain.
 *
 * Intramolecular double rebridging (Karayiannis, Giannousaki, Mavrantzas and Theodorou, J. Chem.
 * Phys. 117, 5465 (2002)) is the single-chain analogue of double bridging: two trimers are excised
 * from the same chain at sites a and b, and the segment between them is reversed by bridging the
 * head of the chain onto the far end of the segment and the near end of the segment onto the tail.
 * Chain length and connectivity are unchanged; the conformation of the whole segment changes at
 * once while every atom outside the two trimers keeps its position. The move is most useful for
 * long chains in dense systems, where intermolecular double bridging cannot find partners for
 * every chain and the chain ends are far apart.
 *
 * The molecule must be acyclic with a backbone of at least twelve units, and the reversed segment
 * must be congruent under reversal (unit k and unit a+4+b-k have the same types, charges and
 * side-group structure; a homopolymer backbone qualifies), see Component::BridgingTopology. Side
 * groups of the segment travel with their backbone atoms, side groups of the re-bridged trimers are
 * carried rigidly with the local frames of their backbone atoms.
 *
 * The site pair is chosen uniformly among the valid pairs, each closure picks one of its discrete
 * solutions uniformly. The acceptance rule is the Metropolis factor times the closure factors of
 * both bridges (see mc_moves_bridging_common) and the ratio of the chart volume elements of the
 * reversed segment, which the two states chart from opposite ends (unity for a chain with uniform
 * fixed bond lengths and bend angles).
 *
 * Polarization is not supported.
 *
 * \param random Random number generator instance.
 * \param system The simulation system.
 * \param selectedComponent Index of the component.
 * \param selectedMolecule Index of the chain within the component.
 *
 * \return Energy difference if the move is accepted; \c std::nullopt otherwise.
 */
std::optional<RunningEnergy> intramolecularDoubleRebridgingMove(RandomNumber& random, System& system,
                                                                std::size_t selectedComponent,
                                                                std::size_t selectedMolecule);
}  // namespace MC_Moves
