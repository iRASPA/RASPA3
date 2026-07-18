module;

export module cbmc_ring_closure;

import std;

import atom;
import double3;
import randomnumbers;
import component;
import intra_molecular_potentials;

export namespace CBMC
{
/**
 * \brief Whether a torsion changes when the ring is rigidly rotated ("spun") about the previous-anchor
 *        bond during the ring-closure orientational bias.
 *
 * Under a rigid rotation of the whole ring about the previous-anchor axis, the anchor 'currentBead'
 * lies on the axis and the ring atoms 'nextBeads' rotate as one body. A torsion transforms rigidly
 * (is spin-invariant) when all four of its atoms belong to the ring body ('currentBead' plus the
 * 'nextBeads'); it is spin-variant when it also involves an atom outside the ring body (the previous
 * bead or any other placed atom). Spin-variant torsions are the junction-crossing torsions that must
 * bias the spin through the coupled-decoupled torsion Rosenbluth step; the internal ring torsions are
 * spin-invariant and are sampled instead by the ring conformational Monte-Carlo.
 *
 * Note this differs from the bend case ('isSpinVariantBend'): a bend on the axis (previous-anchor-next)
 * is spin-invariant, but a torsion involving the previous bead is spin-variant.
 */
bool isSpinVariantRingTorsion(const std::array<std::size_t, 4> &identifiers, std::size_t currentBead,
                              const std::vector<std::size_t> &nextBeads);

/**
 * \brief Potentials restricted to the spin-selected (junction-crossing) terms of a ring-unit segment.
 *
 * Returns a copy of 'intra' whose torsions are only the spin-variant (junction-crossing) torsions and
 * whose bends are kept as-is (the torsion selection filters spin-variant bends itself). All other
 * terms are dropped. Passed to 'selectTorsionOrientation' so the ring spin is biased by exactly the
 * junction-crossing torsions and bends, while the internal ring torsions (sampled by the
 * conformational Monte-Carlo) are excluded to avoid double counting.
 */
Potentials::IntraMolecularPotentials ringSpinPotentials(const Potentials::IntraMolecularPotentials &intra,
                                                        std::size_t currentBead,
                                                        const std::vector<std::size_t> &nextBeads);

/**
 * \brief Samples a base conformation and orientation of a flexible ring hinged on the placed anchor.
 *
 * Ring-unit counterpart of 'generateRigidUnitOrientationMonteCarloScheme'. The ring atoms 'nextBeads'
 * form a closed cycle that bonds back to the anchor 'currentBead'. Because the internal conformation
 * of the ring is independent of the rest of the system, it is sampled from its Boltzmann distribution
 * by an internal Metropolis Monte-Carlo (small Cartesian displacements of the ring atoms accepted on
 * the internal ring bond/bend/torsion energy, which naturally keeps the ring closed through the
 * closure bond), combined with small rigid-body rotations of the whole ring about the anchor that
 * sample the junction bends. The initial conformation is the component reference ring geometry placed
 * with a random orientation, which is a valid closed ring.
 *
 * Like the flexible bond/bend Monte-Carlo and the rigid-unit tilt Monte-Carlo, this routine carries no
 * Rosenbluth weight: the internal energy is absorbed into the sampling and the junction-crossing
 * torsions/bends are weighted afterwards by the shared torsion step ('selectTorsionOrientation') using
 * the spin potentials from 'ringSpinPotentials'. The spin about the junction bond and the
 * junction-crossing terms are excluded here (they are the caller's torsion-selection responsibility).
 *
 * Without a previous bead (the ring is the growth seed) every orientation is equally likely and a
 * uniformly random overall orientation is returned.
 */
/**
 * \brief Rigidly rotates a ring conformation about the anchor with a uniformly random orientation.
 *
 * Used for the seed ring (no previous bead), where every overall orientation is equally likely: one
 * internal conformation is sampled once and the trial directions are cheap random rigid rotations of
 * it (rigid rotation preserves the internal ring geometry). 'anchorPosition' is the fixed anchor.
 */
std::vector<Atom> randomlyOrientRing(RandomNumber &random, const std::vector<Atom> &ringAtoms,
                                     double3 anchorPosition);

std::vector<Atom> generateRingConformationMonteCarloScheme(RandomNumber &random,
                                                           std::size_t numberOfTrialMovesPerOpenBead, double beta,
                                                           const Component &component,
                                                           const std::vector<Atom> &chainAtoms,
                                                           std::optional<std::size_t> previousBead,
                                                           std::size_t currentBead,
                                                           const std::vector<std::size_t> &nextBeads,
                                                           const Potentials::IntraMolecularPotentials &intra);
}  // namespace CBMC
