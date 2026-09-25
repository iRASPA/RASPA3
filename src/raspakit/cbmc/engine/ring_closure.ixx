module;

export module cbmc_ring_closure;

import std;

import atom;
import double3;
import randomnumbers;
import cbmc_grow_context;
import component;
import cbmc_grow_step;

export namespace CBMC
{
/**
 * \brief Samples the internal conformation of a cyclic cluster (a CloseRing step) plus its junction
 * tilt.
 *
 * The cluster is seeded from the component's conformation reservoir (or its reference geometry),
 * rigidly oriented about the anchor with the parity of every determined chiral center preserved, and
 * relaxed by an internal Metropolis Monte-Carlo on the step's ring-internal bonded terms: per-atom
 * displacements and constraint-preserving rotations for flexible ring atoms, whole-unit moves for
 * rigid sub-fragments, whole-ring tilts about the anchor, and conformer-hopping crankshaft rotations.
 * Fixed bonds are kept exactly. Step sizes are adaptive (the anchor bead's CBMC statistics). Carries
 * no Rosenbluth weight: the spin about the junction bond and the junction-crossing terms are weighted
 * afterwards in the torsion selection.
 *
 * Returns the positions of the step's next-beads (in 'step.nextBeads' order). 'chainAtoms' is the chain
 * the step grows in; its next-beads are used as scratch and restored on return (see CBMC::ScratchBeads).
 */
std::vector<Atom> generateRingConformation(RandomNumber &random, const GrowthSettings &settings, double beta,
                                           const Component &component, std::vector<Atom> &chainAtoms,
                                           const GrowStep &step);

/// A uniformly random rigid rotation of 'ringAtoms' about 'anchorPosition' (the other trial
/// directions of a ring seed, which has no orientational reference).
std::vector<Atom> randomlyOrientRing(RandomNumber &random, const std::vector<Atom> &ringAtoms,
                                     double3 anchorPosition);
}  // namespace CBMC
