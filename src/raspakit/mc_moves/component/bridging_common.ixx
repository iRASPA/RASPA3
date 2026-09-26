module;

export module mc_moves_bridging_common;

import std;

import double3;
import atom;
import randomnumbers;
import component;
import mc_moves_concerted_rotation_geometry;

/**
 * Shared kernel of the connectivity-altering bridging moves (double bridging and intramolecular
 * double rebridging; Karayiannis, Mavrantzas and Theodorou, Phys. Rev. Lett. 88, 105503 (2002);
 * Karayiannis, Giannousaki, Mavrantzas and Theodorou, J. Chem. Phys. 117, 5465 (2002)).
 *
 * Both moves excise two trimers and re-bridge them between different anchor dimers, so that chain
 * parts change owner (double bridging) or direction (intramolecular double rebridging) while every
 * atom outside the trimers keeps its position. A closure is the fixed-endpoint problem of the
 * concerted-rotation move (mc_moves_concerted_rotation_geometry): with the anchors a1, a2 upstream
 * and a6, a7 downstream given, a trimer a3, a4, a5 has nine coordinates and nine constraints (four
 * bond lengths, five bend angles), hence a discrete set of solutions. The closure preserves the
 * internal geometry of the trimer being excised from the same site, so that the multiset of bond
 * lengths and bend angles of the system is unchanged and the reverse move (the same sites, the same
 * excisions) finds the old trimer among its solutions.
 *
 * Acceptance rule. Fixing all atoms outside a trimer fixes the six torsion-like coordinates that
 * place its downstream atoms; the constrained Boltzmann density on the discrete solution set is
 * therefore exp(-beta U) / J, with J the closure Jacobian of the concerted rotation (the derivative
 * of the chart 'a6, a7 perpendicular to a6a7, a8 out of the a6a7a8 plane' with respect to the
 * torsions about a1a2 ... a6a7). A solution is picked uniformly, so a move built from closures k
 * carries the factor prod_k (n_forward,k / n_backward,k) (J_old,k / J_new,k).
 */
export namespace Bridging
{
/// The anchors of a closure: the dimer (a1, a2) the trimer grows from, the dimer (a6, a7) it closes
/// onto, and the next fixed atom a8 beyond a7 (absent when a7 is a chain end or the next atom moves).
struct Anchors
{
  double3 a1{}, a2{}, a6{}, a7{};
  std::optional<double3> a8{};
};

struct Closure
{
  ConcertedRotation::Trimer trimer{};  ///< The new trimer.
  double logWeight{};                  ///< log(n_forward / n_backward) + log(J_old / J_new).
};

/// Re-bridges 'oldTrimer' (which sits between 'oldAnchors') between 'newAnchors' with its own internal
/// geometry. Draws one of the closure solutions uniformly. Nullopt when there is no solution, when
/// the backward closure (between the old anchors) does not reproduce the old trimer (a numerical
/// failure of the root search; rejecting keeps the move reversible), or when a closure Jacobian is
/// degenerate.
std::optional<Closure> rebridge(RandomNumber& random, const Anchors& oldAnchors,
                                const ConcertedRotation::Trimer& oldTrimer, const Anchors& newAnchors);

/// The volume element of the closure chart relative to the rigid-motion measure of the downstream
/// body: |a7-a6|^2 |a8-a7| sin(angle a6 a7 a8), or |a7-a6|^2 without a8. The ratio of these constants
/// enters the acceptance rule when the old and the new state chart the same fixed body from different
/// ends (the reversed segment of the intramolecular double rebridging).
double chartVolumeElement(const double3& a6, const double3& a7, const std::optional<double3>& a8);

/// Writes a re-bridged trimer into 'trialAtoms': the backbone atoms of units site+1 ... site+3 take
/// the positions of 'newTrimer', and the side atoms of each of these units (read from 'sourceAtoms')
/// are carried rigidly with the local frame (previous backbone atom, atom, next backbone atom) of
/// their backbone atom, which the closure maps congruently.
void placeTrimer(std::span<Atom> trialAtoms, std::span<const Atom> sourceAtoms,
                 const std::vector<Component::BridgingTopology::Unit>& units, std::size_t site,
                 const Anchors& oldAnchors, const ConcertedRotation::Trimer& oldTrimer, const Anchors& newAnchors,
                 const ConcertedRotation::Trimer& newTrimer);
}  // namespace Bridging
