module;

export module fragment;

import std;

import archive;
import double3;
import double3x3;
import simd_quatd;
import atom;
import molecule;

/**
 * \brief A maximal rigid body inside a molecule.
 *
 * A fragment is the unit the constructive-bias Monte-Carlo (CBMC) growth places at once and the unit
 * molecular dynamics and minimization treat as a rigid body. It is the successor of RASPA2's rigid
 * 'MoleculeGroup': every atom of a molecule belongs to exactly one fragment.
 *
 *  - A single-atom fragment ('atoms.size() == 1') is a flexible bead: it carries no orientation and
 *    is grown/integrated as a point mass.
 *  - A multi-atom fragment ('atoms.size() > 1') is a rigid body: its atoms keep the body-fixed
 *    geometry taken from the molecule's reference atoms, the authoritative dynamical state is a
 *    center of mass plus an orientation quaternion (see 'GroupState'), and the atom positions are
 *    regenerated as 'centerOfMass + q * bodyFixedPositions[k]'. The body frame is the fragment's
 *    principal-axis frame, so the inertia tensor is diagonal ('inertiaVector').
 *
 * Flexible rings (simple, fused, or bridged) are NOT a fragment property: a ring is a set of
 * single-atom fragments connected in a cycle in the fragment graph, closed with ring-closure CBMC.
 */
export struct Fragment
{
  std::uint64_t versionNumber{1};  ///< Version number for serialization.

  std::vector<std::size_t> atoms{};  ///< Indices (into the molecule's atoms) that make up this fragment.
  std::vector<double> atomMasses{};  ///< Mass of every fragment atom, in the same order as 'atoms'.

  double mass{0.0};                           ///< Total mass of the fragment.
  double3 centerOfMassReferencePosition{};    ///< Fragment center of mass in the molecule reference frame.
  double3 inertiaVector{};                    ///< Principal moments of inertia.
  double3 inverseInertiaVector{};             ///< 1/I with (near-)zero principal axes set to zero.
  std::vector<double3> bodyFixedPositions{};  ///< Per fragment-atom offset in the fragment principal-axis frame.
  std::size_t rotationalDegreesOfFreedom{0};  ///< Rotational DOF of the fragment (3 nonlinear, 2 linear, 0 point).
  std::size_t shapeType{2};                   ///< 0 = nonlinear, 1 = linear, 2 = point (single atom).

  Fragment() = default;
  explicit Fragment(std::vector<std::size_t> atoms) : atoms(std::move(atoms)) {}

  /// A fragment representing a rigid body (more than one atom, sharing a single orientation).
  bool isRigidBody() const { return atoms.size() > 1; }

  /**
   * \brief Computes the rigid-body reference data (mass, principal-axis body frame, inertia,
   *        rotational degrees of freedom, shape) from the molecule reference geometry.
   *
   * \param referencePositions Reference position of every atom of the molecule (indexed by atom id).
   * \param masses             Mass of every atom of the molecule (indexed by atom id).
   *
   * A single-atom fragment is left as a point mass (no orientation). A multi-atom fragment gets its
   * inertia tensor diagonalized so molecular dynamics can integrate it as a free rigid rotor and the
   * body-fixed offsets are expressed in that principal-axis frame.
   */
  void computeRigidProperties(std::span<const double3> referencePositions, std::span<const double> masses);

  /**
   * \brief Regenerates the laboratory positions of the atoms of this fragment from its rigid-body
   *        state, writing into 'moleculeAtoms' (local molecule atom order).
   */
  void regenerateAtoms(const GroupState &state, std::span<Atom> moleculeAtoms) const;

  /**
   * \brief Recovers the rigid-body state (center of mass and orientation) of this fragment from the
   *        current laboratory atom positions (orthogonal Procrustes fit; exact for rigid geometry).
   */
  GroupState deriveState(std::span<const Atom> moleculeAtoms) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const Fragment &f);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, Fragment &f);
};
