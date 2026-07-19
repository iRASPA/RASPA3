module;

export module minimization_hessian_scatter;

import std;

import double3;
import generalized_hessian;
import minimization_dof_layout;
import minimization_bend_hessian_geometry;
import minimization_torsion_hessian_geometry;
import minimization_rigid_kinematics;

export namespace Minimization
{
/**
 * Scatter a Cartesian atom-atom Hessian block into the generalized position-position Hessian.
 *
 * Uses RASPA2 conventions with f1 = (1/r) dU/dr and f2 = d²U/dr² terms assembled as
 * f2 * dr_i dr_j + f1 * delta_ij on diagonal blocks and the negative cross block.
 */
void scatterAtomicPositionPosition(GeneralizedHessian& hessian, const MinimizationDofLayout& layout,
                                   std::size_t moleculeA, std::size_t localAtomA, std::size_t moleculeB,
                                   std::size_t localAtomB, double f1, double f2, const double3& dr);
void scatterAtomicPositionPositionByDof(GeneralizedHessian& hessian, std::size_t baseA, std::size_t baseB, double f1,
                                        double f2, const double3& dr);

/** Scatter intermolecular / framework pair Hessian blocks for flexible and rigid molecules. */
void scatterInteractionHessian(GeneralizedHessian& hessian, const MinimizationDofLayout& layout,
                               const RigidDerivativeCache& rigidCache, std::size_t moleculeA, std::size_t localAtomA,
                               bool rigidA, const double3& posA, const double3& comA, std::size_t moleculeB,
                               std::size_t localAtomB, bool rigidB, const double3& posB, const double3& comB, double f1,
                               double f2, const double3& dr);

/**
 * Scatter a framework-molecule pair; the framework atom carries no degrees of freedom, so only the
 * molecule-side diagonal blocks are written. Convention: dr = posMolecule - posFramework.
 */
void scatterFrameworkMoleculeHessian(GeneralizedHessian& hessian, const MinimizationDofLayout& layout,
                                     const RigidDerivativeCache& rigidCache, std::size_t moleculeIndex,
                                     std::size_t localAtom, bool rigid, double f1, double f2, const double3& dr);
void scatterFlexibleFrameworkMoleculeHessian(GeneralizedHessian& hessian, const MinimizationDofLayout& layout,
                                             const RigidDerivativeCache& rigidCache, std::size_t frameworkAtom,
                                             std::size_t moleculeIndex, std::size_t localAtom, bool moleculeRigid,
                                             double f1, double f2, const double3& dr);

/** Scatter a harmonic (and related) bend Hessian for a flexible triplet of atoms. */
void scatterBendHessian(GeneralizedHessian& hessian, const MinimizationDofLayout& layout, std::size_t moleculeIndex,
                        std::size_t localAtomA, std::size_t localAtomB, std::size_t localAtomC,
                        const BendHessianGeometry& geometry);
void scatterBendHessianByDof(GeneralizedHessian& hessian, const std::array<std::size_t, 3>& bases,
                             const BendHessianGeometry& geometry);

/** Scatter a dihedral (torsion / improper torsion) Hessian for a flexible quadruplet of atoms. */
void scatterTorsionHessian(GeneralizedHessian& hessian, const MinimizationDofLayout& layout, std::size_t moleculeIndex,
                           std::size_t localAtomA, std::size_t localAtomB, std::size_t localAtomC,
                           std::size_t localAtomD, const TorsionHessianGeometry& geometry);
void scatterTorsionHessianByDof(GeneralizedHessian& hessian, const std::array<std::size_t, 4>& bases,
                                const TorsionHessianGeometry& geometry);

/**
 * Generalized-DOF descriptor for one atom of a multi-atom intramolecular term. Flexible atoms
 * carry three Cartesian DOFs at 'positionBase'; atoms driven by a rigid body carry the body's
 * center-of-mass DOFs at 'positionBase' plus orientation DOFs at 'orientationBase' with the
 * per-atom orientation Jacobian in 'derivatives'.
 */
struct HessianSite
{
  std::optional<std::size_t> positionBase{};
  std::optional<std::size_t> orientationBase{};
  const RigidAtomDerivatives* derivatives{nullptr};
};

/** Build the site descriptor for atom 'localAtom' of molecule 'moleculeIndex'. */
HessianSite makeHessianSite(const MinimizationDofLayout& layout, const RigidDerivativeCache& rigidCache,
                            std::size_t moleculeIndex, std::size_t localAtom);

/**
 * Build the site descriptor for a (global) framework atom of a mixed framework: flexible atoms
 * carry three Cartesian DOFs, atoms driven by a rigid group carry the group's center-of-mass and
 * orientation DOFs, and fixed atoms carry none (empty site: all contributions drop out).
 */
HessianSite makeFrameworkHessianSite(const MinimizationDofLayout& layout, const RigidDerivativeCache& rigidCache,
                                     std::size_t frameworkAtom);

/**
 * Scatter a radial (two-atom) Hessian where either atom may be driven by a rigid body. 'dr' is
 * posA - posB, and the RASPA2 factors are f1 = (1/r) dU/dr and f2 = d2U/dr2 assembled as
 * f2 dr dr^T + f1 I on the diagonal blocks and the negative cross block; 'gradients' are the
 * Cartesian energy gradients per atom (gradientA = f1 dr, gradientB = -f1 dr).
 */
void scatterRadialHessianSites(GeneralizedHessian& hessian, const std::array<HessianSite, 2>& sites,
                               const std::array<double3, 2>& gradients, double f1, double f2, const double3& dr);

/**
 * Scatter a bend (three-atom) Hessian where any of the atoms may be driven by a rigid body:
 * the Cartesian blocks are projected as J^T H J, plus the gradient-curvature term
 * sum_i g_i . ddVec_i on the orientation blocks. The gradients are the Cartesian energy
 * gradients per atom.
 */
void scatterBendHessianSites(GeneralizedHessian& hessian, const std::array<HessianSite, 3>& sites,
                             const std::array<double3, 3>& gradients, const BendHessianGeometry& geometry);

/** Group-aware torsion (four-atom) analog of scatterBendHessianSites. */
void scatterTorsionHessianSites(GeneralizedHessian& hessian, const std::array<HessianSite, 4>& sites,
                                const std::array<double3, 4>& gradients, const TorsionHessianGeometry& geometry);

/**
 * Group-aware scatter of a dense Cartesian per-term Hessian for a many-body intramolecular term
 * (higher-order cross terms: bond-bond, bond-bend, bend-bend, bond/bend-torsion, inversion and
 * out-of-plane bends). 'localCartesian' is the dense 3M x 3M atom-major Cartesian block of the
 * term, 'sites' describe each atom's generalized degrees of freedom (flexible Cartesian, or a
 * rigid body's center-of-mass plus orientation with the per-atom orientation Jacobian), and
 * 'gradients' are the per-atom Cartesian energy gradients used for the orientation gradient
 * curvature. The Cartesian block is projected as J^T H J plus sum_i g_i . d2p_i on rigid
 * orientation blocks.
 */
void scatterCartesianTermSites(GeneralizedHessian& hessian, std::span<const HessianSite> sites,
                               std::span<const double3> gradients, const GeneralizedHessian& localCartesian);

/**
 * Molecular-strain-convention correction of the strain gradient for the rigid-group atoms of a
 * many-body term: subtracts the non-scaling internal-offset virial g_i (pos_i - bodyCoM_i) of
 * every atom driven by a rigid body (whole-molecule-rigid atoms carry no generalized DOFs and are
 * skipped by the caller). 'localAtoms', 'positions' and 'gradients' are the term's per-atom data.
 */
void removeRigidOffsetStrainGradientSites(GeneralizedHessian& hessian, const MinimizationDofLayout& layout,
                                          const RigidDerivativeCache& rigidCache, std::size_t moleculeIndex,
                                          std::span<const std::size_t> localAtoms, std::span<const double3> positions,
                                          std::span<const double3> gradients);

/** Isotropic strain coupling for a single strain degree of freedom (NPT, nStrain = 1). */
void scatterAtomicPositionStrainIsotropic(GeneralizedHessian& hessian, const MinimizationDofLayout& layout,
                                          std::size_t moleculeA, std::size_t localAtomA, std::size_t moleculeB,
                                          std::size_t localAtomB, double f1, double f2, const double3& dr);

/**
 * One-sided isotropic position-strain scatter for one site of an interaction pair, valid for
 * flexible atoms and rigid molecules (center-of-mass and orientation degrees of freedom).
 *
 * Under the exp(epsilon) strain convention dr(epsilon) = exp(epsilon) c + delta, where
 * c = dr - (dA - dB) with d = pos - com for rigid sites (zero for flexible and framework sites),
 * the mixed second derivatives of a radial pair term are
 *
 *   d2U/d(com_a) dEps   = sign * [ f1 (dr_a + c_a) + f2 (dr.c) dr_a ]
 *   d2U/d(omega_x) dEps = sign * [ f2 (dr.c) (dr.DVec_x) + f1 (c.DVec_x) ]
 *
 * with sign = +1 for the site whose gradient is +f1 dr (atom A of dr = posA - posB) and -1 for
 * the partner. For a flexible-flexible pair (c = dr) this reduces to the classic
 * f2 r^2 dr_a + 2 f1 dr_a form.
 */
void scatterSitePositionStrainIsotropic(GeneralizedHessian& hessian, const MinimizationDofLayout& layout,
                                        const RigidDerivativeCache& rigidCache, std::size_t moleculeIndex,
                                        std::size_t localAtom, bool rigid, double sign, double f1, double f2,
                                        const double3& dr, const double3& drStrainDerivative);

/** Isotropic strain-strain contribution for a single strain degree of freedom. */
void scatterAtomicStrainStrainIsotropic(GeneralizedHessian& hessian, double f1, double f2, const double3& dr,
                                        const double3& posA, const double3& comA, const double3& posB,
                                        const double3& comB, bool rigidA, bool rigidB);

}  // namespace Minimization
