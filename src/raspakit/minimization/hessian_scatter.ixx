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
void scatterAtomicPositionPosition(GeneralizedHessian &hessian, const MinimizationDofLayout &layout,
                                   std::size_t moleculeA, std::size_t localAtomA, std::size_t moleculeB,
                                   std::size_t localAtomB, double f1, double f2, const double3 &dr);

/** Scatter intermolecular / framework pair Hessian blocks for flexible and rigid molecules. */
void scatterInteractionHessian(GeneralizedHessian &hessian, const MinimizationDofLayout &layout,
                               const RigidDerivativeCache &rigidCache, std::size_t moleculeA,
                               std::size_t localAtomA, bool rigidA, const double3 &posA, const double3 &comA,
                               std::size_t moleculeB, std::size_t localAtomB, bool rigidB, const double3 &posB,
                               const double3 &comB, double f1, double f2, const double3 &dr);

/**
 * Scatter a framework-molecule pair; the framework atom carries no degrees of freedom, so only the
 * molecule-side diagonal blocks are written. Convention: dr = posMolecule - posFramework.
 */
void scatterFrameworkMoleculeHessian(GeneralizedHessian &hessian, const MinimizationDofLayout &layout,
                                     const RigidDerivativeCache &rigidCache, std::size_t moleculeIndex,
                                     std::size_t localAtom, bool rigid, double f1, double f2, const double3 &dr);

/** Scatter a harmonic (and related) bend Hessian for a flexible triplet of atoms. */
void scatterBendHessian(GeneralizedHessian &hessian, const MinimizationDofLayout &layout, std::size_t moleculeIndex,
                        std::size_t localAtomA, std::size_t localAtomB, std::size_t localAtomC,
                        const BendHessianGeometry &geometry);

/** Scatter a dihedral (torsion / improper torsion) Hessian for a flexible quadruplet of atoms. */
void scatterTorsionHessian(GeneralizedHessian &hessian, const MinimizationDofLayout &layout, std::size_t moleculeIndex,
                           std::size_t localAtomA, std::size_t localAtomB, std::size_t localAtomC,
                           std::size_t localAtomD, const TorsionHessianGeometry &geometry);

/** Isotropic strain coupling for a single strain degree of freedom (NPT, nStrain = 1). */
void scatterAtomicPositionStrainIsotropic(GeneralizedHessian &hessian, const MinimizationDofLayout &layout,
                                          std::size_t moleculeA, std::size_t localAtomA, std::size_t moleculeB,
                                          std::size_t localAtomB, double f1, double f2, const double3 &dr);

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
void scatterSitePositionStrainIsotropic(GeneralizedHessian &hessian, const MinimizationDofLayout &layout,
                                        const RigidDerivativeCache &rigidCache, std::size_t moleculeIndex,
                                        std::size_t localAtom, bool rigid, double sign, double f1, double f2,
                                        const double3 &dr, const double3 &drStrainDerivative);

/** Isotropic strain-strain contribution for a single strain degree of freedom. */
void scatterAtomicStrainStrainIsotropic(GeneralizedHessian &hessian, double f1, double f2, const double3 &dr,
                                        const double3 &posA, const double3 &comA, const double3 &posB,
                                        const double3 &comB, bool rigidA, bool rigidB);

}  // namespace Minimization
