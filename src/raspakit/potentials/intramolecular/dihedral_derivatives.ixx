module;

export module dihedral_derivatives;

import std;

import double3;

export namespace Potentials
{
/**
 * \brief The signed dihedral angle phi of A-B-C-D (IUPAC sign, as AMBER, CHARMM and LAMMPS: positive when, looking
 *        from B to C, A turns clockwise onto D) with its exact Cartesian derivatives.
 *
 * A potential written in phi takes its gradient U'(phi) grad(phi) and Hessian U''(phi) grad(phi) grad(phi)^T +
 * U'(phi) Hess(phi) from these. Both are finite at every phi, cis and trans included, where the chain rule through
 * cos(phi) divides by sin(phi) = 0 (grad cos(phi) vanishes there): it zeroed the out-of-plane forces of a potential
 * with U'(phi) != 0 and lost the curvature U'' grad(phi) grad(phi)^T at a planar minimum. The derivatives are zero for
 * a degenerate geometry (three collinear atoms), where phi is undefined.
 */
struct DihedralDerivatives
{
  double phi{};
  /// d phi / d x of A, B, C and D (Blondel and Karplus, J. Comput. Chem. 17, 1132 (1996)).
  std::array<double3, 4> gradient{};
  /// d2 phi / dx_i dx_j, coordinates ordered A.x, A.y, A.z, B.x, ..., D.z; filled only when requested.
  std::array<std::array<double, 12>, 12> hessian{};
};

/// The signed (IUPAC) dihedral angle alone, accurate at cis and trans where acos(cos(phi)) loses half the digits.
double dihedralAngle(const double3& posA, const double3& posB, const double3& posC, const double3& posD);

DihedralDerivatives dihedralAngleDerivatives(const double3& posA, const double3& posB, const double3& posC,
                                             const double3& posD, bool computeHessian);
}  // namespace Potentials
