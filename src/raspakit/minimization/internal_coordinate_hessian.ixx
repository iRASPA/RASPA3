module;

export module minimization_internal_coordinate_hessian;

import std;

import double3;

export namespace Minimization
{
/** 3x3 Cartesian second-derivative block; [row][column] with row = first atom's axis. */
using CoordinateBlock = std::array<std::array<double, 3>, 3>;

/**
 * Value and analytic Cartesian derivatives of a single internal coordinate
 * (distance, bend angle theta, or dihedral cos(phi)).
 *
 * gradient[slot] = dq/dx of the atom in that slot, hessian[slotI][slotJ] = d²q/dx_I dx_J.
 * Cross-term potentials U(q_1, ..., q_n) assemble their Hessian by the chain rule:
 *   d²U/dxdy = sum_k U_k d²q_k/dxdy + sum_{k,l} U_{kl} (dq_k/dx)(dq_l/dy).
 */
struct InternalCoordinate
{
  double value{};   // r [Å], theta [rad], or cos(phi)
  double phi{};     // signed dihedral angle (dihedral coordinate only)
  double sinPhi{};  // signed sin(phi), clamped away from zero (dihedral coordinate only)
  std::array<double3, 4> gradient{};
  std::array<std::array<CoordinateBlock, 4>, 4> hessian{};
  std::size_t numAtoms{};
};

/** Distance r = |posI - posJ|; slots {I, J}. */
InternalCoordinate distanceInternalCoordinate(const double3 &posI, const double3 &posJ);

/** Bend angle theta between (A-B) and (C-B); slots {A, B, C}. */
InternalCoordinate angleInternalCoordinate(const double3 &posA, const double3 &posB, const double3 &posC);

/** Dihedral cos(phi) of A-B-C-D (protein convention); slots {A, B, C, D}. */
InternalCoordinate dihedralInternalCoordinate(const double3 &posA, const double3 &posB, const double3 &posC,
                                              const double3 &posD);

/** Signed dihedral angle phi of A-B-C-D (protein convention); slots {A, B, C, D}. */
InternalCoordinate dihedralAngleInternalCoordinate(const double3 &posA, const double3 &posB, const double3 &posC,
                                                   const double3 &posD);
}  // namespace Minimization
