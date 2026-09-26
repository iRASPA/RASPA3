module;

export module mc_moves_concerted_rotation_geometry;

import std;

import double3;
import double3x3;

/**
 * Geometry kernel of the concerted-rotation (ConRot) move of Dodd, Boone and Theodorou,
 * Mol. Phys. 78, 961 (1993).
 *
 * The window is a path of backbone atoms a0-a1-a2-a3-a4-a5-a6-a7(-a8). The driver rotates a2 about
 * the a0-a1 axis; the trimer a3, a4, a5 is then re-bridged between the displaced a2 and the fixed
 * a6, a7 such that the six bond lengths a1a2 ... a6a7 and the seven backbone bend angles at
 * a1 ... a7 keep their values. With a2 and a6, a7 given, the trimer has nine coordinates and nine
 * constraints (bonds a2a3, a3a4, a4a5, a5a6; bends at a2, a3, a4, a5, a6): a discrete set of
 * solutions. It is found as in the original paper, by reducing the closure to a one-dimensional
 * equation in the torsion phi1 that places a3 on its cone about a1-a2: for a given a3 the distance
 * |a3 a5| is fixed by the bend at a4, which puts a5 on at most two points of its circle about the
 * a6-a7 axis; a4 then lies on at most two points of the circle around the a3-a5 axis that satisfy
 * the bend at a3; the bend at a5 is the remaining equation. Its roots are located by a grid scan
 * of the four branches: sign changes are refined by regula falsi, the branch boundaries (where
 * the two solutions of a quadratic merge) are refined so that roots next to a boundary are not
 * lost, and every local minimum of |residual| without a sign change is probed by parabolic
 * interpolation to catch pairs of roots between grid points. With the default one-degree grid,
 * about one root in a thousand is still missed; the move only loses efficiency by that, because
 * the same deterministic search is applied to the forward and the reverse closure.
 *
 * The acceptance rule needs the Jacobian of the closure: the determinant of the derivatives of the
 * fixed downstream positions (a6, a7, a8) with respect to the torsion angles about the bonds a1a2,
 * a2a3, a3a4, a4a5, a5a6 (and a6a7 when a8 exists), in a chart of the six-dimensional manifold of
 * rigid placements of (a6, a7, a8): a6 itself, the two components of a7 perpendicular to the a6a7
 * bond, and the component of a8 along the normal of the a6a7a8 plane. The chart is the same
 * function of the (identical) downstream positions in the old and the new state, so any
 * non-degenerate chart gives the same ratio J_old / J_new.
 */
export namespace ConcertedRotation
{
/// The invariant internal geometry of the window: bond lengths a2a3 ... a5a6 and the cosines of the
/// bends at a2 ... a6 (a1a2a3 ... a5a6a7).
struct BackboneGeometry
{
  double l23{}, l34{}, l45{}, l56{};
  double cos2{}, cos3{}, cos4{}, cos5{}, cos6{};

  /// Measures the geometry from the current positions of a1 ... a7.
  static BackboneGeometry fromPositions(std::span<const double3, 7> a1_to_a7);
};

/// One re-bridged trimer.
struct Trimer
{
  double3 a3{}, a4{}, a5{};
};

/// The default number of grid points of the scan over phi1 (one degree).
inline constexpr std::size_t defaultGridPoints = 360;

/// All trimers (a3, a4, a5) that close the window between the given a1, a2 and a6, a7 with the
/// given internal geometry. Solutions closer than 'duplicateTolerance' (Angstrom) are merged.
std::vector<Trimer> rebridge(const double3 &a1, const double3 &a2, const double3 &a6, const double3 &a7,
                             const BackboneGeometry &geometry, std::size_t gridPoints = defaultGridPoints,
                             double duplicateTolerance = 1e-7);

/// The closure Jacobian |det dP/dphi| (see above) for the backbone a1 ... a7 and the optional a8.
/// Zero for degenerate geometry (collinear a6, a7, a8).
double closureJacobian(std::span<const double3, 7> a1_to_a7, const std::optional<double3> &a8);

/// Rotates 'point' about the axis through 'origin' with unit direction 'axis' by 'angle'.
double3 rotateAboutAxis(const double3 &origin, const double3 &axis, double angle, const double3 &point);

/// The orthonormal local frame of a backbone atom: columns e1 along (previous - atom), e2 the
/// perpendicular part of (next - atom), e3 = e1 x e2.
double3x3 localFrame(const double3 &atom, const double3 &previous, const double3 &next);

/// Maps 'point', given in the old frame at 'oldOrigin', rigidly to the new frame at 'newOrigin'.
double3 transformRigidly(const double3 &oldOrigin, const double3x3 &oldFrame, const double3 &newOrigin,
                         const double3x3 &newFrame, const double3 &point);

/// Determinant of a dense n x n matrix in row-major order (Gaussian elimination with pivoting).
double determinant(std::vector<double> matrix, std::size_t n);
}  // namespace ConcertedRotation
