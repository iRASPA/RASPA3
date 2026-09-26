module;

module mc_moves_concerted_rotation_geometry;

import std;

import double3;
import double3x3;

namespace
{

// An orthonormal basis with 'u' along 'direction'.
struct Basis
{
  double3 u{}, n1{}, n2{};
};

Basis basisAlong(const double3 &direction)
{
  double3 u = double3(direction).normalized();
  double3 helper = std::abs(u.x) < 0.9 ? double3(1.0, 0.0, 0.0) : double3(0.0, 1.0, 0.0);
  double3 n1 = double3::cross(u, helper).normalized();
  double3 n2 = double3::cross(u, n1);
  return {u, n1, n2};
}

// The solutions x = center +/- delta of A cos x + B sin x = C, and the margin R - |C| (negative
// when there is no solution). The degenerate case A = B = 0 (every x or none) is reported as none.
struct CircleSolution
{
  bool exists{false};
  double margin{};
  double center{};
  double delta{};
};

CircleSolution solveCircle(double A, double B, double C)
{
  double R = std::hypot(A, B);
  double margin = R - std::abs(C);
  if (!(R > 0.0) || margin < 0.0) return {false, margin, 0.0, 0.0};
  return {true, margin, std::atan2(B, A), std::acos(std::clamp(C / R, -1.0, 1.0))};
}

// One of the four solution branches (sign of the a5 root x sign of the a4 root) at a given phi1.
struct Branch
{
  bool exists{false};
  double margin{};    // min of the two circle margins: > 0 inside the branch's domain
  double residual{};  // cos(bend at a5) - cos(theta5); zero at a solution
  ConcertedRotation::Trimer trimer{};
};

class Rebridger
{
 public:
  Rebridger(const double3 &a1, const double3 &a2, const double3 &a6, const double3 &a7,
            const ConcertedRotation::BackboneGeometry &geometry)
      : a2(a2), a6(a6), g(geometry), cone2(basisAlong(a2 - a1)), cone6(basisAlong(a7 - a6))
  {
    sin2 = std::sqrt(std::max(0.0, 1.0 - g.cos2 * g.cos2));
    sin6 = std::sqrt(std::max(0.0, 1.0 - g.cos6 * g.cos6));
    // |a3 a5| is fixed by the bend at a4.
    d35Squared = g.l34 * g.l34 + g.l45 * g.l45 - 2.0 * g.l34 * g.l45 * g.cos4;
  }

  // All four branches at phi1.
  std::array<Branch, 4> evaluate(double phi1) const { return evaluate(phi1, std::nullopt); }

  // A single branch at phi1 (the other branches of the returned array are left empty).
  Branch evaluate(double phi1, std::size_t branch) const { return evaluate(phi1, std::optional(branch))[branch]; }

 private:
  std::array<Branch, 4> evaluate(double phi1, std::optional<std::size_t> only) const
  {
    std::array<Branch, 4> branches{};

    // a3 on its cone about the a1-a2 axis: the angle a1-a2-a3 is theta2.
    const double3 a3 =
        a2 + g.l23 * (-g.cos2 * cone2.u + sin2 * (std::cos(phi1) * cone2.n1 + std::sin(phi1) * cone2.n2));

    // a5 on its circle about the a6-a7 axis (angle a5-a6-a7 = theta6) at distance |a3 a5| from a3.
    const double3 v = a6 + g.l56 * g.cos6 * cone6.u - a3;
    const double A = 2.0 * g.l56 * sin6 * double3::dot(v, cone6.n1);
    const double B = 2.0 * g.l56 * sin6 * double3::dot(v, cone6.n2);
    const double C = d35Squared - double3::dot(v, v) - g.l56 * g.l56 * sin6 * sin6;
    const CircleSolution psi = solveCircle(A, B, C);
    if (!psi.exists)
    {
      for (Branch &branch : branches) branch.margin = psi.margin;
      return branches;
    }

    for (std::size_t signPsi = 0; signPsi != 2; ++signPsi)
    {
      if (only.has_value() && *only / 2 != signPsi) continue;
      const double angle = psi.center + (signPsi == 0 ? psi.delta : -psi.delta);
      const double3 a5 =
          a6 + g.l56 * (g.cos6 * cone6.u + sin6 * (std::cos(angle) * cone6.n1 + std::sin(angle) * cone6.n2));

      // a4 on the circle around the a3-a5 axis (|a3 a4| = l34, |a4 a5| = l45) at the bend theta3.
      const double3 w35 = a5 - a3;
      const double d = w35.length();
      const Basis axis35 = basisAlong(w35);
      const double t = (g.l34 * g.l34 - g.l45 * g.l45 + d * d) / (2.0 * d);
      const double rho = std::sqrt(std::max(0.0, g.l34 * g.l34 - t * t));
      const double3 e = a2 - a3;
      const double Ap = rho * double3::dot(axis35.n1, e);
      const double Bp = rho * double3::dot(axis35.n2, e);
      const double Cp = g.l34 * g.l23 * g.cos3 - t * double3::dot(axis35.u, e);
      const CircleSolution chi = solveCircle(Ap, Bp, Cp);

      for (std::size_t signChi = 0; signChi != 2; ++signChi)
      {
        if (only.has_value() && *only % 2 != signChi) continue;
        Branch &branch = branches[2 * signPsi + signChi];
        branch.margin = std::min(psi.margin, chi.margin);
        if (!chi.exists) continue;

        const double chiAngle = chi.center + (signChi == 0 ? chi.delta : -chi.delta);
        const double3 a4 = a3 + t * axis35.u + rho * (std::cos(chiAngle) * axis35.n1 + std::sin(chiAngle) * axis35.n2);

        branch.exists = true;
        branch.residual = double3::dot(a4 - a5, a6 - a5) / (g.l45 * g.l56) - g.cos5;
        branch.trimer = {a3, a4, a5};
      }
    }
    return branches;
  }

  double3 a2, a6;
  ConcertedRotation::BackboneGeometry g;
  Basis cone2, cone6;
  double sin2{}, sin6{}, d35Squared{};
};

bool signChange(const Branch &a, const Branch &b) { return a.residual * b.residual < 0.0; }

// The boundary of branch 'b' between an outside point and an inside point, located as the zero of
// the (continuous) existence margin by regula falsi (Illinois variant) with a bisection fallback.
// Returns the innermost located point at which the branch exists.
double branchBoundary(const Rebridger &solver, std::size_t b, double outside, double inside)
{
  double fOutside = solver.evaluate(outside, b).margin;
  double fInside = solver.evaluate(inside, b).margin;
  if (!(fOutside < 0.0) || !(fInside >= 0.0)) return inside;
  int side = 0;
  for (std::size_t iteration = 0; iteration != 40; ++iteration)
  {
    double trial = (fInside - fOutside) != 0.0 ? inside - fInside * (inside - outside) / (fInside - fOutside)
                                               : 0.5 * (outside + inside);
    if (!(std::min(outside, inside) < trial && trial < std::max(outside, inside))) trial = 0.5 * (outside + inside);
    const Branch atTrial = solver.evaluate(trial, b);
    if (atTrial.exists)
    {
      inside = trial;
      fInside = atTrial.margin;
      if (side == 1) fOutside *= 0.5;
      side = 1;
    }
    else
    {
      outside = trial;
      fOutside = atTrial.margin;
      if (side == -1) fInside *= 0.5;
      side = -1;
    }
    if (std::abs(inside - outside) < 1e-12) break;
  }
  return inside;
}

// The root of the residual of branch 'b' between two inside points with opposite residual signs:
// regula falsi (Illinois variant) with a bisection fallback; the bracket always stays on the branch.
std::optional<ConcertedRotation::Trimer> bisectRoot(const Rebridger &solver, std::size_t b, double lo, Branch atLo,
                                                    double hi, Branch atHi)
{
  int side = 0;
  for (std::size_t iteration = 0; iteration != 60; ++iteration)
  {
    double trial = (atHi.residual - atLo.residual) != 0.0
                       ? hi - atHi.residual * (hi - lo) / (atHi.residual - atLo.residual)
                       : 0.5 * (lo + hi);
    if (!(lo < trial && trial < hi) || iteration % 4 == 3) trial = 0.5 * (lo + hi);  // guard against stalling
    const Branch atTrial = solver.evaluate(trial, b);
    if (!atTrial.exists) return std::nullopt;  // the branch vanished in between: no root on this interval
    if (atTrial.residual == 0.0) return atTrial.trimer;
    if (signChange(atLo, atTrial))
    {
      hi = trial;
      atHi = atTrial;
      if (side == -1) atLo.residual *= 0.5;
      side = -1;
    }
    else
    {
      lo = trial;
      atLo = atTrial;
      if (side == 1) atHi.residual *= 0.5;
      side = 1;
    }
    if (hi - lo < 1e-13) break;
  }
  return std::abs(atLo.residual) < std::abs(atHi.residual) ? atLo.trimer : atHi.trimer;
}

// A pair of roots between grid points shows up as a local minimum of |residual| without a sign
// change. Starting from the three grid points around such a minimum, the minimum of the residual
// is located by successive parabolic interpolation; as soon as an evaluated point has the opposite
// sign the pair is bracketed and both roots are bisected.
template <typename AddRoot>
void refineLocalMinimum(const Rebridger &solver, std::size_t b, std::array<double, 3> phi, std::array<Branch, 3> at,
                        AddRoot &&addRoot)
{
  const double sign = at[1].residual > 0.0 ? 1.0 : -1.0;
  for (std::size_t iteration = 0; iteration != 12; ++iteration)
  {
    // Parabola through the three points; the vertex is the next trial point.
    const double x0 = phi[0], x1 = phi[1], x2 = phi[2];
    const double f0 = at[0].residual, f1 = at[1].residual, f2 = at[2].residual;
    const double denominator = (x1 - x0) * (f1 - f2) - (x1 - x2) * (f1 - f0);
    if (denominator == 0.0) return;
    const double vertex =
        x1 - 0.5 * ((x1 - x0) * (x1 - x0) * (f1 - f2) - (x1 - x2) * (x1 - x2) * (f1 - f0)) / denominator;
    if (!(vertex > x0 && vertex < x2) || vertex == x1) return;

    const Branch atVertex = solver.evaluate(vertex, b);
    if (!atVertex.exists) return;
    if (atVertex.residual * sign <= 0.0)
    {
      // Two roots: one on each side of the vertex.
      const bool left = vertex < x1;
      const double lo = left ? x0 : x1, hi = left ? x1 : x2;
      const Branch &atLo = left ? at[0] : at[1];
      const Branch &atHi = left ? at[1] : at[2];
      if (atVertex.residual == 0.0)
      {
        addRoot(atVertex.trimer);
        return;
      }
      if (std::optional<ConcertedRotation::Trimer> root = bisectRoot(solver, b, lo, atLo, vertex, atVertex))
        addRoot(*root);
      if (std::optional<ConcertedRotation::Trimer> root = bisectRoot(solver, b, vertex, atVertex, hi, atHi))
        addRoot(*root);
      return;
    }

    // Keep the three points with the smallest |residual| that still bracket the minimum.
    if (std::abs(atVertex.residual) >= std::abs(f1))
    {
      // The vertex is not an improvement: shrink the bracket toward the current best point.
      if (vertex < x1)
      {
        phi[0] = vertex;
        at[0] = atVertex;
      }
      else
      {
        phi[2] = vertex;
        at[2] = atVertex;
      }
    }
    else if (vertex < x1)
    {
      phi[2] = x1;
      at[2] = at[1];
      phi[1] = vertex;
      at[1] = atVertex;
    }
    else
    {
      phi[0] = x1;
      at[0] = at[1];
      phi[1] = vertex;
      at[1] = atVertex;
    }
    if (phi[2] - phi[0] < 1e-12) return;
  }
}

}  // namespace

ConcertedRotation::BackboneGeometry ConcertedRotation::BackboneGeometry::fromPositions(
    std::span<const double3, 7> a)
{
  auto cosAngle = [](const double3 &left, const double3 &center, const double3 &right)
  { return std::clamp(double3::dot((left - center).normalized(), (right - center).normalized()), -1.0, 1.0); };
  BackboneGeometry g{};
  g.l23 = (a[2] - a[1]).length();
  g.l34 = (a[3] - a[2]).length();
  g.l45 = (a[4] - a[3]).length();
  g.l56 = (a[5] - a[4]).length();
  g.cos2 = cosAngle(a[0], a[1], a[2]);
  g.cos3 = cosAngle(a[1], a[2], a[3]);
  g.cos4 = cosAngle(a[2], a[3], a[4]);
  g.cos5 = cosAngle(a[3], a[4], a[5]);
  g.cos6 = cosAngle(a[4], a[5], a[6]);
  return g;
}

std::vector<ConcertedRotation::Trimer> ConcertedRotation::rebridge(const double3 &a1, const double3 &a2,
                                                                   const double3 &a6, const double3 &a7,
                                                                   const BackboneGeometry &geometry,
                                                                   std::size_t gridPoints, double duplicateTolerance)
{
  const Rebridger solver(a1, a2, a6, a7, geometry);
  const double h = 2.0 * std::numbers::pi / static_cast<double>(gridPoints);

  std::vector<std::array<Branch, 4>> samples(gridPoints);
  for (std::size_t i = 0; i != gridPoints; ++i)
  {
    samples[i] = solver.evaluate(-std::numbers::pi + static_cast<double>(i) * h);
  }

  std::vector<Trimer> roots{};
  auto addRoot = [&](const Trimer &trimer)
  {
    for (const Trimer &known : roots)
    {
      if ((known.a3 - trimer.a3).length() < duplicateTolerance && (known.a4 - trimer.a4).length() < duplicateTolerance &&
          (known.a5 - trimer.a5).length() < duplicateTolerance)
      {
        return;
      }
    }
    roots.push_back(trimer);
  };

  for (std::size_t b = 0; b != 4; ++b)
  {
    for (std::size_t i = 0; i != gridPoints; ++i)
    {
      const std::size_t j = (i + 1) % gridPoints;
      const double phiI = -std::numbers::pi + static_cast<double>(i) * h;
      const double phiJ = phiI + h;  // the wrap-around interval is continuous in phi1 modulo 2 pi
      const Branch &atI = samples[i][b];
      const Branch &atJ = samples[j][b];

      if (atI.exists && atI.residual == 0.0) addRoot(atI.trimer);

      if (atI.exists && atJ.exists)
      {
        if (signChange(atI, atJ))
        {
          if (std::optional<Trimer> root = bisectRoot(solver, b, phiI, atI, phiJ, atJ)) addRoot(*root);
        }
      }
      else if (!atI.exists && atJ.exists)
      {
        // The branch starts between the two grid points: refine its boundary and look for a root
        // between the boundary and the first inside grid point.
        const double alpha = branchBoundary(solver, b, phiI, phiJ);
        const Branch atAlpha = solver.evaluate(alpha, b);
        if (atAlpha.exists && signChange(atAlpha, atJ))
        {
          if (std::optional<Trimer> root = bisectRoot(solver, b, alpha, atAlpha, phiJ, atJ)) addRoot(*root);
        }
      }
      else if (atI.exists && !atJ.exists)
      {
        const double beta = branchBoundary(solver, b, phiJ, phiI);
        const Branch atBeta = solver.evaluate(beta, b);
        if (atBeta.exists && signChange(atI, atBeta))
        {
          if (std::optional<Trimer> root = bisectRoot(solver, b, phiI, atI, beta, atBeta)) addRoot(*root);
        }
      }
    }

    // Pairs of roots between grid points (no sign change): refine every local minimum of |residual|.
    for (std::size_t i = 0; i != gridPoints; ++i)
    {
      const Branch &atI = samples[i][b];
      if (!atI.exists || atI.residual == 0.0) continue;
      const std::size_t previous = (i + gridPoints - 1) % gridPoints;
      const std::size_t next = (i + 1) % gridPoints;
      const double phiI = -std::numbers::pi + static_cast<double>(i) * h;

      // Neighbouring points on the branch; where the branch ends inside the interval, its boundary
      // is used instead so that root pairs next to a boundary are found as well.
      auto neighbour = [&](std::size_t index, double phiNeighbour) -> std::pair<double, Branch>
      {
        if (samples[index][b].exists) return {phiNeighbour, samples[index][b]};
        const double boundary = branchBoundary(solver, b, phiNeighbour, phiI);
        return {boundary, solver.evaluate(boundary, b)};
      };
      auto [phiPrevious, atPrevious] = neighbour(previous, phiI - h);
      auto [phiNext, atNext] = neighbour(next, phiI + h);
      if (!atPrevious.exists || !atNext.exists) continue;
      if (signChange(atPrevious, atI) || signChange(atI, atNext)) continue;  // handled by the scan
      if (std::abs(atI.residual) > std::abs(atPrevious.residual) || std::abs(atI.residual) >= std::abs(atNext.residual))
        continue;

      refineLocalMinimum(solver, b, {phiPrevious, phiI, phiNext}, {atPrevious, atI, atNext}, addRoot);
    }
  }
  return roots;
}

double ConcertedRotation::determinant(std::vector<double> m, std::size_t n)
{
  double det = 1.0;
  for (std::size_t col = 0; col != n; ++col)
  {
    std::size_t pivot = col;
    for (std::size_t row = col + 1; row < n; ++row)
    {
      if (std::abs(m[row * n + col]) > std::abs(m[pivot * n + col])) pivot = row;
    }
    if (m[pivot * n + col] == 0.0) return 0.0;
    if (pivot != col)
    {
      for (std::size_t k = 0; k != n; ++k) std::swap(m[pivot * n + k], m[col * n + k]);
      det = -det;
    }
    det *= m[col * n + col];
    for (std::size_t row = col + 1; row < n; ++row)
    {
      const double factor = m[row * n + col] / m[col * n + col];
      for (std::size_t k = col; k != n; ++k) m[row * n + k] -= factor * m[col * n + k];
    }
  }
  return det;
}

double ConcertedRotation::closureJacobian(std::span<const double3, 7> a, const std::optional<double3> &a8)
{
  // a[k-1] is atom a_k (k = 1..7). Column k (k = 1..6) is the rotation about the bond a_k a_{k+1};
  // it moves every atom a_m with m >= k+2 by u_k x (a_m - a_k).
  const std::size_t n = a8.has_value() ? 6 : 5;
  const double3 a6 = a[5], a7 = a[6];

  // Chart of the rigid placements of (a6, a7, a8): a6, the two components of a7 perpendicular to
  // the a6a7 bond, and the component of a8 along the normal of the a6a7a8 plane.
  const Basis bond67 = basisAlong(a7 - a6);
  double3 normal{};
  if (a8.has_value())
  {
    normal = double3::cross(bond67.u, a8.value() - a7);
    const double length = normal.length();
    if (length < 1e-12) return 0.0;  // collinear a6, a7, a8: degenerate chart
    normal = normal / length;
  }

  std::vector<double> matrix(n * n, 0.0);
  for (std::size_t k = 1; k <= n; ++k)
  {
    const double3 origin = a[k - 1];
    const double3 u = (a[k] - origin).normalized();
    auto derivative = [&](std::size_t m, const double3 &position) -> double3
    { return m >= k + 2 ? double3::cross(u, position - origin) : double3{}; };

    const double3 d6 = derivative(6, a6);
    const double3 d7 = derivative(7, a7);
    const std::size_t col = k - 1;
    matrix[0 * n + col] = d6.x;
    matrix[1 * n + col] = d6.y;
    matrix[2 * n + col] = d6.z;
    matrix[3 * n + col] = double3::dot(bond67.n1, d7);
    matrix[4 * n + col] = double3::dot(bond67.n2, d7);
    if (a8.has_value())
    {
      matrix[5 * n + col] = double3::dot(normal, derivative(8, a8.value()));
    }
  }
  return std::abs(determinant(std::move(matrix), n));
}

double3 ConcertedRotation::rotateAboutAxis(const double3 &origin, const double3 &axis, double angle,
                                           const double3 &point)
{
  const double3 r = point - origin;
  const double c = std::cos(angle), s = std::sin(angle);
  return origin + r * c + double3::cross(axis, r) * s + axis * (double3::dot(axis, r) * (1.0 - c));
}

double3x3 ConcertedRotation::localFrame(const double3 &atom, const double3 &previous, const double3 &next)
{
  const double3 e1 = (previous - atom).normalized();
  const double3 toNext = next - atom;
  const double3 e2 = (toNext - double3::dot(toNext, e1) * e1).normalized();
  const double3 e3 = double3::cross(e1, e2);
  return double3x3(e1, e2, e3);
}

double3 ConcertedRotation::transformRigidly(const double3 &oldOrigin, const double3x3 &oldFrame,
                                            const double3 &newOrigin, const double3x3 &newFrame, const double3 &point)
{
  return newOrigin + newFrame * transposedMultiply(oldFrame, point - oldOrigin);
}
