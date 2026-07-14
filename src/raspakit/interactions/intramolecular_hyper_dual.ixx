module;

export module interactions_intramolecular_hyper_dual;

import std;

import double3;
import double3x3;
import generalized_hessian;
import minimization_dof_layout;
import minimization_internal_coordinate_hessian;

export namespace Interactions::HyperDual
{
using Minimization::CoordinateBlock;
using Minimization::InternalCoordinate;

inline double component(const double3& v, std::size_t axis) { return (&v.x)[axis]; }

/**
 * Forward-mode hyper-dual number carrying value, gradient, and Hessian with respect to N
 * independent scalar coordinates. All derivatives are exact (analytic chain rule on the
 * elementary operations), so the coordinate-space partials dU/dq_k and d²U/dq_k dq_l are
 * obtained without finite differences.
 */
template <std::size_t N>
struct Dual
{
  double v{};
  std::array<double, N> g{};
  std::array<std::array<double, N>, N> h{};

  static Dual constant(double value)
  {
    Dual result{};
    result.v = value;
    return result;
  }

  static Dual variable(double value, std::size_t index)
  {
    Dual result{};
    result.v = value;
    result.g[index] = 1.0;
    return result;
  }
};

template <std::size_t N>
Dual<N> operator+(const Dual<N>& a, const Dual<N>& b)
{
  Dual<N> result{};
  result.v = a.v + b.v;
  for (std::size_t i = 0; i < N; ++i)
  {
    result.g[i] = a.g[i] + b.g[i];
    for (std::size_t j = 0; j < N; ++j) result.h[i][j] = a.h[i][j] + b.h[i][j];
  }
  return result;
}

template <std::size_t N>
Dual<N> operator-(const Dual<N>& a, const Dual<N>& b)
{
  Dual<N> result{};
  result.v = a.v - b.v;
  for (std::size_t i = 0; i < N; ++i)
  {
    result.g[i] = a.g[i] - b.g[i];
    for (std::size_t j = 0; j < N; ++j) result.h[i][j] = a.h[i][j] - b.h[i][j];
  }
  return result;
}

template <std::size_t N>
Dual<N> operator*(const Dual<N>& a, const Dual<N>& b)
{
  Dual<N> result{};
  result.v = a.v * b.v;
  for (std::size_t i = 0; i < N; ++i)
  {
    result.g[i] = a.g[i] * b.v + a.v * b.g[i];
  }
  for (std::size_t i = 0; i < N; ++i)
  {
    for (std::size_t j = 0; j < N; ++j)
    {
      result.h[i][j] = a.h[i][j] * b.v + a.g[i] * b.g[j] + a.g[j] * b.g[i] + a.v * b.h[i][j];
    }
  }
  return result;
}

template <std::size_t N>
Dual<N> operator*(const Dual<N>& a, double s)
{
  Dual<N> result{};
  result.v = a.v * s;
  for (std::size_t i = 0; i < N; ++i)
  {
    result.g[i] = a.g[i] * s;
    for (std::size_t j = 0; j < N; ++j) result.h[i][j] = a.h[i][j] * s;
  }
  return result;
}

template <std::size_t N>
Dual<N> operator*(double s, const Dual<N>& a)
{
  return a * s;
}

template <std::size_t N>
Dual<N> operator+(const Dual<N>& a, double s)
{
  Dual<N> result = a;
  result.v += s;
  return result;
}

template <std::size_t N>
Dual<N> operator+(double s, const Dual<N>& a)
{
  return a + s;
}

template <std::size_t N>
Dual<N> operator-(const Dual<N>& a, double s)
{
  return a + (-s);
}

template <std::size_t N>
Dual<N> operator-(double s, const Dual<N>& a)
{
  return (a * -1.0) + s;
}

// Apply a smooth unary function with known first and second derivatives.
template <std::size_t N>
Dual<N> applyUnary(const Dual<N>& u, double f, double df, double ddf)
{
  Dual<N> result{};
  result.v = f;
  for (std::size_t i = 0; i < N; ++i)
  {
    result.g[i] = df * u.g[i];
  }
  for (std::size_t i = 0; i < N; ++i)
  {
    for (std::size_t j = 0; j < N; ++j)
    {
      result.h[i][j] = ddf * u.g[i] * u.g[j] + df * u.h[i][j];
    }
  }
  return result;
}

template <std::size_t N>
Dual<N> expDual(const Dual<N>& u)
{
  const double e = std::exp(u.v);
  return applyUnary(u, e, e, e);
}

template <std::size_t N>
Dual<N> sinDual(const Dual<N>& u)
{
  return applyUnary(u, std::sin(u.v), std::cos(u.v), -std::sin(u.v));
}

template <std::size_t N>
Dual<N> cosDual(const Dual<N>& u)
{
  return applyUnary(u, std::cos(u.v), -std::sin(u.v), -std::cos(u.v));
}

template <std::size_t N>
Dual<N> sqrtDual(const Dual<N>& u)
{
  const double f = std::sqrt(u.v);
  return applyUnary(u, f, 0.5 / f, -0.25 / (f * u.v));
}

template <std::size_t N>
Dual<N> acosDual(const Dual<N>& u)
{
  const double x = u.v;
  const double oneMinus = 1.0 - x * x;
  const double root = std::sqrt(oneMinus);
  return applyUnary(u, std::acos(x), -1.0 / root, -x / (oneMinus * root));
}

template <std::size_t N>
Dual<N> powDual(const Dual<N>& u, double p)
{
  const double f = std::pow(u.v, p);
  const double df = p * std::pow(u.v, p - 1.0);
  const double ddf = p * (p - 1.0) * std::pow(u.v, p - 2.0);
  return applyUnary(u, f, df, ddf);
}

template <std::size_t N>
Dual<N> reciprocal(const Dual<N>& u)
{
  const double f = 1.0 / u.v;
  return applyUnary(u, f, -f * f, 2.0 * f * f * f);
}

template <std::size_t N>
Dual<N> operator/(const Dual<N>& a, const Dual<N>& b)
{
  return a * reciprocal(b);
}

// Vector of hyper-dual scalars, used to differentiate a Cartesian energy expression directly.
template <std::size_t N>
using DualVec = std::array<Dual<N>, 3>;

template <std::size_t N>
DualVec<N> operator-(const DualVec<N>& a, const DualVec<N>& b)
{
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

template <std::size_t N>
Dual<N> dotDual(const DualVec<N>& a, const DualVec<N>& b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

template <std::size_t N>
DualVec<N> crossDual(const DualVec<N>& a, const DualVec<N>& b)
{
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

// Torsion smoothing S(theta): 1 below 170 degrees, tapering to 0 at 180 degrees.
template <std::size_t N>
Dual<N> smoothingDual(const Dual<N>& theta)
{
  const double on = 170.0 * (std::numbers::pi / 180.0);
  const double off = std::numbers::pi;

  if (theta.v < on)
  {
    return Dual<N>::constant(1.0);
  }
  const double scale = 1.0 / ((off - on) * (off - on) * (off - on));
  const Dual<N> offMinus = off - theta;
  const Dual<N> bracket = (2.0 * theta + (off - 3.0 * on));
  return (offMinus * offMinus * bracket) * scale;
}

// Multiply a coordinate 3x3 Hessian block by a vector: (block * p).
inline double3 blockTimesVector(const CoordinateBlock& block, const double3& p)
{
  return {block[0][0] * p.x + block[0][1] * p.y + block[0][2] * p.z,
          block[1][0] * p.x + block[1][1] * p.y + block[1][2] * p.z,
          block[2][0] * p.x + block[2][1] * p.y + block[2][2] * p.z};
}

// Extract coordinate-space first/second derivatives from a hyper-dual evaluation.
template <std::size_t N>
void extractDerivatives(const Dual<N>& energy, std::array<double, 3>& dUdq, std::array<std::array<double, 3>, 3>& d2)
{
  dUdq = {};
  d2 = {};
  for (std::size_t i = 0; i < N; ++i)
  {
    dUdq[i] = energy.g[i];
    for (std::size_t j = 0; j < N; ++j)
    {
      d2[i][j] = energy.h[i][j];
    }
  }
}

/**
 * Scatter the strain second-derivative blocks of a cross term into the generalized Hessian.
 *
 * The generalized Hessian uses a single isotropic strain scalar ε under which real positions
 * scale exponentially, r_i(ε) = exp(ε) r_i (matching the pair-distance convention already used by
 * the bond/Urey-Bradley/van der Waals/Coulomb terms). Because d(r_i)/dε = r_i and d²(r_i)/dε² = r_i,
 * the strain derivative of a translation-invariant internal coordinate q is
 *   dq/dε        = sum_i (dq/dr_i) . r_i,
 *   d²q/dε²      = dq/dε + sum_{i,j} r_i . (d²q/dr_i dr_j) . r_j,
 *   d²q/dr_i dε  = dq/dr_i + sum_j (d²q/dr_i dr_j) . r_j.
 * Angles and dihedrals are scale-invariant (degree-0 homogeneous), so only distance coordinates
 * contribute; the remaining terms evaluate to zero automatically. The energy strain derivatives
 * follow from the same chain rule as the position-position block.
 */
inline void scatterCrossStrain(GeneralizedHessian& hessian, const MinimizationDofLayout& layout,
                               std::size_t moleculeIndex, const std::vector<std::size_t>& termAtoms,
                               const std::vector<InternalCoordinate>& coords,
                               const std::vector<std::array<std::size_t, 4>>& slotMaps,
                               const std::array<double3, 4>& termPositions, const std::array<double, 3>& dUdq,
                               const std::array<std::array<double, 3>, 3>& d2Udq2)
{
  const std::size_t M = termAtoms.size();
  const std::size_t K = coords.size();

  std::array<double, 3> dqdEps{};
  std::array<double, 3> d2qdEps2{};
  std::array<std::array<double3, 4>, 3> gradKT{};     // dq_k/dr at term atom t (0 if untouched)
  std::array<std::array<double3, 4>, 3> d2qdxdEps{};  // d²q_k/dr_t dε at term atom t

  for (std::size_t k = 0; k < K; ++k)
  {
    const InternalCoordinate& c = coords[k];

    for (std::size_t si = 0; si < c.numAtoms; ++si)
    {
      const std::size_t t = slotMaps[k][si];
      gradKT[k][t] = c.gradient[si];
      dqdEps[k] += double3::dot(c.gradient[si], termPositions[t]);
    }

    double quadratic = 0.0;
    for (std::size_t si = 0; si < c.numAtoms; ++si)
    {
      const std::size_t t = slotMaps[k][si];
      double3 mixed = c.gradient[si];
      for (std::size_t sj = 0; sj < c.numAtoms; ++sj)
      {
        const double3 contribution = blockTimesVector(c.hessian[si][sj], termPositions[slotMaps[k][sj]]);
        mixed += contribution;
        quadratic += double3::dot(termPositions[t], contribution);
      }
      d2qdxdEps[k][t] = mixed;
    }
    d2qdEps2[k] = dqdEps[k] + quadratic;
  }

  double strainStrain = 0.0;
  for (std::size_t k = 0; k < K; ++k)
  {
    strainStrain += dUdq[k] * d2qdEps2[k];
    for (std::size_t l = 0; l < K; ++l)
    {
      strainStrain += d2Udq2[k][l] * dqdEps[k] * dqdEps[l];
    }
  }
  hessian.addStrainStrain(0, 0, strainStrain);

  for (std::size_t t = 0; t < M; ++t)
  {
    const auto base = layout.flexibleAtomDof(moleculeIndex, termAtoms[t], MinimizationDofAxis::X);
    if (!base)
    {
      continue;
    }
    double3 positionStrain{};
    for (std::size_t k = 0; k < K; ++k)
    {
      positionStrain += dUdq[k] * d2qdxdEps[k][t];
      for (std::size_t l = 0; l < K; ++l)
      {
        positionStrain += (d2Udq2[k][l] * dqdEps[l]) * gradKT[k][t];
      }
    }
    hessian.addPositionStrain(*base + 0, 0, positionStrain.x);
    hessian.addPositionStrain(*base + 1, 0, positionStrain.y);
    hessian.addPositionStrain(*base + 2, 0, positionStrain.z);
  }
}

/**
 * Assemble and scatter the position-position Hessian of a cross term U(q_0, ..., q_{K-1}) from
 * the coordinate-space partials and the analytic Cartesian derivatives of each coordinate:
 *   d²U/dx dy = sum_k (dU/dq_k) d²q_k/dx dy + sum_{k,l} (d²U/dq_k dq_l)(dq_k/dx)(dq_l/dy).
 * When the generalized Hessian carries an isotropic strain degree of freedom, the position-strain
 * and strain-strain blocks are scattered as well.
 */
inline void scatterCrossHessian(GeneralizedHessian& hessian, const MinimizationDofLayout& layout,
                                std::size_t moleculeIndex, const std::vector<std::size_t>& termAtoms,
                                const std::vector<InternalCoordinate>& coords,
                                const std::vector<std::array<std::size_t, 4>>& slotMaps,
                                const std::array<double3, 4>& termPositions, const std::array<double, 3>& dUdq,
                                const std::array<std::array<double, 3>, 3>& d2Udq2)
{
  const std::size_t M = termAtoms.size();
  const std::size_t K = coords.size();

  std::array<std::array<CoordinateBlock, 4>, 4> full{};

  for (std::size_t k = 0; k < K; ++k)
  {
    const InternalCoordinate& c = coords[k];
    for (std::size_t si = 0; si < c.numAtoms; ++si)
    {
      for (std::size_t sj = 0; sj < c.numAtoms; ++sj)
      {
        const std::size_t ti = slotMaps[k][si];
        const std::size_t tj = slotMaps[k][sj];
        for (std::size_t a = 0; a < 3; ++a)
        {
          for (std::size_t b = 0; b < 3; ++b)
          {
            full[ti][tj][a][b] += dUdq[k] * c.hessian[si][sj][a][b];
          }
        }
      }
    }
  }

  for (std::size_t k = 0; k < K; ++k)
  {
    for (std::size_t l = 0; l < K; ++l)
    {
      const double factor = d2Udq2[k][l];
      if (factor == 0.0)
      {
        continue;
      }
      const InternalCoordinate& ck = coords[k];
      const InternalCoordinate& cl = coords[l];
      for (std::size_t si = 0; si < ck.numAtoms; ++si)
      {
        for (std::size_t sj = 0; sj < cl.numAtoms; ++sj)
        {
          const std::size_t ti = slotMaps[k][si];
          const std::size_t tj = slotMaps[l][sj];
          for (std::size_t a = 0; a < 3; ++a)
          {
            for (std::size_t b = 0; b < 3; ++b)
            {
              full[ti][tj][a][b] += factor * component(ck.gradient[si], a) * component(cl.gradient[sj], b);
            }
          }
        }
      }
    }
  }

  for (std::size_t ti = 0; ti < M; ++ti)
  {
    const auto baseI = layout.flexibleAtomDof(moleculeIndex, termAtoms[ti], MinimizationDofAxis::X);
    if (!baseI)
    {
      continue;
    }
    for (std::size_t tj = 0; tj < M; ++tj)
    {
      const auto baseJ = layout.flexibleAtomDof(moleculeIndex, termAtoms[tj], MinimizationDofAxis::X);
      if (!baseJ)
      {
        continue;
      }
      for (std::size_t a = 0; a < 3; ++a)
      {
        for (std::size_t b = 0; b < 3; ++b)
        {
          hessian.add(*baseI + a, *baseJ + b, full[ti][tj][a][b]);
        }
      }
    }
  }

  if (hessian.numStrain() == 1)
  {
    scatterCrossStrain(hessian, layout, moleculeIndex, termAtoms, coords, slotMaps, termPositions, dUdq, d2Udq2);
  }
}

/**
 * Scatter a term whose full Cartesian gradient and Hessian over four atoms are already known
 * (used for the inversion bend and out-of-plane bend, whose angles are not expressible in the
 * distance/angle/dihedral basis). The isotropic strain blocks follow from the same exp-strain
 * convention:
 *   d²E/dε²      = sum_i g_i . x_i + sum_{i,j} x_i . H_{ij} . x_j,
 *   d²E/dx_t dε  = g_t + sum_j H_{tj} . x_j.
 */
inline void scatterCartesianFourBody(GeneralizedHessian& hessian, const MinimizationDofLayout& layout,
                                     std::size_t moleculeIndex, const std::array<std::size_t, 4>& termAtoms,
                                     const std::array<double3, 4>& positions, const std::array<double3, 4>& gradient,
                                     const std::array<std::array<CoordinateBlock, 4>, 4>& H)
{
  std::array<std::optional<std::size_t>, 4> base{};
  for (std::size_t t = 0; t < 4; ++t)
  {
    base[t] = layout.flexibleAtomDof(moleculeIndex, termAtoms[t], MinimizationDofAxis::X);
  }

  for (std::size_t ti = 0; ti < 4; ++ti)
  {
    if (!base[ti]) continue;
    for (std::size_t tj = 0; tj < 4; ++tj)
    {
      if (!base[tj]) continue;
      for (std::size_t a = 0; a < 3; ++a)
      {
        for (std::size_t b = 0; b < 3; ++b)
        {
          hessian.add(*base[ti] + a, *base[tj] + b, H[ti][tj][a][b]);
        }
      }
    }
  }

  if (hessian.numStrain() != 1)
  {
    return;
  }

  double strainStrain = 0.0;
  for (std::size_t i = 0; i < 4; ++i)
  {
    strainStrain += double3::dot(gradient[i], positions[i]);
    for (std::size_t j = 0; j < 4; ++j)
    {
      strainStrain += double3::dot(positions[i], blockTimesVector(H[i][j], positions[j]));
    }
  }
  hessian.addStrainStrain(0, 0, strainStrain);

  for (std::size_t t = 0; t < 4; ++t)
  {
    if (!base[t]) continue;
    double3 positionStrain = gradient[t];
    for (std::size_t j = 0; j < 4; ++j)
    {
      positionStrain += blockTimesVector(H[t][j], positions[j]);
    }
    hessian.addPositionStrain(*base[t] + 0, 0, positionStrain.x);
    hessian.addPositionStrain(*base[t] + 1, 0, positionStrain.y);
    hessian.addPositionStrain(*base[t] + 2, 0, positionStrain.z);
  }
}
}  // namespace Interactions::HyperDual
