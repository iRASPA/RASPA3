module;

module dihedral_derivatives;

import std;

import double3;

namespace
{
// A minimal forward-mode hyper-dual number in the twelve Cartesian coordinates of a dihedral: value, gradient and
// Hessian, exact through the elementary operations.
struct Dual
{
  double v{};
  std::array<double, 12> g{};
  std::array<std::array<double, 12>, 12> h{};
};

Dual operator+(const Dual& a, const Dual& b)
{
  Dual r{};
  r.v = a.v + b.v;
  for (std::size_t i = 0; i < 12; ++i)
  {
    r.g[i] = a.g[i] + b.g[i];
    for (std::size_t j = 0; j < 12; ++j) r.h[i][j] = a.h[i][j] + b.h[i][j];
  }
  return r;
}

Dual operator-(const Dual& a, const Dual& b)
{
  Dual r{};
  r.v = a.v - b.v;
  for (std::size_t i = 0; i < 12; ++i)
  {
    r.g[i] = a.g[i] - b.g[i];
    for (std::size_t j = 0; j < 12; ++j) r.h[i][j] = a.h[i][j] - b.h[i][j];
  }
  return r;
}

Dual operator*(const Dual& a, const Dual& b)
{
  Dual r{};
  r.v = a.v * b.v;
  for (std::size_t i = 0; i < 12; ++i)
  {
    r.g[i] = a.g[i] * b.v + a.v * b.g[i];
    for (std::size_t j = 0; j < 12; ++j)
    {
      r.h[i][j] = a.h[i][j] * b.v + a.g[i] * b.g[j] + a.g[j] * b.g[i] + a.v * b.h[i][j];
    }
  }
  return r;
}

// f(u) with f' and f'' at u
Dual unary(const Dual& u, double f, double df, double ddf)
{
  Dual r{};
  r.v = f;
  for (std::size_t i = 0; i < 12; ++i)
  {
    r.g[i] = df * u.g[i];
    for (std::size_t j = 0; j < 12; ++j) r.h[i][j] = ddf * u.g[i] * u.g[j] + df * u.h[i][j];
  }
  return r;
}

Dual reciprocal(const Dual& u) { return unary(u, 1.0 / u.v, -1.0 / (u.v * u.v), 2.0 / (u.v * u.v * u.v)); }

Dual squareRoot(const Dual& u)
{
  const double f = std::sqrt(u.v);
  return unary(u, f, 0.5 / f, -0.25 / (f * u.v));
}

Dual arcTangent(const Dual& u)
{
  const double denominator = 1.0 + u.v * u.v;
  return unary(u, std::atan(u.v), 1.0 / denominator, -2.0 * u.v / (denominator * denominator));
}

// atan2(y, x) as atan of the ratio whose magnitude is at most one: smooth everywhere but at x = y = 0
Dual arcTangent2(const Dual& y, const Dual& x)
{
  if (std::fabs(x.v) >= std::fabs(y.v))
  {
    Dual angle = arcTangent(y * reciprocal(x));
    if (x.v < 0.0) angle.v += y.v >= 0.0 ? std::numbers::pi : -std::numbers::pi;
    return angle;
  }
  Dual angle = arcTangent(x * reciprocal(y));
  for (std::size_t i = 0; i < 12; ++i)
  {
    angle.g[i] = -angle.g[i];
    for (std::size_t j = 0; j < 12; ++j) angle.h[i][j] = -angle.h[i][j];
  }
  angle.v = (y.v > 0.0 ? 0.5 * std::numbers::pi : -0.5 * std::numbers::pi) - angle.v;
  return angle;
}

using DualVector = std::array<Dual, 3>;

DualVector difference(const DualVector& a, const DualVector& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }

Dual dot(const DualVector& a, const DualVector& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

DualVector cross(const DualVector& a, const DualVector& b)
{
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
}  // namespace

double Potentials::dihedralAngle(const double3& posA, const double3& posB, const double3& posC, const double3& posD)
{
  const double3 b1 = posB - posA;
  const double3 b2 = posC - posB;
  const double3 b3 = posD - posC;
  return std::atan2(std::sqrt(b2.length_squared()) * double3::dot(b1, double3::cross(b2, b3)),
                    double3::dot(double3::cross(b1, b2), double3::cross(b2, b3)));
}

Potentials::DihedralDerivatives Potentials::dihedralAngleDerivatives(const double3& posA, const double3& posB,
                                                                     const double3& posC, const double3& posD,
                                                                     bool computeHessian)
{
  DihedralDerivatives result{};

  // b1 = B - A, b2 = C - B, b3 = D - C; phi = atan2(|b2| b1.(b2 x b3), (b1 x b2).(b2 x b3))
  const double3 b1 = posB - posA;
  const double3 b2 = posC - posB;
  const double3 b3 = posD - posC;
  const double3 n1 = double3::cross(b1, b2);
  const double3 n2 = double3::cross(b2, b3);
  const double n1Squared = double3::dot(n1, n1);
  const double n2Squared = double3::dot(n2, n2);
  const double b2Squared = double3::dot(b2, b2);
  const double b2Length = std::sqrt(b2Squared);
  result.phi = std::atan2(b2Length * double3::dot(b1, n2), double3::dot(n1, n2));

  // collinear arms: phi is undefined and so are its derivatives
  constexpr double degenerate = 1.0e-20;
  if (!(n1Squared > degenerate) || !(n2Squared > degenerate) || !(b2Squared > degenerate)) return result;

  // Blondel-Karplus, with F = A - B = -b1, G = B - C = -b2, H = D - C = b3, so that their A = F x G = n1 and
  // B = H x G = n2, F.G = b1.b2 and H.G = -b3.b2
  const double fDotG = double3::dot(b1, b2);
  const double hDotG = -double3::dot(b3, b2);
  const double3 gradientA = (-b2Length / n1Squared) * n1;
  const double3 gradientD = (b2Length / n2Squared) * n2;
  const double3 cross1 = (fDotG / (n1Squared * b2Length)) * n1;
  const double3 cross2 = (hDotG / (n2Squared * b2Length)) * n2;
  result.gradient[0] = gradientA;
  result.gradient[1] = -1.0 * gradientA + cross1 - cross2;
  result.gradient[2] = -1.0 * gradientD - cross1 + cross2;
  result.gradient[3] = gradientD;

  if (!computeHessian) return result;

  // the Hessian of the same expression through hyper-dual numbers
  std::array<DualVector, 4> positions{};
  const std::array<const double3*, 4> sources = {&posA, &posB, &posC, &posD};
  for (std::size_t atom = 0; atom < 4; ++atom)
  {
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      Dual& coordinate = positions[atom][axis];
      coordinate.v = (&sources[atom]->x)[axis];
      coordinate.g[3 * atom + axis] = 1.0;
    }
  }
  const DualVector B1 = difference(positions[1], positions[0]);
  const DualVector B2 = difference(positions[2], positions[1]);
  const DualVector B3 = difference(positions[3], positions[2]);
  const DualVector N2 = cross(B2, B3);
  const Dual y = squareRoot(dot(B2, B2)) * dot(B1, N2);
  const Dual x = dot(cross(B1, B2), N2);
  const Dual phi = arcTangent2(y, x);
  result.hessian = phi.h;
  return result;
}
