module;

module interactions_hessian_intramolecular;

import std;

import double3;
import double3x3;
import units;
import intra_molecular_potentials;
import bond_bond_potential;
import bond_bend_potential;
import bend_bend_potential;
import bond_torsion_potential;
import bend_torsion_potential;
import inversion_bend_potential;
import minimization_dof_layout;
import minimization_internal_coordinate_hessian;

namespace
{
using Minimization::CoordinateBlock;
using Minimization::InternalCoordinate;

double component(const double3 &v, std::size_t axis) { return (&v.x)[axis]; }

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
Dual<N> operator+(const Dual<N> &a, const Dual<N> &b)
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
Dual<N> operator-(const Dual<N> &a, const Dual<N> &b)
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
Dual<N> operator*(const Dual<N> &a, const Dual<N> &b)
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
Dual<N> operator*(const Dual<N> &a, double s)
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
Dual<N> operator*(double s, const Dual<N> &a)
{
  return a * s;
}

template <std::size_t N>
Dual<N> operator+(const Dual<N> &a, double s)
{
  Dual<N> result = a;
  result.v += s;
  return result;
}

template <std::size_t N>
Dual<N> operator+(double s, const Dual<N> &a)
{
  return a + s;
}

template <std::size_t N>
Dual<N> operator-(const Dual<N> &a, double s)
{
  return a + (-s);
}

template <std::size_t N>
Dual<N> operator-(double s, const Dual<N> &a)
{
  return (a * -1.0) + s;
}

// Apply a smooth unary function with known first and second derivatives.
template <std::size_t N>
Dual<N> applyUnary(const Dual<N> &u, double f, double df, double ddf)
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
Dual<N> expDual(const Dual<N> &u)
{
  const double e = std::exp(u.v);
  return applyUnary(u, e, e, e);
}

template <std::size_t N>
Dual<N> sinDual(const Dual<N> &u)
{
  return applyUnary(u, std::sin(u.v), std::cos(u.v), -std::sin(u.v));
}

template <std::size_t N>
Dual<N> cosDual(const Dual<N> &u)
{
  return applyUnary(u, std::cos(u.v), -std::sin(u.v), -std::cos(u.v));
}

template <std::size_t N>
Dual<N> sqrtDual(const Dual<N> &u)
{
  const double f = std::sqrt(u.v);
  return applyUnary(u, f, 0.5 / f, -0.25 / (f * u.v));
}

template <std::size_t N>
Dual<N> acosDual(const Dual<N> &u)
{
  const double x = u.v;
  const double oneMinus = 1.0 - x * x;
  const double root = std::sqrt(oneMinus);
  return applyUnary(u, std::acos(x), -1.0 / root, -x / (oneMinus * root));
}

template <std::size_t N>
Dual<N> powDual(const Dual<N> &u, double p)
{
  const double f = std::pow(u.v, p);
  const double df = p * std::pow(u.v, p - 1.0);
  const double ddf = p * (p - 1.0) * std::pow(u.v, p - 2.0);
  return applyUnary(u, f, df, ddf);
}

template <std::size_t N>
Dual<N> reciprocal(const Dual<N> &u)
{
  const double f = 1.0 / u.v;
  return applyUnary(u, f, -f * f, 2.0 * f * f * f);
}

template <std::size_t N>
Dual<N> operator/(const Dual<N> &a, const Dual<N> &b)
{
  return a * reciprocal(b);
}

// Vector of hyper-dual scalars, used to differentiate a Cartesian energy expression directly.
template <std::size_t N>
using DualVec = std::array<Dual<N>, 3>;

template <std::size_t N>
DualVec<N> operator-(const DualVec<N> &a, const DualVec<N> &b)
{
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

template <std::size_t N>
Dual<N> dotDual(const DualVec<N> &a, const DualVec<N> &b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

template <std::size_t N>
DualVec<N> crossDual(const DualVec<N> &a, const DualVec<N> &b)
{
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

// Torsion smoothing S(theta): 1 below 170 degrees, tapering to 0 at 180 degrees.
template <std::size_t N>
Dual<N> smoothingDual(const Dual<N> &theta)
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
double3 blockTimesVector(const CoordinateBlock &block, const double3 &p)
{
  return {block[0][0] * p.x + block[0][1] * p.y + block[0][2] * p.z,
          block[1][0] * p.x + block[1][1] * p.y + block[1][2] * p.z,
          block[2][0] * p.x + block[2][1] * p.y + block[2][2] * p.z};
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
void scatterCrossStrain(GeneralizedHessian &hessian, const MinimizationDofLayout &layout, std::size_t moleculeIndex,
                        const std::vector<std::size_t> &termAtoms, const std::vector<InternalCoordinate> &coords,
                        const std::vector<std::array<std::size_t, 4>> &slotMaps,
                        const std::array<double3, 4> &termPositions, const std::array<double, 3> &dUdq,
                        const std::array<std::array<double, 3>, 3> &d2Udq2)
{
  const std::size_t M = termAtoms.size();
  const std::size_t K = coords.size();

  std::array<double, 3> dqdEps{};
  std::array<double, 3> d2qdEps2{};
  std::array<std::array<double3, 4>, 3> gradKT{};      // dq_k/dr at term atom t (0 if untouched)
  std::array<std::array<double3, 4>, 3> d2qdxdEps{};   // d²q_k/dr_t dε at term atom t

  for (std::size_t k = 0; k < K; ++k)
  {
    const InternalCoordinate &c = coords[k];

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
void scatterCrossHessian(GeneralizedHessian &hessian, const MinimizationDofLayout &layout, std::size_t moleculeIndex,
                         const std::vector<std::size_t> &termAtoms, const std::vector<InternalCoordinate> &coords,
                         const std::vector<std::array<std::size_t, 4>> &slotMaps,
                         const std::array<double3, 4> &termPositions, const std::array<double, 3> &dUdq,
                         const std::array<std::array<double, 3>, 3> &d2Udq2)
{
  const std::size_t M = termAtoms.size();
  const std::size_t K = coords.size();

  std::array<std::array<CoordinateBlock, 4>, 4> full{};

  for (std::size_t k = 0; k < K; ++k)
  {
    const InternalCoordinate &c = coords[k];
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
      const InternalCoordinate &ck = coords[k];
      const InternalCoordinate &cl = coords[l];
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
 * (used for the inversion bend, whose inversion angle is not expressible in the distance/angle/
 * dihedral basis). The isotropic strain blocks follow from the same exp-strain convention:
 *   d²E/dε²      = sum_i g_i . x_i + sum_{i,j} x_i . H_{ij} . x_j,
 *   d²E/dx_t dε  = g_t + sum_j H_{tj} . x_j.
 */
void scatterCartesianFourBody(GeneralizedHessian &hessian, const MinimizationDofLayout &layout,
                              std::size_t moleculeIndex, const std::array<std::size_t, 4> &termAtoms,
                              const std::array<double3, 4> &positions, const std::array<double3, 4> &gradient,
                              const std::array<std::array<CoordinateBlock, 4>, 4> &H)
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

// Extract coordinate-space first/second derivatives from a hyper-dual evaluation.
template <std::size_t N>
void extractDerivatives(const Dual<N> &energy, std::array<double, 3> &dUdq, std::array<std::array<double, 3>, 3> &d2)
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
}  // namespace

RunningEnergy Interactions::computeIntraMolecularBondBondHessian(std::span<const Molecule> moleculeData,
                                                                 std::span<const Atom> atoms,
                                                                 std::span<const Component> components,
                                                                 const MinimizationDofLayout &layout,
                                                                 GeneralizedHessian &hessian,
                                                                 std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule &molecule = moleculeData[moleculeIndex];
    if (components[molecule.componentId].rigid)
    {
      continue;
    }

    const Potentials::IntraMolecularPotentials &potentials = components[molecule.componentId].intraMolecularPotentials;
    std::span<const Atom> atom_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};

    for (const BondBondPotential &term : potentials.bondBonds)
    {
      const std::size_t A = term.identifiers[0];
      const std::size_t B = term.identifiers[1];
      const std::size_t C = term.identifiers[2];

      auto [energy, gradient, strain] =
          term.potentialEnergyGradientStrain(atom_span[A].position, atom_span[B].position, atom_span[C].position);
      energies.bondBond += energy;
      dynamics_span[A].gradient += gradient[0];
      dynamics_span[B].gradient += gradient[1];
      dynamics_span[C].gradient += gradient[2];
      hessian.strainGradient() += strain;

      const InternalCoordinate rab = Minimization::distanceInternalCoordinate(atom_span[A].position, atom_span[B].position);
      const InternalCoordinate rbc = Minimization::distanceInternalCoordinate(atom_span[C].position, atom_span[B].position);

      const Dual<2> qab = Dual<2>::variable(rab.value, 0);
      const Dual<2> qbc = Dual<2>::variable(rbc.value, 1);
      Dual<2> U{};
      switch (term.type)
      {
        case BondBondType::CVFF:
        case BondBondType::CFF:
          U = term.parameters[0] * (qab - term.parameters[1]) * (qbc - term.parameters[2]);
          break;
      }

      std::array<double, 3> dUdq{};
      std::array<std::array<double, 3>, 3> d2{};
      extractDerivatives(U, dUdq, d2);

      scatterCrossHessian(hessian, layout, moleculeIndex, {A, B, C}, {rab, rbc}, {{{0, 1}}, {{2, 1}}},
                          {atom_span[A].position, atom_span[B].position, atom_span[C].position, double3{}}, dUdq, d2);
    }
  }

  return energies;
}

RunningEnergy Interactions::computeIntraMolecularBondBendHessian(std::span<const Molecule> moleculeData,
                                                                 std::span<const Atom> atoms,
                                                                 std::span<const Component> components,
                                                                 const MinimizationDofLayout &layout,
                                                                 GeneralizedHessian &hessian,
                                                                 std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule &molecule = moleculeData[moleculeIndex];
    if (components[molecule.componentId].rigid)
    {
      continue;
    }

    const Potentials::IntraMolecularPotentials &potentials = components[molecule.componentId].intraMolecularPotentials;
    std::span<const Atom> atom_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};

    for (const BondBendPotential &term : potentials.bondBends)
    {
      const std::size_t A = term.identifiers[0];
      const std::size_t B = term.identifiers[1];
      const std::size_t C = term.identifiers[2];

      auto [energy, gradient, strain] =
          term.potentialEnergyGradientStrain(atom_span[A].position, atom_span[B].position, atom_span[C].position);
      energies.bondBend += energy;
      dynamics_span[A].gradient += gradient[0];
      dynamics_span[B].gradient += gradient[1];
      dynamics_span[C].gradient += gradient[2];
      hessian.strainGradient() += strain;

      const InternalCoordinate rab = Minimization::distanceInternalCoordinate(atom_span[A].position, atom_span[B].position);
      const InternalCoordinate rbc = Minimization::distanceInternalCoordinate(atom_span[C].position, atom_span[B].position);
      const InternalCoordinate theta =
          Minimization::angleInternalCoordinate(atom_span[A].position, atom_span[B].position, atom_span[C].position);

      const Dual<3> qab = Dual<3>::variable(rab.value, 0);
      const Dual<3> qbc = Dual<3>::variable(rbc.value, 1);
      const Dual<3> qth = Dual<3>::variable(theta.value, 2);
      const auto &p = term.parameters;
      Dual<3> U{};
      switch (term.type)
      {
        case BondBendType::CVFF:
        case BondBendType::CFF:
          U = (qth - p[0]) * (p[1] * (qab - p[2]) + p[3] * (qbc - p[4]));
          break;
        case BondBendType::MM3:
          U = p[0] * ((qab - p[1]) + (qbc - p[2])) * Units::RadiansToDegrees * (qth - p[3]);
          break;
        case BondBendType::TruncatedHarmonic:
        {
          const Dual<3> dth = qth - p[1];
          const Dual<3> exponent =
              (powDual(qab, 8.0) + powDual(qbc, 8.0)) * (-1.0 / std::pow(p[2], 8));
          U = 0.5 * p[0] * dth * dth * expDual(exponent);
          break;
        }
        case BondBendType::ScreenedHarmonic:
        {
          const Dual<3> dth = qth - p[1];
          const Dual<3> exponent = (qab * (-1.0 / p[2])) + (qbc * (-1.0 / p[3]));
          U = 0.5 * p[0] * dth * dth * expDual(exponent);
          break;
        }
        case BondBendType::ScreenedVessal:
        {
          const double pi = std::numbers::pi;
          const Dual<3> thetaMinusPi = qth - pi;
          const Dual<3> bracket = ((p[1] - pi) * (p[1] - pi)) - (thetaMinusPi * thetaMinusPi);
          const Dual<3> exponent = (qab * (-1.0 / p[2])) + (qbc * (-1.0 / p[3]));
          U = (Dual<3>::constant(p[0] / 8.0) * reciprocal(thetaMinusPi * thetaMinusPi)) * (bracket * bracket) *
              expDual(exponent);
          break;
        }
        case BondBendType::TruncatedVessal:
        {
          const double pi = std::numbers::pi;
          const Dual<3> dth = qth - p[1];
          const Dual<3> dth2 = qth + (p[1] - 2.0 * pi);
          const Dual<3> theTerm = powDual(qth, p[2]) * (dth * dth) * (dth2 * dth2);
          const double constTerm = 0.5 * p[2] * std::pow(pi, p[2] - 1.0) * std::pow(pi - p[1], 3.0);
          const Dual<3> theBracket = theTerm - (dth * dth) * constTerm;
          const Dual<3> exponent =
              (powDual(qab, 8.0) + powDual(qbc, 8.0)) * (-1.0 / std::pow(p[3], 8));
          U = p[0] * theBracket * expDual(exponent);
          break;
        }
      }

      std::array<double, 3> dUdq{};
      std::array<std::array<double, 3>, 3> d2{};
      extractDerivatives(U, dUdq, d2);

      scatterCrossHessian(hessian, layout, moleculeIndex, {A, B, C}, {rab, rbc, theta},
                          {{{0, 1}}, {{2, 1}}, {{0, 1, 2}}},
                          {atom_span[A].position, atom_span[B].position, atom_span[C].position, double3{}}, dUdq, d2);
    }
  }

  return energies;
}

RunningEnergy Interactions::computeIntraMolecularBendBendHessian(std::span<const Molecule> moleculeData,
                                                                 std::span<const Atom> atoms,
                                                                 std::span<const Component> components,
                                                                 const MinimizationDofLayout &layout,
                                                                 GeneralizedHessian &hessian,
                                                                 std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule &molecule = moleculeData[moleculeIndex];
    if (components[molecule.componentId].rigid)
    {
      continue;
    }

    const Potentials::IntraMolecularPotentials &potentials = components[molecule.componentId].intraMolecularPotentials;
    std::span<const Atom> atom_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};

    for (const BendBendPotential &term : potentials.bendBends)
    {
      const std::size_t A = term.identifiers[0];
      const std::size_t B = term.identifiers[1];
      const std::size_t C = term.identifiers[2];
      const std::size_t D = term.identifiers[3];

      auto [energy, gradient, strain] = term.potentialEnergyGradientStrain(
          atom_span[A].position, atom_span[B].position, atom_span[C].position, atom_span[D].position);
      energies.bendBend += energy;
      dynamics_span[A].gradient += gradient[0];
      dynamics_span[B].gradient += gradient[1];
      dynamics_span[C].gradient += gradient[2];
      dynamics_span[D].gradient += gradient[3];
      hessian.strainGradient() += strain;

      const InternalCoordinate theta1 =
          Minimization::angleInternalCoordinate(atom_span[A].position, atom_span[B].position, atom_span[C].position);
      const InternalCoordinate theta2 =
          Minimization::angleInternalCoordinate(atom_span[A].position, atom_span[B].position, atom_span[D].position);

      const Dual<2> q1 = Dual<2>::variable(theta1.value, 0);
      const Dual<2> q2 = Dual<2>::variable(theta2.value, 1);
      const auto &p = term.parameters;
      Dual<2> U{};
      switch (term.type)
      {
        case BendBendType::CVFF:
        case BendBendType::CFF:
          U = p[0] * (q1 - p[1]) * (q2 - p[2]);
          break;
        case BendBendType::MM3:
          U = (-p[0] * Units::RadiansToDegrees * Units::RadiansToDegrees) * (q1 - p[1]) * (q2 - p[2]);
          break;
      }

      std::array<double, 3> dUdq{};
      std::array<std::array<double, 3>, 3> d2{};
      extractDerivatives(U, dUdq, d2);

      scatterCrossHessian(
          hessian, layout, moleculeIndex, {A, B, C, D}, {theta1, theta2}, {{{0, 1, 2}}, {{0, 1, 3}}},
          {atom_span[A].position, atom_span[B].position, atom_span[C].position, atom_span[D].position}, dUdq, d2);
    }
  }

  return energies;
}

RunningEnergy Interactions::computeIntraMolecularBondTorsionHessian(std::span<const Molecule> moleculeData,
                                                                    std::span<const Atom> atoms,
                                                                    std::span<const Component> components,
                                                                    const MinimizationDofLayout &layout,
                                                                    GeneralizedHessian &hessian,
                                                                    std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule &molecule = moleculeData[moleculeIndex];
    if (components[molecule.componentId].rigid)
    {
      continue;
    }

    const Potentials::IntraMolecularPotentials &potentials = components[molecule.componentId].intraMolecularPotentials;
    std::span<const Atom> atom_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};

    for (const BondTorsionPotential &term : potentials.bondTorsions)
    {
      const std::size_t A = term.identifiers[0];
      const std::size_t B = term.identifiers[1];
      const std::size_t C = term.identifiers[2];
      const std::size_t D = term.identifiers[3];

      auto [energy, gradient, strain] = term.potentialEnergyGradientStrain(
          atom_span[A].position, atom_span[B].position, atom_span[C].position, atom_span[D].position);
      energies.bondTorsion += energy;
      dynamics_span[A].gradient += gradient[0];
      dynamics_span[B].gradient += gradient[1];
      dynamics_span[C].gradient += gradient[2];
      dynamics_span[D].gradient += gradient[3];
      hessian.strainGradient() += strain;

      const InternalCoordinate rbc = Minimization::distanceInternalCoordinate(atom_span[B].position, atom_span[C].position);
      const InternalCoordinate cphi = Minimization::dihedralInternalCoordinate(
          atom_span[A].position, atom_span[B].position, atom_span[C].position, atom_span[D].position);

      const Dual<2> qr = Dual<2>::variable(rbc.value, 0);
      const Dual<2> qc = Dual<2>::variable(cphi.value, 1);
      const auto &p = term.parameters;
      Dual<2> U{};
      switch (term.type)
      {
        case BondTorsionType::MM3:
        {
          const Dual<2> temp = qr - p[3];
          const Dual<2> c2 = qc * qc;
          U = p[0] * temp * qc + p[1] * temp * (2.0 * c2 - 1.0) +
              p[2] * temp * ((4.0 * c2 * qc) - 3.0 * qc);
          break;
        }
      }

      std::array<double, 3> dUdq{};
      std::array<std::array<double, 3>, 3> d2{};
      extractDerivatives(U, dUdq, d2);

      scatterCrossHessian(
          hessian, layout, moleculeIndex, {A, B, C, D}, {rbc, cphi}, {{{1, 2}}, {{0, 1, 2, 3}}},
          {atom_span[A].position, atom_span[B].position, atom_span[C].position, atom_span[D].position}, dUdq, d2);
    }
  }

  return energies;
}

RunningEnergy Interactions::computeIntraMolecularBendTorsionHessian(std::span<const Molecule> moleculeData,
                                                                    std::span<const Atom> atoms,
                                                                    std::span<const Component> components,
                                                                    const MinimizationDofLayout &layout,
                                                                    GeneralizedHessian &hessian,
                                                                    std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule &molecule = moleculeData[moleculeIndex];
    if (components[molecule.componentId].rigid)
    {
      continue;
    }

    const Potentials::IntraMolecularPotentials &potentials = components[molecule.componentId].intraMolecularPotentials;
    std::span<const Atom> atom_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};

    for (const BendTorsionPotential &term : potentials.bendTorsions)
    {
      const std::size_t A = term.identifiers[0];
      const std::size_t B = term.identifiers[1];
      const std::size_t C = term.identifiers[2];
      const std::size_t D = term.identifiers[3];

      auto [energy, gradient, strain] = term.potentialEnergyGradientStrain(
          atom_span[A].position, atom_span[B].position, atom_span[C].position, atom_span[D].position);
      energies.bendTorsion += energy;
      dynamics_span[A].gradient += gradient[0];
      dynamics_span[B].gradient += gradient[1];
      dynamics_span[C].gradient += gradient[2];
      dynamics_span[D].gradient += gradient[3];
      hessian.strainGradient() += strain;

      const InternalCoordinate theta1 =
          Minimization::angleInternalCoordinate(atom_span[A].position, atom_span[B].position, atom_span[C].position);
      const InternalCoordinate theta2 =
          Minimization::angleInternalCoordinate(atom_span[B].position, atom_span[C].position, atom_span[D].position);

      // Smoothed uses the signed angle phi; all other types are polynomials in cos(phi).
      const bool usesPhi = (term.type == BendTorsionType::Smoothed);
      const InternalCoordinate torsion =
          usesPhi ? Minimization::dihedralAngleInternalCoordinate(atom_span[A].position, atom_span[B].position,
                                                                  atom_span[C].position, atom_span[D].position)
                  : Minimization::dihedralInternalCoordinate(atom_span[A].position, atom_span[B].position,
                                                             atom_span[C].position, atom_span[D].position);

      const Dual<3> q1 = Dual<3>::variable(theta1.value, 0);
      const Dual<3> q2 = Dual<3>::variable(theta2.value, 1);
      const Dual<3> qt = Dual<3>::variable(torsion.value, 2);  // cos(phi) or phi
      const auto &p = term.parameters;
      Dual<3> U{};
      switch (term.type)
      {
        case BendTorsionType::CVFF:
        case BendTorsionType::CFF:
          U = p[0] * (q1 - p[1]) * (q2 - p[2]) * qt;
          break;
        case BendTorsionType::Smoothed:
          U = p[0] * (1.0 + cosDual(p[1] * qt - p[2])) * smoothingDual(q1) * smoothingDual(q2);
          break;
        case BendTorsionType::SmoothedThreeCosine:
        {
          const Dual<3> c2 = qt * qt;
          U = (0.5 * p[0] * (1.0 + qt) + p[1] * (1.0 - c2) + 0.5 * p[2] * (1.0 - 3.0 * qt + 4.0 * qt * c2)) *
              smoothingDual(q1) * smoothingDual(q2);
          break;
        }
        case BendTorsionType::Nicholas:
        {
          const Dual<3> c2 = qt * qt;
          U = (0.5 * p[0] * (1.0 + qt) + p[1] * (1.0 - c2) + 0.5 * p[2] * (1.0 - 3.0 * qt + 4.0 * qt * c2)) *
              smoothingDual(q1);
          break;
        }
        case BendTorsionType::SmoothedCFF:
        {
          const Dual<3> c2 = qt * qt;
          U = (p[0] * (1.0 - qt) + 2.0 * p[1] * (1.0 - c2) + p[2] * (1.0 + 3.0 * qt - 4.0 * qt * c2)) *
              smoothingDual(q1) * smoothingDual(q2);
          break;
        }
        case BendTorsionType::SmoothedCFF2:
          U = (p[0] * (1.0 + qt) + p[2] + qt * (-3.0 * p[2] + 2.0 * qt * (p[1] + 2.0 * p[2] * qt))) *
              smoothingDual(q1) * smoothingDual(q2);
          break;
        case BendTorsionType::SmoothedCFF3:
          U = p[0] * (q1 - p[1]) * (q2 - p[2]) * qt * smoothingDual(q1) * smoothingDual(q2);
          break;
      }

      std::array<double, 3> dUdq{};
      std::array<std::array<double, 3>, 3> d2{};
      extractDerivatives(U, dUdq, d2);

      scatterCrossHessian(
          hessian, layout, moleculeIndex, {A, B, C, D}, {theta1, theta2, torsion},
          {{{0, 1, 2}}, {{1, 2, 3}}, {{0, 1, 2, 3}}},
          {atom_span[A].position, atom_span[B].position, atom_span[C].position, atom_span[D].position}, dUdq, d2);
    }
  }

  return energies;
}

RunningEnergy Interactions::computeIntraMolecularInversionBendHessian(std::span<const Molecule> moleculeData,
                                                                      std::span<const Atom> atoms,
                                                                      std::span<const Component> components,
                                                                      const MinimizationDofLayout &layout,
                                                                      GeneralizedHessian &hessian,
                                                                      std::span<AtomDynamics> dynamics)
{
  RunningEnergy energies{};

  for (std::size_t moleculeIndex = 0; moleculeIndex < moleculeData.size(); ++moleculeIndex)
  {
    const Molecule &molecule = moleculeData[moleculeIndex];
    if (components[molecule.componentId].rigid)
    {
      continue;
    }

    const Potentials::IntraMolecularPotentials &potentials = components[molecule.componentId].intraMolecularPotentials;
    std::span<const Atom> atom_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};

    for (const InversionBendPotential &term : potentials.inversionBends)
    {
      const std::size_t A = term.identifiers[0];
      const std::size_t B = term.identifiers[1];
      const std::size_t C = term.identifiers[2];
      const std::size_t D = term.identifiers[3];

      auto [energy, gradient, strain] = term.potentialEnergyGradientStrain(
          atom_span[A].position, atom_span[B].position, atom_span[C].position, atom_span[D].position);
      energies.inversionBend += energy;
      dynamics_span[A].gradient += gradient[0];
      dynamics_span[B].gradient += gradient[1];
      dynamics_span[C].gradient += gradient[2];
      dynamics_span[D].gradient += gradient[3];
      hessian.strainGradient() += strain;

      // Differentiate the full Cartesian inversion-bend energy through 12 hyper-dual variables
      // (three per atom); the inversion angle chi is not a distance/angle/dihedral coordinate.
      const std::array<double3, 4> positions = {atom_span[A].position, atom_span[B].position, atom_span[C].position,
                                                atom_span[D].position};
      auto seed = [&](std::size_t atomSlot) -> DualVec<12>
      {
        return {Dual<12>::variable((&positions[atomSlot].x)[0], atomSlot * 3 + 0),
                Dual<12>::variable((&positions[atomSlot].x)[1], atomSlot * 3 + 1),
                Dual<12>::variable((&positions[atomSlot].x)[2], atomSlot * 3 + 2)};
      };
      const DualVec<12> posA = seed(0);
      const DualVec<12> posB = seed(1);
      const DualVec<12> posC = seed(2);
      const DualVec<12> posD = seed(3);

      const DualVec<12> Rab = posA - posB;
      const DualVec<12> Rbc = posC - posB;
      const DualVec<12> Rbd = posD - posB;
      const DualVec<12> Rcd = posD - posC;
      const DualVec<12> Rad = posD - posA;

      const bool useAcdPlane = term.type == InversionBendType::Harmonic2 ||
                               term.type == InversionBendType::HarmonicCosine2 ||
                               term.type == InversionBendType::Planar2 || term.type == InversionBendType::MM3;

      Dual<12> c{};
      if (useAcdPlane)
      {
        const Dual<12> temp = dotDual(Rad, Rcd);
        c = dotDual(Rcd, Rcd) * dotDual(Rad, Rad) - temp * temp;
      }
      else
      {
        const Dual<12> temp = dotDual(Rbc, Rbd);
        c = dotDual(Rbc, Rbc) * dotDual(Rbd, Rbd) - temp * temp;
      }

      const Dual<12> rab2 = dotDual(Rab, Rab);
      const Dual<12> e = dotDual(Rab, crossDual(Rbd, Rbc));
      const Dual<12> cosChi = sqrtDual(rab2 - (e * e) / c) / sqrtDual(rab2);

      const auto &p = term.parameters;
      Dual<12> U{};
      switch (term.type)
      {
        case InversionBendType::Harmonic:
        case InversionBendType::Harmonic2:
        {
          const Dual<12> t = acosDual(cosChi) - p[1];
          U = 0.5 * p[0] * t * t;
          break;
        }
        case InversionBendType::HarmonicCosine:
        case InversionBendType::HarmonicCosine2:
        {
          const Dual<12> t = cosChi - p[1];
          U = 0.5 * p[0] * t * t;
          break;
        }
        case InversionBendType::Planar:
        case InversionBendType::Planar2:
          U = p[0] * (1.0 - cosChi);
          break;
        case InversionBendType::MM3:
        {
          const Dual<12> t = (acosDual(cosChi) - p[1]) * Units::RadiansToDegrees;
          const Dual<12> t2 = t * t;
          U = p[0] * t2 *
              (1.0 - 0.014 * t + 5.6e-5 * t2 - 7.0e-7 * (t * t2) + 2.2e-8 * (t2 * t2));
          break;
        }
      }

      std::array<double3, 4> analyticGradient{};
      std::array<std::array<CoordinateBlock, 4>, 4> analyticHessian{};
      for (std::size_t i = 0; i < 4; ++i)
      {
        analyticGradient[i] = {U.g[i * 3 + 0], U.g[i * 3 + 1], U.g[i * 3 + 2]};
        for (std::size_t j = 0; j < 4; ++j)
        {
          for (std::size_t a = 0; a < 3; ++a)
          {
            for (std::size_t b = 0; b < 3; ++b)
            {
              analyticHessian[i][j][a][b] = U.h[i * 3 + a][j * 3 + b];
            }
          }
        }
      }

      scatterCartesianFourBody(hessian, layout, moleculeIndex, {A, B, C, D}, positions, analyticGradient,
                               analyticHessian);
    }
  }

  return energies;
}
