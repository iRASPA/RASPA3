module;

export module potential_pair_vdw;

import std;

import double4;

import vdwparameters;
import forcefield;
import potential_pair_derivatives;
import potential_vdw_rare_derivatives;

namespace Potentials
{
namespace Detail
{
/**
 * \brief Converts radial derivatives dU/dr and d2U/dr2 to the caller-facing factor convention.
 *
 * The pair-loop callers consume (1/r) dU/dr and (U'' - U'/r) / r^2 so that forces and Hessian
 * blocks can be assembled from dr without recomputing the square root.
 */
template <std::size_t Order>
[[clang::always_inline]] inline PairDerivatives<Order> vdwFromRareDerivatives(const RareVDWDerivatives& derivatives,
                                                                              double rr)
{
  if constexpr (Order == 0)
  {
    return PairDerivatives<0>{derivatives.energy, derivatives.dUdlambda};
  }
  else
  {
    double firstDerivativeFactor = rr > 0.0 ? derivatives.radialFirstDerivative / std::sqrt(rr) : 0.0;
    if constexpr (Order == 1)
    {
      return PairDerivatives<1>{derivatives.energy, derivatives.dUdlambda, firstDerivativeFactor};
    }
    else
    {
      double secondDerivativeFactor =
          rr > 0.0 ? (derivatives.radialSecondDerivative - firstDerivativeFactor) / rr : 0.0;
      return PairDerivatives<2>{derivatives.energy, derivatives.dUdlambda, firstDerivativeFactor,
                                secondDerivativeFactor};
    }
  }
}

/**
 * \brief Builds the caller-facing derivative struct from the energy and the radial derivatives at full coupling.
 *
 * dU/dlambda of a fully coupled pair is the pair energy itself (the soft-core term vanishes at lambda = 1).
 */
template <std::size_t Order>
[[clang::always_inline]] inline PairDerivatives<Order> fromRadial(double energy, double firstRadial,
                                                                  double secondRadial, double r, double rr)
{
  if constexpr (Order == 0)
  {
    return PairDerivatives<0>{energy, energy};
  }
  else
  {
    double firstDerivativeFactor = firstRadial / r;
    if constexpr (Order == 1)
    {
      return PairDerivatives<1>{energy, energy, firstDerivativeFactor};
    }
    else
    {
      return PairDerivatives<2>{energy, energy, firstDerivativeFactor, (secondRadial - firstDerivativeFactor) / rr};
    }
  }
}

/**
 * \brief The switched Lennard-Jones forms at full coupling (VDWParameters::Type::LennardJonesSwitched and
 * LennardJonesForceSwitched); see VDWParameters::potentialEnergyAtFullCoupling for the formulas.
 *
 * Below the switching distance the potential is plain Lennard-Jones (minus a constant for the force switch) and
 * the square root is avoided; the switching region [r_s, rc] needs r.
 */
template <std::size_t Order>
[[clang::always_inline]] inline PairDerivatives<Order> switchedLennardJones(const VDWParameters& p, const double rr,
                                                                            bool forceSwitch)
{
  double eps4 = 4.0 * p.parameters.x;
  double sigma2 = p.parameters.y * p.parameters.y;
  double temp = rr / sigma2;
  double rri3 = 1.0 / (temp * temp * temp);
  double rri6 = rri3 * rri3;
  double rc = p.parameters2.x;
  double rs = p.parameters2.y;
  double energy = eps4 * (rri6 - rri3);
  double firstFactor = 12.0 * eps4 * rri3 * (0.5 - rri3) / rr;  // (1/r) dU/dr of Lennard-Jones

  if (rr <= rs * rs) [[likely]]
  {
    if (forceSwitch)
    {
      // constant offset that makes the energy continuous with the switching region
      double q = p.parameters2.z;
      double invRc3 = 1.0 / (rc * rc * rc);
      double invRs3 = invRc3 / q;
      double sigma6 = sigma2 * sigma2 * sigma2;
      double c6 = eps4 * sigma6;
      energy -= c6 * (sigma6 * invRc3 * invRc3 * invRs3 * invRs3 - invRc3 * invRs3);
    }
    if constexpr (Order == 0)
    {
      return PairDerivatives<0>{energy, energy};
    }
    else if constexpr (Order == 1)
    {
      return PairDerivatives<1>{energy, energy, firstFactor};
    }
    else
    {
      return PairDerivatives<2>{energy, energy, firstFactor, 24.0 * eps4 * rri3 * (7.0 * rri3 - 2.0) / (rr * rr)};
    }
  }
  if (rr >= rc * rc) return {};

  double r = std::sqrt(rr);
  if (!forceSwitch)
  {
    double inverseWidth = p.parameters2.z;
    double x = (r - rs) * inverseWidth;
    double x2 = x * x;
    double oneMinusX = 1.0 - x;
    double s = 1.0 + x2 * x * (-10.0 + x * (15.0 - 6.0 * x));
    double ds = -30.0 * x2 * oneMinusX * oneMinusX * inverseWidth;
    double dds = -60.0 * x * oneMinusX * (1.0 - 2.0 * x) * inverseWidth * inverseWidth;
    double first = firstFactor * r;                                              // dU/dr
    double second = 24.0 * eps4 * rri3 * (7.0 * rri3 - 2.0) / rr + firstFactor;  // d2U/dr2
    return fromRadial<Order>(energy * s, first * s + energy * ds, second * s + 2.0 * first * ds + energy * dds, r,
                             rr);
  }

  double q = p.parameters2.z;
  double a12 = 1.0 / (1.0 - q * q);
  double a6 = 1.0 / (1.0 - q);
  double sigma6 = sigma2 * sigma2 * sigma2;
  double c6 = eps4 * sigma6;
  double c12 = c6 * sigma6;
  double invRc3 = 1.0 / (rc * rc * rc);
  double invR3 = 1.0 / (r * rr);
  double invR4 = invR3 / r;
  double invR5 = invR3 / rr;
  double d3 = invR3 - invRc3;
  double d6 = invR3 * invR3 - invRc3 * invRc3;
  double switchedEnergy = c12 * a12 * d6 * d6 - c6 * a6 * d3 * d3;
  double first = -12.0 * c12 * a12 * d6 * invR3 * invR4 + 6.0 * c6 * a6 * d3 * invR4;
  double second = 2.0 * c12 * a12 * (36.0 * invR3 * invR3 * invR4 * invR4 + 42.0 * d6 * invR4 * invR4) -
                  2.0 * c6 * a6 * (9.0 * invR4 * invR4 + 12.0 * d3 * invR5);
  return fromRadial<Order>(switchedEnergy, first, second, r, rr);
}
}  // namespace Detail

/**
 * \brief Computes the van der Waals pair potential and its radial derivatives up to 'Order'.
 *
 * Single implementation for energy (Order 0), gradient (Order 1), and Hessian (Order 2)
 * evaluation. Only the squared distance (rr) is required; the square root is avoided on the
 * Lennard-Jones hot path.
 *
 * The Lennard-Jones potential is dispatched with a single compare-and-branch, and a pair without
 * van der Waals interaction (Type::None) returns zeros immediately. For gradient and
 * Hessian evaluation the two shifted Lennard-Jones variants keep an inlined fast path (used in
 * minimization hot loops). All other potential types are handled by the non-inlined
 * evaluateRareVDWDerivatives so that the hot code path stays free of code bloat; that function
 * computes energy, dU/dlambda, dU/dr, and d2U/dr2 for all remaining types in one pass.
 *
 * See PairDerivatives for the field conventions of the returned struct.
 *
 * \param p The pair parameters (an entry of the force-field pair table, or of its 1-4 table).
 * \param scalingA Scaling factor for atom A.
 * \param scalingB Scaling factor for atom B.
 * \param rr The squared distance between the two atoms.
 *
 * \return A PairDerivatives<Order> object with the energy and requested derivative factors.
 */
export template <std::size_t Order>
[[clang::always_inline]] inline PairDerivatives<Order> potentialVDW(const VDWParameters& p, const double scalingA,
                                                                    const double scalingB, const double rr)
{
  static_assert(Order <= 2, "potentialVDW supports derivative orders 0, 1, and 2");

  VDWParameters::Type potentialType = p.type;

  double scaling = scalingA * scalingB;

  if (potentialType == VDWParameters::Type::LennardJones) [[likely]]
  {
    double arg1 = 4.0 * p.parameters.x;
    double arg2 = p.parameters.y * p.parameters.y;
    double arg3 = p.shift;
    double temp = (rr / arg2);          // (r/sigma)^2
    double temp3 = temp * temp * temp;  // (r/sigma)^6
    double inv_scaling = 1.0 - scaling;
    double rri3 = 1.0 / (temp3 + 0.5 * inv_scaling * inv_scaling);  // 1.0 / [0.5 (1-l)^2 + (r/sigma)^6]
    double rri6 = rri3 * rri3;
    double term = arg1 * (rri3 * (rri3 - 1.0)) - arg3;
    double dlambda_term = arg1 * scaling * inv_scaling * (2.0 * rri6 * rri3 - rri6);

    if constexpr (Order == 0)
    {
      return {scaling * term, term + dlambda_term};
    }
    else if constexpr (Order == 1)
    {
      return {scaling * term, term + dlambda_term, 12.0 * scaling * arg1 * (rri6 * temp3 * (0.5 - rri3)) / rr};
    }
    else
    {
      return {scaling * term, term + dlambda_term, 12.0 * scaling * arg1 * (rri6 * temp3 * (0.5 - rri3)) / rr,
              24.0 * arg1 * scaling * rri6 * temp3 * (1.0 + rri3 * (temp3 * (-3.0 + 9.0 * rri3) - 2.0)) / (rr * rr)};
    }
  }

  // No van der Waals interaction between these pseudo-atom types (for example the charged, massless sites of
  // multi-site water models against everything): the pair loops still reach this point for every pair within
  // the cutoff, so return before the non-inlined dispatch of the rare potentials.
  if (potentialType == VDWParameters::Type::None)
  {
    return {};
  }

  if constexpr (Order >= 1)
  {
    if (potentialType == VDWParameters::Type::LennardJonesShiftedForce ||
        potentialType == VDWParameters::Type::LennardJonesSecondOrderTaylorShifted)
    {
      double eps4 = 4.0 * p.parameters.x;
      double sigma2 = p.parameters.y * p.parameters.y;
      double sigma6 = sigma2 * sigma2 * sigma2;
      double c6 = p.parameters2.x;
      double rc = p.parameters2.y;
      double linearCoefficient = 12.0 * c6 * c6 - 6.0 * c6;
      double quadraticCoefficient = potentialType == VDWParameters::Type::LennardJonesSecondOrderTaylorShifted
                                        ? 156.0 * c6 * c6 - 42.0 * c6
                                        : 0.0;
      double invScaling = 1.0 - scaling;
      double x6 = rr * rr * rr + 0.5 * invScaling * invScaling * p.parameters2.w;
      double rs = std::sqrt(std::cbrt(x6));
      double displacement = (rs - rc) / rc;
      double u6 = sigma6 / x6;
      double u12 = u6 * u6;
      double term = eps4 * (u12 - u6 - c6 * (c6 - 1.0) + linearCoefficient * displacement -
                            0.5 * quadraticCoefficient * displacement * displacement);
      double deriv = eps4 * (12.0 * u12 - 6.0 * u6 - linearCoefficient * rs / rc +
                             quadraticCoefficient * rs * (rs - rc) / (rc * rc));
      double dlambdaTerm = scaling * invScaling * p.parameters2.w * deriv / (6.0 * x6);
      double firstDerivativeFactor = -scaling * deriv * rr * rr / x6;

      if constexpr (Order == 1)
      {
        return {scaling * term, term + dlambdaTerm, firstDerivativeFactor};
      }
      else
      {
        double radialFirstFactor =
            -12.0 * u12 + 6.0 * u6 + linearCoefficient * rs / rc - quadraticCoefficient * rs * (rs - rc) / (rc * rc);
        double radialSecondFactor = 24.0 * u12 - 6.0 * u6 + linearCoefficient * rs / (6.0 * rc) -
                                    quadraticCoefficient * (2.0 * rs * rs - rc * rs) / (6.0 * rc * rc);
        double secondDerivativeFactor =
            2.0 * scaling * eps4 * rr / (x6 * x6) *
            ((2.0 * x6 - 3.0 * rr * rr * rr) * radialFirstFactor + 3.0 * rr * rr * rr * radialSecondFactor);

        return {scaling * term, term + dlambdaTerm, firstDerivativeFactor, secondDerivativeFactor};
      }
    }
  }

  // The switched Lennard-Jones forms at full coupling (the biomolecular force fields: no lambda-scaling); the
  // soft-core version goes through the rare path.
  if ((potentialType == VDWParameters::Type::LennardJonesSwitched ||
       potentialType == VDWParameters::Type::LennardJonesForceSwitched) &&
      scaling == 1.0)
  {
    return Detail::switchedLennardJones<Order>(p, rr, potentialType == VDWParameters::Type::LennardJonesForceSwitched);
  }

  Detail::RareVDWDerivatives derivatives =
      Detail::evaluateRareVDWDerivatives(p, scalingA, scalingB, rr, potentialType);
  return Detail::vdwFromRareDerivatives<Order>(derivatives, rr);
}

/**
 * \brief The van der Waals pair potential of the pseudo-atom types 'typeA' and 'typeB' of the force field (the
 * regular pair table); see the VDWParameters overload.
 */
export template <std::size_t Order>
[[clang::always_inline]] inline PairDerivatives<Order> potentialVDW(const ForceField& forcefield, const double scalingA,
                                                                    const double scalingB, const double rr,
                                                                    const std::size_t typeA, const std::size_t typeB)
{
  return potentialVDW<Order>(forcefield(typeA, typeB), scalingA, scalingB, rr);
}
}  // namespace Potentials
