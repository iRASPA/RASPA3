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
}  // namespace Detail

/**
 * \brief Computes the van der Waals pair potential and its radial derivatives up to 'Order'.
 *
 * Single implementation for energy (Order 0), gradient (Order 1), and Hessian (Order 2)
 * evaluation. Only the squared distance (rr) is required; the square root is avoided on the
 * Lennard-Jones hot path.
 *
 * The Lennard-Jones potential is dispatched with a single compare-and-branch. For gradient and
 * Hessian evaluation the two shifted Lennard-Jones variants keep an inlined fast path (used in
 * minimization hot loops). All other potential types are handled by the non-inlined
 * evaluateRareVDWDerivatives so that the hot code path stays free of code bloat; that function
 * computes energy, dU/dlambda, dU/dr, and d2U/dr2 for all remaining types in one pass.
 *
 * See PairDerivatives for the field conventions of the returned struct.
 *
 * \param forcefield The force field parameters defining the interaction.
 * \param scalingA Scaling factor for atom A.
 * \param scalingB Scaling factor for atom B.
 * \param rr The squared distance between the two atoms.
 * \param typeA The type identifier for atom A.
 * \param typeB The type identifier for atom B.
 *
 * \return A PairDerivatives<Order> object with the energy and requested derivative factors.
 */
export template <std::size_t Order>
[[clang::always_inline]] inline PairDerivatives<Order> potentialVDW(const ForceField& forcefield, const double scalingA,
                                                                    const double scalingB, const double rr,
                                                                    const std::size_t typeA, const std::size_t typeB)
{
  static_assert(Order <= 2, "potentialVDW supports derivative orders 0, 1, and 2");

  VDWParameters::Type potentialType = forcefield(typeA, typeB).type;

  double scaling = scalingA * scalingB;

  if (potentialType == VDWParameters::Type::LennardJones) [[likely]]
  {
    double arg1 = 4.0 * forcefield(typeA, typeB).parameters.x;
    double arg2 = forcefield(typeA, typeB).parameters.y * forcefield(typeA, typeB).parameters.y;
    double arg3 = forcefield(typeA, typeB).shift;
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

  if constexpr (Order >= 1)
  {
    if (potentialType == VDWParameters::Type::LennardJonesShiftedForce ||
        potentialType == VDWParameters::Type::LennardJonesSecondOrderTaylorShifted)
    {
      const VDWParameters& p = forcefield(typeA, typeB);
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

  Detail::RareVDWDerivatives derivatives =
      Detail::evaluateRareVDWDerivatives(forcefield, scalingA, scalingB, rr, typeA, typeB, potentialType);
  return Detail::vdwFromRareDerivatives<Order>(derivatives, rr);
}
}  // namespace Potentials
