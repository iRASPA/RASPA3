module;

export module potential_pair_coulomb;

import std;

import double4;

import units;
import forcefield;
import potential_pair_derivatives;
import potential_coulomb_real_space;

namespace Potentials
{
/**
 * \brief Computes the Coulomb pair potential and its radial derivatives up to 'Order'.
 *
 * Single implementation for energy (Order 0), gradient (Order 1), and Hessian (Order 2)
 * evaluation, handling all supported charge methods: Ewald real-space, plain Coulomb, Wolf,
 * damped/modified shifted force, and zero-dipole summation. The Ewald real-space erfc terms
 * are computed inline (only up to the requested order); the remaining charge methods share
 * coulombRealSpaceFactors.
 *
 * The scaling is linear: it first switches the VDW interaction on in the range 0-0.5 and
 * then the electrostatics from 0.5 to 1.0.
 *
 * See PairDerivatives for the field conventions of the returned struct.
 *
 * \param forcefield The force field parameters, including charge method and Ewald alpha.
 * \param scalingA Scaling factor for atom A.
 * \param scalingB Scaling factor for atom B.
 * \param r Distance between the two atoms.
 * \param chargeA Electric charge of atom A.
 * \param chargeB Electric charge of atom B.
 *
 * \return A PairDerivatives<Order> object with the energy and requested derivative factors.
 */
export template <std::size_t Order>
[[clang::always_inline]] inline PairDerivatives<Order> potentialCoulomb(const ForceField& forcefield,
                                                                        const double scalingA, const double scalingB,
                                                                        const double r, const double chargeA,
                                                                        const double chargeB)
{
  static_assert(Order <= 3, "potentialCoulomb supports derivative orders 0, 1, 2, and 3");

  double scaling = scalingA * scalingB;

  switch (forcefield.chargeMethod)
  {
    [[likely]] case ForceField::ChargeMethod::Ewald:
    {
      // Offset form of the scaled Ewald real-space term (Brick-CFCMC, SI Eqs. (S13) and (S19)-(S20)):
      // U = lambda_t * C q_A q_B * erfc(alpha s)/s with s = r + Q, Q = delta * (1 - lambda_t) and
      // lambda_t = scalingA * scalingB. For full interactions Q = 0 and plain Ewald is recovered; for
      // scaled (fractional) interactions the potential stays finite even at r = 0, so no divergence
      // and no 0 * inf = NaN for zero Coulomb scaling. Spatial derivatives are taken at fixed lambda
      // (Q constant), so d^n U/dr^n = lambda_t * C q q * phi^(n)(s) with phi(s) = erfc(alpha s)/s.
      double alpha = forcefield.EwaldAlpha;
      double prefactor = Units::CoulombicConversionFactor * chargeA * chargeB;
      double s = r + EwaldChargeOffsetDelta * (1.0 - scaling);
      double inverseS = 1.0 / s;
      double inverseSS = inverseS * inverseS;
      double erfcTerm = std::erfc(alpha * s);
      double gaussian = std::exp(-alpha * alpha * s * s) * std::numbers::inv_sqrtpi_v<double>;

      double phi = erfcTerm * inverseS;
      double phiPrime = -(erfcTerm * inverseSS + 2.0 * alpha * gaussian * inverseS);

      // dU/d(scalingA) = scalingB * dUdlambda with the chain rule through Q (Eq. (S19)):
      // dUdlambda = C q q [ phi(s) - lambda_t * delta * phi'(s) ].
      double dUdlambda = prefactor * (phi - scaling * EwaldChargeOffsetDelta * phiPrime);

      if constexpr (Order == 0)
      {
        return {scaling * prefactor * phi, dUdlambda};
      }
      else
      {
        double inverseR = 1.0 / r;
        double firstDerivativeFactor = scaling * prefactor * phiPrime * inverseR;

        if constexpr (Order == 1)
        {
          return {scaling * prefactor * phi, dUdlambda, firstDerivativeFactor};
        }
        else
        {
          double alphaSquared = alpha * alpha;
          double phiSecond = 2.0 * erfcTerm * inverseSS * inverseS +
                             4.0 * alpha * gaussian * (inverseSS + alphaSquared);
          double secondDerivativeFactor =
              scaling * prefactor * (phiSecond - phiPrime * inverseR) * inverseR * inverseR;

          if constexpr (Order == 2)
          {
            return {scaling * prefactor * phi, dUdlambda, firstDerivativeFactor, secondDerivativeFactor};
          }
          else
          {
            double phiThird = -6.0 * erfcTerm * inverseSS * inverseSS -
                              2.0 * alpha * gaussian *
                                  (4.0 * alphaSquared * alphaSquared * s + 4.0 * alphaSquared * inverseS +
                                   6.0 * inverseS * inverseSS);
            double thirdDerivativeFactor = scaling * prefactor *
                                           (phiThird - 3.0 * phiSecond * inverseR + 3.0 * phiPrime * inverseR * inverseR) *
                                           inverseR * inverseR * inverseR;

            return {scaling * prefactor * phi, dUdlambda, firstDerivativeFactor, secondDerivativeFactor,
                    thirdDerivativeFactor};
          }
        }
      }
    }
    case ForceField::ChargeMethod::Coulomb:
    case ForceField::ChargeMethod::Wolf:
    case ForceField::ChargeMethod::DampedShiftedForce:
    case ForceField::ChargeMethod::ModifiedShiftedForce:
    case ForceField::ChargeMethod::ZeroDipole:
    {
      const CoulombRealSpaceFactors factors = coulombRealSpaceFactors(forcefield, r);
      const double prefactor = Units::CoulombicConversionFactor * chargeA * chargeB;

      if constexpr (Order == 0)
      {
        return {scaling * prefactor * factors.potential, prefactor * factors.potential};
      }
      else if constexpr (Order == 1)
      {
        return {scaling * prefactor * factors.potential, prefactor * factors.potential,
                scaling * prefactor * factors.firstDerivativeFactor};
      }
      else if constexpr (Order == 2)
      {
        return {scaling * prefactor * factors.potential, prefactor * factors.potential,
                scaling * prefactor * factors.firstDerivativeFactor,
                scaling * prefactor * factors.secondDerivativeFactor};
      }
      else
      {
        return {scaling * prefactor * factors.potential, prefactor * factors.potential,
                scaling * prefactor * factors.firstDerivativeFactor,
                scaling * prefactor * factors.secondDerivativeFactor,
                scaling * prefactor * factors.thirdDerivativeFactor};
      }
    }
  }

  std::unreachable();
}
}  // namespace Potentials
