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
      double alpha = forcefield.EwaldAlpha;
      double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erfc(alpha * r) / r;

      if constexpr (Order == 0)
      {
        return {scaling * temp, temp};
      }
      else
      {
        double firstDerivativeFactor = -Units::CoulombicConversionFactor * scaling * chargeA * chargeB *
                                       ((std::erfc(alpha * r) + 2.0 * alpha * r * std::exp(-alpha * alpha * r * r) *
                                                                    std::numbers::inv_sqrtpi_v<double>) /
                                        (r * r * r));

        if constexpr (Order == 1)
        {
          return {scaling * temp, temp, firstDerivativeFactor};
        }
        else
        {
          double rr = r * r;
          double gaussian = std::exp(-alpha * alpha * rr) * std::numbers::inv_sqrtpi_v<double>;
          double secondDerivativeFactor =
              Units::CoulombicConversionFactor * scaling * chargeA * chargeB *
              (3.0 * std::erfc(alpha * r) / (rr * rr * r) + 4.0 * alpha * alpha * alpha * gaussian / rr +
               6.0 * alpha * gaussian / (rr * rr));

          if constexpr (Order == 2)
          {
            return {scaling * temp, temp, firstDerivativeFactor, secondDerivativeFactor};
          }
          else
          {
            // thirdDerivativeFactor = (1/r) d(secondDerivativeFactor)/dr for U = C q_A q_B erfc(alpha r)/r,
            // giving -C q_A q_B [ 15 erfc/r^7 + (2 alpha/sqrt(pi)) e^{-alpha^2 r^2} (4 alpha^4/r^2 + 10 alpha^2/r^4
            //                                                                      + 15/r^6) ].
            double thirdDerivativeFactor =
                -Units::CoulombicConversionFactor * scaling * chargeA * chargeB *
                (15.0 * std::erfc(alpha * r) / (rr * rr * rr * r) +
                 2.0 * alpha *
                     (4.0 * alpha * alpha * alpha * alpha / rr + 10.0 * alpha * alpha / (rr * rr) +
                      15.0 / (rr * rr * rr)) *
                     gaussian);

            return {scaling * temp, temp, firstDerivativeFactor, secondDerivativeFactor, thirdDerivativeFactor};
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
