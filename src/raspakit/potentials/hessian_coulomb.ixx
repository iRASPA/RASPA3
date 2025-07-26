module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <iostream>
#include <numbers>
#endif

export module potential_hessian_coulomb;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double4;

import units;
import forcefield;
import hessian_factor;

export namespace Potentials
{
/**
 * \brief Computes the gradient of the Coulomb potential.
 *
 * This function calculates the energy, gradient, and hessian of the Coulomb potential based
 * on the specified charge method in the provided force field. It handles different charge
 * method in the provided force field. It handles different charge calculation
 * methods such as Ewald, Coulomb, Wolf, and ModifiedWolf.
 *
 * \param forcefield The force field configuration containing charge method and parameters.
 * \param groupIdA Boolean indicating if the first group ID is active.
 * \param groupIdB Boolean indicating if the second group ID is active.
 * \param scalingA Scaling factor for the first interaction.
 * \param scalingB Scaling factor for the second interaction.
 * \param rr The distance squared between the two charges.
 * \param r The distance between the two charges.
 * \param chargeA The charge of the first particle.
 * \param chargeB The charge of the second particle.
 * \return A ForceFactor object representing the computed force factors.
 *
 * \note This function returns derivates as D[U[r], r] / r, where U[r] is the potential energy.
 */
[[clang::always_inline]] inline HessianFactor potentialCoulombHessian(const ForceField& forcefield,
                                                                      const bool& groupIdA, const bool& groupIdB,
                                                                      const double& scalingA, const double& scalingB,
                                                                      const double& rr, const double& r,
                                                                      const double& chargeA, const double& chargeB)
{
  double scaling = scalingA * scalingB;  ///< Combined scaling factor for interactions.

  switch (forcefield.chargeMethod)
  {
    [[likely]] case ForceField::ChargeMethod::Ewald:
    {
      double alpha = forcefield.EwaldAlpha;
      double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erfc(alpha * r) / r;
      return HessianFactor(
          scaling * temp, (groupIdA ? scalingB * temp : 0.0) + (groupIdB ? scalingA * temp : 0.0),
          -Units::CoulombicConversionFactor * scaling * chargeA * chargeB *
              ((std::erfc(alpha * r) +
                2.0 * alpha * r * std::exp(-alpha * alpha * rr) * std::numbers::inv_sqrtpi_v<double>) /
               (rr * r)),
          Units::CoulombicConversionFactor * scaling * chargeA * chargeB *
              (3.0 * std::erfc(alpha * r) / (rr * rr * r) +
               4.0 * alpha * alpha * alpha * std::exp(-alpha * alpha * rr) * std::numbers::inv_sqrtpi_v<double> / rr +
               6.0 * alpha * std::exp(-alpha * alpha * rr) * std::numbers::inv_sqrtpi_v<double> / (rr * rr)));
    }
    case ForceField::ChargeMethod::Coulomb:
    {
      return HessianFactor(Units::CoulombicConversionFactor * scaling * chargeA * chargeB / r,
                           -Units::CoulombicConversionFactor * scaling * chargeA * chargeB / (rr * r),
                           Units::CoulombicConversionFactor * 3.0 * scaling * chargeA * chargeB / (rr * rr * r),
                           -Units::CoulombicConversionFactor * 15.0 * scaling * chargeA * chargeB / (rr * rr * rr * r));
    }
    default:
      break;
  }

  // In case of an unsupported charge method, return a default ForceFactor.
  return HessianFactor(0.0, 0.0, 0.0, 0.0);
};
}  // namespace Potentials
