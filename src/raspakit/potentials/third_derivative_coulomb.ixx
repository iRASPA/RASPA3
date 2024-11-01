module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <iostream>
#include <numbers>
#endif

export module potential_third_derivative_coulomb;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <iostream>;
import <numbers>;
#endif

import double4;

import units;
import vdwparameters;
import forcefield;
import third_derivative_factor;

/**
 * \brief Computes the first, second, and third-derivatives of the Coulombic potential.
 *
 * This function calculates the derivatives of the VDW potential between two atom types.
 * It returns D[U[r], r] / r, D[D[U[r], r] / r, r] / r, 
 * D[D[D[U[r], r] / r, r] / r, r] / r, and D[D[D[D[U[r], r] / r, r] / r, r] / r, r] / r.
 *
 * \param forcefield The force field parameters defining the interaction.
 * \param groupIdA The group identifier for atom A.
 * \param groupIdB The group identifier for atom B.
 * \param scalingA Scaling factor for atom A.
 * \param scalingB Scaling factor for atom B.
 * \param rr The squared distance between the two atoms.
 * \param typeA The type identifier for atom A.
 * \param typeB The type identifier for atom B.
 *
 * \return A ThirdDerivativeFactor object containing the computed derivatives.
 */
export [[clang::always_inline]] inline ThirdDerivativeFactor potentialCoulombThirdDerivative(
    const ForceField& forcefield, const bool& groupIdA, const bool& groupIdB, const double& scalingA,
    const double& scalingB, const double &rr, const double& r, const double& chargeA, const double& chargeB)
{
  double scaling = scalingA * scalingB;  ///< Combined scaling factor for interactions.

  switch (forcefield.chargeMethod)
  {
    [[likely]] case ForceField::ChargeMethod::Ewald:
    {
      double alpha = forcefield.EwaldAlpha;
      double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erfc(alpha * r) / r;
      return ThirdDerivativeFactor(scaling * temp,
                                   (groupIdA ? scalingB * temp : 0.0) + (groupIdB ? scalingA * temp : 0.0),
                                   -Units::CoulombicConversionFactor * scaling * chargeA * chargeB *
                                   ((std::erfc(alpha * r) + 2.0 * alpha * r * std::exp(-alpha * alpha * rr) *
                                   std::numbers::inv_sqrtpi_v<double>) / (rr * r)),
                                   Units::CoulombicConversionFactor * scaling * chargeA * chargeB *
                                   (3.0 * std::erfc(alpha * r) / (rr * rr * r) +
                                    4.0 * alpha * alpha * alpha * std::exp(-alpha * alpha * rr) * std::numbers::inv_sqrtpi_v<double> / rr +
                                    6.0 * alpha * std::exp(-alpha * alpha * rr) * std::numbers::inv_sqrtpi_v<double> / (rr * rr )),
                                   -Units::CoulombicConversionFactor * scaling * chargeA * chargeB *
                                   (15.0 * std::erfc(alpha * r) / (rr * rr * rr * r) +
                                    30.0 * alpha * std::exp(-alpha * alpha * rr) * std::numbers::inv_sqrtpi_v<double> / (rr * rr  * rr) +
                                    20.0 * alpha * alpha * alpha * std::exp(-alpha * alpha * rr) * std::numbers::inv_sqrtpi_v<double> / (rr * rr ) +
                                    8.0 * alpha * alpha * alpha * alpha * alpha * std::exp(-alpha * alpha * rr) * std::numbers::inv_sqrtpi_v<double> / rr)
                                   );
    }
    case ForceField::ChargeMethod::Coulomb:
    {
      return ThirdDerivativeFactor(Units::CoulombicConversionFactor * scaling * chargeA * chargeB / r,
                                   -Units::CoulombicConversionFactor * scaling * chargeA * chargeB / (rr *r),
                                   Units::CoulombicConversionFactor * 3.0 * scaling * chargeA * chargeB / (rr * rr * r),
                                   -Units::CoulombicConversionFactor * 15.0 * scaling * chargeA * chargeB / (rr * rr * rr * r),
                                   Units::CoulombicConversionFactor * 105.0 * scaling * chargeA * chargeB / (rr * rr * rr * rr * r));
    }
    default:
      break;
  }

  // In case of an unsupported charge method, return a default ForceFactor.
  return ThirdDerivativeFactor(0.0, 0.0, 0.0, 0.0, 0.0);
};

