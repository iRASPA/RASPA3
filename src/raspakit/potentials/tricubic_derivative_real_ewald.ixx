module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <iostream>
#include <memory>
#include <numbers>
#endif

export module potential_tricubic_derivative_real_ewald;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <iostream>;
import <numbers>;
import <memory>;
#endif

import double4;

import units;
import vdwparameters;
import forcefield;
import tricubic_derivative_factor;

export namespace Potentials
{
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
 * \param rr The squared distance between the two atoms.
 * \param typeA The type identifier for atom A.
 * \param typeB The type identifier for atom B.
 *
 * \return A ThirdDerivativeFactor object containing the computed derivatives.
 */
[[clang::always_inline]] inline TricubicDerivativeFactor potentialRealEwaldTricubicDerivative(
    const ForceField& forceField, const double& rr, const double& r, const double& chargeA,
    const double& chargeB)
{
  double alpha = forceField.EwaldAlpha;
  double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erfc(alpha * r) / r;
  return TricubicDerivativeFactor(
      temp,
      -Units::CoulombicConversionFactor * chargeA * chargeB *
          ((std::erfc(alpha * r) +
            2.0 * alpha * r * std::exp(-alpha * alpha * rr) * std::numbers::inv_sqrtpi_v<double>) /
           (rr * r)),
      Units::CoulombicConversionFactor * chargeA * chargeB *
          (3.0 * std::erfc(alpha * r) / (rr * rr * r) +
           4.0 * alpha * alpha * alpha * std::exp(-alpha * alpha * rr) * std::numbers::inv_sqrtpi_v<double> / rr +
           6.0 * alpha * std::exp(-alpha * alpha * rr) * std::numbers::inv_sqrtpi_v<double> / (rr * rr)),
      -Units::CoulombicConversionFactor * chargeA * chargeB *
          (15.0 * std::erfc(alpha * r) / (rr * rr * rr * r) +
           30.0 * alpha * std::exp(-alpha * alpha * rr) * std::numbers::inv_sqrtpi_v<double> / (rr * rr * rr) +
           20.0 * alpha * alpha * alpha * std::exp(-alpha * alpha * rr) * std::numbers::inv_sqrtpi_v<double> /
               (rr * rr) +
           8.0 * alpha * alpha * alpha * alpha * alpha * std::exp(-alpha * alpha * rr) *
               std::numbers::inv_sqrtpi_v<double> / rr));
};
}  // namespace Potentials
