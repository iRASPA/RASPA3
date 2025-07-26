module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <iostream>
#include <numbers>
#endif

export module potential_tricubic_derivative_real_ewald;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <iostream>;
import <numbers>;
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
    const ForceField& forcefield, const double& rr, const double& r, const double& chargeA, const double& chargeB)
{
  double alpha = forcefield.EwaldAlpha;
  double exp_factor = std::erfc(alpha * r);
  double temp = Units::CoulombicConversionFactor * chargeA * chargeB * exp_factor / r;
  double exp_term_rr = 2.0 * alpha * r * std::exp(-alpha * alpha * rr) * std::numbers::inv_sqrtpi_v<double>;
  double alpha_r = alpha * r;
  double alpha_r_2 = alpha_r * alpha_r;
  double alpha_r_4 = alpha_r_2 * alpha_r_2;
  double r3 = rr * r;
  double r5 = r3 * rr;
  double r7 = r5 * rr;

  return TricubicDerivativeFactor(
      temp,

      -Units::CoulombicConversionFactor * chargeA * chargeB * ((exp_term_rr + exp_factor) / r3),

      Units::CoulombicConversionFactor * chargeA * chargeB *
          ((exp_term_rr * (2.0 * alpha_r_2 + 3.0) + 3.0 * exp_factor) / r5),

      -Units::CoulombicConversionFactor * chargeA * chargeB *
          ((exp_term_rr * (4.0 * alpha_r_4 + 10 * alpha_r_2 + 15.0) + 15.0 * exp_factor) / r7));
};
}  // namespace Potentials
