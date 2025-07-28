module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <iostream>
#endif

export module potential_tricubic_derivative_lj;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double4;

import vdwparameters;
import forcefield;
import tricubic_derivative_factor;

export namespace Potentials
{
/**
 * \brief Computes the first, second, and third-derivatives of the van der Waals (VDW) potential.
 *
 * This function calculates the derivatives of the VDW potential between two atom types.
 * It returns D[U[r], r] / r to avoid computing the square root for Lennard-Jones (LJ) potential,
 * as only the squared distance (rr) is required.
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
[[clang::always_inline]] inline TricubicDerivativeFactor potentialLennardJonesTricubicDerivative(
    const ForceField& forcefield, const double& r2, const std::size_t& typeA, const std::size_t& typeB)
{
  double arg1 = 4.0 * forcefield(typeA, typeB).parameters.x;
  double arg2 = forcefield(typeA, typeB).parameters.y * forcefield(typeA, typeB).parameters.y;
  double arg3 = forcefield(typeA, typeB).shift;
  double temp3 = (arg2 / r2) * (arg2 / r2) * (arg2 / r2);
  double r4 = r2 * r2;
  double r6 = r4 * r2;

  return TricubicDerivativeFactor(
      arg1 * (temp3 * temp3 - temp3) - arg3, -6.0 * arg1 * (2.0 * temp3 * temp3 - temp3) / r2,
      24.0 * arg1 * (7.0 * temp3 * temp3 - 2.0 * temp3) / r4, -96.0 * arg1 * (28.0 * temp3 * temp3 - 5.0 * temp3) / r6);
};

[[clang::always_inline]] inline TricubicDerivativeFactor potentialLennardJonesRepulsionTricubicDerivative(
    const ForceField& forcefield, const double& r2, const double& cutoff2, const std::size_t& typeA,
    const std::size_t& typeB)
{
  double arg1 = 4.0 * forcefield(typeA, typeB).parameters.x;
  double arg2 = forcefield(typeA, typeB).parameters.y * forcefield(typeA, typeB).parameters.y;
  double temp3 = (arg2 / r2) * (arg2 / r2) * (arg2 / r2);
  double temp3_rc = (arg2 / cutoff2) * (arg2 / cutoff2) * (arg2 / cutoff2);
  temp3_rc = 0.0;
  double r4 = r2 * r2;
  double r6 = r4 * r2;

  return TricubicDerivativeFactor(arg1 * temp3 * temp3 - temp3_rc * temp3_rc, -12.0 * arg1 * temp3 * temp3 / r2,
                                  168.0 * arg1 * temp3 * temp3 / r4, -2688.0 * arg1 * temp3 * temp3 / r6);
};

[[clang::always_inline]] inline TricubicDerivativeFactor potentialLennardJonesAttractionTricubicDerivative(
    const ForceField& forcefield, const double& r2, const double& cutoff2, const std::size_t& typeA,
    const std::size_t& typeB)
{
  double arg1 = 4.0 * forcefield(typeA, typeB).parameters.x;
  double arg2 = forcefield(typeA, typeB).parameters.y * forcefield(typeA, typeB).parameters.y;
  double temp3 = (arg2 / r2) * (arg2 / r2) * (arg2 / r2);
  double temp3_rc = (arg2 / cutoff2) * (arg2 / cutoff2) * (arg2 / cutoff2);
  temp3_rc = 0.0;
  double r4 = r2 * r2;
  double r6 = r4 * r2;

  return TricubicDerivativeFactor(arg1 * temp3 - temp3_rc, -6.0 * arg1 * temp3 / r2, 48.0 * arg1 * temp3 / r4,
                                  -480.0 * arg1 * temp3 / r6);
};

}  // namespace Potentials
