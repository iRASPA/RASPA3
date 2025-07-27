module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <iostream>
#endif

export module potential_gradient_vdw;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double4;

import vdwparameters;
import forcefield;
import gradient_factor;

export namespace Potentials
{
/**
 * \brief Computes the gradient of the van der Waals (VDW) potential.
 *
 * This function calculates the gradient of the VDW potential between two atom types.
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
 * \return A ForceFactor object containing the computed forces.
 */
[[clang::always_inline]] inline GradientFactor potentialVDWGradient(const ForceField& forcefield, const bool& groupIdA,
                                                                    const bool& groupIdB, const double& scalingA,
                                                                    const double& scalingB, const double& rr,
                                                                    const std::size_t& typeA, const std::size_t& typeB)
{
  VDWParameters::Type potentialType = forcefield(typeA, typeB).type;

  double scaling = scalingA * scalingB;
  switch (potentialType)
  {
    [[likely]] case VDWParameters::Type::LennardJones:
    {
      double arg1 = 4.0 * forcefield(typeA, typeB).parameters.x;
      double arg2 = forcefield(typeA, typeB).parameters.y * forcefield(typeA, typeB).parameters.y;
      double arg3 = forcefield(typeA, typeB).shift;
      double temp = (rr / arg2);
      double temp3 = temp * temp * temp;
      double inv_scaling = 1.0 - scaling;
      double rri3 = 1.0 / (temp3 + 0.5 * inv_scaling * inv_scaling);
      double rri6 = rri3 * rri3;
      double term = arg1 * (rri3 * (rri3 - 1.0)) - arg3;
      double dlambda_term = arg1 * scaling * inv_scaling * (2.0 * rri6 * rri3 - rri6);
      return GradientFactor(
          scaling * term,
          (groupIdA ? scalingB * (term + dlambda_term) : 0.0) + (groupIdB ? scalingA * (term + dlambda_term) : 0.0),
          12.0 * scaling * arg1 * (rri6 * temp3 * (0.5 - rri3)) / rr);
    }
    case VDWParameters::Type::RepulsiveHarmonic:
    {
      return GradientFactor(0.0, 0.0, 0.0);
    }
    default:
      return GradientFactor(0.0, 0.0, 0.0);
  }
};
}  // namespace Potentials
