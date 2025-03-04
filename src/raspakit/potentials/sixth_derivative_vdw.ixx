module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <iostream>
#endif

export module potential_sixth_derivative_vdw;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <iostream>;
#endif

import double4;

import vdwparameters;
import forcefield;
import sixth_derivative_factor;

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
export [[clang::always_inline]] inline SixthDerivativeFactor 
  potentialVDWSixthDerivative(const ForceField& forcefield,
                              const bool& groupIdA, const bool& groupIdB,
                              const double& scalingA, const double& scalingB,
                              const double& rr, const size_t& typeA,
                              const size_t& typeB)
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
      double temp = (rr / arg2);              // (r/sigma)^2
      double temp3 = temp * temp * temp;      // (r/sigma)^6
      double inv_scaling = 1.0 - scaling;
      double rri3 = 1.0 / (temp3 + 0.5 * inv_scaling * inv_scaling);   // 1.0 / [0.5 (1-l)^2 + (r/sigma)^6]
      double rri6 = rri3 * rri3;
      double term = arg1 * (rri3 * (rri3 - 1.0)) - arg3;
      double dlambda_term = arg1 * scaling * inv_scaling * (2.0 * rri6 * rri3 - rri6);

      /*
      temp3 = (arg2 / rr) * (arg2 / rr) * (arg2 / rr);
      return SixthDerivativeFactor(
          arg1 * (temp3 * temp3 - temp3) - arg3,
          0.0,
          -6.0 * arg1 * (2.0 * temp3 * temp3 - temp3) / rr,
          24.0 * arg1 * (7.0 * temp3 * temp3 - 2.0 * temp3) / (rr * rr),
          -96.0 * arg1 * (28.0 * temp3 * temp3 - 5.0 * temp3) / (rr * rr * rr),
          1152.0 * arg1 * (42.0 * temp3 * temp3 - 5.0 * temp3) / (rr * rr * rr * rr),
          -80640.0 * arg1 * (12.0 * temp3 * temp3 - temp3) / (rr * rr * rr * rr * rr),
          645120.0 * arg1 * (33.0 * temp3 * temp3 - 2.0 * temp3) / (rr * rr * rr * rr * rr * rr));
          */


      return SixthDerivativeFactor(
          scaling * term, 
          (groupIdA ? scalingB * (term + dlambda_term) : 0.0) + (groupIdB ? scalingA * (term + dlambda_term) : 0.0),
          6.0 * arg1 * scaling * rri6 * temp3 * (1.0 - 2.0 * rri3) / rr,
          24.0 * arg1 * scaling * rri6 * temp3 * (1.0 + rri3 * (temp3 * (-3.0 + 9.0 * rri3) - 2.0)) / (rr * rr),
          48.0 * arg1 * scaling * rri6 * temp3 * (1.0 - rri3 * (2.0 + 9.0 * temp3 * (2.0 + 12.0 * rri6 * temp3 - 3.0 * rri3 * (2.0 + temp3)))) / (rr * rr * rr),
          1152.0 * arg1 * scaling * rri6 * rri3 * temp3 * temp3 * (-5.0 + 3.0 * rri3 * (5.0 + 9.0 * temp3 * (1.0 + 5.0 * rri6 * temp3 - rri3 * (4.0 + temp3)))) / (rr * rr * rr * rr),
          11520.0 * arg1 * scaling * rri6 * rri3 * temp3 * temp3 * (-2.0 + 6.0 * rri3 - 9.0 * rri3 * temp3 * (-2.0 + 3.0 * rri3 * temp3) * 
            (2.0 + rri3 * (-8.0 + 3.0 * (-1.0 + 6.0 * rri3) * temp3))) / (rr * rr * rr * rr * rr),
          46080.0 * arg1 * scaling * rri6 * rri3 * temp3 * temp3 * (-1.0 + 3.0 * rri3 * (1.0 + 9.0 * temp3 * (3.0 + rri3 * (-12.0 + temp3 * (-22.0 + rri3 * (110.0 + 9.0 * temp3 * (5.0 + 21.0 * rri6 * temp3 - 3.0 * rri3 * (10.0 + temp3)))))))) / (rr * rr * rr * rr * rr * rr)
          );
    }
    case VDWParameters::Type::RepulsiveHarmonic:
    {
      double r = std::sqrt(rr);
      double arg1 = forcefield(typeA, typeB).parameters.x;
      double arg2 = forcefield(typeA, typeB).parameters.y;
      double temp = (1.0 - r / arg2);
      return (r >= arg2) ? SixthDerivativeFactor(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) :
                            SixthDerivativeFactor(0.5 * arg1 * temp  * temp,
                                   0.0, 
                                  -arg1 * temp / (arg2 * r),
                                   arg1 / (rr * r * arg2), 
                                   -3.0 * arg1 / (rr * rr * r * arg2),
                                   15.0 * arg1 / (rr * rr * rr * r * arg2),
                                  -105.0 * arg1 / (rr * rr * rr * rr * r * arg2),
                                   945.0 * arg1 / (rr * rr * rr * rr * rr * r * arg2));
    }
    default:
      return SixthDerivativeFactor(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  }
};
