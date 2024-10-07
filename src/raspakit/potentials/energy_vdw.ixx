module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <iostream>
#endif

export module potential_energy_vdw;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <iostream>;
#endif

import vdwparameters;
import forcefield;
import energy_factor;

import double4;

export [[clang::always_inline]] inline EnergyFactor potentialVDWEnergy(const ForceField& forcefield,
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
      double temp = (rr / arg2);
      double temp3 = temp * temp * temp;
      double inv_scaling = 1.0 - scaling;
      double rri3 = 1.0 / (temp3 + 0.5 * inv_scaling * inv_scaling);
      double rri6 = rri3 * rri3;
      double term = arg1 * (rri3 * (rri3 - 1.0)) - arg3;
      double dlambda_term = arg1 * scaling * inv_scaling * (2.0 * rri6 * rri3 - rri6);
      return EnergyFactor(scaling * term, (groupIdA ? scalingB * (term + dlambda_term) : 0.0) +
                                              (groupIdB ? scalingA * (term + dlambda_term) : 0.0));
    }
    case VDWParameters::Type::BuckingHam:
    {
      double arg1 = forcefield(typeA, typeB).parameters.x;
      double arg2 = forcefield(typeA, typeB).parameters.y * forcefield(typeA, typeB).parameters.y;
      double arg3 = forcefield(typeA, typeB).shift;
      double temp = (rr / arg2);
      double rri3 = 1.0 / ((temp * temp * temp) + 0.5 * (1.0 - scaling) * (1.0 - scaling));
      return EnergyFactor(scaling * (4.0 * arg1 * (rri3 * (rri3 - 1.0)) - arg3), 0.0);
    }
    case VDWParameters::Type::Morse:
    {
      double arg1 = forcefield(typeA, typeB).parameters.x;
      double arg2 = forcefield(typeA, typeB).parameters.y * forcefield(typeA, typeB).parameters.y;
      double arg3 = forcefield(typeA, typeB).shift;
      double temp = (rr / arg2);
      double rri3 = 1.0 / ((temp * temp * temp) + 0.5 * (1.0 - scaling) * (1.0 - scaling));
      return EnergyFactor(scaling * (4.0 * arg1 * (rri3 * (rri3 - 1.0)) - arg3), 0.0);
    }
    case VDWParameters::Type::FeynmannHibbs:
    {
      double arg1 = forcefield(typeA, typeB).parameters.x;
      double arg2 = forcefield(typeA, typeB).parameters.y * forcefield(typeA, typeB).parameters.y;
      double arg3 = forcefield(typeA, typeB).shift;
      double temp = (rr / arg2);
      double rri3 = 1.0 / ((temp * temp * temp) + 0.5 * (1.0 - scaling) * (1.0 - scaling));
      return EnergyFactor(scaling * (4.0 * arg1 * (rri3 * (rri3 - 1.0)) - arg3), 0.0);
    }
    case VDWParameters::Type::MM3:
    {
      double arg1 = forcefield(typeA, typeB).parameters.x;
      double arg2 = forcefield(typeA, typeB).parameters.y * forcefield(typeA, typeB).parameters.y;
      double arg3 = forcefield(typeA, typeB).shift;
      double temp = (rr / arg2);
      double rri3 = 1.0 / ((temp * temp * temp) + 0.5 * (1.0 - scaling) * (1.0 - scaling));
      return EnergyFactor(scaling * (4.0 * arg1 * (rri3 * (rri3 - 1.0)) - arg3), 0.0);
    }
    case VDWParameters::Type::BornHugginsMeyer:
    {
      double arg1 = forcefield(typeA, typeB).parameters.x;
      double arg2 = forcefield(typeA, typeB).parameters.y * forcefield(typeA, typeB).parameters.y;
      double arg3 = forcefield(typeA, typeB).shift;
      double temp = (rr / arg2);
      double rri3 = 1.0 / ((temp * temp * temp) + 0.5 * (1.0 - scaling) * (1.0 - scaling));
      return EnergyFactor(scaling * (4.0 * arg1 * (rri3 * (rri3 - 1.0)) - arg3), 0.0);
    }
    default:
      break;
  }
};
