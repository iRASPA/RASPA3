module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <iostream>
#include <numbers>
#endif

export module potential_electrostatics;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <numbers>;
import <iostream>;
#endif

import double4;

import units;
import forcefield;
import energy_factor;

export [[clang::always_inline]] inline double potentialElectrostatics(const ForceField& forcefield,
                                                                      const double& scalingB, const double& r,
                                                                      const double& chargeB)
{
  switch (forcefield.chargeMethod)
  {
    [[likely]] case ForceField::ChargeMethod::Ewald:
    {
      double alpha = forcefield.EwaldAlpha;
      double temp = Units::CoulombicConversionFactor * scalingB * chargeB * std::erfc(alpha * r) / r;
      return temp;
    }
    case ForceField::ChargeMethod::Coulomb:
    {
      return 0.0;
    }
    case ForceField::ChargeMethod::Wolf:
    {
      return 0.0;
    }
    case ForceField::ChargeMethod::ModifiedWolf:
    {
      return 0.0;
    }
    default:
      break;
  }
};
