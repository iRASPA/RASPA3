module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <numbers>
#include <iostream>
#endif

export module potential_gradient_coulomb;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <numbers>;
import <iostream>;
#endif

import double4;

import units;
import forcefield;
import force_factor;

// return D[U[r],r] / r

export [[clang::always_inline]] inline ForceFactor 
potentialCoulombGradient(const ForceField& forcefield, const bool &groupIdA, const bool &groupIdB, 
                         const double &scalingA, const double &scalingB, const double& r, 
                         const double& chargeA, const double& chargeB)
{
  double scaling = scalingA * scalingB;
  switch(forcefield.chargeMethod)
  {
    [[likely]] case ForceField::ChargeMethod::Ewald:
    {
      double alpha = forcefield.EwaldAlpha;
      double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erfc(alpha * r) / r;
      ForceFactor result =  ForceFactor(scaling * temp,
            -Units::CoulombicConversionFactor *  scaling * chargeA * chargeB * 
            ((std::erfc(alpha * r) + 2.0 * alpha * r * std::exp(-alpha * alpha * r * r) * 
            std::numbers::inv_sqrtpi_v<double>) / (r * r * r)), 
            (groupIdA ? scalingB * temp : 0.0) + (groupIdB ? scalingA * temp : 0.0));
      return result;
    }
    case ForceField::ChargeMethod::Coulomb:
    {
      return ForceFactor(scaling * chargeA * chargeB / r, 0.0, 0.0);
    }
    case ForceField::ChargeMethod::Wolf:
    {
      return ForceFactor(scaling * chargeA * chargeB * std::erfc(forcefield.EwaldAlpha * r) / r, 0.0, 0.0);
    }
    case ForceField::ChargeMethod::ModifiedWolf:
    {
      return ForceFactor(scaling * chargeA * chargeB * std::erfc(forcefield.EwaldAlpha * r) / r, 0.0, 0.0);
    }
    default:
      break;
  }
};
