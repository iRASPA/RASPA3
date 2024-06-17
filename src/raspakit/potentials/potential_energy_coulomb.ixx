module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <numbers>
#endif

export module potential_energy_coulomb;

#ifndef USE_LEGACY_HEADERS
import <cmath>;
import <numbers>;
#endif

import forcefield;
import energy_factor;
import units;

import double4;

export [[clang::always_inline]] inline EnergyFactor potentialCoulombEnergy(const ForceField& forcefield,
                                                                           const bool& groupIdA, const bool& groupIdB,
                                                                           const double& scalingA,
                                                                           const double& scalingB, const double& r,
                                                                           const double& chargeA, const double& chargeB)
{
  double scaling = scalingA * scalingB;
  switch (forcefield.chargeMethod)
  {
    [[likely]] case ForceField::ChargeMethod::Ewald:
    {
      double alpha = forcefield.EwaldAlpha;
      double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erfc(alpha * r) / r;
      EnergyFactor result =
          EnergyFactor(scaling * temp, (groupIdA ? scalingB * temp : 0.0) + (groupIdB ? scalingA * temp : 0.0));
      return result;
    }
    case ForceField::ChargeMethod::Coulomb:
    {
      return EnergyFactor(scaling * chargeA * chargeB / r, 0.0);
    }
    case ForceField::ChargeMethod::Wolf:
    {
      return EnergyFactor(scaling * chargeA * chargeB * std::erfc(forcefield.EwaldAlpha * r) / r, 0.0);
    }
    case ForceField::ChargeMethod::ModifiedWolf:
    {
      return EnergyFactor(scaling * chargeA * chargeB * std::erfc(forcefield.EwaldAlpha * r) / r, 0.0);
    }
    default:
      break;
  }
};
