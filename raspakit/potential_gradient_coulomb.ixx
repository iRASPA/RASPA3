export module potential_gradient_coulomb;

import <cmath>;
import <numbers>;

import double4;

import units;
import forcefield;
import force_factor;

// return D[U[r],r] / r

export inline ForceFactor potentialCoulombGradient(const ForceField& forcefield, const double& scaling, const double& r, const double& chargeA, const double& chargeB)
{
    switch(forcefield.chargeMethod)
    {
      default:
      case ForceField::ChargeMethod::Ewald:
      {
        double alpha = forcefield.EwaldAlpha;
        ForceFactor result =  ForceFactor(Units::CoulombicConversionFactor * scaling * chargeA * chargeB * std::erfc(alpha * r) / r,
              -Units::CoulombicConversionFactor *  scaling * chargeA * chargeB * ((std::erfc(alpha * r) + 2.0 * alpha * r * std::exp(-alpha * alpha * r * r) * std::numbers::inv_sqrtpi_v<double>) / (r * r * r)), 0.0);
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
    }
};
