export module potential_energy_coulomb;

import forcefield;
import energy_factor;
import units;

import double4;
import <cmath>;
import <numbers>;

export inline EnergyFactor potentialCoulombEnergy(const ForceField& forcefield, const double& scalingA, const double& scalingB, const double& r, const double& chargeA, const double& chargeB)
{
  double scaling = scalingA * scalingB;
    switch(forcefield.chargeMethod)
    {
      default:
      case ForceField::ChargeMethod::Ewald:
      {
         double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erfc(forcefield.EwaldAlpha * r) / r;
         EnergyFactor result(scaling * temp, 
                             temp);
        return  result;
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
    }
};
