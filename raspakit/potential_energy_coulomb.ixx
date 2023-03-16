export module potential_energy_coulomb;

import forcefield;
import energy_factor;
import units;

import double4;
import <cmath>;
import <numbers>;

export inline EnergyFactor potentialCoulombEnergy(const ForceField& forcefield, const double& scaling, const double& r, const double& chargeA, const double& chargeB)
{
    switch(forcefield.chargeMethod)
    {
      default:
      case ForceField::ChargeMethod::Ewald:
      {
        EnergyFactor result(Units::CoulombicConversionFactor *  scaling * chargeA * chargeB * std::erfc(forcefield.EwaldAlpha * r) / r, 0.0);
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
