export module potential_energy_coulomb;

import forcefield;

import double4;
import <cmath>;

export inline double potentialCoulombEnergy(const ForceField& forcefield, const double& scaling, const double& r, const double& chargeA, const double& chargeB)
{
    switch(forcefield.chargeMethod)
    {
      case ForceField::ChargeMethod::Ewald:
      {
		return scaling * chargeA * chargeB * std::erfc(forcefield.alpha * r) / r;
      }
      case ForceField::ChargeMethod::Coulomb:
      {
		return scaling * chargeA * chargeB / r;
      }
      case ForceField::ChargeMethod::Wolf:
      {
		return scaling * chargeA * chargeB * std::erfc(forcefield.alpha * r) / r;
      }
      case ForceField::ChargeMethod::ModifiedWolf:
      {
		return scaling * chargeA * chargeB * std::erfc(forcefield.alpha * r) / r;
      }
    }
};
