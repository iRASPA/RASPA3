export module potential_gradient_coulomb;

import forcefield;
import force_factor;

import double4;
import <cmath>;
import <numbers>;

export inline ForceFactor potentialCoulombGradient(const ForceField& forcefield, const double& scaling, const double& r, const double& chargeA, const double& chargeB)
{
    switch(forcefield.chargeMethod)
    {
      case ForceField::ChargeMethod::Ewald:
      {
        double alpha = forcefield.alpha;
		return ForceFactor(scaling * chargeA * chargeB * std::erfc(alpha * r) / r, 
              scaling * chargeA * chargeB * ((std::erfc(alpha * r) + 2.0 * alpha * r * std::exp(-alpha * alpha * r * r) * std::numbers::inv_sqrtpi_v<double>) / (r * r * r)), 0.0);

        //-COULOMBIC_CONVERSION_FACTOR*chargeA*chargeB*
        //     (erfc(Alpha[CurrentSystem]*r)+2.0*Alpha[CurrentSystem]*r*exp(-SQR(Alpha[CurrentSystem])*rr)/sqrt(M_PI))/
        //     (r*rr);
      }
      case ForceField::ChargeMethod::Coulomb:
      {
		return ForceFactor(scaling * chargeA * chargeB / r, 0.0, 0.0);
      }
      case ForceField::ChargeMethod::Wolf:
      {
		return ForceFactor(scaling * chargeA * chargeB * std::erfc(forcefield.alpha * r) / r, 0.0, 0.0);
      }
      case ForceField::ChargeMethod::ModifiedWolf:
      {
		return ForceFactor(scaling * chargeA * chargeB * std::erfc(forcefield.alpha * r) / r, 0.0, 0.0);
      }
    }
};
