module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
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

export namespace Potentials
{
/**
 * \brief Calculates the Coulomb potential energy factor between two atoms.
 *
 * This function computes the Coulombic energy between two atoms based on the specified
 * force field's charge method. It accounts for scaling factors and the distance between
 * atoms. Supported charge methods include Ewald, Coulomb, Wolf, and ModifiedWolf.
 *
 * \param forcefield The force field parameters, including charge method and Ewald alpha.
 * \param groupIdA Indicates if atom A belongs to a group affecting the scaling.
 * \param groupIdB Indicates if atom B belongs to a group affecting the scaling.
 * \param scalingA Scaling factor for atom A.
 * \param scalingB Scaling factor for atom B.
 * \param r Distance between the two atoms.
 * \param chargeA Electric charge of atom A.
 * \param chargeB Electric charge of atom B.
 * \return An EnergyFactor object containing the computed Coulomb energy and any group-related scaling.
 */
[[clang::always_inline]] inline EnergyFactor potentialCoulombEnergy(const ForceField& forcefield, const bool& groupIdA,
                                                                    const bool& groupIdB, const double& scalingA,
                                                                    const double& scalingB, const double& r,
                                                                    const double& chargeA, const double& chargeB)
{
  double scaling = scalingA * scalingB;
  // scaling is linear and first switch LJ on in 0-0.5, then the electrostatics from 0.5 to 1.0
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
    default:
      break;
  }

  return EnergyFactor(0.0, 0.0);
};
}  // namespace Potentials
