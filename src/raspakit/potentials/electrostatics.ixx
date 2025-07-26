module;

#ifdef USE_LEGACY_HEADERS
#include <cmath>
#include <cstddef>
#include <iostream>
#include <numbers>
#endif

export module potential_electrostatics;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double4;

import units;
import forcefield;
import energy_factor;

export namespace Potentials
{
/**
 * \brief Calculates the electrostatic potential based on the provided force field parameters.
 *
 * This function computes the electrostatic potential using different charge methods
 * defined in the forcefield. For the Ewald charge method, it calculates the potential
 * using the Ewald summation technique. For other charge methods like Coulomb, Wolf,
 * and ModifiedWolf, it currently returns 0.0.
 *
 * \param forcefield The force field parameters, including charge method and Ewald alpha.
 * \param scalingB Scaling factor for the electrostatic interaction.
 * \param r The distance between interacting particles.
 * \param chargeB The charge of the second particle.
 * \return The computed electrostatic potential.
 */
[[clang::always_inline]] inline double potentialElectrostatics(const ForceField& forcefield, const double& scalingB,
                                                               const double& r, const double& chargeB)
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
}  // namespace Potentials
