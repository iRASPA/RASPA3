module;

export module potential_electrostatics;

import std;

import double4;

import units;
import forcefield;
import potential_coulomb_real_space;

export namespace Potentials
{
/**
 * \brief Calculates the electrostatic potential based on the provided force field parameters.
 *
 * This function computes the electrostatic potential using different charge methods
 * defined in the forcefield. Ewald uses its real-space erfc term; finite-cutoff
 * methods use their corresponding shifted or multipole-neutralized potentials.
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
    case ForceField::ChargeMethod::Wolf:
    case ForceField::ChargeMethod::DampedShiftedForce:
    case ForceField::ChargeMethod::ModifiedShiftedForce:
    case ForceField::ChargeMethod::ZeroDipole:
    {
      const CoulombRealSpaceFactors factors = coulombRealSpaceFactors(forcefield, r);
      return Units::CoulombicConversionFactor * scalingB * chargeB * factors.potential;
    }
  }
  std::unreachable();
};
}  // namespace Potentials
