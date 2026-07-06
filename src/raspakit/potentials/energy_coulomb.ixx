module;

export module potential_energy_coulomb;

import std;

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
 * The returned EnergyFactor.dUdlambda holds the symmetric derivative factor X such that
 *   dU/d(scalingA) = scalingB * X   and   dU/d(scalingB) = scalingA * X.
 * The caller routes these per-atom derivatives into the dU/dlambda accumulator of the
 * thermodynamic-integration group (Atom::groupId) each atom belongs to.
 *
 * \param forcefield The force field parameters, including charge method and Ewald alpha.
 * \param scalingA Scaling factor for atom A.
 * \param scalingB Scaling factor for atom B.
 * \param r Distance between the two atoms.
 * \param chargeA Electric charge of atom A.
 * \param chargeB Electric charge of atom B.
 * \return An EnergyFactor object containing the computed Coulomb energy and the derivative factor.
 */
[[clang::always_inline]] inline EnergyFactor potentialCoulombEnergy(const ForceField& forcefield,
                                                                    const double& scalingA, const double& scalingB,
                                                                    const double& r, const double& chargeA,
                                                                    const double& chargeB)
{
  double scaling = scalingA * scalingB;
  // scaling is linear and first switch LJ on in 0-0.5, then the electrostatics from 0.5 to 1.0
  switch (forcefield.chargeMethod)
  {
    [[likely]] case ForceField::ChargeMethod::Ewald:
    {
      double alpha = forcefield.EwaldAlpha;
      double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erfc(alpha * r) / r;
      return EnergyFactor(scaling * temp, temp);
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
