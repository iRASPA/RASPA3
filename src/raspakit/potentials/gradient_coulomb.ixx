module;

export module potential_gradient_coulomb;

import std;

import double4;

import units;
import forcefield;
import gradient_factor;

export namespace Potentials
{
/**
 * \brief Computes the gradient of the Coulomb potential.
 *
 * This function calculates the gradient of the Coulomb potential based on the specified
 * charge method in the provided force field. It handles different charge calculation
 * methods such as Ewald, Coulomb, Wolf, and ModifiedWolf.
 *
 * The returned GradientFactor.dUdlambda holds the symmetric derivative factor X such that
 *   dU/d(scalingA) = scalingB * X   and   dU/d(scalingB) = scalingA * X.
 *
 * \param forcefield The force field configuration containing charge method and parameters.
 * \param scalingA Scaling factor for the first interaction.
 * \param scalingB Scaling factor for the second interaction.
 * \param r The distance between the two charges.
 * \param chargeA The charge of the first particle.
 * \param chargeB The charge of the second particle.
 * \return A ForceFactor object representing the computed force factors.
 *
 * \note This function returns D[U[r], r] / r, where U[r] is the potential energy.
 */
[[clang::always_inline]] inline GradientFactor potentialCoulombGradient(const ForceField& forcefield,
                                                                        const double& scalingA, const double& scalingB,
                                                                        const double& r, const double& chargeA,
                                                                        const double& chargeB)
{
  double scaling = scalingA * scalingB;  ///< Combined scaling factor for interactions.

  switch (forcefield.chargeMethod)
  {
    [[likely]] case ForceField::ChargeMethod::Ewald:
    {
      double alpha = forcefield.EwaldAlpha;
      double temp = Units::CoulombicConversionFactor * chargeA * chargeB * std::erfc(alpha * r) / r;
      return GradientFactor(scaling * temp, temp,
                            -Units::CoulombicConversionFactor * scaling * chargeA * chargeB *
                                ((std::erfc(alpha * r) + 2.0 * alpha * r * std::exp(-alpha * alpha * r * r) *
                                                             std::numbers::inv_sqrtpi_v<double>) /
                                 (r * r * r)));
    }
    case ForceField::ChargeMethod::Coulomb:
    {
      return GradientFactor(scaling * chargeA * chargeB / r, 0.0, 0.0);
    }
    default:
      break;
  }

  // In case of an unsupported charge method, return a default ForceFactor.
  return GradientFactor(0.0, 0.0, 0.0);
};
}  // namespace Potentials
