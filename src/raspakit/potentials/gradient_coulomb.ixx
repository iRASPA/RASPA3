module;

export module potential_gradient_coulomb;

import std;

import double4;

import units;
import forcefield;
import gradient_factor;
import potential_pair_derivatives;
import potential_pair_coulomb;

export namespace Potentials
{
/**
 * \brief Computes the energy and gradient of the Coulomb potential.
 *
 * Thin wrapper around the unified potentialCoulomb<1> evaluator; see potential_pair_coulomb
 * for the implementation shared with the energy and Hessian entry points. It handles different
 * charge calculation methods such as Ewald, Coulomb, Wolf, shifted-force, and zero-dipole.
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
 * \return A GradientFactor object representing the computed force factors.
 *
 * \note This function returns D[U[r], r] / r, where U[r] is the potential energy.
 */
[[clang::always_inline]] inline GradientFactor potentialCoulombGradient(const ForceField& forcefield,
                                                                        const double& scalingA, const double& scalingB,
                                                                        const double& r, const double& chargeA,
                                                                        const double& chargeB)
{
  const PairDerivatives<1> derivatives = potentialCoulomb<1>(forcefield, scalingA, scalingB, r, chargeA, chargeB);
  return GradientFactor(derivatives.energy, derivatives.dUdlambda, derivatives.firstDerivativeFactor);
};
}  // namespace Potentials
