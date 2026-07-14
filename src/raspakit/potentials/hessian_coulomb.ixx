module;

export module potential_hessian_coulomb;

import std;

import double4;

import units;
import forcefield;
import hessian_factor;
import potential_pair_derivatives;
import potential_pair_coulomb;

export namespace Potentials
{
/**
 * \brief Computes the energy, gradient, and Hessian of the Coulomb potential.
 *
 * Thin wrapper around the unified potentialCoulomb<2> evaluator; see potential_pair_coulomb
 * for the implementation shared with the energy and gradient entry points. It handles different
 * charge calculation methods such as Ewald, Coulomb, Wolf, shifted-force, and zero-dipole.
 *
 * The returned HessianFactor.dUdlambda holds the symmetric derivative factor X such that
 *   dU/d(scalingA) = scalingB * X   and   dU/d(scalingB) = scalingA * X.
 *
 * \param forcefield The force field configuration containing charge method and parameters.
 * \param scalingA Scaling factor for the first interaction.
 * \param scalingB Scaling factor for the second interaction.
 * \param rr The distance squared between the two charges (unused; kept for API compatibility).
 * \param r The distance between the two charges.
 * \param chargeA The charge of the first particle.
 * \param chargeB The charge of the second particle.
 * \return A HessianFactor object representing the computed force factors.
 *
 * \note This function returns derivatives as D[U[r], r] / r, where U[r] is the potential energy.
 */
[[clang::always_inline]] inline HessianFactor potentialCoulombHessian(const ForceField& forcefield,
                                                                      const double& scalingA, const double& scalingB,
                                                                      [[maybe_unused]] const double& rr,
                                                                      const double& r, const double& chargeA,
                                                                      const double& chargeB)
{
  const PairDerivatives<2> derivatives = potentialCoulomb<2>(forcefield, scalingA, scalingB, r, chargeA, chargeB);
  return HessianFactor(derivatives.energy, derivatives.dUdlambda, derivatives.firstDerivativeFactor,
                       derivatives.secondDerivativeFactor);
};
}  // namespace Potentials
