module;

export module potential_hessian_vdw;

import std;

import double4;

import vdwparameters;
import forcefield;
import hessian_factor;
import potential_pair_derivatives;
import potential_pair_vdw;

export namespace Potentials
{
/**
 * \brief Computes the energy, gradient, and Hessian of the van der Waals (VDW) potential.
 *
 * Thin wrapper around the unified potentialVDW<2> evaluator; see potential_pair_vdw for the
 * implementation shared with the energy and gradient entry points.
 *
 * It returns derivatives as D[U[r], r] / r to avoid computing the square root for
 * Lennard-Jones (LJ) potential, as only the squared distance (rr) is required.
 *
 * The returned HessianFactor.dUdlambda holds the symmetric derivative factor X such that
 *   dU/d(scalingA) = scalingB * X   and   dU/d(scalingB) = scalingA * X.
 *
 * \param forcefield The force field parameters defining the interaction.
 * \param scalingA Scaling factor for atom A.
 * \param scalingB Scaling factor for atom B.
 * \param rr The squared distance between the two atoms.
 * \param typeA The type identifier for atom A.
 * \param typeB The type identifier for atom B.
 *
 * \return A HessianFactor object containing the computed energy, gradient, and Hessian.
 */
[[clang::always_inline]] inline HessianFactor potentialVDWHessian(const ForceField& forcefield, const double& scalingA,
                                                                  const double& scalingB, const double& rr,
                                                                  const std::size_t& typeA, const std::size_t& typeB)
{
  const PairDerivatives<2> derivatives = potentialVDW<2>(forcefield, scalingA, scalingB, rr, typeA, typeB);
  return HessianFactor(derivatives.energy, derivatives.dUdlambda, derivatives.firstDerivativeFactor,
                       derivatives.secondDerivativeFactor);
};
}  // namespace Potentials
