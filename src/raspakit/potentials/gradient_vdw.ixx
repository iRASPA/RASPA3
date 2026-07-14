module;

export module potential_gradient_vdw;

import std;

import double4;

import vdwparameters;
import forcefield;
import gradient_factor;
import potential_pair_derivatives;
import potential_pair_vdw;

export namespace Potentials
{
/**
 * \brief Computes the energy and gradient of the van der Waals (VDW) potential.
 *
 * Thin wrapper around the unified potentialVDW<1> evaluator; see potential_pair_vdw for the
 * implementation shared with the energy and Hessian entry points.
 *
 * It returns D[U[r], r] / r to avoid computing the square root for Lennard-Jones (LJ)
 * potential, as only the squared distance (rr) is required.
 *
 * The returned GradientFactor.dUdlambda holds the symmetric derivative factor X such that
 *   dU/d(scalingA) = scalingB * X   and   dU/d(scalingB) = scalingA * X.
 *
 * \param forcefield The force field parameters defining the interaction.
 * \param scalingA Scaling factor for atom A.
 * \param scalingB Scaling factor for atom B.
 * \param rr The squared distance between the two atoms.
 * \param typeA The type identifier for atom A.
 * \param typeB The type identifier for atom B.
 *
 * \return A GradientFactor object containing the computed energy and gradient factor.
 */
[[clang::always_inline]] inline GradientFactor potentialVDWGradient(const ForceField& forcefield,
                                                                    const double& scalingA, const double& scalingB,
                                                                    const double& rr, const std::size_t& typeA,
                                                                    const std::size_t& typeB)
{
  const PairDerivatives<1> derivatives = potentialVDW<1>(forcefield, scalingA, scalingB, rr, typeA, typeB);
  return GradientFactor(derivatives.energy, derivatives.dUdlambda, derivatives.firstDerivativeFactor);
};
}  // namespace Potentials
