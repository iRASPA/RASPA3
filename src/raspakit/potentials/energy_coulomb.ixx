module;

export module potential_energy_coulomb;

import std;

import forcefield;
import energy_factor;
import potential_pair_derivatives;
import potential_pair_coulomb;

import double4;

export namespace Potentials
{
/**
 * \brief Calculates the Coulomb potential energy factor between two atoms.
 *
 * Thin wrapper around the unified potentialCoulomb<0> evaluator; see potential_pair_coulomb
 * for the implementation shared with the gradient and Hessian entry points. Supported methods
 * include Ewald, Coulomb, Wolf, damped/modified shifted force, and zero-dipole summation.
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
  const PairDerivatives<0> derivatives = potentialCoulomb<0>(forcefield, scalingA, scalingB, r, chargeA, chargeB);
  return EnergyFactor(derivatives.energy, derivatives.dUdlambda);
};
}  // namespace Potentials
