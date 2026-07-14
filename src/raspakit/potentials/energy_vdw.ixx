module;

export module potential_energy_vdw;

import std;

import double4;

import vdwparameters;
import forcefield;
import energy_factor;
import potential_pair_derivatives;
import potential_pair_vdw;

export namespace Potentials
{
/**
 * \brief Computes the van der Waals (VDW) energy between two atoms.
 *
 * Thin wrapper around the unified potentialVDW<0> evaluator; see potential_pair_vdw for the
 * implementation shared with the gradient and Hessian entry points.
 *
 * The returned EnergyFactor.dUdlambda holds the symmetric derivative factor X such that
 *   dU/d(scalingA) = scalingB * X   and   dU/d(scalingB) = scalingA * X.
 * The caller routes these per-atom derivatives into the dU/dlambda accumulator of the
 * thermodynamic-integration group (Atom::groupId) each atom belongs to.
 *
 * \param forcefield The force field parameters defining the interaction.
 * \param scalingA Scaling factor for atom A.
 * \param scalingB Scaling factor for atom B.
 * \param rr The squared distance between the two atoms.
 * \param typeA The type identifier for atom A.
 * \param typeB The type identifier for atom B.
 *
 * \return An EnergyFactor object containing the calculated potential energy and the derivative factor.
 */
[[clang::always_inline]] inline EnergyFactor potentialVDWEnergy(const ForceField& forcefield, const double& scalingA,
                                                                const double& scalingB, const double& rr,
                                                                const std::size_t& typeA, const std::size_t& typeB)
{
  const PairDerivatives<0> derivatives = potentialVDW<0>(forcefield, scalingA, scalingB, rr, typeA, typeB);
  return EnergyFactor(derivatives.energy, derivatives.dUdlambda);
};
}  // namespace Potentials
