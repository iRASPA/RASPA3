module;

export module distance_potential_gradient_strain;

import std;

import double3;
import double3x3;
import bond_potential;

export namespace Potentials::Internal
{
/**
 * \brief Energy, atomic gradients, and strain derivative for bond-like distance potentials.
 *
 * Shared by bond and Urey-Bradley potentials. \p DF is dU/dr divided by r so that the
 * gradient contribution on atom A is DF * (posA - posB).
 */
std::tuple<double, std::array<double3, 2>, double3x3> distancePotentialEnergyGradientStrain(
    BondType type, const std::array<double, maximumNumberOfBondParameters> &parameters, const double3 &posA,
    const double3 &posB, bool zero_gradient = false);
}  // namespace Potentials::Internal
