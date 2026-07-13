module;

export module distance_potential_gradient_hessian_strain;

import std;

import double3;
import double3x3;
import bond_potential;

export namespace Potentials::Internal
{
/**
 * Energy, gradients, strain derivative, and distance Hessian factors for bond-like potentials.
 *
 * Returns f1 = (1/r) dU/dr and f2 = U''/r^2 - U'/r^3 used by RASPA2 Hessian assembly.
 */
std::tuple<double, std::array<double3, 2>, double3x3, double, double> distancePotentialEnergyGradientHessianStrain(
    BondType type, const std::array<double, maximumNumberOfBondParameters> &parameters, const double3 &posA,
    const double3 &posB, bool zero_gradient = false);
}  // namespace Potentials::Internal
