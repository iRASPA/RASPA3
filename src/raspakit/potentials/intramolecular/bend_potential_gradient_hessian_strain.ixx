module;

export module bend_potential_gradient_hessian_strain;

import std;

import double3;
import double3x3;
import bend_potential;
import minimization_bend_hessian_geometry;

export namespace Potentials::Internal
{
std::tuple<double, std::array<double3, 3>, double3x3, Minimization::BendHessianGeometry>
bendPotentialEnergyGradientHessianStrain(BendType type,
                                         const std::array<double, maximumNumberOfBendParameters> &parameters,
                                         const double3 &posA, const double3 &posB, const double3 &posC);
}  // namespace Potentials::Internal
