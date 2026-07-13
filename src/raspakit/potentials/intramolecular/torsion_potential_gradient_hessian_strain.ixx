module;

export module torsion_potential_gradient_hessian_strain;

import std;

import double3;
import double3x3;
import torsion_potential;
import minimization_torsion_hessian_geometry;

export namespace Potentials::Internal
{
/**
 * Energy, gradients, strain derivative, and Hessian geometry for dihedral potentials.
 *
 * Returns DF = dU/dCosPhi and DDF = d²U/dCosPhi² (RASPA2 conventions) inside the geometry,
 * together with the first-derivative vectors dtA..dtD needed to assemble the Hessian.
 */
std::tuple<double, std::array<double3, 4>, double3x3, Minimization::TorsionHessianGeometry>
torsionPotentialEnergyGradientHessianStrain(TorsionType type,
                                            const std::array<double, maximumNumberOfTorsionParameters> &parameters,
                                            const double3 &posA, const double3 &posB, const double3 &posC,
                                            const double3 &posD);
}  // namespace Potentials::Internal
