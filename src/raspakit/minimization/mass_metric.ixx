module;

export module minimization_mass_metric;

import std;

import system;
import minimization_dof_layout;

/**
 * Inverse-square-root mass metric in the generalized minimization DOF basis.
 *
 * Cartesian atom DOFs and rigid center-of-mass DOFs carry scalar weights 1/sqrt(m); rigid orientation
 * DOFs carry a 3x3 block, the pseudo-inverse square root of the space-frame inertia tensor. The metric is
 * used both by the Gamma-point normal-mode analysis and by the (k-dependent) generalized dynamical matrix,
 * so that at k = 0 the phonon eigenvalues reproduce the normal-mode frequencies exactly.
 */
export struct MassMetric
{
  // Zero entries mark DOFs that belong to a 3x3 orientation block.
  std::vector<double> scalarInverseSqrt;
  // Base DOF index plus row-major inverse-square-root inertia block.
  std::vector<std::pair<std::size_t, std::array<double, 9>>> blocks;
  std::size_t discardedRotationalDofs{};
};

/** Build the mass metric for `system` in the given DOF layout. */
export MassMetric buildMassMetric(const System& system, const MinimizationDofLayout& layout);

/** Apply the symmetric block-diagonal metric to a dense row-major matrix: D = W H W. */
export std::vector<double> applyMassMetric(std::span<const double> hessian, std::size_t size, const MassMetric& metric);

/** Apply the symmetric block-diagonal metric to a single vector: q = W v. */
export std::vector<double> metricTimesVector(const MassMetric& metric, std::span<const double> vector, std::size_t size);
