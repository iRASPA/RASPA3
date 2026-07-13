module;

export module minimization_generalized_coordinates;

import std;

import system;
import minimization_dof_layout;

/** Apply one local Cartesian/quaternion-tangent displacement to a System. */
export void applyGeneralizedDisplacement(System &system, const MinimizationDofLayout &layout,
                                         std::span<const double> displacement);
