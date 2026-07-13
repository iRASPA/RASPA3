module;

export module interactions_hessian_intermolecular;

import std;

import atom_dynamics;
import generalized_hessian;
import minimization_dof_layout;
import running_energy;
import system;

export namespace Interactions
{
RunningEnergy computeInterMolecularHessian(const System &system, const MinimizationDofLayout &layout,
                                           GeneralizedHessian &hessian, std::span<AtomDynamics> dynamics = {});

}  // namespace Interactions
