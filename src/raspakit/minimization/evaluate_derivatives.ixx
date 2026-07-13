module;

export module minimization_evaluate_derivatives;

import std;

import atom;
import atom_dynamics;
import molecule;
import component;
import system;
import running_energy;
import generalized_hessian;
import minimization_dof_layout;
import minimization_rigid_kinematics;
import interactions_hessian_intramolecular;
import interactions_hessian_intermolecular;
import interactions_hessian_framework_molecule;
import interactions_hessian_ewald;

export struct DerivativeCapabilities
{
  bool energy{true};
  bool gradient{false};
  bool hessianPositionPosition{false};
  bool hessianPositionStrain{false};
  bool hessianStrainStrain{false};
};

export struct DerivativeResults
{
  double energy{};
  std::span<double> gradient;
  GeneralizedHessian &hessian;
};

export void evaluateDerivatives(System &system, const MinimizationDofLayout &layout, DerivativeCapabilities capabilities,
                                DerivativeResults &results);
