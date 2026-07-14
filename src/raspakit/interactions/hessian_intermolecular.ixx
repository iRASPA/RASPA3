module;

export module interactions_hessian_intermolecular;

import std;

import atom;
import atom_dynamics;
import molecule;
import component;
import forcefield;
import simulationbox;
import generalized_hessian;
import minimization_dof_layout;
import minimization_cell_layout;
import running_energy;

export namespace Interactions
{
RunningEnergy computeInterMolecularHessian(const ForceField& forceField, const SimulationBox& simulationBox,
                                           std::span<const Molecule> moleculeData,
                                           std::span<const Component> components, std::span<const Atom> moleculeAtoms,
                                           const MinimizationDofLayout& layout, GeneralizedHessian& hessian,
                                           std::span<AtomDynamics> dynamics = {},
                                           const CellMinimizationLayout& cellLayout = {});

}  // namespace Interactions
