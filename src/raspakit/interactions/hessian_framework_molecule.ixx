module;

export module interactions_hessian_framework_molecule;

import std;

import atom;
import atom_dynamics;
import molecule;
import component;
import framework;
import forcefield;
import simulationbox;
import generalized_hessian;
import minimization_dof_layout;
import minimization_cell_layout;
import running_energy;

export namespace Interactions
{
RunningEnergy computeFrameworkMoleculeHessian(
    const ForceField& forceField, const SimulationBox& simulationBox, std::span<const Molecule> moleculeData,
    std::span<const Component> components, std::span<const Atom> frameworkAtoms, std::span<const Atom> moleculeAtoms,
    const MinimizationDofLayout& layout, GeneralizedHessian& hessian, std::span<AtomDynamics> moleculeDynamics = {},
    std::span<AtomDynamics> frameworkDynamics = {}, const CellMinimizationLayout& cellLayout = {},
    const Framework* framework = nullptr);

}  // namespace Interactions
