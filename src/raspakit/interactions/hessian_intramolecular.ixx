module;

export module interactions_hessian_intramolecular;

import std;

import atom;
import atom_dynamics;
import molecule;
import component;
import framework;
import simulationbox;
import forcefield;
import running_energy;
import generalized_hessian;
import minimization_dof_layout;
import minimization_cell_layout;

export namespace Interactions
{
RunningEnergy computeFrameworkIntraMolecularHessian(const ForceField& forceField, const Framework& framework,
                                                    const SimulationBox& simulationBox, std::span<const Atom> atoms,
                                                    const MinimizationDofLayout& layout, GeneralizedHessian& hessian,
                                                    std::span<AtomDynamics> dynamics,
                                                    const CellMinimizationLayout& cellLayout = {});

RunningEnergy computeIntraMolecularBondHessian(std::span<const Molecule> moleculeData, std::span<const Atom> atoms,
                                               std::span<const Component> components,
                                               const MinimizationDofLayout& layout, GeneralizedHessian& hessian,
                                               std::span<AtomDynamics> dynamics);

RunningEnergy computeIntraMolecularBendHessian(std::span<const Molecule> moleculeData, std::span<const Atom> atoms,
                                               std::span<const Component> components,
                                               const MinimizationDofLayout& layout, GeneralizedHessian& hessian,
                                               std::span<AtomDynamics> dynamics);

RunningEnergy computeIntraMolecularUreyBradleyHessian(std::span<const Molecule> moleculeData,
                                                      std::span<const Atom> atoms,
                                                      std::span<const Component> components,
                                                      const MinimizationDofLayout& layout, GeneralizedHessian& hessian,
                                                      std::span<AtomDynamics> dynamics);

/// The non-excluded intramolecular van der Waals pairs (Potentials::intraMolecularVDW: the regular force-field
/// pair potential inside its cutoff times the pair scaling, minimum image).
RunningEnergy computeIntraMolecularVanDerWaalsHessian(const ForceField& forceField, const SimulationBox& simulationBox,
                                                      std::span<const Molecule> moleculeData,
                                                      std::span<const Atom> atoms,
                                                      std::span<const Component> components,
                                                      const MinimizationDofLayout& layout, GeneralizedHessian& hessian,
                                                      std::span<AtomDynamics> dynamics);

/// The non-excluded intramolecular Coulomb pairs (Potentials::intraMolecularCoulomb, including the Ewald or
/// shifted-potential completion of the pair).
RunningEnergy computeIntraMolecularCoulombHessian(const ForceField& forceField, const SimulationBox& simulationBox,
                                                  std::span<const Molecule> moleculeData, std::span<const Atom> atoms,
                                                  std::span<const Component> components,
                                                  const MinimizationDofLayout& layout, GeneralizedHessian& hessian,
                                                  std::span<AtomDynamics> dynamics);

RunningEnergy computeIntraMolecularTorsionHessian(std::span<const Molecule> moleculeData, std::span<const Atom> atoms,
                                                  std::span<const Component> components,
                                                  const MinimizationDofLayout& layout, GeneralizedHessian& hessian,
                                                  std::span<AtomDynamics> dynamics);

RunningEnergy computeIntraMolecularBondBondHessian(std::span<const Molecule> moleculeData, std::span<const Atom> atoms,
                                                   std::span<const Component> components,
                                                   const MinimizationDofLayout& layout, GeneralizedHessian& hessian,
                                                   std::span<AtomDynamics> dynamics);

RunningEnergy computeIntraMolecularBondBendHessian(std::span<const Molecule> moleculeData, std::span<const Atom> atoms,
                                                   std::span<const Component> components,
                                                   const MinimizationDofLayout& layout, GeneralizedHessian& hessian,
                                                   std::span<AtomDynamics> dynamics);

RunningEnergy computeIntraMolecularBendBendHessian(std::span<const Molecule> moleculeData, std::span<const Atom> atoms,
                                                   std::span<const Component> components,
                                                   const MinimizationDofLayout& layout, GeneralizedHessian& hessian,
                                                   std::span<AtomDynamics> dynamics);

RunningEnergy computeIntraMolecularBondTorsionHessian(std::span<const Molecule> moleculeData,
                                                      std::span<const Atom> atoms,
                                                      std::span<const Component> components,
                                                      const MinimizationDofLayout& layout, GeneralizedHessian& hessian,
                                                      std::span<AtomDynamics> dynamics);

RunningEnergy computeIntraMolecularBendTorsionHessian(std::span<const Molecule> moleculeData,
                                                      std::span<const Atom> atoms,
                                                      std::span<const Component> components,
                                                      const MinimizationDofLayout& layout, GeneralizedHessian& hessian,
                                                      std::span<AtomDynamics> dynamics);

RunningEnergy computeIntraMolecularInversionBendHessian(std::span<const Molecule> moleculeData,
                                                        std::span<const Atom> atoms,
                                                        std::span<const Component> components,
                                                        const MinimizationDofLayout& layout,
                                                        GeneralizedHessian& hessian, std::span<AtomDynamics> dynamics);

RunningEnergy computeIntraMolecularOutOfPlaneBendHessian(std::span<const Molecule> moleculeData,
                                                         std::span<const Atom> atoms,
                                                         std::span<const Component> components,
                                                         const MinimizationDofLayout& layout,
                                                         GeneralizedHessian& hessian, std::span<AtomDynamics> dynamics);

}  // namespace Interactions
