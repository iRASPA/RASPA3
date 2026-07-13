module;

export module interactions_hessian_intramolecular;

import std;

import atom;
import atom_dynamics;
import molecule;
import component;
import running_energy;
import generalized_hessian;
import minimization_dof_layout;

export namespace Interactions
{
RunningEnergy computeIntraMolecularBondHessian(std::span<const Molecule> moleculeData, std::span<const Atom> atoms,
                                               std::span<const Component> components,
                                               const MinimizationDofLayout &layout, GeneralizedHessian &hessian,
                                               std::span<AtomDynamics> dynamics);

RunningEnergy computeIntraMolecularBendHessian(std::span<const Molecule> moleculeData, std::span<const Atom> atoms,
                                               std::span<const Component> components,
                                               const MinimizationDofLayout &layout, GeneralizedHessian &hessian,
                                               std::span<AtomDynamics> dynamics);

RunningEnergy computeIntraMolecularUreyBradleyHessian(std::span<const Molecule> moleculeData,
                                                      std::span<const Atom> atoms,
                                                      std::span<const Component> components,
                                                      const MinimizationDofLayout &layout, GeneralizedHessian &hessian,
                                                      std::span<AtomDynamics> dynamics);

RunningEnergy computeIntraMolecularVanDerWaalsHessian(std::span<const Molecule> moleculeData,
                                                      std::span<const Atom> atoms,
                                                      std::span<const Component> components,
                                                      const MinimizationDofLayout &layout, GeneralizedHessian &hessian,
                                                      std::span<AtomDynamics> dynamics);

RunningEnergy computeIntraMolecularCoulombHessian(std::span<const Molecule> moleculeData, std::span<const Atom> atoms,
                                                  std::span<const Component> components,
                                                  const MinimizationDofLayout &layout, GeneralizedHessian &hessian,
                                                  std::span<AtomDynamics> dynamics);

RunningEnergy computeIntraMolecularTorsionHessian(std::span<const Molecule> moleculeData, std::span<const Atom> atoms,
                                                  std::span<const Component> components,
                                                  const MinimizationDofLayout &layout, GeneralizedHessian &hessian,
                                                  std::span<AtomDynamics> dynamics);

RunningEnergy computeIntraMolecularBondBondHessian(std::span<const Molecule> moleculeData, std::span<const Atom> atoms,
                                                   std::span<const Component> components,
                                                   const MinimizationDofLayout &layout, GeneralizedHessian &hessian,
                                                   std::span<AtomDynamics> dynamics);

RunningEnergy computeIntraMolecularBondBendHessian(std::span<const Molecule> moleculeData, std::span<const Atom> atoms,
                                                   std::span<const Component> components,
                                                   const MinimizationDofLayout &layout, GeneralizedHessian &hessian,
                                                   std::span<AtomDynamics> dynamics);

RunningEnergy computeIntraMolecularBendBendHessian(std::span<const Molecule> moleculeData, std::span<const Atom> atoms,
                                                   std::span<const Component> components,
                                                   const MinimizationDofLayout &layout, GeneralizedHessian &hessian,
                                                   std::span<AtomDynamics> dynamics);

RunningEnergy computeIntraMolecularBondTorsionHessian(std::span<const Molecule> moleculeData,
                                                      std::span<const Atom> atoms,
                                                      std::span<const Component> components,
                                                      const MinimizationDofLayout &layout, GeneralizedHessian &hessian,
                                                      std::span<AtomDynamics> dynamics);

RunningEnergy computeIntraMolecularBendTorsionHessian(std::span<const Molecule> moleculeData,
                                                      std::span<const Atom> atoms,
                                                      std::span<const Component> components,
                                                      const MinimizationDofLayout &layout, GeneralizedHessian &hessian,
                                                      std::span<AtomDynamics> dynamics);

RunningEnergy computeIntraMolecularInversionBendHessian(std::span<const Molecule> moleculeData,
                                                        std::span<const Atom> atoms,
                                                        std::span<const Component> components,
                                                        const MinimizationDofLayout &layout, GeneralizedHessian &hessian,
                                                        std::span<AtomDynamics> dynamics);

}  // namespace Interactions
