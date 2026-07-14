module;

export module interactions_internal;

import std;

import double3;
import double3x3;
import atom;
import atom_dynamics;
import molecule;
import running_energy;
import energy_status;
import simulationbox;
import energy_factor;
import intra_molecular_potentials;
import bond_potential;
import gradient_factor;
import forcefield;
import component;
import framework;

export namespace Interactions
{
RunningEnergy computeFrameworkIntraMolecularEnergy(const Framework& framework, const SimulationBox& simulationBox,
                                                   std::span<const Atom> frameworkAtoms);

RunningEnergy computeIntraMolecularEnergy(const Potentials::IntraMolecularPotentials& intraMolecularPotentials,
                                          std::span<const Molecule> moleculeData,
                                          std::span<const Atom> moleculeAtoms) noexcept;

RunningEnergy computeIntraMolecularBondEnergy(const Potentials::IntraMolecularPotentials& intraMolecularPotentials,
                                              std::span<const Molecule> moleculeData,
                                              std::span<const Atom> moleculeAtoms) noexcept;

RunningEnergy computeIntraMolecularBendEnergy(const Potentials::IntraMolecularPotentials& intraMolecularPotentials,
                                              std::span<const Molecule> moleculeData,
                                              std::span<const Atom> moleculeAtoms) noexcept;

std::pair<double, double3x3> computeIntraMolecularBondStrainDerivative(
    const Potentials::IntraMolecularPotentials& intraMolecularPotentials, std::span<const Molecule> moleculeData,
    std::span<const Atom> atoms, std::span<AtomDynamics> dynamics);

std::pair<double, double3x3> computeIntraMolecularBondStrainDerivative(
    const Potentials::IntraMolecularPotentials& intraMolecularPotentials, std::span<const Atom> atoms,
    std::span<AtomDynamics> dynamics);

std::pair<double, double3x3> computeIntraMolecularBendStrainDerivative(
    const Potentials::IntraMolecularPotentials& intraMolecularPotentials, std::span<const Molecule> moleculeData,
    std::span<const Atom> atoms, std::span<AtomDynamics> dynamics);

std::pair<double, double3x3> computeIntraMolecularBendStrainDerivative(
    const Potentials::IntraMolecularPotentials& intraMolecularPotentials, std::span<const Atom> atoms,
    std::span<AtomDynamics> dynamics);

RunningEnergy computeIntraMolecularGradient(const Potentials::IntraMolecularPotentials& intraMolecularPotentials,
                                            std::span<const Molecule> moleculeData, std::span<const Atom> moleculeAtoms,
                                            std::span<AtomDynamics> moleculeDynamics) noexcept;

std::pair<double, double3x3> computeIntraMolecularStrainDerivative(
    const Potentials::IntraMolecularPotentials& intraMolecularPotentials, std::span<const Molecule> moleculeData,
    std::span<const Atom> atoms, std::span<AtomDynamics> dynamics);

std::pair<double, double3x3> computeIntraMolecularStrainDerivative(
    const Potentials::IntraMolecularPotentials& intraMolecularPotentials, std::span<const Atom> atoms,
    std::span<AtomDynamics> dynamics);

};  // namespace Interactions
