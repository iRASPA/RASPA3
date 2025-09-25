module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#endif

export module interactions_internal;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double3;
import double3x3;
import atom;
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

export namespace Interactions
{
RunningEnergy computeIntraMolecularEnergy(const Potentials::IntraMolecularPotentials &intraMolecularPotentials,
                                          std::span<const Molecule> moleculeData,
                                          std::span<const Atom> moleculeAtoms) noexcept;

RunningEnergy computeIntraMolecularBondEnergy(const Potentials::IntraMolecularPotentials &intraMolecularPotentials,
                                              std::span<const Molecule> moleculeData,
                                              std::span<const Atom> moleculeAtoms) noexcept;

RunningEnergy computeIntraMolecularBendEnergy(const Potentials::IntraMolecularPotentials &intraMolecularPotentials,
                                              std::span<const Molecule> moleculeData,
                                              std::span<const Atom> moleculeAtoms) noexcept;

std::pair<double, double3x3> computeIntraMolecularBondStrainDerivative(
    const Potentials::IntraMolecularPotentials &intraMolecularPotentials, std::span<const Molecule> moleculeData,
    const std::span<Atom> atoms);

std::pair<double, double3x3> computeIntraMolecularBondStrainDerivative(
    const Potentials::IntraMolecularPotentials &intraMolecularPotentials, const std::span<Atom> atoms);

std::pair<double, double3x3> computeIntraMolecularBendStrainDerivative(
    const Potentials::IntraMolecularPotentials &intraMolecularPotentials, std::span<const Molecule> moleculeData,
    const std::span<Atom> atoms);

std::pair<double, double3x3> computeIntraMolecularBendStrainDerivative(
    const Potentials::IntraMolecularPotentials &intraMolecularPotentials, const std::span<Atom> atoms);

};  // namespace Interactions
