module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <future>
#include <iostream>
#include <numbers>
#include <optional>
#include <span>
#include <thread>
#include <tuple>
#include <vector>
#endif

module interactions_internal;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import energy_status;
import potential_energy_vdw;
import potential_gradient_vdw;
import potential_energy_coulomb;
import potential_gradient_coulomb;
import potential_correction_vdw;
import potential_electrostatics;
import simulationbox;
import double3;
import double3x3;
import forcefield;
import atom;
import molecule;
import energy_factor;
import gradient_factor;
import energy_status_inter;
import running_energy;
import intra_molecular_potentials;
import bond_potential;
import bend_potential;
import component;
import units;
import threadpool;

RunningEnergy Interactions::computeIntraMolecularEnergy(
    const Potentials::IntraMolecularPotentials &intraMolecularPotentials, std::span<const Molecule> moleculeData,
    std::span<const Atom> moleculeAtoms) noexcept
{
  RunningEnergy energy{};

  for (const Molecule &molecule : moleculeData)
  {
    std::span<const Atom> atom_molecule_span = {&moleculeAtoms[molecule.atomIndex], molecule.numberOfAtoms};
    energy += intraMolecularPotentials.computeInternalEnergies(atom_molecule_span);
  }

  return energy;
}

RunningEnergy Interactions::computeIntraMolecularBondEnergy(
    const Potentials::IntraMolecularPotentials &intraMolecularPotentials, std::span<const Molecule> moleculeData,
    std::span<const Atom> moleculeAtoms) noexcept
{
  RunningEnergy energy{};

  for (const Molecule &molecule : moleculeData)
  {
    std::span<const Atom> atom_molecule_span = {&moleculeAtoms[molecule.atomIndex], molecule.numberOfAtoms};
    energy += intraMolecularPotentials.computeInternalBondEnergies(atom_molecule_span);
  }

  return energy;
}

RunningEnergy Interactions::computeIntraMolecularBendEnergy(
    const Potentials::IntraMolecularPotentials &intraMolecularPotentials, std::span<const Molecule> moleculeData,
    std::span<const Atom> moleculeAtoms) noexcept
{
  RunningEnergy energy{};

  for (const Molecule &molecule : moleculeData)
  {
    std::span<const Atom> atom_molecule_span = {&moleculeAtoms[molecule.atomIndex], molecule.numberOfAtoms};
    energy += intraMolecularPotentials.computeInternalBendEnergies(atom_molecule_span);
  }

  return energy;
}

std::pair<double, double3x3> Interactions::computeIntraMolecularBondStrainDerivative(
    const Potentials::IntraMolecularPotentials &intraMolecularPotentials, std::span<const Molecule> moleculeData,
    const std::span<Atom> atoms)
{
  double energy{};
  double3x3 strain_dervative{};

  for (const Molecule &molecule : moleculeData)
  {
    std::span<Atom> atom_molecule_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    auto [intra_energy, intra_strain_dervative] =
        Interactions::computeIntraMolecularBondStrainDerivative(intraMolecularPotentials, atom_molecule_span);

    energy += intra_energy;
    strain_dervative += intra_strain_dervative;
  }

  return {energy, strain_dervative};
}

std::pair<double, double3x3> Interactions::computeIntraMolecularBondStrainDerivative(
    const Potentials::IntraMolecularPotentials &intraMolecularPotentials, const std::span<Atom> atoms)
{
  double bond_energy{};
  double3x3 bond_strain_derivative_tensor{};

  for (const BondPotential &bond_potential : intraMolecularPotentials.bonds)
  {
    std::size_t A = bond_potential.identifiers[0];
    std::size_t B = bond_potential.identifiers[1];

    double3 posA = atoms[A].position;
    double3 posB = atoms[B].position;

    auto [energy, energy_derivative, strain_derivative] = bond_potential.potentialEnergyGradientStrain(posA, posB);

    bond_energy += energy;

    atoms[A].gradient += energy_derivative[0];
    atoms[B].gradient += energy_derivative[1];

    bond_strain_derivative_tensor += strain_derivative;
  }

  return {bond_energy, bond_strain_derivative_tensor};
}

std::pair<double, double3x3> Interactions::computeIntraMolecularBendStrainDerivative(
    const Potentials::IntraMolecularPotentials &intraMolecularPotentials, std::span<const Molecule> moleculeData,
    const std::span<Atom> atoms)
{
  double energy{};
  double3x3 strain_dervative{};

  for (const Molecule &molecule : moleculeData)
  {
    std::span<Atom> atom_molecule_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    auto [intra_energy, intra_strain_dervative] =
        Interactions::computeIntraMolecularBendStrainDerivative(intraMolecularPotentials, atom_molecule_span);

    energy += intra_energy;
    strain_dervative += intra_strain_dervative;
  }

  return {energy, strain_dervative};
}

std::pair<double, double3x3> Interactions::computeIntraMolecularBendStrainDerivative(
    const Potentials::IntraMolecularPotentials &intraMolecularPotentials, const std::span<Atom> atoms)
{
  double bend_energy{};
  double3x3 bend_strain_derivative_tensor{};

  for (const BendPotential &bend_potential : intraMolecularPotentials.bends)
  {
    std::size_t A = bend_potential.identifiers[0];
    std::size_t B = bend_potential.identifiers[1];
    std::size_t C = bend_potential.identifiers[2];

    double3 posA = atoms[A].position;
    double3 posB = atoms[B].position;
    double3 posC = atoms[C].position;

    auto [energy, energy_derivative, strain_derivative] =
        bend_potential.potentialEnergyGradientStrain(posA, posB, posC);

    bend_energy += energy;

    atoms[A].gradient += energy_derivative[0];
    atoms[B].gradient += energy_derivative[1];
    atoms[C].gradient += energy_derivative[2];

    bend_strain_derivative_tensor += strain_derivative;
  }

  return {bend_energy, bend_strain_derivative_tensor};
}
