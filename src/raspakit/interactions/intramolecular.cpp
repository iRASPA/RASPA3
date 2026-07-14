module;

module interactions_internal;

import std;

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
import int3;
import forcefield;
import framework;
import atom;
import atom_dynamics;
import molecule;
import energy_factor;
import gradient_factor;
import energy_status_inter;
import running_energy;
import intra_molecular_potentials;
import bond_potential;
import bend_potential;
import torsion_potential;
import van_der_waals_potential;
import coulomb_potential;
import component;
import units;
import threadpool;

RunningEnergy Interactions::computeFrameworkIntraMolecularEnergy(const Framework& framework,
                                                                 const SimulationBox& simulationBox,
                                                                 std::span<const Atom> atoms)
{
  RunningEnergy energy{};
  if (framework.rigid) return energy;

  const auto shiftedPosition = [&](std::size_t atom, const int3& shift)
  {
    return atoms[atom].position + simulationBox.cell * double3(static_cast<double>(shift.x),
                                                               static_cast<double>(shift.y),
                                                               static_cast<double>(shift.z));
  };
  const auto& potentials = framework.intraMolecularPotentials;
  const auto& images = framework.intraMolecularImageShifts;

  for (std::size_t index = 0; index < potentials.bonds.size(); ++index)
  {
    const BondPotential& potential = potentials.bonds[index];
    energy.bond += potential.calculateEnergy(shiftedPosition(potential.identifiers[0], images.bonds[index][0]),
                                             shiftedPosition(potential.identifiers[1], images.bonds[index][1]));
  }
  for (std::size_t index = 0; index < potentials.bends.size(); ++index)
  {
    const BendPotential& potential = potentials.bends[index];
    const auto& ids = potential.identifiers;
    const auto& shifts = images.bends[index];
    energy.bend += potential.calculateEnergy(shiftedPosition(ids[0], shifts[0]), shiftedPosition(ids[1], shifts[1]),
                                             shiftedPosition(ids[2], shifts[2]), std::nullopt);
  }
  const auto accumulateTorsions = [&](const std::vector<TorsionPotential>& torsions,
                                      const std::vector<std::array<int3, 4>>& termImages,
                                      double RunningEnergy::* member)
  {
    for (std::size_t index = 0; index < torsions.size(); ++index)
    {
      const TorsionPotential& potential = torsions[index];
      const auto& ids = potential.identifiers;
      const auto& shifts = termImages[index];
      energy.*member +=
          potential.calculateEnergy(shiftedPosition(ids[0], shifts[0]), shiftedPosition(ids[1], shifts[1]),
                                    shiftedPosition(ids[2], shifts[2]), shiftedPosition(ids[3], shifts[3]));
    }
  };
  accumulateTorsions(potentials.torsions, images.torsions, &RunningEnergy::torsion);
  accumulateTorsions(potentials.improperTorsions, images.improperTorsions, &RunningEnergy::improperTorsion);

  for (std::size_t index = 0; index < potentials.vanDerWaals.size(); ++index)
  {
    const VanDerWaalsPotential& potential = potentials.vanDerWaals[index];
    energy.intraVDW +=
        potential.calculateEnergy(shiftedPosition(potential.identifiers[0], images.vanDerWaals[index][0]),
                                  shiftedPosition(potential.identifiers[1], images.vanDerWaals[index][1]));
  }
  for (std::size_t index = 0; index < potentials.coulombs.size(); ++index)
  {
    const CoulombPotential& potential = potentials.coulombs[index];
    energy.intraCoul += potential.calculateEnergy(shiftedPosition(potential.identifiers[0], images.coulombs[index][0]),
                                                  shiftedPosition(potential.identifiers[1], images.coulombs[index][1]));
  }
  return energy;
}

RunningEnergy Interactions::computeIntraMolecularEnergy(
    const Potentials::IntraMolecularPotentials& intraMolecularPotentials, std::span<const Molecule> moleculeData,
    std::span<const Atom> moleculeAtoms) noexcept
{
  RunningEnergy energy{};

  for (const Molecule& molecule : moleculeData)
  {
    std::span<const Atom> atom_molecule_span = {&moleculeAtoms[molecule.atomIndex], molecule.numberOfAtoms};
    energy += intraMolecularPotentials.computeInternalEnergies(atom_molecule_span);
  }

  return energy;
}

RunningEnergy Interactions::computeIntraMolecularBondEnergy(
    const Potentials::IntraMolecularPotentials& intraMolecularPotentials, std::span<const Molecule> moleculeData,
    std::span<const Atom> moleculeAtoms) noexcept
{
  RunningEnergy energy{};

  for (const Molecule& molecule : moleculeData)
  {
    std::span<const Atom> atom_molecule_span = {&moleculeAtoms[molecule.atomIndex], molecule.numberOfAtoms};
    energy += intraMolecularPotentials.computeInternalBondEnergies(atom_molecule_span);
  }

  return energy;
}

RunningEnergy Interactions::computeIntraMolecularBendEnergy(
    const Potentials::IntraMolecularPotentials& intraMolecularPotentials, std::span<const Molecule> moleculeData,
    std::span<const Atom> moleculeAtoms) noexcept
{
  RunningEnergy energy{};

  for (const Molecule& molecule : moleculeData)
  {
    std::span<const Atom> atom_molecule_span = {&moleculeAtoms[molecule.atomIndex], molecule.numberOfAtoms};
    energy += intraMolecularPotentials.computeInternalBendEnergies(atom_molecule_span);
  }

  return energy;
}

std::pair<double, double3x3> Interactions::computeIntraMolecularBondStrainDerivative(
    const Potentials::IntraMolecularPotentials& intraMolecularPotentials, std::span<const Molecule> moleculeData,
    std::span<const Atom> atoms, std::span<AtomDynamics> dynamics)
{
  double energy{};
  double3x3 strain_dervative{};

  for (const Molecule& molecule : moleculeData)
  {
    std::span<const Atom> atom_molecule_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_molecule_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};
    auto [intra_energy, intra_strain_dervative] = Interactions::computeIntraMolecularBondStrainDerivative(
        intraMolecularPotentials, atom_molecule_span, dynamics_molecule_span);

    energy += intra_energy;
    strain_dervative += intra_strain_dervative;
  }

  return {energy, strain_dervative};
}

std::pair<double, double3x3> Interactions::computeIntraMolecularBondStrainDerivative(
    const Potentials::IntraMolecularPotentials& intraMolecularPotentials, std::span<const Atom> atoms,
    std::span<AtomDynamics> dynamics)
{
  double bond_energy{};
  double3x3 bond_strain_derivative_tensor{};

  for (const BondPotential& bond_potential : intraMolecularPotentials.bonds)
  {
    std::size_t A = bond_potential.identifiers[0];
    std::size_t B = bond_potential.identifiers[1];

    double3 posA = atoms[A].position;
    double3 posB = atoms[B].position;

    auto [energy, energy_derivative, strain_derivative] = bond_potential.potentialEnergyGradientStrain(posA, posB);

    bond_energy += energy;

    dynamics[A].gradient += energy_derivative[0];
    dynamics[B].gradient += energy_derivative[1];

    bond_strain_derivative_tensor += strain_derivative;
  }

  return {bond_energy, bond_strain_derivative_tensor};
}

std::pair<double, double3x3> Interactions::computeIntraMolecularBendStrainDerivative(
    const Potentials::IntraMolecularPotentials& intraMolecularPotentials, std::span<const Molecule> moleculeData,
    std::span<const Atom> atoms, std::span<AtomDynamics> dynamics)
{
  double energy{};
  double3x3 strain_dervative{};

  for (const Molecule& molecule : moleculeData)
  {
    std::span<const Atom> atom_molecule_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_molecule_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};
    auto [intra_energy, intra_strain_dervative] = Interactions::computeIntraMolecularBendStrainDerivative(
        intraMolecularPotentials, atom_molecule_span, dynamics_molecule_span);

    energy += intra_energy;
    strain_dervative += intra_strain_dervative;
  }

  return {energy, strain_dervative};
}

std::pair<double, double3x3> Interactions::computeIntraMolecularBendStrainDerivative(
    const Potentials::IntraMolecularPotentials& intraMolecularPotentials, std::span<const Atom> atoms,
    std::span<AtomDynamics> dynamics)
{
  double bend_energy{};
  double3x3 bend_strain_derivative_tensor{};

  for (const BendPotential& bend_potential : intraMolecularPotentials.bends)
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

    dynamics[A].gradient += energy_derivative[0];
    dynamics[B].gradient += energy_derivative[1];
    dynamics[C].gradient += energy_derivative[2];

    bend_strain_derivative_tensor += strain_derivative;
  }

  return {bend_energy, bend_strain_derivative_tensor};
}

RunningEnergy Interactions::computeIntraMolecularGradient(
    const Potentials::IntraMolecularPotentials& intraMolecularPotentials, std::span<const Molecule> moleculeData,
    std::span<const Atom> moleculeAtoms, std::span<AtomDynamics> moleculeDynamics) noexcept
{
  RunningEnergy energy{};

  for (const Molecule& molecule : moleculeData)
  {
    std::span<const Atom> atom_molecule_span = {&moleculeAtoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_molecule_span = {&moleculeDynamics[molecule.atomIndex], molecule.numberOfAtoms};
    energy += intraMolecularPotentials.computeInternalGradient(atom_molecule_span, dynamics_molecule_span);
  }

  return energy;
}

std::pair<double, double3x3> Interactions::computeIntraMolecularStrainDerivative(
    const Potentials::IntraMolecularPotentials& intraMolecularPotentials, std::span<const Molecule> moleculeData,
    std::span<const Atom> atoms, std::span<AtomDynamics> dynamics)
{
  double energy{};
  double3x3 strain_derivative{};

  for (const Molecule& molecule : moleculeData)
  {
    std::span<const Atom> atom_molecule_span = {&atoms[molecule.atomIndex], molecule.numberOfAtoms};
    std::span<AtomDynamics> dynamics_molecule_span = {&dynamics[molecule.atomIndex], molecule.numberOfAtoms};
    auto [intra_energy, intra_strain_derivative] =
        intraMolecularPotentials.computeInternalStrainDerivative(atom_molecule_span, dynamics_molecule_span);

    energy += intra_energy.potentialEnergy();
    strain_derivative += intra_strain_derivative;
  }

  return {energy, strain_derivative};
}

std::pair<double, double3x3> Interactions::computeIntraMolecularStrainDerivative(
    const Potentials::IntraMolecularPotentials& intraMolecularPotentials, std::span<const Atom> atoms,
    std::span<AtomDynamics> dynamics)
{
  auto [running_energy, strain_derivative] = intraMolecularPotentials.computeInternalStrainDerivative(atoms, dynamics);
  return {running_energy.potentialEnergy(), strain_derivative};
}
