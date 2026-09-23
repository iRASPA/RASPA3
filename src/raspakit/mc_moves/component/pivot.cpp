module;

module mc_moves_pivot;

import std;

import component;
import atom;
import molecule;
import double3;
import double3x3;
import simd_quatd;
import simulationbox;
import randomnumbers;
import system;
import running_energy;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import mc_moves_move_types;

std::optional<RunningEnergy> MC_Moves::pivotMove(RandomNumber &random, System &system, std::size_t selectedComponent,
                                                 std::size_t selectedMolecule)
{
  std::chrono::steady_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::Pivot;
  Component &component = system.components[selectedComponent];

  // Restriction: the polarization energy update is not implemented. Molecules without a valid
  // pivot axis (fully rigid molecules, pure rings, diatomics) cannot use this move. Neither
  // condition changes during a simulation, so the attempts are not counted as trials.
  const std::vector<Component::PivotBond> &pivotBonds = component.pivotBonds();
  if (pivotBonds.empty() || system.forceField.computePolarization)
  {
    return std::nullopt;
  }

  // Mixed step sizes: most attempts perturb the angle within the adaptive window (statistics
  // channel 0, efficient in dense states), a fraction of the attempts randomizes the angle
  // completely (channel 1, whose window is pinned at pi; efficient in expanded states and able to
  // escape local minima). The channels have separate statistics so that the full randomizations do
  // not bias the adaptive maximum angle of the small-step channel.
  std::size_t channel = (random.uniform() < component.pivotRandomizationFraction) ? 1uz : 0uz;

  component.mc_moves_statistics.addTrial(move, channel);

  std::span<Atom> molecule_atoms = system.spanOfMolecule(selectedComponent, selectedMolecule);
  Molecule &molecule = system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, selectedMolecule)];

  // Select a pivot axis uniformly from the precomputed valid bonds (not part of a ring, not
  // interior to a rigid fragment, with a non-empty part to rotate). The rotated part is the
  // smaller side of the molecule, chosen deterministically, so the proposal stays symmetric (the
  // reverse move selects the same bond, the same part and the negated angle).
  const Component::PivotBond &pivotBond =
      pivotBonds[static_cast<std::size_t>(random.uniform() * static_cast<double>(pivotBonds.size()))];
  const std::array<std::size_t, 2> &bond = pivotBond.bond;
  const std::vector<std::size_t> &rotatedAtoms = pivotBond.rotatedAtoms;

  // Construct the trial positions: rigid rotation of the selected part about the bond axis.
  double maxAngle = component.mc_moves_statistics.getMaxChange(move, channel);
  double rotationAngle = maxAngle * 2.0 * (random.uniform() - 0.5);
  double3 axisOrigin = molecule_atoms[bond[0]].position;
  double3 rotationAxis = (molecule_atoms[bond[1]].position - axisOrigin).normalized();
  simd_quatd q = simd_quatd::fromAxisAngle(rotationAngle, rotationAxis);
  double3x3 rotationMatrix = double3x3::buildRotationMatrixInverse(q);

  std::vector<Atom> trialAtoms(molecule_atoms.begin(), molecule_atoms.end());
  for (std::size_t index : rotatedAtoms)
  {
    trialAtoms[index].position = axisOrigin + rotationMatrix * (trialAtoms[index].position - axisOrigin);
  }

  if (system.insideBlockedPockets(component, trialAtoms))
  {
    return std::nullopt;
  }

  // Compute external field energy contribution
  time_begin = std::chrono::steady_clock::now();
  std::optional<RunningEnergy> externalFieldMolecule = Interactions::computeExternalFieldEnergyDifference(
      system.hasExternalField, system.forceField, system.simulationBox, system.externalFieldInterpolationGrid,
      trialAtoms, molecule_atoms);
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::ExternalFieldMolecule] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::ExternalFieldMolecule] += (time_end - time_begin);
  if (!externalFieldMolecule.has_value()) return std::nullopt;

  // Compute framework-molecule energy contribution
  time_begin = std::chrono::steady_clock::now();
  std::optional<RunningEnergy> frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
      system.spanOfFrameworkAtoms(), trialAtoms, molecule_atoms);
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::FrameworkMolecule] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::FrameworkMolecule] += (time_end - time_begin);
  if (!frameworkMolecule.has_value()) return std::nullopt;

  // Compute molecule-molecule energy contribution
  time_begin = std::chrono::steady_clock::now();
  std::optional<RunningEnergy> interMolecule = Interactions::computeInterMolecularEnergyDifference(
      system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), trialAtoms, molecule_atoms);
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::MoleculeMolecule] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::MoleculeMolecule] += (time_end - time_begin);
  if (!interMolecule.has_value()) return std::nullopt;

  // Compute Ewald energy contribution
  time_begin = std::chrono::steady_clock::now();
  RunningEnergy ewaldFourierEnergy = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
      system.simulationBox, trialAtoms, molecule_atoms);
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);

  // Intramolecular energy contribution: bond lengths and bend angles are invariant under the pivot
  // rotation, but the torsions through the pivot bond and the intramolecular non-bonded energy
  // change. Recomputing all internal terms keeps the bookkeeping exact for any force field.
  RunningEnergy internalDifference = component.intraMolecularPotentials.computeInternalEnergies(trialAtoms) -
                                     component.intraMolecularPotentials.computeInternalEnergies(molecule_atoms);

  // Calculate the total energy difference
  RunningEnergy energyDifference = externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() +
                                   ewaldFourierEnergy + internalDifference;

  component.mc_moves_statistics.addConstructed(move, channel);

  // Apply acceptance/rejection rule based on Metropolis criterion
  if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy()))
  {
    component.mc_moves_statistics.addAccepted(move, channel);

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);

    std::copy(trialAtoms.cbegin(), trialAtoms.cend(), molecule_atoms.begin());
    molecule.centerOfMassPosition = component.computeCenterOfMass(molecule_atoms);

    return energyDifference;
  }
  return std::nullopt;
}
