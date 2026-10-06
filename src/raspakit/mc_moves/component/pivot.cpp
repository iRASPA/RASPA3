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
import mc_moves_cputime;

std::optional<RunningEnergy> MC_Moves::pivotMove(RandomNumber &random, System &system, std::size_t selectedComponent,
                                                 std::size_t selectedMolecule)
{
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
  const std::size_t bondIndex = static_cast<std::size_t>(random.uniform() * static_cast<double>(pivotBonds.size()));
  component.mc_moves_statistics.addSubTrial(move, bondIndex, pivotBonds.size());
  const Component::PivotBond &pivotBond = pivotBonds[bondIndex];
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

  // Only the rotated part interacts differently with its surroundings: the external-field, framework,
  // inter-molecular and Ewald differences are evaluated for those atoms alone (every other atom contributes
  // identically to the new and the old configuration, so its terms cancel exactly). The intramolecular
  // terms below still use the whole molecule.
  const std::span<const std::size_t> movedIndices(rotatedAtoms);
  std::vector<Atom> movedNew;
  std::vector<Atom> movedOld;
  movedNew.reserve(movedIndices.size());
  movedOld.reserve(movedIndices.size());
  for (std::size_t index : movedIndices)
  {
    movedNew.push_back(trialAtoms[index]);
    movedOld.push_back(molecule_atoms[index]);
  }

  // Compute external field energy contribution
  std::optional<RunningEnergy> externalFieldMolecule =
      timed(system, component, move, Move::Timing::ExternalFieldMolecule,
            [&]
            {
              return Interactions::computeExternalFieldEnergyDifference(
                  system.hasExternalField, system.forceField, system.simulationBox,
                  system.externalFieldInterpolationGrid, movedNew, movedOld);
            });
  if (!externalFieldMolecule.has_value()) return std::nullopt;

  // Compute framework-molecule energy contribution
  std::optional<RunningEnergy> frameworkMolecule =
      timed(system, component, move, Move::Timing::FrameworkMolecule,
            [&]
            {
              return Interactions::computeFrameworkMoleculeEnergyDifference(
                  system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
                  system.spanOfFrameworkAtoms(), movedNew, movedOld);
            });
  if (!frameworkMolecule.has_value()) return std::nullopt;

  // Compute molecule-molecule energy contribution
  std::optional<RunningEnergy> interMolecule =
      timed(system, component, move, Move::Timing::MoleculeMolecule,
            [&]
            {
              return Interactions::computeInterMolecularEnergyDifference(
                  system.forceField, system.simulationBox, system.cellList(), system.spanOfMoleculeAtoms(), movedNew,
                  movedOld);
            });
  if (!interMolecule.has_value()) return std::nullopt;

  // Compute Ewald energy contribution
  RunningEnergy ewaldFourierEnergy =
      timed(system, component, move, Move::Timing::Ewald,
            [&]
            {
              return Interactions::energyDifferenceEwaldFourierMovedAtoms(
                  system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik,
                  system.forceField, system.simulationBox, system.components, trialAtoms, molecule_atoms, movedIndices);
            });

  // Intramolecular energy contribution: bond lengths and bend angles are invariant under the pivot
  // rotation, but the torsions through the pivot bond and the intramolecular non-bonded energy
  // change. Recomputing all internal terms keeps the bookkeeping exact for any force field.
  RunningEnergy internalDifference = component.intraMolecularPotentials.computeInternalEnergies(system.forceField, system.simulationBox, trialAtoms) -
                                     component.intraMolecularPotentials.computeInternalEnergies(system.forceField, system.simulationBox, molecule_atoms);

  // Calculate the total energy difference
  RunningEnergy energyDifference =
      externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() + ewaldFourierEnergy +
      internalDifference + system.crossLinkEnergyDifference(selectedComponent, selectedMolecule, trialAtoms, molecule_atoms);

  component.mc_moves_statistics.addConstructed(move, channel);
  component.mc_moves_statistics.addSubConstructed(move, bondIndex);

  // Apply acceptance/rejection rule based on Metropolis criterion
  if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy()))
  {
    component.mc_moves_statistics.addAccepted(move, channel);
    component.mc_moves_statistics.addSubAccepted(move, bondIndex);

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);

    std::copy(trialAtoms.cbegin(), trialAtoms.cend(), molecule_atoms.begin());
    system.cellListAtomsMoved(molecule_atoms);
    molecule.centerOfMassPosition = component.computeCenterOfMass(molecule_atoms);

    return energyDifference;
  }
  return std::nullopt;
}
