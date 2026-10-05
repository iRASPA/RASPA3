module;

module mc_moves_bead_flip;

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

std::optional<RunningEnergy> MC_Moves::beadFlipMove(RandomNumber &random, System &system,
                                                    std::size_t selectedComponent, std::size_t selectedMolecule)
{
  Move::Types move = Move::Types::BeadFlip;
  Component &component = system.components[selectedComponent];

  // Restriction: the polarization energy update is not implemented. Molecules without a valid flip
  // bead cannot use this move. Neither condition changes during a simulation, so the attempts are
  // not counted as trials.
  const std::vector<Component::FlipBead> &beads = component.flipBeads();
  if (beads.empty() || system.forceField.computePolarization)
  {
    return std::nullopt;
  }

  // Mixed step sizes (see the pivot move): channel 0 is the adaptive small-step channel, channel 1
  // the full randomization pinned at pi.
  std::size_t channel = (random.uniform() < component.beadFlipRandomizationFraction) ? 1uz : 0uz;

  component.mc_moves_statistics.addTrial(move, channel);

  std::span<Atom> molecule_atoms = system.spanOfMolecule(selectedComponent, selectedMolecule);
  Molecule &molecule = system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, selectedMolecule)];

  // The bead is chosen uniformly from the fixed list, so the proposal is symmetric.
  const Component::FlipBead &bead = beads[static_cast<std::size_t>(random.uniform() * static_cast<double>(beads.size()))];

  const double3 anchor = molecule_atoms[bead.neighbours[0]].position;
  const double3 oldPosition = molecule_atoms[bead.atom].position;
  double3 newPosition{};
  if (bead.neighbours.size() == 2)
  {
    // Kink jump: rotation about the axis through the two neighbours by a random angle. Both bond
    // lengths (distances to points on the axis) are preserved.
    const double maxAngle = component.mc_moves_statistics.getMaxChange(move, channel);
    const double angle = maxAngle * 2.0 * (random.uniform() - 0.5);
    const double3 axis = (molecule_atoms[bead.neighbours[1]].position - anchor).normalized();
    const double3x3 rotation = double3x3::buildRotationMatrixInverse(simd_quatd::fromAxisAngle(angle, axis));
    newPosition = anchor + rotation * (oldPosition - anchor);
  }
  else if (channel == 1)
  {
    // End rotation, full randomization: a uniformly random bond direction at the old bond length.
    // The proposal density does not depend on the old direction, hence it is symmetric.
    newPosition = anchor + (oldPosition - anchor).length() * random.randomVectorOnUnitSphere();
  }
  else
  {
    // End rotation, small step: rotation about a uniformly random axis through the neighbour by a
    // random angle in [-maxAngle, maxAngle]. The reverse move uses the same axis with the negated
    // angle, which has the same proposal density, hence the proposal is symmetric.
    const double maxAngle = component.mc_moves_statistics.getMaxChange(move, channel);
    const double angle = maxAngle * 2.0 * (random.uniform() - 0.5);
    const double3 axis = random.randomVectorOnUnitSphere();
    const double3x3 rotation = double3x3::buildRotationMatrixInverse(simd_quatd::fromAxisAngle(angle, axis));
    newPosition = anchor + rotation * (oldPosition - anchor);
  }

  std::vector<Atom> trialAtoms(molecule_atoms.begin(), molecule_atoms.end());
  trialAtoms[bead.atom].position = newPosition;

  if (system.insideBlockedPockets(component, trialAtoms))
  {
    return std::nullopt;
  }

  // Only the flipped bead interacts differently with its surroundings: the external-field, framework,
  // inter-molecular and Ewald differences are evaluated for that atom alone (every other atom contributes
  // identically to the new and the old configuration, so its terms cancel exactly). The intramolecular
  // terms below still use the whole molecule.
  const std::array<std::size_t, 1> movedIndices{bead.atom};
  const std::span<const Atom> movedNew(&trialAtoms[bead.atom], 1);
  const std::span<const Atom> movedOld(&molecule_atoms[bead.atom], 1);

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
                  system.forceField, system.simulationBox, trialAtoms, molecule_atoms, movedIndices);
            });

  // Intramolecular energy contribution: the bond lengths of the bead are invariant, the bends and
  // torsions at the junctions and the intramolecular non-bonded energy change. Recomputing all
  // internal terms keeps the bookkeeping exact for any force field.
  RunningEnergy internalDifference = component.intraMolecularPotentials.computeInternalEnergies(trialAtoms) -
                                     component.intraMolecularPotentials.computeInternalEnergies(molecule_atoms);

  RunningEnergy energyDifference =
      externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() + ewaldFourierEnergy +
      internalDifference + system.crossLinkEnergyDifference(selectedComponent, selectedMolecule, trialAtoms, molecule_atoms);

  component.mc_moves_statistics.addConstructed(move, channel);

  // Metropolis criterion (symmetric proposal).
  if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy()))
  {
    component.mc_moves_statistics.addAccepted(move, channel);

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);

    std::copy(trialAtoms.cbegin(), trialAtoms.cend(), molecule_atoms.begin());
    system.cellListAtomsMoved(molecule_atoms);
    molecule.centerOfMassPosition = component.computeCenterOfMass(molecule_atoms);

    return energyDifference;
  }
  return std::nullopt;
}
