module;

module mc_moves_crankshaft;

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

std::optional<RunningEnergy> MC_Moves::crankshaftMove(RandomNumber &random, System &system,
                                                      std::size_t selectedComponent, std::size_t selectedMolecule)
{
  Move::Types move = Move::Types::Crankshaft;
  Component &component = system.components[selectedComponent];

  // Restriction: the polarization energy update is not implemented. Molecules without a valid
  // crankshaft unit (fully rigid molecules, diatomics, too-small chains) cannot use this move.
  // Neither condition changes during a simulation, so the attempts are not counted as trials.
  const std::vector<Component::CrankshaftUnit> &crankshaftUnits = component.crankshaftUnits();
  if (crankshaftUnits.empty() || system.forceField.computePolarization)
  {
    return std::nullopt;
  }

  // Mixed step sizes: most attempts perturb the angle within the adaptive window (statistics
  // channel 0, efficient in dense states), a fraction of the attempts randomizes the angle
  // completely (channel 1, whose window is pinned at pi; efficient in expanded states and able to
  // escape local minima). The channels have separate statistics so that the full randomizations do
  // not bias the adaptive maximum angle of the small-step channel.
  std::size_t channel = (random.uniform() < component.crankshaftRandomizationFraction) ? 1uz : 0uz;

  component.mc_moves_statistics.addTrial(move, channel);

  std::span<Atom> molecule_atoms = system.spanOfMolecule(selectedComponent, selectedMolecule);
  Molecule &molecule = system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, selectedMolecule)];

  // Select a crankshaft unit uniformly from the precomputed valid units (a segment attached to the
  // rest of the molecule only through the two anchor atoms, with rigid fragments entirely inside or
  // outside the segment). The unit list is fixed, so the proposal stays symmetric (the reverse move
  // selects the same unit and the negated angle).
  const Component::CrankshaftUnit &unit =
      crankshaftUnits[static_cast<std::size_t>(random.uniform() * static_cast<double>(crankshaftUnits.size()))];
  const std::array<std::size_t, 2> &axis = unit.axis;
  const std::vector<std::size_t> &rotatedAtoms = unit.rotatedAtoms;

  // Construct the trial positions: rigid rotation of the segment about the anchor-anchor axis. The
  // anchors themselves do not move, so every rotated atom keeps its distance to both anchors and
  // all bond lengths are preserved.
  double maxAngle = component.mc_moves_statistics.getMaxChange(move, channel);
  double rotationAngle = maxAngle * 2.0 * (random.uniform() - 0.5);
  double3 axisOrigin = molecule_atoms[axis[0]].position;
  double3 rotationAxis = (molecule_atoms[axis[1]].position - axisOrigin).normalized();
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
  std::optional<RunningEnergy> externalFieldMolecule =
      timed(system, component, move, Move::Timing::ExternalFieldMolecule,
            [&]
            {
              return Interactions::computeExternalFieldEnergyDifference(
                  system.hasExternalField, system.forceField, system.simulationBox,
                  system.externalFieldInterpolationGrid, trialAtoms, molecule_atoms);
            });
  if (!externalFieldMolecule.has_value()) return std::nullopt;

  // Compute framework-molecule energy contribution
  std::optional<RunningEnergy> frameworkMolecule =
      timed(system, component, move, Move::Timing::FrameworkMolecule,
            [&]
            {
              return Interactions::computeFrameworkMoleculeEnergyDifference(
                  system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
                  system.spanOfFrameworkAtoms(), trialAtoms, molecule_atoms);
            });
  if (!frameworkMolecule.has_value()) return std::nullopt;

  // Compute molecule-molecule energy contribution
  std::optional<RunningEnergy> interMolecule =
      timed(system, component, move, Move::Timing::MoleculeMolecule,
            [&]
            {
              return Interactions::computeInterMolecularEnergyDifference(
                  system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), trialAtoms, molecule_atoms);
            });
  if (!interMolecule.has_value()) return std::nullopt;

  // Compute Ewald energy contribution
  RunningEnergy ewaldFourierEnergy =
      timed(system, component, move, Move::Timing::Ewald,
            [&]
            {
              return Interactions::energyDifferenceEwaldFourier(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                                system.storedEik, system.trialEik, system.forceField,
                                                                system.simulationBox, trialAtoms, molecule_atoms);
            });

  // Intramolecular energy contribution: bond lengths are invariant under the crankshaft rotation,
  // but the bend angles and torsions at the junctions and the intramolecular non-bonded energy
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
