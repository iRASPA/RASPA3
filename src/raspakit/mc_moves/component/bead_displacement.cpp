module;

module mc_moves_bead_displacement;

import std;

import component;
import atom;
import molecule;
import double3;
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

std::optional<RunningEnergy> MC_Moves::beadDisplacementMove(RandomNumber &random, System &system,
                                                            std::size_t selectedComponent,
                                                            std::size_t selectedMolecule)
{
  Move::Types move = Move::Types::BeadDisplacement;
  Component &component = system.components[selectedComponent];

  // Restriction: the polarization energy update is not implemented. Molecules without a
  // displaceable bead (rigid molecules, molecules whose beads all carry FIXED bonds) cannot use
  // this move. Neither condition changes during a simulation, so the attempts are not counted.
  const std::vector<std::size_t> &beads = component.displaceableBeads();
  if (beads.empty() || system.forceField.computePolarization)
  {
    return std::nullopt;
  }

  // One random Cartesian direction per attempt, with its own adaptive maximum displacement (the
  // layout of the translation move).
  std::size_t selectedDirection = static_cast<std::size_t>(3.0 * random.uniform());
  double maxDisplacement = component.mc_moves_statistics.getMaxChange(move, selectedDirection);
  double3 displacement{};
  displacement[selectedDirection] = maxDisplacement * 2.0 * (random.uniform() - 0.5);

  component.mc_moves_statistics.addTrial(move, selectedDirection);

  std::span<Atom> molecule_atoms = system.spanOfMolecule(selectedComponent, selectedMolecule);
  Molecule &molecule = system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, selectedMolecule)];

  // The bead is chosen uniformly from the fixed list, so the proposal is symmetric (the reverse
  // move selects the same bead, direction, and the negated displacement).
  const std::size_t bead = beads[static_cast<std::size_t>(random.uniform() * static_cast<double>(beads.size()))];

  std::vector<Atom> trialAtoms(molecule_atoms.begin(), molecule_atoms.end());
  trialAtoms[bead].position += displacement;

  if (system.insideBlockedPockets(component, trialAtoms))
  {
    return std::nullopt;
  }

  // Only the displaced bead interacts differently with its surroundings: the external-field, framework,
  // inter-molecular and Ewald differences are evaluated for that atom alone (every other atom contributes
  // identically to the new and the old configuration, so its terms cancel exactly). The cost of the move
  // is then independent of the chain length. The intramolecular terms below still use the whole molecule.
  const std::array<std::size_t, 1> movedIndices{bead};
  const std::span<const Atom> movedNew(&trialAtoms[bead], 1);
  const std::span<const Atom> movedOld(&molecule_atoms[bead], 1);

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

  // Intramolecular energy contribution: every bonded term containing the bead and the
  // intramolecular non-bonded energy change. Recomputing all internal terms keeps the bookkeeping
  // exact for any force field.
  RunningEnergy internalDifference = component.intraMolecularPotentials.computeInternalEnergies(trialAtoms) -
                                     component.intraMolecularPotentials.computeInternalEnergies(molecule_atoms);

  RunningEnergy energyDifference =
      externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() + ewaldFourierEnergy +
      internalDifference + system.crossLinkEnergyDifference(selectedComponent, selectedMolecule, trialAtoms, molecule_atoms);

  component.mc_moves_statistics.addConstructed(move, selectedDirection);

  // Metropolis criterion (symmetric proposal).
  if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy()))
  {
    component.mc_moves_statistics.addAccepted(move, selectedDirection);

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);

    std::copy(trialAtoms.cbegin(), trialAtoms.cend(), molecule_atoms.begin());
    system.cellListAtomsMoved(molecule_atoms);
    molecule.centerOfMassPosition = component.computeCenterOfMass(molecule_atoms);

    return energyDifference;
  }
  return std::nullopt;
}
