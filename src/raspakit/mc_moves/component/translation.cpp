module;

module mc_moves_translation;

import std;

import component;
import molecule;
import atom;
import double3;
import double3x3;
import simd_quatd;
import simulationbox;
import cbmc;
import randomnumbers;
import system;
import energy_status;
import energy_status_inter;
import running_energy;
import property_lambda_probability_histogram;
import property_widom;
import averages;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import interactions_polarization;
import mc_moves_move_types;

std::optional<RunningEnergy> MC_Moves::translationMove(RandomNumber &random, System &system,
                                                       std::size_t selectedComponent, std::size_t selectedMolecule)
{
  std::span<Atom> molecule_atoms = system.spanOfMolecule(selectedComponent, selectedMolecule);
  Molecule &molecule = system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, selectedMolecule)];

  double3 displacement{};
  std::chrono::steady_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::Translation;
  Component &component = system.components[selectedComponent];

  // Randomly select a direction (0 for x, 1 for y, 2 for z)
  std::size_t selectedDirection = std::size_t(3.0 * random.uniform());

  // Get the maximum displacement allowed for the selected component
  double maxDisplacement = component.mc_moves_statistics.getMaxChange(move, selectedDirection);

  // Compute the displacement along the selected direction
  displacement[selectedDirection] = maxDisplacement * 2.0 * (random.uniform() - 0.5);

  // Update move statistics for the selected direction
  component.mc_moves_statistics.addTrial(move, selectedDirection);

  // Construct the trial molecule with the new displacement
  std::pair<Molecule, std::vector<Atom>> trialMolecule =
      component.translate(molecule, molecule_atoms, displacement);

  // Check if the trial molecule is inside blocked pockets
  if (system.insideBlockedPockets(component, trialMolecule.second))
  {
    return std::nullopt;
  }

  std::vector<double3> electricFieldMoleculeNew(molecule_atoms.size());
  std::vector<double3> electricFieldMoleculeOld(molecule_atoms.size());

  // Compute external field energy contribution
  time_begin = std::chrono::steady_clock::now();
  std::optional<RunningEnergy> externalFieldMolecule = Interactions::computeExternalFieldEnergyDifference(
      system.hasExternalField, system.forceField, system.simulationBox, 
      system.externalFieldInterpolationGrid, trialMolecule.second, molecule_atoms);
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::ExternalFieldMolecule] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::ExternalFieldMolecule] += (time_end - time_begin);
  if (!externalFieldMolecule.has_value()) return std::nullopt;

  // Compute framework-molecule energy contribution
  time_begin = std::chrono::steady_clock::now();
  std::optional<RunningEnergy> frameworkMolecule;
  if (system.forceField.computePolarization)
  {
    frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), electricFieldMoleculeNew, electricFieldMoleculeOld, trialMolecule.second,
        molecule_atoms);
  }
  else
  {
    frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
        system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
        system.spanOfFrameworkAtoms(), trialMolecule.second, molecule_atoms);
  }
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::FrameworkMolecule] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::FrameworkMolecule] += (time_end - time_begin);
  if (!frameworkMolecule.has_value()) return std::nullopt;

  // Compute molecule-molecule energy contribution. When molecule-molecule polarization is enabled the same
  // neighbor loop also produces (i) the inter-molecular electric field on the moved molecule (new and old) and
  // (ii) the change of the electric field on every other atom, so that the neighbor polarization energy can be
  // updated incrementally.
  time_begin = std::chrono::steady_clock::now();
  std::vector<double3> electricFieldNeighborDelta;
  std::optional<RunningEnergy> interMolecule;
  if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
  {
    electricFieldNeighborDelta.assign(system.spanOfMoleculeAtoms().size(), double3(0.0, 0.0, 0.0));
    interMolecule = Interactions::computeInterMolecularPolarizationElectricFieldDifference(
        system.forceField, system.simulationBox, electricFieldNeighborDelta, electricFieldMoleculeNew,
        electricFieldMoleculeOld, system.spanOfMoleculeAtoms(), trialMolecule.second, molecule_atoms);
  }
  else
  {
    interMolecule = Interactions::computeInterMolecularEnergyDifference(
        system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), trialMolecule.second, molecule_atoms);
  }
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::MoleculeMolecule] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::MoleculeMolecule] += (time_end - time_begin);
  if (!interMolecule.has_value()) return std::nullopt;

  // Compute Ewald energy contribution
  time_begin = std::chrono::steady_clock::now();
  RunningEnergy ewaldFourierEnergy;
  if (system.forceField.computePolarization)
  {
    ewaldFourierEnergy = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.trialEik, system.forceField, system.simulationBox, electricFieldMoleculeNew, electricFieldMoleculeOld,
        trialMolecule.second, molecule_atoms);
  }
  else
  {
    ewaldFourierEnergy = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
        system.simulationBox, trialMolecule.second, molecule_atoms);
  }
  time_end = std::chrono::steady_clock::now();
  component.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);

  RunningEnergy polarizationDifference;
  if (system.forceField.computePolarization)
  {
    // Polarization energy change of the moved molecule (framework + reciprocal [+ inter-molecular] field).
    polarizationDifference = Interactions::computePolarizationEnergyDifference(
        system.forceField, electricFieldMoleculeNew, electricFieldMoleculeOld, trialMolecule.second, molecule_atoms);

    // Polarization energy change of all other molecules whose field changed because this molecule moved.
    if (!system.forceField.omitInterPolarization)
    {
      polarizationDifference += Interactions::computePolarizationEnergyNeighborDifference(
          system.forceField, system.spanOfMoleculeElectricField(), electricFieldNeighborDelta,
          system.spanOfMoleculeAtoms());
    }
  }

  // Calculate the total energy difference
  RunningEnergy energyDifference = externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() +
                                   ewaldFourierEnergy + polarizationDifference;

  // Update move construction statistics
  component.mc_moves_statistics.addConstructed(move, selectedDirection);

  // Apply acceptance/rejection rule based on Metropolis criterion
  if (random.uniform() < std::exp(-system.beta * energyDifference.potentialEnergy()))
  {
    // Update acceptance statistics
    component.mc_moves_statistics.addAccepted(move, selectedDirection);

    // Accept the Ewald move
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);

    // Update the molecule atoms with the trial positions
    std::copy(trialMolecule.second.cbegin(), trialMolecule.second.cend(), molecule_atoms.begin());
    molecule = trialMolecule.first;

    // Commit the electric field to the stored (committed) field so that the running polarization energy stays
    // consistent with a full recomputation. In the framework-only model moving one molecule only changes its own
    // field; with molecule-molecule polarization the neighbor field changes are applied as well.
    if (system.forceField.computePolarization)
    {
      if (!system.forceField.omitInterPolarization)
      {
        std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
        for (std::size_t i = 0; i < storedElectricField.size(); ++i)
        {
          storedElectricField[i] += electricFieldNeighborDelta[i];
        }
      }
      std::span<double3> electricFieldMolecule = system.spanElectricFieldOld(selectedComponent, selectedMolecule);
      std::copy(electricFieldMoleculeNew.begin(), electricFieldMoleculeNew.end(), electricFieldMolecule.begin());
    }

    return energyDifference;
  };
  return std::nullopt;
}
