module;

module mc_moves_deletion;

import std;

import double3;
import double3x3;
import simd_quatd;
import component;
import atom;
import simulationbox;
import cbmc;
import cbmc_chain_data;
import randomnumbers;
import system;
import energy_status;
import energy_status_inter;
import property_lambda_probability_histogram;
import property_widom;
import averages;
import running_energy;
import forcefield;
import transition_matrix;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import interactions_polarization;
import mc_moves_move_types;

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::deletionMove(RandomNumber& random, System& system,
                                                                        std::size_t selectedComponent,
                                                                        std::size_t selectedMolecule)
{
  std::chrono::steady_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::Swap;
  Component& component = system.components[selectedComponent];

  // Increment swap deletion move counts for the selected component
  component.mc_moves_statistics.addTrial(move, 1);

  if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] > 0)
  {
    std::span<Atom> molecule = system.spanOfMolecule(selectedComponent, selectedMolecule);

    // Copy the current electric field if polarization is computed
    std::vector<double3> electricFieldMoleculeOld(molecule.size());
    std::vector<double3> electricFieldNeighborDelta;

    // Compute external field energy contribution
    std::optional<RunningEnergy> externalFieldMolecule = Interactions::computeExternalFieldEnergyDifference(
        system.hasExternalField, system.forceField, system.simulationBox, system.externalFieldInterpolationGrid, {},
        molecule);
    if (!externalFieldMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    // Compute framework-molecule energy contribution
    std::optional<RunningEnergy> frameworkMolecule;
    if (system.forceField.computePolarization)
    {
      frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
          system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
          system.spanOfFrameworkAtoms(), {}, electricFieldMoleculeOld, {}, molecule);
    }
    else
    {
      frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
          system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
          system.spanOfFrameworkAtoms(), {}, molecule);
    }
    if (!frameworkMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    // Compute molecule-molecule energy contribution (and, for molecule-molecule polarization, the electric field
    // on the deleted molecule as well as the change of the field on every remaining molecule)
    std::optional<RunningEnergy> interMolecule;
    if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
    {
      electricFieldNeighborDelta.assign(system.spanOfMoleculeAtoms().size(), double3(0.0, 0.0, 0.0));
      interMolecule = Interactions::computeInterMolecularPolarizationElectricFieldDifference(
          system.forceField, system.simulationBox, electricFieldNeighborDelta, std::span<double3>{},
          electricFieldMoleculeOld, system.spanOfMoleculeAtoms(), {}, molecule);
    }
    else
    {
      interMolecule = Interactions::computeInterMolecularEnergyDifference(system.forceField, system.simulationBox,
                                                                          system.spanOfMoleculeAtoms(), {}, molecule);
    }
    if (!interMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    // Compute Ewald Fourier energy difference
    time_begin = std::chrono::steady_clock::now();
    RunningEnergy energyFourierDifference;
    if (system.forceField.computePolarization)
    {
      energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
          system.totalEik, system.forceField, system.simulationBox, {}, electricFieldMoleculeOld, {}, molecule,
          system.netCharge);
    }
    else
    {
      energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
          system.simulationBox, {}, molecule, system.netCharge);
    }
    time_end = std::chrono::steady_clock::now();

    // Update CPU time statistics for Ewald calculations
    component.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);

    // Compute tail correction energy difference
    time_begin = std::chrono::steady_clock::now();
    [[maybe_unused]] RunningEnergy tailEnergyDifference =
        Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                                system.spanOfMoleculeAtoms(), {}, molecule) +
        Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                   system.spanOfFrameworkAtoms(), {}, molecule);
    time_end = std::chrono::steady_clock::now();

    // Update CPU time statistics for tail corrections
    component.mc_moves_cputime[move][Move::Timing::Tail] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::Tail] += (time_end - time_begin);

    RunningEnergy polarizationDifference;
    if (system.forceField.computePolarization)
    {
      // Polarization energy of the deleted molecule
      polarizationDifference = Interactions::computePolarizationEnergyDifference(
          system.forceField, {}, electricFieldMoleculeOld, {}, molecule);

      // Polarization energy change of all remaining molecules whose field changes due to the deletion
      if (!system.forceField.omitInterPolarization)
      {
        polarizationDifference += Interactions::computePolarizationEnergyNeighborDifference(
            system.forceField, system.spanOfMoleculeElectricField(), electricFieldNeighborDelta,
            system.spanOfMoleculeAtoms());
      }
    }

    // Get the total difference in energy
    RunningEnergy energyDifference = externalFieldMolecule.value() + frameworkMolecule.value() + interMolecule.value() +
                                     energyFourierDifference + tailEnergyDifference + polarizationDifference;

    // Increment constructed swap deletion move counts
    component.mc_moves_statistics.addConstructed(move, 1);

    // Calculate the acceptance probability
    double fugacity = component.molFraction * component.fugacityCoefficient.value_or(1.0) * system.pressure;
    double preFactor = double(system.numberOfIntegerMoleculesPerComponent[selectedComponent]) /
                       (system.beta * fugacity * system.simulationBox.volume);
    double Pacc = preFactor * std::exp(-system.beta * energyDifference.potentialEnergy());
    std::size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];

    // Check if TMMC is enabled and if new state is below minimum macrostate
    if (system.tmmc.doTMMC && system.tmmc.rejectOutOfBound && oldN <= system.tmmc.minMacrostate)
    {
      return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
    }

    const std::size_t newN = oldN == 0 ? 0 : oldN - 1;
    double biasTransitionMatrix = system.tmmc.biasFactor(newN, oldN);

    // Apply acceptance/rejection rule
    if (random.uniform() < biasTransitionMatrix * Pacc)
    {
      // Move accepted; update acceptance statistics
      component.mc_moves_statistics.addAccepted(move, 1);

      // Accept Ewald move and delete molecule from system
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);

      // Apply the field changes on the remaining molecules before the deleted molecule (and its field) is removed.
      if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
      {
        std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
        for (std::size_t i = 0; i < storedElectricField.size(); ++i)
        {
          storedElectricField[i] += electricFieldNeighborDelta[i];
        }
      }

      system.deleteMolecule(selectedComponent, selectedMolecule, molecule);

      return {-energyDifference, double3(Pacc, 1.0 - Pacc, 0.0)};
    };
    return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
  }

  // No molecules to delete; return default values
  return {std::nullopt, double3(0.0, 1.0, 0.0)};
}
