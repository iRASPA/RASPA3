module;

module mc_moves_deletion_cbmc;

import std;

import double3;
import double3x3;
import simd_quatd;
import component;
import atom;
import simulationbox;
import cbmc;
import cbmc_chain_data;
import cbmc_interactions;
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

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::deletionMoveCBMC(RandomNumber& random, System& system,
                                                                            std::size_t selectedComponent,
                                                                            std::size_t selectedMolecule)
{
  std::chrono::steady_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::SwapCBMC;
  Component& component = system.components[selectedComponent];

  // Increment the count of swap deletion moves for the selected component
  component.mc_moves_statistics.addTrial(move, 1);

  // Proceed only if there is at least one molecule of the selected component
  if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] > 0)
  {
    // Get a reference to the molecule being deleted
    std::span<Atom> molecule = system.spanOfMolecule(selectedComponent, selectedMolecule);
    std::copy(system.electricField.begin(), system.electricField.end(), system.electricFieldNew.begin());
    // std::span<double3> electricFieldMoleculeNew = system.spanElectricFieldNew(selectedComponent, selectedMolecule);

    // Determine cutoff distances based on whether dual cutoff is used.
    double cutOffFrameworkVDW =
        system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffFrameworkVDW;
    double cutOffMoleculeVDW =
        system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffMoleculeVDW;
    double cutOffCoulomb =
        system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffCoulomb;
    Component::GrowType growType = component.growType;

    const CBMC::GrowContext retraceContext{system.hasExternalField, system.forceField, system.simulationBox,
                                           system.interpolationGrids, system.externalFieldInterpolationGrid,
                                           system.framework, system.spanOfFrameworkAtoms(),
                                           system.spanOfMoleculeAtoms(), system.beta, cutOffFrameworkVDW,
                                           cutOffMoleculeVDW, cutOffCoulomb};

    // Retrace the molecule for the swap deletion using CBMC algorithm
    time_begin = std::chrono::steady_clock::now();
    ChainRetraceData retraceData = CBMC::retraceMoleculeSwapDeletion(
        random, retraceContext, system.components[selectedComponent], growType, molecule);
    time_end = std::chrono::steady_clock::now();

    // Update the CPU time statistics for the non-Ewald part of the move
    component.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);

    if (system.forceField.useDualCutOff)
    {
      // Dual cut-off scheme: correct the retraced configuration from the inner cut-off to the full
      // cut-offs, so that Rosenbluth weight and energies behave as if retraced at the full cut-offs.
      std::vector<Atom> oldMolecule = std::vector(molecule.begin(), molecule.end());
      std::optional<RunningEnergy> correctionOld =
          CBMC::computeDualCutOffCorrection(retraceContext, component, oldMolecule);
      if (!correctionOld.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

      retraceData.energies += correctionOld.value();
      retraceData.RosenbluthWeight *= std::exp(-system.beta * correctionOld->potentialEnergy());
    }

    // Compute the energy difference in Fourier space due to the deletion
    time_begin = std::chrono::steady_clock::now();
    RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
        system.simulationBox, {}, molecule, system.netCharge);
    time_end = std::chrono::steady_clock::now();
    // Update the CPU time statistics for the Ewald part of the move
    component.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);

    // Compute the tail energy difference due to the deletion
    time_begin = std::chrono::steady_clock::now();
    [[maybe_unused]] RunningEnergy tailEnergyDifference =
        Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                                system.spanOfMoleculeAtoms(), {}, molecule) +
        Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                   system.spanOfFrameworkAtoms(), {}, molecule);
    time_end = std::chrono::steady_clock::now();
    // Update the CPU time statistics for the tail corrections
    component.mc_moves_cputime[move][Move::Timing::Tail] += (time_end - time_begin);
    system.mc_moves_cputime[move][Move::Timing::Tail] += (time_end - time_begin);

    // Update the constructed count for the move statistics
    component.mc_moves_statistics.addConstructed(move, 1);

    std::vector<double3> electricFieldNeighborDelta;
    RunningEnergy polarizationDifference;
    if (system.forceField.computePolarization)
    {
      std::vector<Atom> old_molecule = std::vector(molecule.begin(), molecule.end());
      std::vector<double3> old_electric_field = std::vector<double3>(old_molecule.size());

      Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                    system.spanOfFrameworkAtoms(), {},
                                                                    old_electric_field, {}, old_molecule);

      Interactions::computeEwaldFourierElectricFieldDifference(
          system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
          system.trialEik, system.forceField, system.simulationBox, {}, old_electric_field, {}, old_molecule);

      // Molecule-molecule polarization: field on the deleted molecule plus the change of the field on every
      // remaining molecule (inter-molecular energy is already accounted for through CBMC, discard returned energy).
      if (!system.forceField.omitInterPolarization)
      {
        electricFieldNeighborDelta.assign(system.spanOfMoleculeAtoms().size(), double3(0.0, 0.0, 0.0));
        [[maybe_unused]] std::optional<RunningEnergy> interPolarizationEnergy =
            Interactions::computeInterMolecularPolarizationElectricFieldDifference(
                system.forceField, system.simulationBox, electricFieldNeighborDelta, std::span<double3>{},
                old_electric_field, system.spanOfMoleculeAtoms(), {}, old_molecule);
      }

      // Compute polarization energy difference
      polarizationDifference = Interactions::computePolarizationEnergyDifference(system.forceField, {},
                                                                                 old_electric_field, {}, old_molecule);

      if (!system.forceField.omitInterPolarization)
      {
        polarizationDifference += Interactions::computePolarizationEnergyNeighborDifference(
            system.forceField, system.spanOfMoleculeElectricField(), electricFieldNeighborDelta,
            system.spanOfMoleculeAtoms());
      }
    }

    // Calculate the correction factor for Ewald summation
    double correctionFactorEwald =
        std::exp(-system.beta * (energyFourierDifference.potentialEnergy() + tailEnergyDifference.potentialEnergy() +
                                 polarizationDifference.potentialEnergy()));

    // Compute acceptance probability factors
    double fugacity = component.molFraction * component.fugacityCoefficient.value_or(1.0) * system.pressure;
    double idealGasRosenbluthWeight = component.idealGasRosenbluthWeight.value_or(1.0);
    double preFactor = correctionFactorEwald * double(system.numberOfIntegerMoleculesPerComponent[selectedComponent]) /
                       (system.beta * fugacity * system.simulationBox.volume);
    double Pacc = preFactor * idealGasRosenbluthWeight / retraceData.RosenbluthWeight;
    std::size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];

    // Check if the new macrostate is within the allowed TMMC range
    if (system.tmmc.doTMMC && system.tmmc.rejectOutOfBound && oldN <= system.tmmc.minMacrostate)
    {
      return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
    }

    const std::size_t newN = oldN == 0 ? 0 : oldN - 1;
    double biasTransitionMatrix = system.tmmc.biasFactor(newN, oldN);

    // Apply acceptance/rejection rule
    if (random.uniform() < biasTransitionMatrix * Pacc)
    {
      component.mc_moves_statistics.addAccepted(move, 1);

      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);

      if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
      {
        std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
        for (std::size_t i = 0; i < storedElectricField.size(); ++i)
        {
          storedElectricField[i] += electricFieldNeighborDelta[i];
        }
      }

      system.deleteMolecule(selectedComponent, selectedMolecule, molecule);

      return {retraceData.energies - energyFourierDifference - tailEnergyDifference - polarizationDifference,
              double3(Pacc, 1.0 - Pacc, 0.0)};
    };
    return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
  }

  // Return default values if no molecules are available for deletion
  return {std::nullopt, double3(0.0, 1.0, 0.0)};
}
