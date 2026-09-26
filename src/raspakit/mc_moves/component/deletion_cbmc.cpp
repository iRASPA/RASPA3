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
import mc_moves_cputime;

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::deletionMoveCBMC(RandomNumber& random, System& system,
                                                                            std::size_t selectedComponent,
                                                                            std::size_t selectedMolecule)
{
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

    const CBMC::GrowContext retraceContext = system.makeGrowContext();

    // Retrace the molecule for the swap deletion using CBMC algorithm
    CBMC::RetraceResult retraceData =
        timed(system, component, move, Move::Timing::NonEwald,
              [&] { return CBMC::retraceMolecule(random, retraceContext, component, molecule); });

    // Compute the energy difference in Fourier space due to the deletion
    RunningEnergy energyFourierDifference =
        timed(system, component, move, Move::Timing::Ewald,
              [&]
              {
                return Interactions::energyDifferenceEwaldFourier(
                    system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik,
                    system.forceField, system.simulationBox, {}, molecule, system.netCharge);
              });

    // Compute the tail energy difference due to the deletion
    [[maybe_unused]] RunningEnergy tailEnergyDifference =
        timed(system, component, move, Move::Timing::Tail,
              [&]
              {
                return Interactions::computeInterMolecularTailEnergyDifference(
                           system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), {}, molecule) +
                       Interactions::computeFrameworkMoleculeTailEnergyDifference(
                           system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), {}, molecule);
              });

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
    // Rosenbluth weight through its exact logarithm: the raw weight of a long chain underflows to zero,
    // which would turn this quotient into +inf and unconditionally accept every deletion.
    double Pacc =
        std::exp(std::log(preFactor) + std::log(idealGasRosenbluthWeight) - retraceData.logRosenbluthWeight);
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
