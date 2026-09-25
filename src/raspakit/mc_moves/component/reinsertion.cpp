module;

module mc_moves_reinsertion;

import std;

import component;
import atom;
import molecule;
import double3;
import double3x3;
import simd_quatd;
import simulationbox;
import cbmc;
import cbmc_results;
import cbmc_external_energy;
import cbmc_grow_context;
import randomnumbers;
import system;
import energy_status;
import energy_status_inter;
import running_energy;
import property_lambda_probability_histogram;
import property_widom;
import averages;
import forcefield;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import interactions_polarization;
import mc_moves_move_types;

std::optional<RunningEnergy> MC_Moves::reinsertionMove(RandomNumber &random, System &system,
                                                       std::size_t selectedComponent, std::size_t selectedMolecule)
{
  std::span<Atom> molecule_atoms = system.spanOfMolecule(selectedComponent, selectedMolecule);
  Molecule &molecule = system.moleculeData[system.moleculeIndexOfComponent(selectedComponent, selectedMolecule)];

  // Variables to record timing for performance measurement.
  std::chrono::steady_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::ReinsertionCBMC;
  Component &component = system.components[selectedComponent];

  // Increment move counts for reinsertion CBMC statistics.
  component.mc_moves_statistics.addTrial(move);

  // If no molecules of selected component are present, exit the move.
  if (system.numberOfMoleculesPerComponent[selectedComponent] == 0)
  {
    return std::nullopt;
  }

  const CBMC::GrowContext context = system.makeGrowContext();

  time_begin = std::chrono::steady_clock::now();
  // Attempt to grow the molecule using CBMC reinsertion.
  std::optional<CBMC::GrowResult> growData = CBMC::regrowMolecule(
      random, context, component, molecule, molecule_atoms, {.firstBead = CBMC::FirstBeadScheme::Reinsertion});
  time_end = std::chrono::steady_clock::now();
  // Record CPU time taken for the non-Ewald part of the move.
  component.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);

  // If growth was unsuccessful, exit the move.
  if (!growData) return std::nullopt;

  // Get the new molecule configuration.
  std::span<const Atom> newMolecule = std::span(growData->atoms.begin(), growData->atoms.end());

  std::vector<Atom> old_molecule = std::vector(molecule_atoms.begin(), molecule_atoms.end());
  std::vector<double3> old_electric_field = std::vector<double3>(old_molecule.size());
  std::vector<double3> new_electric_field = std::vector<double3>(old_molecule.size());

  // Increment the constructed moves count.
  component.mc_moves_statistics.addConstructed(move);

  // Retrace the old molecule configuration using CBMC retracing.
  time_begin = std::chrono::steady_clock::now();
  CBMC::RetraceResult retraceData = CBMC::retraceMolecule(
      random, context, component, molecule_atoms,
      {.firstBead = CBMC::FirstBeadScheme::Reinsertion,
       .storedR = growData->firstBeadStoredR,
       .skipBackgroundMolecule = molecule_atoms[component.startingBead].moleculeId});
  time_end = std::chrono::steady_clock::now();

  // Record CPU time taken for the retracing step.
  component.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);

  // Compute the energy difference in the Fourier space due to Ewald summation.
  time_begin = std::chrono::steady_clock::now();
  RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
      system.simulationBox, newMolecule, molecule_atoms);
  time_end = std::chrono::steady_clock::now();
  // Record CPU time taken for the Ewald Fourier part of the move.
  component.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);

  // Dual cut-off scheme: correct the grown and retraced configurations from the inner cut-off to the
  // full cut-offs.
  if (!CBMC::applyDualCutOffCorrection(context, component, *growData) ||
      !CBMC::applyDualCutOffCorrection(context, component, old_molecule, retraceData))
  {
    return std::nullopt;
  }

  std::vector<double3> electricFieldNeighborDelta;
  RunningEnergy polarizationDifference;
  if (system.forceField.computePolarization)
  {
    Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                  system.spanOfFrameworkAtoms(), new_electric_field,
                                                                  old_electric_field, growData->atoms, old_molecule);

    Interactions::computeEwaldFourierElectricFieldDifference(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.trialEik, system.forceField, system.simulationBox, new_electric_field, old_electric_field,
        growData->atoms, old_molecule);

    // Molecule-molecule polarization: add the inter-molecular field on the reinserted molecule and the change of
    // the field on every other molecule (the inter-molecular energy is already accounted for through CBMC, so the
    // returned energy is intentionally discarded here).
    if (!system.forceField.omitInterPolarization)
    {
      electricFieldNeighborDelta.assign(system.spanOfMoleculeAtoms().size(), double3(0.0, 0.0, 0.0));
      [[maybe_unused]] std::optional<RunningEnergy> interPolarizationEnergy =
          Interactions::computeInterMolecularPolarizationElectricFieldDifference(
              system.forceField, system.simulationBox, electricFieldNeighborDelta, new_electric_field,
              old_electric_field, system.spanOfMoleculeAtoms(), growData->atoms, old_molecule);
    }

    // Compute polarization energy difference
    polarizationDifference = Interactions::computePolarizationEnergyDifference(
        system.forceField, new_electric_field, old_electric_field, growData->atoms, old_molecule);

    if (!system.forceField.omitInterPolarization)
    {
      polarizationDifference += Interactions::computePolarizationEnergyNeighborDifference(
          system.forceField, system.spanOfMoleculeElectricField(), electricFieldNeighborDelta,
          system.spanOfMoleculeAtoms());
    }
  }

  // Metropolis acceptance with the weight ratio evaluated in log space: the raw Rosenbluth weights of a
  // long chain underflow to zero (their product spans hundreds of exp(-beta u) factors), which turns the
  // raw quotient into 0/0 = NaN (silent reject) or x/0 = inf (unconditional accept). The log difference
  // is exact for any chain length.
  double logAcceptance =
      -system.beta * (energyFourierDifference.potentialEnergy() + polarizationDifference.potentialEnergy()) +
      growData->logRosenbluthWeight - retraceData.logRosenbluthWeight;

  // Apply Metropolis acceptance criterion.
  if (random.uniform() < std::exp(logAcceptance))
  {
    // Move is accepted; update statistics and state.
    component.mc_moves_statistics.addAccepted(move);

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);
    std::copy(newMolecule.begin(), newMolecule.end(), molecule_atoms.begin());

    if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
    {
      std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
      for (std::size_t i = 0; i < storedElectricField.size(); ++i)
      {
        storedElectricField[i] += electricFieldNeighborDelta[i];
      }
    }

    std::span<double3> electricFieldMolecule = system.spanElectricFieldOld(selectedComponent, selectedMolecule);
    std::copy(new_electric_field.begin(), new_electric_field.end(), electricFieldMolecule.begin());

    molecule = growData->molecule;

    return (growData->energies - retraceData.energies) + energyFourierDifference + polarizationDifference;
  };

  // Move is rejected.
  return std::nullopt;
}
