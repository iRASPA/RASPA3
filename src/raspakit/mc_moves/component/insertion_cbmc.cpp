module;

module mc_moves_insertion_cbmc;

import std;

import component;
import molecule;
import atom;
import double3;
import double3x3;
import simd_quatd;
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

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::insertionMoveCBMC(RandomNumber& random, System& system,
                                                                             std::size_t selectedComponent)
{
  std::chrono::steady_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::SwapCBMC;
  Component& component = system.components[selectedComponent];

  // Set trial moleculeId to something that does not overlap with the current molecules
  // The 'insertMolecule' routine will place it after the last molecule of the component
  std::size_t selectedMolecule = system.numberOfMolecules();

  // Update move counts statistics for swap insertion move
  component.mc_moves_statistics.addTrial(move, 0);

  // Determine cutoff distances based on whether dual cutoff is used.
  double cutOffFrameworkVDW =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffFrameworkVDW;
  double cutOffMoleculeVDW =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffMoleculeVDW;
  double cutOffCoulomb =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffCoulomb;
  Component::GrowType growType = component.growType;

  const CBMC::GrowContext growContext{system.hasExternalField, system.forceField, system.simulationBox,
                                      system.interpolationGrids, system.externalFieldInterpolationGrid,
                                      system.framework, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
                                      system.beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb};

  // Attempt to grow a new molecule using CBMC
  time_begin = std::chrono::steady_clock::now();
  std::optional<ChainGrowData> growData = CBMC::growMoleculeSwapInsertion(
      random, growContext, component, selectedComponent, growType, selectedMolecule, 1.0, false, false);
  time_end = std::chrono::steady_clock::now();

  // Update CPU time statistics for the non-Ewald part of the move
  component.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);

  // If growth failed, reject the move
  if (!growData) return {std::nullopt, double3(0.0, 1.0, 0.0)};

  if (system.forceField.useDualCutOff)
  {
    // Dual cut-off scheme: correct the grown configuration from the inner cut-off to the full
    // cut-offs, so that Rosenbluth weight and energies behave as if grown at the full cut-offs.
    std::optional<RunningEnergy> correctionNew =
        CBMC::computeDualCutOffCorrection(growContext, component, growData->atoms);
    if (!correctionNew.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    growData->energies += correctionNew.value();
    growData->RosenbluthWeight *= std::exp(-system.beta * correctionNew->potentialEnergy());
  }

  std::span<const Atom> newMolecule = std::span(growData->atoms.begin(), growData->atoms.end());
  std::vector<double3> new_electric_field = std::vector<double3>(newMolecule.size());

  // Check if the new molecule is inside blocked pockets
  if (system.insideBlockedPockets(system.components[selectedComponent], newMolecule))
  {
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  // Update statistics for successfully constructed molecules
  system.components[selectedComponent].mc_moves_statistics.addConstructed(move, 0);

  // Compute energy difference due to Ewald Fourier components
  time_begin = std::chrono::steady_clock::now();
  RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, newMolecule, {}, system.netCharge);
  time_end = std::chrono::steady_clock::now();

  // Update CPU time statistics for the Ewald part of the move
  component.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);

  // Compute tail energy difference due to long-range corrections
  time_begin = std::chrono::steady_clock::now();
  RunningEnergy tailEnergyDifference =
      Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                              system.spanOfMoleculeAtoms(), newMolecule, {}) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.spanOfFrameworkAtoms(), newMolecule, {});
  time_end = std::chrono::steady_clock::now();

  // Update CPU time statistics for the tail corrections
  component.mc_moves_cputime[move][Move::Timing::Tail] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::Tail] += (time_end - time_begin);

  std::vector<double3> electricFieldNeighborDelta;
  RunningEnergy polarizationDifference;
  if (system.forceField.computePolarization)
  {
    Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                  system.spanOfFrameworkAtoms(), new_electric_field, {},
                                                                  growData->atoms, {});

    Interactions::computeEwaldFourierElectricFieldDifference(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.totalEik, system.forceField, system.simulationBox, new_electric_field, {}, growData->atoms, {});

    // Molecule-molecule polarization: field on the inserted molecule plus the change of the field on every
    // existing molecule (inter-molecular energy is already accounted for through CBMC, discard returned energy).
    if (!system.forceField.omitInterPolarization)
    {
      electricFieldNeighborDelta.assign(system.spanOfMoleculeAtoms().size(), double3(0.0, 0.0, 0.0));
      [[maybe_unused]] std::optional<RunningEnergy> interPolarizationEnergy =
          Interactions::computeInterMolecularPolarizationElectricFieldDifference(
              system.forceField, system.simulationBox, electricFieldNeighborDelta, new_electric_field,
              std::span<double3>{}, system.spanOfMoleculeAtoms(), growData->atoms, {});
    }

    // Compute polarization energy difference
    polarizationDifference = Interactions::computePolarizationEnergyDifference(system.forceField, new_electric_field,
                                                                               {}, growData->atoms, {});

    if (!system.forceField.omitInterPolarization)
    {
      polarizationDifference += Interactions::computePolarizationEnergyNeighborDifference(
          system.forceField, system.spanOfMoleculeElectricField(), electricFieldNeighborDelta,
          system.spanOfMoleculeAtoms());
    }
  }

  // Calculate correction factor for Ewald energy difference
  double correctionFactorEwald =
      std::exp(-system.beta * (energyFourierDifference.potentialEnergy() + tailEnergyDifference.potentialEnergy() +
                               polarizationDifference.potentialEnergy()));

  // Compute the acceptance probability pre-factor
  double fugacity = component.molFraction * component.fugacityCoefficient.value_or(1.0) * system.pressure;
  double idealGasRosenbluthWeight = component.idealGasRosenbluthWeight.value_or(1.0);
  double preFactor = correctionFactorEwald * system.beta * fugacity * system.simulationBox.volume /
                     double(1 + system.numberOfIntegerMoleculesPerComponent[selectedComponent]);

  // Calculate the acceptance probability Pacc
  double Pacc = preFactor * growData->RosenbluthWeight / idealGasRosenbluthWeight;

  std::size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];

  // Check if TMMC is enabled and macrostate limit is not exceeded
  if (system.tmmc.doTMMC && system.tmmc.rejectOutOfBound && oldN >= system.tmmc.maxMacrostate)
  {
    return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
  }

  const std::size_t newN = oldN == std::numeric_limits<std::size_t>::max() ? oldN : oldN + 1;
  double biasTransitionMatrix = system.tmmc.biasFactor(newN, oldN);

  // Apply acceptance/rejection criterion
  if (random.uniform() < biasTransitionMatrix * Pacc)
  {
    // Move accepted; update acceptance statistics
    component.mc_moves_statistics.addAccepted(move, 0);

    // Accept Ewald move and insert the new molecule into the system
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);

    if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
    {
      std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
      for (std::size_t i = 0; i < storedElectricField.size(); ++i)
      {
        storedElectricField[i] += electricFieldNeighborDelta[i];
      }
    }

    system.insertMoleculePolarization(selectedComponent, growData->molecule, growData->atoms, new_electric_field);

    return {growData->energies + energyFourierDifference + tailEnergyDifference + polarizationDifference,
            double3(0.0, 1.0 - Pacc, Pacc)};
  };

  // Move rejected
  return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
}
