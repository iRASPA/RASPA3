module;

module mc_moves_identity_change;

import std;

import component;
import atom;
import molecule;
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

std::optional<RunningEnergy> MC_Moves::identityChangeMove(RandomNumber &random, System &system,
                                                          std::size_t selectedComponent)
{
  std::chrono::steady_clock::time_point time_begin, time_end;
  Move::Types move = Move::Types::IdentityChangeCBMC;
  Component &startComponent = system.components[selectedComponent];

  if (startComponent.identityChanges.empty())
  {
    return std::nullopt;
  }

  std::size_t oldComponent = selectedComponent;
  std::size_t newComponent =
      startComponent.identityChanges[random.uniform_integer(0, startComponent.identityChanges.size() - 1)];

  if (newComponent >= system.components.size())
  {
    return std::nullopt;
  }

  if (newComponent == selectedComponent)
  {
    return std::nullopt;
  }

  if (system.components[newComponent].type != system.components[oldComponent].type)
  {
    return std::nullopt;
  }

  // Guarantee detailed balance by randomly swapping old and new components.
  if (random.uniform() < 0.5)
  {
    std::swap(oldComponent, newComponent);
  }

  Component &oldComponentData = system.components[oldComponent];
  Component &newComponentData = system.components[newComponent];

  if (system.numberOfIntegerMoleculesPerComponent[oldComponent] == 0)
  {
    return std::nullopt;
  }

  oldComponentData.mc_moves_statistics.addTrial(move);

  const std::size_t selectedMoleculeOld = system.randomIntegerMoleculeOfComponent(random, oldComponent);
  std::span<Atom> oldMoleculeAtoms = system.spanOfMolecule(oldComponent, selectedMoleculeOld);
  const Atom &oldStartingBead = oldMoleculeAtoms[oldComponentData.startingBead];

  // Determine cutoff distances based on whether dual cutoff is used.
  double cutOffFrameworkVDW =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffFrameworkVDW;
  double cutOffMoleculeVDW =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffMoleculeVDW;
  double cutOffCoulomb =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffCoulomb;
  Component::GrowType newGrowType = newComponentData.growType;
  Component::GrowType oldGrowType = oldComponentData.growType;

  const std::size_t oldGlobalMoleculeId =
      system.moleculeIndexOfComponent(oldComponent, selectedMoleculeOld);
  const std::size_t trialMoleculeId = system.numberOfMolecules();
  const std::make_signed_t<std::size_t> skipBackgroundMolecule =
      static_cast<std::make_signed_t<std::size_t>>(oldGlobalMoleculeId);

  const CBMC::GrowContext growContext{system.hasExternalField, system.forceField, system.simulationBox,
                                      system.interpolationGrids, system.externalFieldInterpolationGrid,
                                      system.framework, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
                                      system.beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb};

  time_begin = std::chrono::steady_clock::now();
  std::optional<ChainGrowData> growData = CBMC::growMoleculeIdentityChangeInsertion(
      random, growContext, newComponentData, newComponent, newGrowType, trialMoleculeId, oldStartingBead, 1.0, false,
      false, skipBackgroundMolecule);
  time_end = std::chrono::steady_clock::now();
  oldComponentData.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);

  if (!growData)
  {
    return std::nullopt;
  }

  if (system.forceField.useDualCutOff)
  {
    // Dual cut-off scheme: correct the grown configuration from the inner cut-off to the full
    // cut-offs, using the same background (the old molecule excluded) as the growth.
    std::optional<RunningEnergy> correctionNew = CBMC::computeDualCutOffCorrection(
        growContext, newComponentData, growData->atoms, skipBackgroundMolecule);
    if (!correctionNew.has_value())
    {
      return std::nullopt;
    }

    growData->energies += correctionNew.value();
    growData->RosenbluthWeight *= std::exp(-system.beta * correctionNew->potentialEnergy());
  }

  std::span<const Atom> newMolecule = std::span(growData->atoms.begin(), growData->atoms.end());
  std::vector<Atom> old_molecule(oldMoleculeAtoms.begin(), oldMoleculeAtoms.end());
  std::vector<double3> old_electric_field(old_molecule.size());
  std::vector<double3> new_electric_field(newMolecule.size());

  if (system.insideBlockedPockets(newComponentData, newMolecule))
  {
    return std::nullopt;
  }

  oldComponentData.mc_moves_statistics.addConstructed(move);

  time_begin = std::chrono::steady_clock::now();
  ChainRetraceData retraceData = CBMC::retraceMoleculeIdentityChangeDeletion(
      random, growContext, oldComponentData, oldGrowType, oldMoleculeAtoms);
  time_end = std::chrono::steady_clock::now();
  oldComponentData.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);

  if (system.forceField.useDualCutOff)
  {
    // Dual cut-off scheme: correct the retraced configuration from the inner cut-off to the full
    // cut-offs (the old molecule excludes itself from the background through its molecule id).
    std::optional<RunningEnergy> correctionOld =
        CBMC::computeDualCutOffCorrection(growContext, oldComponentData, old_molecule);
    if (!correctionOld.has_value())
    {
      return std::nullopt;
    }

    retraceData.energies += correctionOld.value();
    retraceData.RosenbluthWeight *= std::exp(-system.beta * correctionOld->potentialEnergy());
  }

  time_begin = std::chrono::steady_clock::now();
  RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, newMolecule, oldMoleculeAtoms, system.netCharge);
  time_end = std::chrono::steady_clock::now();
  oldComponentData.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);

  time_begin = std::chrono::steady_clock::now();
  RunningEnergy tailEnergyDifference =
      Interactions::computeInterMolecularTailEnergyDifferenceAddRemove(
          system.forceField, system.simulationBox, system.totalNumberOfPseudoAtoms, newComponentData,
          oldComponentData) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.spanOfFrameworkAtoms(), newMolecule,
                                                                 oldMoleculeAtoms);
  time_end = std::chrono::steady_clock::now();
  oldComponentData.mc_moves_cputime[move][Move::Timing::Tail] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::Tail] += (time_end - time_begin);

  std::vector<double3> electricFieldNeighborDelta;
  RunningEnergy polarizationDifference;
  if (system.forceField.computePolarization)
  {
    Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                  system.spanOfFrameworkAtoms(), new_electric_field,
                                                                  old_electric_field, growData->atoms, old_molecule);

    Interactions::computeEwaldFourierElectricFieldDifference(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.totalEik, system.forceField, system.simulationBox, new_electric_field, old_electric_field,
        growData->atoms, old_molecule);

    // Molecule-molecule polarization: the old and new molecule occupy the same "slot", so the grown atoms are given
    // the old molecule's id for the field-difference so that the transformed molecule's own atoms are skipped (only
    // the field change on the surrounding molecules is retained). The inter-molecular energy is already accounted for
    // through CBMC, so the returned energy is discarded.
    if (!system.forceField.omitInterPolarization)
    {
      std::vector<Atom> newAtomsForField(growData->atoms.begin(), growData->atoms.end());
      if (!old_molecule.empty())
      {
        for (Atom &atom : newAtomsForField) atom.moleculeId = old_molecule.front().moleculeId;
      }

      electricFieldNeighborDelta.assign(system.spanOfMoleculeAtoms().size(), double3(0.0, 0.0, 0.0));
      [[maybe_unused]] std::optional<RunningEnergy> interPolarizationEnergy =
          Interactions::computeInterMolecularPolarizationElectricFieldDifference(
              system.forceField, system.simulationBox, electricFieldNeighborDelta, new_electric_field,
              old_electric_field, system.spanOfMoleculeAtoms(), newAtomsForField, old_molecule);
    }

    polarizationDifference = Interactions::computePolarizationEnergyDifference(
        system.forceField, new_electric_field, old_electric_field, growData->atoms, old_molecule);

    if (!system.forceField.omitInterPolarization)
    {
      polarizationDifference += Interactions::computePolarizationEnergyNeighborDifference(
          system.forceField, system.spanOfMoleculeElectricField(), electricFieldNeighborDelta,
          system.spanOfMoleculeAtoms());
    }
  }

  const double correctionFactorEwald =
      std::exp(-system.beta * (energyFourierDifference.potentialEnergy() + polarizationDifference.potentialEnergy()));

  const double idealGasRosenbluthWeightNew = newComponentData.idealGasRosenbluthWeight.value_or(1.0);
  const double idealGasRosenbluthWeightOld = oldComponentData.idealGasRosenbluthWeight.value_or(1.0);
  const double fugacityNew =
      newComponentData.molFraction * newComponentData.fugacityCoefficient.value_or(1.0) * system.pressure;
  const double fugacityOld =
      oldComponentData.molFraction * oldComponentData.fugacityCoefficient.value_or(1.0) * system.pressure;
  const double numberOfMoleculesOld =
      static_cast<double>(system.numberOfIntegerMoleculesPerComponent[oldComponent]);
  const double numberOfMoleculesNew =
      static_cast<double>(system.numberOfIntegerMoleculesPerComponent[newComponent]);

  const double rosenbluthNew =
      growData->RosenbluthWeight * std::exp(-system.beta * tailEnergyDifference.potentialEnergy());
  const double acceptanceProbability =
      correctionFactorEwald * (rosenbluthNew / idealGasRosenbluthWeightNew) * fugacityNew * numberOfMoleculesOld /
      ((retraceData.RosenbluthWeight / idealGasRosenbluthWeightOld) * fugacityOld * (numberOfMoleculesNew + 1.0));

  if (random.uniform() < acceptanceProbability)
  {
    oldComponentData.mc_moves_statistics.addAccepted(move);

    const RunningEnergy energyDifference = (growData->energies - retraceData.energies) + energyFourierDifference +
                                           tailEnergyDifference + polarizationDifference;

    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);

    // Apply the field changes on the surrounding molecules before the old molecule is removed and the new one added.
    if (system.forceField.computePolarization && !system.forceField.omitInterPolarization)
    {
      std::span<double3> storedElectricField = system.spanOfMoleculeElectricField();
      for (std::size_t i = 0; i < storedElectricField.size(); ++i)
      {
        storedElectricField[i] += electricFieldNeighborDelta[i];
      }
    }

    std::vector<Atom> acceptedAtoms(growData->atoms.begin(), growData->atoms.end());
    for (Atom &atom : acceptedAtoms)
    {
      atom.componentId = static_cast<std::uint8_t>(newComponent);
    }

    Molecule acceptedMolecule = growData->molecule;
    acceptedMolecule.componentId = newComponent;

    system.deleteMolecule(oldComponent, selectedMoleculeOld, oldMoleculeAtoms);
    if (system.forceField.computePolarization)
    {
      system.insertMoleculePolarization(newComponent, acceptedMolecule, acceptedAtoms, new_electric_field);
    }
    else
    {
      system.insertMolecule(newComponent, acceptedMolecule, acceptedAtoms);
    }

    return energyDifference;
  }

  return std::nullopt;
}
