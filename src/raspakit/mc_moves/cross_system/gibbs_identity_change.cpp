module;

module mc_moves_gibbs_identity_change;

import std;

import randomnumbers;
import running_energy;
import system;
import atom;
import molecule;
import component;
import double3;
import cbmc;
import cbmc_chain_data;
import cbmc_interactions;
import forcefield;
import simulationbox;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_polarization;
import mc_moves_move_types;

namespace
{

struct BoxIdentityChangeData
{
  ChainGrowData growData;
  ChainRetraceData retraceData;
  RunningEnergy energyFourierDifference;
  RunningEnergy tailEnergyDifference;
  RunningEnergy polarizationDifference;
  double correctionFactorEwald{1.0};
  double rosenbluthNew{1.0};
  double rosenbluthOld{1.0};
  std::size_t oldComponent{};
  std::size_t newComponent{};
  std::size_t selectedMoleculeOld{};
  std::span<Atom> oldMoleculeAtoms;
  std::vector<Atom> oldMoleculeCopy;
  std::vector<double3> oldElectricField;
  std::vector<double3> newElectricField;
};

bool performBoxIdentityChange(RandomNumber& random, System& system, Move::Types move, std::size_t oldComponent,
                              std::size_t newComponent, BoxIdentityChangeData& data)
{
  Component& oldComponentData = system.components[oldComponent];
  Component& newComponentData = system.components[newComponent];

  data.oldComponent = oldComponent;
  data.newComponent = newComponent;
  data.selectedMoleculeOld = system.randomIntegerMoleculeOfComponent(random, oldComponent);
  data.oldMoleculeAtoms = system.spanOfMolecule(oldComponent, data.selectedMoleculeOld);
  const Atom& oldStartingBead = data.oldMoleculeAtoms[oldComponentData.startingBead];

  // Determine cutoff distances based on whether dual cutoff is used.
  const double cutOffFrameworkVDW =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffFrameworkVDW;
  const double cutOffMoleculeVDW =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffMoleculeVDW;
  const double cutOffCoulomb =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffCoulomb;
  const Component::GrowType newGrowType = newComponentData.growType;
  const Component::GrowType oldGrowType = oldComponentData.growType;

  const std::size_t oldGlobalMoleculeId = system.moleculeIndexOfComponent(oldComponent, data.selectedMoleculeOld);
  const std::size_t trialMoleculeId = system.numberOfMolecules();
  const std::make_signed_t<std::size_t> skipBackgroundMolecule =
      static_cast<std::make_signed_t<std::size_t>>(oldGlobalMoleculeId);

  const CBMC::GrowContext growContext{system.hasExternalField, system.forceField, system.simulationBox,
                                      system.interpolationGrids, system.externalFieldInterpolationGrid,
                                      system.framework, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
                                      system.beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb};

  std::chrono::steady_clock::time_point time_begin = std::chrono::steady_clock::now();
  std::optional<ChainGrowData> growData = CBMC::growMoleculeIdentityChangeInsertion(
      random, growContext, newComponentData, newComponent, newGrowType, trialMoleculeId, oldStartingBead, 1.0, false,
      false, skipBackgroundMolecule);
  std::chrono::steady_clock::time_point time_end = std::chrono::steady_clock::now();
  oldComponentData.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);

  if (!growData)
  {
    return false;
  }

  if (system.forceField.useDualCutOff)
  {
    // Dual cut-off scheme: correct the grown configuration from the inner cut-off to the full
    // cut-offs, using the same background (the old molecule excluded) as the growth.
    std::optional<RunningEnergy> correctionNew = CBMC::computeDualCutOffCorrection(
        growContext, newComponentData, growData->atoms, skipBackgroundMolecule);
    if (!correctionNew.has_value())
    {
      return false;
    }

    growData->energies += correctionNew.value();
    growData->RosenbluthWeight *= std::exp(-system.beta * correctionNew->potentialEnergy());
  }

  data.growData = std::move(*growData);
  const std::span<const Atom> newMolecule(data.growData.atoms.begin(), data.growData.atoms.end());

  data.oldMoleculeCopy.assign(data.oldMoleculeAtoms.begin(), data.oldMoleculeAtoms.end());
  data.oldElectricField.resize(data.oldMoleculeCopy.size());
  data.newElectricField.resize(newMolecule.size());

  if (system.insideBlockedPockets(newComponentData, newMolecule))
  {
    return false;
  }

  oldComponentData.mc_moves_statistics.addConstructed(move);

  time_begin = std::chrono::steady_clock::now();
  data.retraceData = CBMC::retraceMoleculeIdentityChangeDeletion(random, growContext, oldComponentData, oldGrowType,
                                                                 data.oldMoleculeAtoms);
  time_end = std::chrono::steady_clock::now();
  oldComponentData.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::NonEwald] += (time_end - time_begin);

  if (system.forceField.useDualCutOff)
  {
    // Dual cut-off scheme: correct the retraced configuration from the inner cut-off to the full
    // cut-offs (the old molecule excludes itself from the background through its molecule id).
    std::optional<RunningEnergy> correctionOld =
        CBMC::computeDualCutOffCorrection(growContext, oldComponentData, data.oldMoleculeCopy);
    if (!correctionOld.has_value())
    {
      return false;
    }

    data.retraceData.energies += correctionOld.value();
    data.retraceData.RosenbluthWeight *= std::exp(-system.beta * correctionOld->potentialEnergy());
  }

  time_begin = std::chrono::steady_clock::now();
  data.energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
      system.simulationBox, newMolecule, data.oldMoleculeAtoms, system.netCharge);
  time_end = std::chrono::steady_clock::now();
  oldComponentData.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::Ewald] += (time_end - time_begin);

  time_begin = std::chrono::steady_clock::now();
  data.tailEnergyDifference =
      Interactions::computeInterMolecularTailEnergyDifferenceAddRemove(
          system.forceField, system.simulationBox, system.totalNumberOfPseudoAtoms,
          system.components[newComponent], system.components[oldComponent]) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.spanOfFrameworkAtoms(), newMolecule,
                                                                 data.oldMoleculeAtoms);
  time_end = std::chrono::steady_clock::now();
  oldComponentData.mc_moves_cputime[move][Move::Timing::Tail] += (time_end - time_begin);
  system.mc_moves_cputime[move][Move::Timing::Tail] += (time_end - time_begin);

  if (system.forceField.computePolarization)
  {
    Interactions::computeFrameworkMoleculeElectricFieldDifference(system.forceField, system.simulationBox,
                                                                  system.spanOfFrameworkAtoms(), data.newElectricField,
                                                                  data.oldElectricField, data.growData.atoms,
                                                                  data.oldMoleculeCopy);

    Interactions::computeEwaldFourierElectricFieldDifference(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
        system.trialEik, system.forceField, system.simulationBox, data.newElectricField, data.oldElectricField,
        data.growData.atoms, data.oldMoleculeCopy);

    data.polarizationDifference = Interactions::computePolarizationEnergyDifference(
        system.forceField, data.newElectricField, data.oldElectricField, data.growData.atoms, data.oldMoleculeCopy);
  }

  data.correctionFactorEwald = std::exp(
      -system.beta * (data.energyFourierDifference.potentialEnergy() + data.polarizationDifference.potentialEnergy()));
  data.rosenbluthNew =
      data.growData.RosenbluthWeight * std::exp(-system.beta * data.tailEnergyDifference.potentialEnergy());
  data.rosenbluthOld = data.retraceData.RosenbluthWeight;

  return true;
}

void acceptBoxIdentityChange(System& system, BoxIdentityChangeData& data)
{
  Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);

  std::vector<Atom> acceptedAtoms(data.growData.atoms.begin(), data.growData.atoms.end());
  for (Atom& atom : acceptedAtoms)
  {
    atom.componentId = static_cast<std::uint8_t>(data.newComponent);
  }

  Molecule acceptedMolecule = data.growData.molecule;
  acceptedMolecule.componentId = data.newComponent;

  system.deleteMolecule(data.oldComponent, data.selectedMoleculeOld, data.oldMoleculeAtoms);
  if (system.forceField.computePolarization)
  {
    system.insertMoleculePolarization(data.newComponent, acceptedMolecule, acceptedAtoms, data.newElectricField);
  }
  else
  {
    system.insertMolecule(data.newComponent, acceptedMolecule, acceptedAtoms);
  }
}

RunningEnergy boxEnergyDifference(const BoxIdentityChangeData& data)
{
  return (data.growData.energies - data.retraceData.energies) + data.energyFourierDifference +
         data.tailEnergyDifference + data.polarizationDifference;
}

}  // namespace

std::optional<std::pair<RunningEnergy, RunningEnergy>> MC_Moves::GibbsIdentityChangeMove_CBMC(
    RandomNumber& random, System& systemI, System& systemII, std::size_t componentA)
{
  const Move::Types move = Move::Types::GibbsIdentityChangeCBMC;
  Component& startComponent = systemI.components[componentA];

  if (startComponent.gibbsIdentityChanges.empty())
  {
    return std::nullopt;
  }

  const std::size_t componentB =
      startComponent.gibbsIdentityChanges[random.uniform_integer(0, startComponent.gibbsIdentityChanges.size() - 1)];

  if (componentB >= systemI.components.size() || componentB >= systemII.components.size())
  {
    return std::nullopt;
  }

  if (componentB == componentA)
  {
    return std::nullopt;
  }

  if (systemI.components[componentB].type != startComponent.type)
  {
    return std::nullopt;
  }

  System& boxI = random.uniform() < 0.5 ? systemI : systemII;
  System& boxII = (&boxI == &systemI) ? systemII : systemI;

  if (boxI.numberOfIntegerMoleculesPerComponent[componentA] == 0)
  {
    return std::nullopt;
  }
  if (boxII.numberOfIntegerMoleculesPerComponent[componentB] == 0)
  {
    return std::nullopt;
  }

  boxI.components[componentA].mc_moves_statistics.addTrial(move);
  boxII.components[componentB].mc_moves_statistics.addTrial(move);

  const double numberOfMoleculesA_boxI =
      static_cast<double>(boxI.numberOfIntegerMoleculesPerComponent[componentA]);
  const double numberOfMoleculesB_boxII =
      static_cast<double>(boxII.numberOfIntegerMoleculesPerComponent[componentB]);
  const double numberOfMoleculesA_boxII =
      static_cast<double>(boxII.numberOfIntegerMoleculesPerComponent[componentA]);
  const double numberOfMoleculesB_boxI =
      static_cast<double>(boxI.numberOfIntegerMoleculesPerComponent[componentB]);

  BoxIdentityChangeData boxIData{};
  if (!performBoxIdentityChange(random, boxI, move, componentA, componentB, boxIData))
  {
    return std::nullopt;
  }

  BoxIdentityChangeData boxIIData{};
  if (!performBoxIdentityChange(random, boxII, move, componentB, componentA, boxIIData))
  {
    return std::nullopt;
  }

  const double acceptanceProbability =
      boxIData.correctionFactorEwald * boxIIData.correctionFactorEwald * boxIData.rosenbluthNew *
      boxIIData.rosenbluthNew / (boxIData.rosenbluthOld * boxIIData.rosenbluthOld) * numberOfMoleculesA_boxI *
      numberOfMoleculesB_boxII / ((numberOfMoleculesA_boxII + 1.0) * (numberOfMoleculesB_boxI + 1.0));

  if (random.uniform() < acceptanceProbability)
  {
    boxI.components[componentA].mc_moves_statistics.addAccepted(move);
    boxII.components[componentB].mc_moves_statistics.addAccepted(move);

    acceptBoxIdentityChange(boxI, boxIData);
    acceptBoxIdentityChange(boxII, boxIIData);

    const RunningEnergy energySystemI =
        (&boxI == &systemI) ? boxEnergyDifference(boxIData) : boxEnergyDifference(boxIIData);
    const RunningEnergy energySystemII =
        (&boxI == &systemI) ? boxEnergyDifference(boxIIData) : boxEnergyDifference(boxIData);
    return std::make_pair(energySystemI, energySystemII);
  }

  return std::nullopt;
}
