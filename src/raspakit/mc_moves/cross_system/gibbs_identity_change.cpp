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
import forcefield;
import simulationbox;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_polarization;
import mc_moves_move_types;
import mc_moves_cputime;

namespace
{

struct BoxIdentityChangeData
{
  CBMC::GrowResult growData;
  CBMC::RetraceResult retraceData;
  RunningEnergy energyFourierDifference;
  RunningEnergy tailEnergyDifference;
  RunningEnergy polarizationDifference;
  double correctionFactorEwald{1.0};
  double logRosenbluthNew{0.0};
  double logRosenbluthOld{0.0};
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

  const std::size_t oldGlobalMoleculeId = system.moleculeIndexOfComponent(oldComponent, data.selectedMoleculeOld);
  const std::size_t trialMoleculeId = system.numberOfMolecules();

  // The new molecule carries a fresh id while the old one is still in the background: exclude it.
  const CBMC::GrowContext growContext = system.makeGrowContext().withSkippedMolecule(oldGlobalMoleculeId);

  std::optional<CBMC::GrowResult> growData =
      MC_Moves::timed(system, oldComponentData, move, Move::Timing::NonEwald,
            [&]
            {
              return CBMC::growNewMolecule(
                  random, growContext, newComponentData, {.componentId = newComponent, .moleculeId = trialMoleculeId},
                  {.firstBead = CBMC::FirstBeadScheme::Pinned, .firstBeadPosition = oldStartingBead.position});
            });

  // Dual cut-off scheme: correct the grown configuration from the inner cut-off to the full
  // cut-offs, using the same background (the old molecule excluded) as the growth.
  if (!growData || !CBMC::applyDualCutOffCorrection(growContext, newComponentData, *growData))
  {
    return false;
  }

  data.growData = std::move(*growData);
  const std::span<const Atom> newMolecule(data.growData.atoms.begin(), data.growData.atoms.end());

  data.oldMoleculeCopy.assign(data.oldMoleculeAtoms.begin(), data.oldMoleculeAtoms.end());
  data.oldElectricField.resize(data.oldMoleculeCopy.size());
  data.newElectricField.resize(newMolecule.size());

  oldComponentData.mc_moves_statistics.addConstructed(move);

  data.retraceData = MC_Moves::timed(system, oldComponentData, move, Move::Timing::NonEwald,
                           [&]
                           {
                             return CBMC::retraceMolecule(random, growContext, oldComponentData, data.oldMoleculeAtoms,
                                                          {.firstBead = CBMC::FirstBeadScheme::Pinned});
                           });

  // Dual cut-off scheme: correct the retraced configuration from the inner cut-off to the full
  // cut-offs (the old molecule excludes itself from the background through its molecule id).
  if (!CBMC::applyDualCutOffCorrection(growContext, oldComponentData, data.oldMoleculeCopy, data.retraceData))
  {
    return false;
  }

  data.energyFourierDifference =
      MC_Moves::timed(system, oldComponentData, move, Move::Timing::Ewald,
            [&]
            {
              return Interactions::energyDifferenceEwaldFourier(
                  system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik,
                  system.forceField, system.simulationBox, newMolecule, data.oldMoleculeAtoms, system.netCharge);
            });

  data.tailEnergyDifference =
      MC_Moves::timed(system, oldComponentData, move, Move::Timing::Tail,
            [&]
            {
              return Interactions::computeInterMolecularTailEnergyDifferenceAddRemove(
                         system.forceField, system.simulationBox, system.totalNumberOfPseudoAtoms,
                         system.components[newComponent], system.components[oldComponent]) +
                     Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                                system.spanOfFrameworkAtoms(),
                                                                                newMolecule, data.oldMoleculeAtoms);
            });

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
  // Rosenbluth weights are carried as exact logarithms: the raw weights of long chains underflow to zero.
  data.logRosenbluthNew =
      data.growData.logRosenbluthWeight - system.beta * data.tailEnergyDifference.potentialEnergy();
  data.logRosenbluthOld = data.retraceData.logRosenbluthWeight;

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
      std::exp(std::log(boxIData.correctionFactorEwald) + std::log(boxIIData.correctionFactorEwald) +
               boxIData.logRosenbluthNew + boxIIData.logRosenbluthNew - boxIData.logRosenbluthOld -
               boxIIData.logRosenbluthOld +
               std::log(numberOfMoleculesA_boxI * numberOfMoleculesB_boxII /
                        ((numberOfMoleculesA_boxII + 1.0) * (numberOfMoleculesB_boxI + 1.0))));

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
