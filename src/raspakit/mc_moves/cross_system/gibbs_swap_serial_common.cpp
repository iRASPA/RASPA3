module;

module mc_moves_gibbs_swap_serial_common;

import std;

import running_energy;
import atom;
import system;
import component;
import property_lambda_probability_histogram;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import mc_moves_move_types;

namespace MC_Moves::GibbsSwapSerialCommon
{

FractionalSnapshot saveFractionalMolecules(const System& system)
{
  FractionalSnapshot snapshot;
  for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
  {
    if (system.numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC[componentId] == 0)
    {
      continue;
    }
    const std::size_t index = system.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCFCMC, componentId);
    std::span<const Atom> fractionalMolecule = system.spanOfMolecule(componentId, index);
    snapshot.emplace_back(std::make_pair(componentId, index),
                          std::vector<Atom>(fractionalMolecule.begin(), fractionalMolecule.end()));
  }
  return snapshot;
}

void restoreFractionalMolecules(System& system, const FractionalSnapshot& snapshot)
{
  for (const auto& [indices, moleculeAtoms] : snapshot)
  {
    const auto& [componentId, index] = indices;
    std::span<Atom> fractionalMolecule = system.spanOfMolecule(componentId, index);
    std::copy(moleculeAtoms.begin(), moleculeAtoms.end(), fractionalMolecule.begin());
    system.cellListAtomsMoved(fractionalMolecule);
  }
}

RunningEnergy computeTailEnergyDifference(
    System& system, std::vector<double>& tailEffectiveCounts,
    std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups>& tailGroupCounts, std::span<const Atom> newAtoms,
    std::span<const Atom> oldAtoms)
{
  RunningEnergy result =
      Interactions::computeInterMolecularTailEnergyDifferenceAggregated(
          system.forceField, system.simulationBox, tailEffectiveCounts, tailGroupCounts, newAtoms, oldAtoms) +
      Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.spanOfFrameworkAtoms(), newAtoms, oldAtoms);
  Interactions::updateEffectiveTypeCounts(tailEffectiveCounts, tailGroupCounts, newAtoms, oldAtoms);
  return result;
}

namespace
{
// 'chainedOnTrialEik': this is not the first Ewald sub-step of the move in 'system'; the difference then starts
// from the structure factor the preceding sub-steps left in 'trialEik' (otherwise their update would be discarded).
std::optional<RunningEnergy> computeFractionalMoleculesEnergyDifference(
    System& system, std::vector<double>& tailEffectiveCounts,
    std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups>& tailGroupCounts,
    std::span<const Atom> newAtomsCombined, std::span<const Atom> oldAtomsCombined, bool chainedOnTrialEik) noexcept
{
  std::optional<RunningEnergy> frameworkDifference = Interactions::computeFrameworkMoleculeEnergyDifference(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework, system.frameworkCellList(),
      system.spanOfFrameworkAtoms(), newAtomsCombined, oldAtomsCombined);
  if (!frameworkDifference.has_value())
  {
    return std::nullopt;
  }

  std::optional<RunningEnergy> moleculeDifference = Interactions::computeInterMolecularEnergyDifference(
      system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), newAtomsCombined, oldAtomsCombined);
  if (!moleculeDifference.has_value())
  {
    return std::nullopt;
  }

  RunningEnergy ewaldDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, chainedOnTrialEik ? system.trialEik : system.storedEik,
      system.trialEik, system.forceField, system.simulationBox, system.components, newAtomsCombined, oldAtomsCombined,
      system.netCharge);

  RunningEnergy tailDifference =
      computeTailEnergyDifference(system, tailEffectiveCounts, tailGroupCounts, newAtomsCombined, oldAtomsCombined);

  return frameworkDifference.value() + moleculeDifference.value() + ewaldDifference + tailDifference;
}

void setTargetScaling(Atom& atom, const System& system, std::size_t componentId, std::size_t selectedComponent,
                      bool systemBecomesActive, double selectedActiveLambda, const System& lambdaSourceSystem)
{
  if (systemBecomesActive)
  {
    const double lambda = (componentId == selectedComponent)
                              ? selectedActiveLambda
                              : lambdaSourceSystem.components[componentId].lambdaGC.lambdaValue();
    atom.setScalingToFractional(lambda, system.components[componentId].lambdaGC.dUdlambdaGroupId);
  }
  else
  {
    // inactive fractional molecules must not contribute dUdlambda
    atom.setScalingToInactiveFractional();
  }
}
}  // namespace

std::optional<RunningEnergy> computeFlagSwapFractionalEnergyDifference(
    System& system, std::vector<double>& tailEffectiveCounts,
    std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups>& tailGroupCounts,
    const FractionalSnapshot& beforeSnapshot, std::size_t selectedComponent, bool systemBecomesActive,
    double selectedActiveLambda, const System& lambdaSourceSystem) noexcept
{
  std::vector<Atom> oldAtomsCombined;
  for (const auto& [indices, moleculeAtoms] : beforeSnapshot)
  {
    (void)indices;
    oldAtomsCombined.insert(oldAtomsCombined.end(), moleculeAtoms.begin(), moleculeAtoms.end());
  }

  // trial state: all fractional molecules take the scaling of the target activity state
  std::vector<Atom> newAtomsCombined;
  for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
  {
    if (system.numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC[componentId] == 0)
    {
      continue;
    }

    const std::size_t indexFractionalMolecule =
        system.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCFCMC, componentId);
    std::span<const Atom> fractionalMolecule = system.spanOfMolecule(componentId, indexFractionalMolecule);
    for (Atom atom : fractionalMolecule)
    {
      setTargetScaling(atom, system, componentId, selectedComponent, systemBecomesActive, selectedActiveLambda,
                       lambdaSourceSystem);
      newAtomsCombined.push_back(atom);
    }
  }

  return computeFractionalMoleculesEnergyDifference(system, tailEffectiveCounts, tailGroupCounts, newAtomsCombined,
                                                    oldAtomsCombined, false);
}

std::optional<RunningEnergy> computeFlagSwapOtherComponentsEnergyDifference(
    System& system, std::vector<double>& tailEffectiveCounts,
    std::array<std::vector<double>, maximumNumberOfDUDlambdaGroups>& tailGroupCounts,
    const FractionalSnapshot& beforeSnapshot, std::size_t selectedComponent, bool systemBecomesActive,
    const System& lambdaSourceSystem) noexcept
{
  std::vector<Atom> oldAtomsCombined;
  std::vector<Atom> newAtomsCombined;

  for (const auto& [indices, moleculeAtoms] : beforeSnapshot)
  {
    const auto& [componentId, index] = indices;
    if (componentId == selectedComponent)
    {
      continue;
    }
    oldAtomsCombined.insert(oldAtomsCombined.end(), moleculeAtoms.begin(), moleculeAtoms.end());

    std::span<const Atom> fractionalMolecule = system.spanOfMolecule(componentId, index);
    for (Atom atom : fractionalMolecule)
    {
      setTargetScaling(atom, system, componentId, selectedComponent, systemBecomesActive, 0.0, lambdaSourceSystem);
      newAtomsCombined.push_back(atom);
    }
  }

  if (oldAtomsCombined.empty())
  {
    return RunningEnergy{};
  }

  return computeFractionalMoleculesEnergyDifference(system, tailEffectiveCounts, tailGroupCounts, newAtomsCombined,
                                                    oldAtomsCombined, true);
}

void applyFlagSwapScaling(System& system, std::size_t selectedComponent, bool systemBecomesActive,
                          double selectedActiveLambda, const System& lambdaSourceSystem)
{
  for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
  {
    if (system.numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC[componentId] == 0)
    {
      continue;
    }

    const std::size_t indexFractionalMolecule =
        system.indexOfFractionalMoleculeForMove(Move::Types::GibbsSwapCFCMC, componentId);
    std::span<Atom> fractionalMolecule = system.spanOfMolecule(componentId, indexFractionalMolecule);
    for (Atom& atom : fractionalMolecule)
    {
      setTargetScaling(atom, system, componentId, selectedComponent, systemBecomesActive, selectedActiveLambda,
                       lambdaSourceSystem);
    }
    // The cell list caches the position, scaling, and group id of every atom: refresh the edited molecule, or a
    // CBMC grow/retrace later in the move sees the fractional molecules at their old position and scaling.
    system.cellListAtomsMoved(fractionalMolecule);
  }
}

void syncOtherComponentLambdaBins(System& systemBecomingInactive, System& systemBecomingActive,
                                  std::size_t selectedComponent)
{
  for (std::size_t componentId = 0; componentId < systemBecomingInactive.components.size(); ++componentId)
  {
    if (componentId == selectedComponent ||
        systemBecomingInactive.numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC[componentId] == 0)
    {
      continue;
    }

    systemBecomingActive.components[componentId].lambdaGC.setCurrentBin(
        systemBecomingInactive.components[componentId].lambdaGC.currentBin);
  }
}

}  // namespace MC_Moves::GibbsSwapSerialCommon
