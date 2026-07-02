module;

module mc_moves_reaction_common;

import std;

import component;
import atom;
import molecule;
import double3;
import simulationbox;
import cbmc;
import molecule;
import cbmc_chain_data;
import randomnumbers;
import system;
import forcefield;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import property_lambda_probability_histogram;
import interactions_external_field;
import reaction;
import mc_moves_move_types;

namespace MC_Moves::ReactionCommon
{

[[nodiscard]] std::optional<RunningEnergy> computeScaledMoleculeEnergyDifference(
    System& system, std::span<const Atom> newAtoms, std::span<const Atom> oldAtoms,
    bool includeTailCorrections) noexcept;

[[nodiscard]] std::vector<std::vector<Atom>> splitAtomsByMoleculeId(std::span<const Atom> atoms)
{
  std::unordered_map<std::uint32_t, std::vector<Atom>> groups;
  groups.reserve(atoms.size());
  for (const Atom& atom : atoms)
  {
    groups[atom.moleculeId].push_back(atom);
  }

  std::vector<std::vector<Atom>> molecules;
  molecules.reserve(groups.size());
  for (auto& [moleculeId, moleculeAtoms] : groups)
  {
    (void)moleculeId;
    molecules.push_back(std::move(moleculeAtoms));
  }
  return molecules;
}

[[nodiscard]] RunningEnergy computePerMoleculeEwaldExclusionDifference(
    System& system, std::span<const std::vector<Atom>> newMolecules,
    std::span<const std::vector<Atom>> oldMolecules) noexcept
{
  RunningEnergy exclusion;
  const std::vector<std::pair<std::complex<double>, std::complex<double>>> totalEikSnapshot = system.totalEik;

  for (const std::vector<Atom>& oldMolecule : oldMolecules)
  {
    const RunningEnergy ewald = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
        system.simulationBox, std::span<const Atom>{}, oldMolecule);
    exclusion.ewald_exclusion += ewald.ewald_exclusion;
    system.totalEik = totalEikSnapshot;
  }
  for (const std::vector<Atom>& newMolecule : newMolecules)
  {
    const RunningEnergy ewald = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
        system.simulationBox, newMolecule, std::span<const Atom>{});
    exclusion.ewald_exclusion += ewald.ewald_exclusion;
    system.totalEik = totalEikSnapshot;
  }
  return exclusion;
}

void updateTotalEikByMoleculeReplacements(System& system, std::span<const Atom> newAtoms,
                                          std::span<const Atom> oldAtoms) noexcept
{
  const std::vector<std::pair<std::complex<double>, std::complex<double>>> storedEikSnapshot = system.storedEik;
  std::vector<std::pair<std::complex<double>, std::complex<double>>> workingStoredEik = storedEikSnapshot;

  for (const std::vector<Atom>& oldMolecule : splitAtomsByMoleculeId(oldAtoms))
  {
    (void)Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, workingStoredEik, system.totalEik, system.forceField,
        system.simulationBox, std::span<const Atom>{}, oldMolecule);
    workingStoredEik = system.totalEik;
  }
  for (const std::vector<Atom>& newMolecule : splitAtomsByMoleculeId(newAtoms))
  {
    (void)Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, workingStoredEik, system.totalEik, system.forceField,
        system.simulationBox, newMolecule, std::span<const Atom>{});
    workingStoredEik = system.totalEik;
  }

  system.storedEik = storedEikSnapshot;
}

void acceptChargedEwaldMove(System& system) noexcept
{
  if (!system.forceField.useCharge || system.forceField.omitEwaldFourier)
  {
    return;
  }

  (void)Interactions::computeEwaldFourierEnergy(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.fixedFrameworkStoredEik, system.storedEik,
      system.forceField, system.simulationBox, system.components, system.numberOfMoleculesPerComponent,
      system.spanOfMoleculeAtoms(), system.netChargeFramework);
}

[[nodiscard]] RunningEnergy computeChargedGroupEwaldDifference(System& system, std::span<const Atom> newAtoms,
                                                                 std::span<const Atom> oldAtoms) noexcept
{
  RunningEnergy ewaldCombined = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, newAtoms, oldAtoms, system.netCharge);

  const RunningEnergy correctedExclusion = computePerMoleculeEwaldExclusionDifference(
      system, splitAtomsByMoleculeId(newAtoms), splitAtomsByMoleculeId(oldAtoms));
  ewaldCombined.ewald_exclusion = correctedExclusion.ewald_exclusion;

  updateTotalEikByMoleculeReplacements(system, newAtoms, oldAtoms);

  return ewaldCombined;
}

[[nodiscard]] RunningEnergy computeChargedWholeMoleculeEwaldDifference(
    System& system, std::span<const Atom> newAtoms, std::span<const Atom> oldAtoms,
    std::span<const Atom> newIntegerAtoms, std::span<const Atom> oldIntegerAtoms,
    std::span<const Atom> newFractionalAtoms, std::span<const Atom> oldFractionalAtoms) noexcept
{
  (void)newIntegerAtoms;
  (void)oldIntegerAtoms;
  (void)newFractionalAtoms;
  (void)oldFractionalAtoms;
  return computeChargedGroupEwaldDifference(system, newAtoms, oldAtoms);
}

[[nodiscard]] std::optional<RunningEnergy> computeGroupSwapEnergyDifference(
    System& system, std::span<const Atom> newAtoms, std::span<const Atom> oldAtoms, bool includeTailCorrections,
    bool includeEwaldCorrections) noexcept;

void applySerialFractionalScaling(std::span<Atom> atoms, double lambda) noexcept;

[[nodiscard]] RunningEnergy internalEnergyFromGrowMolecules(System& system,
                                                            std::span<const ChainGrowData> molecules) noexcept;

void appendAllReactionFractionalMoleculeExclusions(
    const System& system, std::vector<std::pair<std::size_t, std::size_t>>& exclude) noexcept
{
  for (const Reaction& reaction : system.reactions.list)
  {
    if (reaction.reactantFractionalMoleculeIds.size() != system.components.size())
    {
      continue;
    }
    for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
    {
      for (const std::size_t moleculeId : reaction.reactantFractionalMoleculeIds[componentId])
      {
        exclude.emplace_back(componentId, moleculeId);
      }
      for (const std::size_t moleculeId : reaction.productFractionalMoleculeIds[componentId])
      {
        exclude.emplace_back(componentId, moleculeId);
      }
    }
  }
}

void applyLinearReactionScaling(std::span<Atom> atoms, bool isReactant, double lambda) noexcept
{
  const double vdwScale = isReactant ? (1.0 - lambda) : lambda;
  const double coulombScale = isReactant ? std::pow(1.0 - lambda, 5) : std::pow(lambda, 5);
  for (Atom& atom : atoms)
  {
    atom.scalingVDW = vdwScale;
    atom.scalingCoulomb = coulombScale;
    atom.isFractional = true;
  }
}

[[nodiscard]] std::size_t totalStoichiometry(std::span<const std::size_t> stoichiometry) noexcept
{
  std::size_t total{0};
  for (const std::size_t count : stoichiometry)
  {
    total += count;
  }
  return total;
}

[[nodiscard]] bool selectRandomIntegerMolecules(RandomNumber& random, System& system,
                                                std::span<const std::size_t> stoichiometry,
                                                std::vector<std::pair<std::size_t, std::size_t>>& selected) noexcept
{
  selected.clear();
  for (std::size_t componentId = 0; componentId < stoichiometry.size(); ++componentId)
  {
    const std::size_t count = stoichiometry[componentId];
    if (count == 0)
    {
      continue;
    }

    if (count > system.numberOfIntegerMoleculesPerComponent[componentId])
    {
      return false;
    }

    std::vector<std::size_t> available;
    const std::size_t fractionalOffset = system.numberOfFractionalMoleculesPerComponent[componentId];
    for (std::size_t moleculeId = 0; moleculeId < system.numberOfIntegerMoleculesPerComponent[componentId]; ++moleculeId)
    {
      available.push_back(fractionalOffset + moleculeId);
    }

    for (std::size_t k = 0; k < count; ++k)
    {
      const std::size_t pick = random.uniform_integer(0, available.size() - 1);
      selected.emplace_back(componentId, available[pick]);
      available.erase(available.begin() + static_cast<std::ptrdiff_t>(pick));
    }
  }
  return true;
}

[[nodiscard]] std::optional<MoleculeGroupGrowData> growMoleculeGroupInsertion(
    RandomNumber& random, System& system, std::span<const std::size_t> stoichiometry,
    std::span<const std::pair<std::size_t, std::size_t>> excludeMolecules,
    std::optional<double> serialFractionalLambda) noexcept
{
  MoleculeGroupGrowData result;

  std::unordered_set<std::size_t> excludedGlobalMoleculeIds;
  excludedGlobalMoleculeIds.reserve(excludeMolecules.size());
  for (const auto& [componentId, moleculeId] : excludeMolecules)
  {
    excludedGlobalMoleculeIds.insert(system.moleculeIndexOfComponent(componentId, moleculeId));
  }

  std::vector<Atom> background;
  background.reserve(system.spanOfMoleculeAtoms().size());
  for (const Atom& atom : system.spanOfMoleculeAtoms())
  {
    if (!excludedGlobalMoleculeIds.contains(static_cast<std::size_t>(atom.moleculeId)))
    {
      background.push_back(atom);
    }
  }

  const double cutOffFrameworkVDW = system.forceField.cutOffFrameworkVDW;
  const double cutOffMoleculeVDW = system.forceField.cutOffMoleculeVDW;
  const double cutOffCoulomb = system.forceField.cutOffCoulomb;
  const bool useChargedEwald = system.forceField.useCharge;

  for (std::size_t componentId = 0; componentId < stoichiometry.size(); ++componentId)
  {
    for (std::size_t n = 0; n < stoichiometry[componentId]; ++n)
    {
      Component& component = system.components[componentId];
      const std::size_t selectedMolecule = system.numberOfMolecules() + result.molecules.size();

      std::optional<ChainGrowData> growData = CBMC::growMoleculeSwapInsertion(
          random, component, componentId, system.hasExternalField, system.forceField, system.simulationBox,
          system.interpolationGrids, system.externalFieldInterpolationGrid, system.framework,
          system.spanOfFrameworkAtoms(), background, system.beta, component.growType, cutOffFrameworkVDW,
          cutOffMoleculeVDW, cutOffCoulomb, selectedMolecule, 1.0, false, false);

      if (!growData)
      {
        return std::nullopt;
      }

      result.RosenbluthWeight *= growData->RosenbluthWeight;
      if (serialFractionalLambda)
      {
        applySerialFractionalScaling(growData->atom, serialFractionalLambda.value());
      }
      else
      {
        result.energies += growData->energies;
      }

      for (const Atom& atom : growData->atom)
      {
        background.push_back(atom);
      }
      result.molecules.push_back(std::move(*growData));
    }
  }

  if (serialFractionalLambda)
  {
    std::vector<Atom> scaledNewAtoms;
    scaledNewAtoms.reserve(result.molecules.size() * system.components.front().atoms.size());
    for (const ChainGrowData& data : result.molecules)
    {
      scaledNewAtoms.insert(scaledNewAtoms.end(), data.atom.begin(), data.atom.end());
    }

    std::optional<RunningEnergy> insertDifference =
        computeGroupSwapEnergyDifference(system, scaledNewAtoms, {}, true, !useChargedEwald);
    if (!insertDifference)
    {
      return std::nullopt;
    }
    result.energies = insertDifference.value() + internalEnergyFromGrowMolecules(system, result.molecules);
  }

  return result;
}

[[nodiscard]] std::optional<MoleculeGroupRetraceData> retraceMoleculeGroupDeletion(
    RandomNumber& random, System& system,
    std::span<const std::pair<std::size_t, std::size_t>> selectedMolecules) noexcept
{
  MoleculeGroupRetraceData result;

  std::unordered_set<std::size_t> selectedGlobalMoleculeIds;
  selectedGlobalMoleculeIds.reserve(selectedMolecules.size());
  for (const auto& [componentId, moleculeId] : selectedMolecules)
  {
    selectedGlobalMoleculeIds.insert(system.moleculeIndexOfComponent(componentId, moleculeId));
  }

  const double cutOffFrameworkVDW = system.forceField.cutOffFrameworkVDW;
  const double cutOffMoleculeVDW = system.forceField.cutOffMoleculeVDW;
  const double cutOffCoulomb = system.forceField.cutOffCoulomb;

  for (const auto& [componentId, moleculeId] : selectedMolecules)
  {
    const std::size_t currentGlobalMoleculeId = system.moleculeIndexOfComponent(componentId, moleculeId);

    std::vector<Atom> background;
    background.reserve(system.spanOfMoleculeAtoms().size());
    for (const Atom& atom : system.spanOfMoleculeAtoms())
    {
      const std::size_t globalMoleculeId = static_cast<std::size_t>(atom.moleculeId);
      if (selectedGlobalMoleculeIds.contains(globalMoleculeId) && globalMoleculeId != currentGlobalMoleculeId)
      {
        continue;
      }
      background.push_back(atom);
    }

    const Component& component = system.components[componentId];
    std::span<Atom> moleculeAtoms = system.spanOfMolecule(componentId, moleculeId);

    ChainRetraceData retraceData;
    try
    {
      retraceData = CBMC::retraceMoleculeSwapDeletion(
          random, component, system.hasExternalField, system.forceField, system.simulationBox,
          system.interpolationGrids, system.externalFieldInterpolationGrid, system.framework,
          system.spanOfFrameworkAtoms(), background, system.beta, component.growType, cutOffFrameworkVDW,
          cutOffMoleculeVDW, cutOffCoulomb, moleculeAtoms);
    }
    catch (const std::runtime_error&)
    {
      return std::nullopt;
    }

    result.RosenbluthWeight *= retraceData.RosenbluthWeight;
    result.energies += retraceData.energies;
    result.molecules.push_back(std::move(retraceData));
  }

  if (selectedMolecules.size() > 1)
  {
    std::vector<Atom> selectedAtoms;
    for (const auto& [componentId, moleculeId] : selectedMolecules)
    {
      std::span<Atom> molecule = system.spanOfMolecule(componentId, moleculeId);
      selectedAtoms.insert(selectedAtoms.end(), molecule.begin(), molecule.end());
    }
    result.energies += Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox,
                                                                 selectedAtoms);
  }

  return result;
}

[[nodiscard]] std::vector<Atom> collectReactionFractionalAtoms(System& system, Reaction& reaction) noexcept
{
  std::vector<Atom> atoms;
  for (std::size_t componentId = 0; componentId < reaction.reactantFractionalMoleculeIds.size(); ++componentId)
  {
    for (const std::size_t moleculeId : reaction.reactantFractionalMoleculeIds[componentId])
    {
      std::span<Atom> molecule = system.spanOfMolecule(componentId, moleculeId);
      atoms.insert(atoms.end(), molecule.begin(), molecule.end());
    }
    for (const std::size_t moleculeId : reaction.productFractionalMoleculeIds[componentId])
    {
      std::span<Atom> molecule = system.spanOfMolecule(componentId, moleculeId);
      atoms.insert(atoms.end(), molecule.begin(), molecule.end());
    }
  }
  return atoms;
}

[[nodiscard]] std::vector<Atom> collectReactionFractionalAtomsAtLambda(System& system, Reaction& reaction,
                                                                       double lambda) noexcept
{
  std::vector<Atom> atoms;
  for (std::size_t componentId = 0; componentId < reaction.reactantFractionalMoleculeIds.size(); ++componentId)
  {
    for (const std::size_t moleculeId : reaction.reactantFractionalMoleculeIds[componentId])
    {
      std::span<Atom> molecule = system.spanOfMolecule(componentId, moleculeId);
      std::vector<Atom> trialMolecule(molecule.begin(), molecule.end());
      applyLinearReactionScaling(trialMolecule, true, lambda);
      atoms.insert(atoms.end(), trialMolecule.begin(), trialMolecule.end());
    }
    for (const std::size_t moleculeId : reaction.productFractionalMoleculeIds[componentId])
    {
      std::span<Atom> molecule = system.spanOfMolecule(componentId, moleculeId);
      std::vector<Atom> trialMolecule(molecule.begin(), molecule.end());
      applyLinearReactionScaling(trialMolecule, false, lambda);
      atoms.insert(atoms.end(), trialMolecule.begin(), trialMolecule.end());
    }
  }
  return atoms;
}

void setReactionFractionalScaling(System& system, Reaction& reaction, double lambda) noexcept
{
  for (std::size_t componentId = 0; componentId < reaction.reactantFractionalMoleculeIds.size(); ++componentId)
  {
    for (const std::size_t moleculeId : reaction.reactantFractionalMoleculeIds[componentId])
    {
      applyLinearReactionScaling(system.spanOfMolecule(componentId, moleculeId), true, lambda);
    }
    for (const std::size_t moleculeId : reaction.productFractionalMoleculeIds[componentId])
    {
      applyLinearReactionScaling(system.spanOfMolecule(componentId, moleculeId), false, lambda);
    }
  }
}

[[nodiscard]] RunningEnergy computeIntraEnergyDifference(const MoleculeGroupGrowData& growData,
                                                         const MoleculeGroupRetraceData& retraceData) noexcept
{
  RunningEnergy intraDifference;
  for (const ChainGrowData& data : growData.molecules)
  {
    intraDifference.intraVDW += data.energies.intraVDW;
    intraDifference.intraCoul += data.energies.intraCoul;
  }
  for (const ChainRetraceData& data : retraceData.molecules)
  {
    intraDifference.intraVDW -= data.energies.intraVDW;
    intraDifference.intraCoul -= data.energies.intraCoul;
  }
  return intraDifference;
}

[[nodiscard]] RunningEnergy internalEnergyFromGrowMolecules(System& system,
                                                            std::span<const ChainGrowData> molecules) noexcept
{
  RunningEnergy internal;
  for (const ChainGrowData& data : molecules)
  {
    if (data.atom.empty())
    {
      continue;
    }
    const std::size_t componentId = static_cast<std::size_t>(data.atom.front().componentId);
    internal += system.components[componentId].intraMolecularPotentials.computeInternalEnergies(data.atom);
  }
  return internal;
}

[[nodiscard]] RunningEnergy internalEnergyFromScaledGrowMolecules(
    System& system, std::span<const ChainGrowData> molecules, double serialLambda) noexcept
{
  RunningEnergy internal;
  for (const ChainGrowData& data : molecules)
  {
    if (data.atom.empty())
    {
      continue;
    }
    std::vector<Atom> scaledAtoms(data.atom.begin(), data.atom.end());
    applySerialFractionalScaling(scaledAtoms, serialLambda);
    const std::size_t componentId = static_cast<std::size_t>(data.atom.front().componentId);
    internal += system.components[componentId].intraMolecularPotentials.computeInternalEnergies(scaledAtoms);
  }
  return internal;
}

[[nodiscard]] RunningEnergy internalEnergyFromSelectedMolecules(
    System& system, std::span<const std::pair<std::size_t, std::size_t>> selectedMolecules) noexcept
{
  RunningEnergy internal;
  for (const auto& [componentId, moleculeId] : selectedMolecules)
  {
    internal += system.components[componentId].intraMolecularPotentials.computeInternalEnergies(
        system.spanOfMolecule(componentId, moleculeId));
  }
  return internal;
}

[[nodiscard]] RunningEnergy internalEnergyFromActiveSerialFractional(System& system, const Reaction& reaction) noexcept
{
  RunningEnergy internal;
  const std::vector<std::vector<std::size_t>>& activeIds =
      reaction.fractionalSideIsReactants ? reaction.reactantFractionalMoleculeIds
                                         : reaction.productFractionalMoleculeIds;
  for (std::size_t componentId = 0; componentId < activeIds.size(); ++componentId)
  {
    for (const std::size_t moleculeId : activeIds[componentId])
    {
      internal += system.components[componentId].intraMolecularPotentials.computeInternalEnergies(
          system.spanOfMolecule(componentId, moleculeId));
    }
  }
  return internal;
}

[[nodiscard]] std::optional<RunningEnergy> computeReactionFractionalScalingEnergyDifference(
    System& system, Reaction& reaction, double lambdaOld, double lambdaNew, bool includeEwaldCorrections) noexcept
{
  (void)lambdaOld;
  std::vector<Atom> oldAtoms = collectReactionFractionalAtoms(system, reaction);
  if (oldAtoms.empty())
  {
    return RunningEnergy{};
  }

  std::vector<Atom> newAtoms;
  newAtoms.reserve(oldAtoms.size());
  for (std::size_t componentId = 0; componentId < reaction.reactantFractionalMoleculeIds.size(); ++componentId)
  {
    for (const std::size_t moleculeId : reaction.reactantFractionalMoleculeIds[componentId])
    {
      std::span<Atom> molecule = system.spanOfMolecule(componentId, moleculeId);
      std::vector<Atom> trialMolecule(molecule.begin(), molecule.end());
      applyLinearReactionScaling(trialMolecule, true, lambdaNew);
      newAtoms.insert(newAtoms.end(), trialMolecule.begin(), trialMolecule.end());
    }
    for (const std::size_t moleculeId : reaction.productFractionalMoleculeIds[componentId])
    {
      std::span<Atom> molecule = system.spanOfMolecule(componentId, moleculeId);
      std::vector<Atom> trialMolecule(molecule.begin(), molecule.end());
      applyLinearReactionScaling(trialMolecule, false, lambdaNew);
      newAtoms.insert(newAtoms.end(), trialMolecule.begin(), trialMolecule.end());
    }
  }

  return computeGroupSwapEnergyDifference(system, newAtoms, oldAtoms, true, includeEwaldCorrections);
}

[[nodiscard]] std::size_t numberOfReactionMoleculesForComponent(const Reaction& reaction, std::size_t componentId,
                                                                ReactionMoveKind moveKind) noexcept
{
  switch (moveKind)
  {
    case ReactionMoveKind::ForwardInsert:
      return reaction.reactantStoichiometry[componentId];
    case ReactionMoveKind::BackwardDelete:
      return reaction.productStoichiometry[componentId];
    case ReactionMoveKind::LambdaChange:
    default:
      return 0;
  }
}

[[nodiscard]] double computeReactionEquilibriumLogTerm(const System& system, const Reaction& reaction,
                                                       ReactionMoveKind moveKind) noexcept
{
  double term = 0.0;
  const double logVolume = std::log(system.simulationBox.volume);

  for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
  {
    const Component& component = system.components[componentId];
    const std::size_t reactionMolecules =
        numberOfReactionMoleculesForComponent(reaction, componentId, moveKind);
    const std::size_t N =
        system.numberOfIntegerMoleculesPerComponent[componentId] - reactionMolecules;

    if (moveKind == ReactionMoveKind::ForwardInsert)
    {
      if (reaction.reactantStoichiometry[componentId] > 0)
      {
        for (std::size_t j = 0; j < reaction.reactantStoichiometry[componentId]; ++j)
        {
          term += std::log(static_cast<double>(N - j));
        }
        term -= static_cast<double>(reaction.reactantStoichiometry[componentId]) *
                (component.lnPartitionFunction + logVolume);
      }
      if (reaction.productStoichiometry[componentId] > 0)
      {
        for (std::size_t j = 0; j < reaction.productStoichiometry[componentId]; ++j)
        {
          term -= std::log(static_cast<double>(N + j + 1));
        }
        term += static_cast<double>(reaction.productStoichiometry[componentId]) *
                (component.lnPartitionFunction + logVolume);
      }
    }
    else if (moveKind == ReactionMoveKind::BackwardDelete)
    {
      if (reaction.productStoichiometry[componentId] > 0)
      {
        for (std::size_t j = 0; j < reaction.productStoichiometry[componentId]; ++j)
        {
          term += std::log(static_cast<double>(N - j));
        }
        term -= static_cast<double>(reaction.productStoichiometry[componentId]) *
                (component.lnPartitionFunction + logVolume);
      }
      if (reaction.reactantStoichiometry[componentId] > 0)
      {
        for (std::size_t j = 0; j < reaction.reactantStoichiometry[componentId]; ++j)
        {
          term -= std::log(static_cast<double>(N + j + 1));
        }
        term += static_cast<double>(reaction.reactantStoichiometry[componentId]) *
                (component.lnPartitionFunction + logVolume);
      }
    }
  }

  return term;
}

[[nodiscard]] std::optional<RunningEnergy> computeScaledMoleculeEnergyDifference(
    System& system, std::span<const Atom> newAtoms, std::span<const Atom> oldAtoms,
    bool includeTailCorrections) noexcept
{
  if (oldAtoms.empty())
  {
    return RunningEnergy{};
  }

  RunningEnergy energyDifference;

  if (system.hasExternalField)
  {
    std::optional<RunningEnergy> externalFieldDifference = Interactions::computeExternalFieldEnergyDifference(
        system.hasExternalField, system.forceField, system.simulationBox, system.externalFieldInterpolationGrid,
        newAtoms, oldAtoms);
    if (!externalFieldDifference)
    {
      return std::nullopt;
    }
    energyDifference += externalFieldDifference.value();
  }

  std::optional<RunningEnergy> frameworkDifference = Interactions::computeFrameworkMoleculeEnergyDifference(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
      system.spanOfFrameworkAtoms(), newAtoms, oldAtoms);
  if (!frameworkDifference)
  {
    return std::nullopt;
  }
  energyDifference += frameworkDifference.value();

  std::optional<RunningEnergy> moleculeDifference = Interactions::computeInterMolecularEnergyDifference(
      system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), newAtoms, oldAtoms);
  if (!moleculeDifference)
  {
    return std::nullopt;
  }
  energyDifference += moleculeDifference.value();

  RunningEnergy ewaldDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, newAtoms, oldAtoms, system.netCharge);
  energyDifference += ewaldDifference;

  if (includeTailCorrections)
  {
    energyDifference +=
        Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                                system.spanOfMoleculeAtoms(), newAtoms, oldAtoms) +
        Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.spanOfFrameworkAtoms(), newAtoms, oldAtoms);
  }

  return energyDifference;
}

[[nodiscard]] std::optional<RunningEnergy> computeGroupSwapEnergyDifference(
    System& system, std::span<const Atom> newAtoms, std::span<const Atom> oldAtoms, bool includeTailCorrections,
    bool includeEwaldCorrections) noexcept
{
  RunningEnergy energyDifference;

  std::unordered_set<std::size_t> oldGlobalMoleculeIds;
  oldGlobalMoleculeIds.reserve(oldAtoms.size());
  for (const Atom& atom : oldAtoms)
  {
    oldGlobalMoleculeIds.insert(static_cast<std::size_t>(atom.moleculeId));
  }

  std::vector<Atom> background;
  background.reserve(system.spanOfMoleculeAtoms().size());
  for (const Atom& atom : system.spanOfMoleculeAtoms())
  {
    if (!oldGlobalMoleculeIds.contains(static_cast<std::size_t>(atom.moleculeId)))
    {
      background.push_back(atom);
    }
  }

  if (system.hasExternalField)
  {
    std::optional<RunningEnergy> externalFieldDifference = Interactions::computeExternalFieldEnergyDifference(
        system.hasExternalField, system.forceField, system.simulationBox, system.externalFieldInterpolationGrid,
        newAtoms, oldAtoms);
    if (!externalFieldDifference)
    {
      return std::nullopt;
    }
    energyDifference += externalFieldDifference.value();
  }

  std::optional<RunningEnergy> frameworkDifference = Interactions::computeFrameworkMoleculeEnergyDifference(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
      system.spanOfFrameworkAtoms(), newAtoms, oldAtoms);
  if (!frameworkDifference)
  {
    return std::nullopt;
  }
  energyDifference += frameworkDifference.value();

  std::optional<RunningEnergy> moleculeDifference = Interactions::computeInterMolecularEnergyDifference(
      system.forceField, system.simulationBox, background, newAtoms, oldAtoms);
  if (!moleculeDifference)
  {
    return std::nullopt;
  }
  energyDifference += moleculeDifference.value();
  energyDifference += Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, newAtoms) -
                      Interactions::computeInterMolecularEnergy(system.forceField, system.simulationBox, oldAtoms);

  if (includeEwaldCorrections)
  {
    RunningEnergy ewaldDifference = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
        system.simulationBox, newAtoms, oldAtoms, system.netCharge);
    energyDifference += ewaldDifference;
  }

  if (includeTailCorrections)
  {
    energyDifference +=
        Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox, background,
                                                                newAtoms, oldAtoms) +
        Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                 system.spanOfFrameworkAtoms(), newAtoms, oldAtoms);
  }

  return energyDifference;
}

[[nodiscard]] RunningEnergy computeGroupSwapTailEnergyDifference(System& system, std::span<const Atom> newAtoms,
                                                                 std::span<const Atom> oldAtoms) noexcept
{
  std::unordered_set<std::size_t> oldGlobalMoleculeIds;
  oldGlobalMoleculeIds.reserve(oldAtoms.size());
  for (const Atom& atom : oldAtoms)
  {
    oldGlobalMoleculeIds.insert(static_cast<std::size_t>(atom.moleculeId));
  }

  std::vector<Atom> background;
  background.reserve(system.spanOfMoleculeAtoms().size());
  for (const Atom& atom : system.spanOfMoleculeAtoms())
  {
    if (!oldGlobalMoleculeIds.contains(static_cast<std::size_t>(atom.moleculeId)))
    {
      background.push_back(atom);
    }
  }

  return Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox, background,
                                                                 newAtoms, oldAtoms) +
         Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                    system.spanOfFrameworkAtoms(), newAtoms, oldAtoms);
}

void swapMoleculesInComponent(System& system, std::size_t componentId, std::size_t moleculeIndexA,
                              std::size_t moleculeIndexB) noexcept
{
  if (moleculeIndexA == moleculeIndexB)
  {
    return;
  }

  std::span<Atom> moleculeA = system.spanOfMolecule(componentId, moleculeIndexA);
  std::span<Atom> moleculeB = system.spanOfMolecule(componentId, moleculeIndexB);
  std::swap_ranges(moleculeA.begin(), moleculeA.end(), moleculeB.begin());
  std::swap(system.moleculeData[system.moleculeIndexOfComponent(componentId, moleculeIndexA)],
            system.moleculeData[system.moleculeIndexOfComponent(componentId, moleculeIndexB)]);
}

void deleteSelectedMolecules(System& system,
                             std::span<const std::pair<std::size_t, std::size_t>> selectedMolecules) noexcept
{
  std::vector<std::pair<std::size_t, std::size_t>> sorted(selectedMolecules.begin(), selectedMolecules.end());
  std::ranges::sort(sorted, [](const auto& a, const auto& b)
                    {
                      if (a.first != b.first)
                      {
                        return a.first < b.first;
                      }
                      return a.second > b.second;
                    });

  for (const auto& [componentId, moleculeId] : sorted)
  {
    const std::size_t lastIntegerIndex = system.numberOfMoleculesPerComponent[componentId] - 1;
    if (moleculeId != lastIntegerIndex)
    {
      swapMoleculesInComponent(system, componentId, moleculeId, lastIntegerIndex);
    }

    std::span<Atom> molecule = system.spanOfMolecule(componentId, lastIntegerIndex);
    system.deleteMolecule(componentId, lastIntegerIndex, molecule);
  }
}

void insertGrownMolecules(System& system, std::span<const ChainGrowData> growData,
                          std::span<const std::size_t> productStoichiometry) noexcept
{
  std::size_t growIndex = 0;
  for (std::size_t componentId = 0; componentId < productStoichiometry.size(); ++componentId)
  {
    for (std::size_t n = 0; n < productStoichiometry[componentId]; ++n)
    {
      const ChainGrowData& data = growData[growIndex++];
      std::vector<Atom> acceptedAtoms(data.atom.begin(), data.atom.end());
      for (Atom& atom : acceptedAtoms)
      {
        atom.componentId = static_cast<std::uint8_t>(componentId);
      }
      Molecule acceptedMolecule = data.molecule;
      acceptedMolecule.componentId = componentId;
      system.insertMolecule(componentId, acceptedMolecule, acceptedAtoms);
    }
  }
}

void overwriteReactionFractionalMolecules(System& system, Reaction& reaction,
                                         std::span<const ChainGrowData> growData,
                                         std::span<const std::size_t> stoichiometry,
                                         const std::vector<std::vector<std::size_t>>& fractionalMoleculeIds,
                                         bool isReactant, double lambda) noexcept
{
  (void)reaction;
  std::size_t growIndex = 0;
  for (std::size_t componentId = 0; componentId < stoichiometry.size(); ++componentId)
  {
    for (std::size_t k = 0; k < stoichiometry[componentId]; ++k)
    {
      const std::size_t moleculeId = fractionalMoleculeIds[componentId][k];
      std::span<Atom> moleculeAtoms = system.spanOfMolecule(componentId, moleculeId);
      const ChainGrowData& data = growData[growIndex++];

      for (std::size_t atomIndex = 0; atomIndex < moleculeAtoms.size(); ++atomIndex)
      {
        moleculeAtoms[atomIndex].position = data.atom[atomIndex].position;
        moleculeAtoms[atomIndex].charge = data.atom[atomIndex].charge;
      }

      std::vector<Molecule>::iterator moleculeIterator = system.indexForMolecule(componentId, moleculeId);
      moleculeIterator->centerOfMassPosition = data.molecule.centerOfMassPosition;
      moleculeIterator->orientation = data.molecule.orientation;

      applyLinearReactionScaling(moleculeAtoms, isReactant, lambda);
    }
  }
}

void acceptReactionForwardInsert(System& system, Reaction& reaction, double lambdaNew,
                                 std::span<const std::pair<std::size_t, std::size_t>> selectedMolecules,
                                 std::span<const ChainGrowData> growData) noexcept
{
  (void)lambdaNew;
  deleteSelectedMolecules(system, selectedMolecules);
  insertGrownMolecules(system, growData, reaction.productStoichiometry);
}

void acceptReactionBackwardDelete(System& system, Reaction& reaction, double lambdaNew,
                                  std::span<const std::pair<std::size_t, std::size_t>> selectedMolecules,
                                  std::span<const ChainGrowData> growData) noexcept
{
  (void)lambdaNew;
  deleteSelectedMolecules(system, selectedMolecules);
  insertGrownMolecules(system, growData, reaction.reactantStoichiometry);
}

void applySerialFractionalScaling(std::span<Atom> atoms, double lambda) noexcept
{
  const double coulombScale = std::pow(lambda, 5);
  for (Atom& atom : atoms)
  {
    atom.scalingVDW = lambda;
    atom.scalingCoulomb = coulombScale;
    atom.isFractional = true;
  }
}

void setSerialReactionFractionalScaling(System& system, Reaction& reaction, double lambda) noexcept
{
  std::vector<std::vector<std::size_t>>& activeIds =
      reaction.fractionalSideIsReactants ? reaction.reactantFractionalMoleculeIds
                                         : reaction.productFractionalMoleculeIds;
  for (std::size_t componentId = 0; componentId < activeIds.size(); ++componentId)
  {
    for (const std::size_t moleculeId : activeIds[componentId])
    {
      applySerialFractionalScaling(system.spanOfMolecule(componentId, moleculeId), lambda);
    }
  }
}

void shiftStoredFractionalIds(System& system, std::size_t componentId, std::size_t removedMoleculeId) noexcept
{
  for (Reaction& reaction : system.reactions.list)
  {
    if (reaction.reactantFractionalMoleculeIds.size() != system.components.size())
    {
      continue;
    }
    for (std::size_t& moleculeId : reaction.reactantFractionalMoleculeIds[componentId])
    {
      if (moleculeId > removedMoleculeId)
      {
        --moleculeId;
      }
    }
    for (std::size_t& moleculeId : reaction.productFractionalMoleculeIds[componentId])
    {
      if (moleculeId > removedMoleculeId)
      {
        --moleculeId;
      }
    }
  }
}

void deleteReactionSideFractionalMolecules(System& system, Reaction& reaction, bool reactantSide) noexcept
{
  std::vector<std::vector<std::size_t>>& sideIds =
      reactantSide ? reaction.reactantFractionalMoleculeIds : reaction.productFractionalMoleculeIds;

  std::vector<std::pair<std::size_t, std::size_t>> toDelete;
  for (std::size_t componentId = 0; componentId < sideIds.size(); ++componentId)
  {
    for (const std::size_t moleculeId : sideIds[componentId])
    {
      toDelete.emplace_back(componentId, moleculeId);
    }
  }

  std::ranges::sort(toDelete, [](const auto& a, const auto& b)
                    {
                      if (a.first != b.first)
                      {
                        return a.first < b.first;
                      }
                      return a.second > b.second;
                    });

  for (const auto& [componentId, moleculeId] : toDelete)
  {
    std::span<Atom> molecule = system.spanOfMolecule(componentId, moleculeId);
    system.deleteFractionalMolecule(componentId, moleculeId, molecule);
    shiftStoredFractionalIds(system, componentId, moleculeId);
  }

  for (std::size_t componentId = 0; componentId < sideIds.size(); ++componentId)
  {
    sideIds[componentId].clear();
  }
  system.syncReactionFractionalMoleculeIndices();
}

void insertSerialSideFractionalMolecules(System& system, Reaction& reaction,
                                         std::span<const ChainGrowData> growData,
                                         std::span<const std::size_t> stoichiometry, double lambda,
                                         std::vector<std::vector<std::size_t>>& targetIds) noexcept
{
  for (std::size_t componentId = 0; componentId < targetIds.size(); ++componentId)
  {
    targetIds[componentId].clear();
  }

  std::size_t growIndex = 0;
  for (std::size_t componentId = 0; componentId < stoichiometry.size(); ++componentId)
  {
    for (std::size_t n = 0; n < stoichiometry[componentId]; ++n)
    {
      const ChainGrowData& data = growData[growIndex++];
      const std::size_t moleculeIndex =
          system.serialReactionFractionalMoleculeIndex(reaction.id, componentId, n);
      system.insertSerialReactionFractionalMolecule(componentId, moleculeIndex, data.molecule, data.atom, lambda);
      targetIds[componentId].push_back(moleculeIndex);
    }
  }
  system.syncReactionFractionalMoleculeIndices();
}

[[nodiscard]] std::vector<Atom> collectActiveSerialFractionalAtoms(System& system, Reaction& reaction) noexcept
{
  std::vector<Atom> atoms;
  std::vector<std::vector<std::size_t>>& activeIds =
      reaction.fractionalSideIsReactants ? reaction.reactantFractionalMoleculeIds
                                         : reaction.productFractionalMoleculeIds;
  for (std::size_t componentId = 0; componentId < activeIds.size(); ++componentId)
  {
    for (const std::size_t moleculeId : activeIds[componentId])
    {
      std::span<Atom> molecule = system.spanOfMolecule(componentId, moleculeId);
      atoms.insert(atoms.end(), molecule.begin(), molecule.end());
    }
  }
  return atoms;
}

[[nodiscard]] std::optional<RunningEnergy> computeSerialFractionalScalingEnergyDifference(
    System& system, Reaction& reaction, double lambdaNew) noexcept
{
  std::vector<Atom> oldAtoms = collectActiveSerialFractionalAtoms(system, reaction);
  if (oldAtoms.empty())
  {
    return RunningEnergy{};
  }

  std::vector<Atom> newAtoms = oldAtoms;
  applySerialFractionalScaling(newAtoms, lambdaNew);

  const bool useChargedEwald = system.forceField.useCharge;
  std::optional<RunningEnergy> swapDifference =
      computeGroupSwapEnergyDifference(system, newAtoms, oldAtoms, true, !useChargedEwald);
  if (!swapDifference)
  {
    return std::nullopt;
  }
  if (!useChargedEwald)
  {
    return swapDifference.value();
  }
  return swapDifference.value() + computeChargedGroupEwaldDifference(system, newAtoms, oldAtoms);
}

[[nodiscard]] double computeSerialFractionalReactionEquilibriumTerm(const System& system, const Reaction& reaction,
                                                                    bool reactantsToProducts) noexcept
{
  double term = 0.0;
  const double logVolume = std::log(system.simulationBox.volume);

  for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
  {
    const Component& component = system.components[componentId];
    if (reactantsToProducts)
    {
      if (reaction.reactantStoichiometry[componentId] > 0)
      {
        term -= static_cast<double>(reaction.reactantStoichiometry[componentId]) *
                (component.lnPartitionFunction + logVolume);
      }
      if (reaction.productStoichiometry[componentId] > 0)
      {
        term += static_cast<double>(reaction.productStoichiometry[componentId]) *
                (component.lnPartitionFunction + logVolume);
      }
    }
    else
    {
      if (reaction.productStoichiometry[componentId] > 0)
      {
        term -= static_cast<double>(reaction.productStoichiometry[componentId]) *
                (component.lnPartitionFunction + logVolume);
      }
      if (reaction.reactantStoichiometry[componentId] > 0)
      {
        term += static_cast<double>(reaction.reactantStoichiometry[componentId]) *
                (component.lnPartitionFunction + logVolume);
      }
    }
  }

  return term;
}

[[nodiscard]] PropertyLambdaProbabilityHistogram& activeSerialLambdaHistogram(Reaction& reaction) noexcept
{
  return reaction.fractionalSideIsReactants ? reaction.lambda : reaction.lambdaProductSide;
}

[[nodiscard]] std::optional<RunningEnergy> serialLambdaChangeMove(RandomNumber& random, System& system,
                                                                   Reaction& reaction, Move::Types move) noexcept
{
  PropertyLambdaProbabilityHistogram& lambdaHistogram = activeSerialLambdaHistogram(reaction);
  const double lambdaOld = reaction.currentLambda;
  const std::size_t oldBin = std::min(
      static_cast<std::size_t>(lambdaHistogram.numberOfSamplePoints * lambdaOld),
      lambdaHistogram.numberOfSamplePoints > 0 ? lambdaHistogram.numberOfSamplePoints - 1 : 0);
  const double biasOld = lambdaHistogram.biasFactor[oldBin];

  const double maximumChange = reaction.fractionalSideIsReactants ? reaction.maximumLambdaChange
                                                                  : reaction.maximumLambdaChangeProducts;
  const double lambdaNew = lambdaOld + (2.0 * random.uniform() - 1.0) * maximumChange;
  if (lambdaNew < 0.0 || lambdaNew > 1.0)
  {
    return std::nullopt;
  }

  const std::size_t newBin = std::min(
      static_cast<std::size_t>(lambdaHistogram.numberOfSamplePoints * lambdaNew),
      lambdaHistogram.numberOfSamplePoints > 0 ? lambdaHistogram.numberOfSamplePoints - 1 : 0);
  const double biasNew = lambdaHistogram.biasFactor[newBin];

  std::optional<RunningEnergy> scalingDifference =
      computeSerialFractionalScalingEnergyDifference(system, reaction, lambdaNew);
  if (!scalingDifference)
  {
    return std::nullopt;
  }

  const double acceptanceProbability =
      std::exp(-system.beta * scalingDifference->potentialEnergy() + biasNew - biasOld);
  if (random.uniform() >= acceptanceProbability)
  {
    return std::nullopt;
  }

  system.mc_moves_statistics.addConstructed(move);
  system.mc_moves_statistics.addAccepted(move);
  setSerialReactionFractionalScaling(system, reaction, lambdaNew);
  reaction.currentLambda = lambdaNew;
  lambdaHistogram.setCurrentBin(newBin);
  acceptChargedEwaldMove(system);

  return scalingDifference.value();
}

[[nodiscard]] std::optional<RunningEnergy> serialFractionalReactionMove(RandomNumber& random, System& system,
                                                                        Reaction& reaction, Move::Types move,
                                                                        bool useCBMC) noexcept
{
  const bool reactantsToProducts = reaction.fractionalSideIsReactants;
  const std::vector<std::size_t>& newStoichiometry =
      reactantsToProducts ? reaction.productStoichiometry : reaction.reactantStoichiometry;

  if (totalStoichiometry(newStoichiometry) == 0)
  {
    return std::nullopt;
  }

  std::vector<Atom> oldAtoms = collectActiveSerialFractionalAtoms(system, reaction);

  std::vector<std::pair<std::size_t, std::size_t>> excludeFractionals;
  appendAllReactionFractionalMoleculeExclusions(system, excludeFractionals);

  std::optional<MoleculeGroupGrowData> growData =
      growMoleculeGroupInsertion(random, system, newStoichiometry, excludeFractionals);
  if (!growData)
  {
    return std::nullopt;
  }

  std::vector<Atom> newAtoms;
  std::uint32_t trialMoleculeId = 1'000'000;
  for (const ChainGrowData& data : growData->molecules)
  {
    std::vector<Atom> scaledAtoms(data.atom.begin(), data.atom.end());
    applySerialFractionalScaling(scaledAtoms, reaction.currentLambda);
    for (Atom& atom : scaledAtoms)
    {
      atom.moleculeId = trialMoleculeId;
    }
    ++trialMoleculeId;
    newAtoms.insert(newAtoms.end(), scaledAtoms.begin(), scaledAtoms.end());
  }

  const bool useChargedEwald = system.forceField.useCharge;
  std::optional<RunningEnergy> swapDifference =
      computeGroupSwapEnergyDifference(system, newAtoms, oldAtoms, true, !useChargedEwald);
  if (!swapDifference)
  {
    return std::nullopt;
  }

  const RunningEnergy internalDifference =
      internalEnergyFromScaledGrowMolecules(system, growData->molecules, reaction.currentLambda) -
      internalEnergyFromActiveSerialFractional(system, reaction);
  RunningEnergy energyDifference = swapDifference.value() + internalDifference;

  if (useChargedEwald)
  {
    energyDifference += computeChargedGroupEwaldDifference(system, newAtoms, oldAtoms);
  }

  const double equilibriumTerm = computeSerialFractionalReactionEquilibriumTerm(system, reaction, reactantsToProducts);
  double rosenbluthRatio = useCBMC ? growData->RosenbluthWeight : 1.0;
  const double acceptanceProbability = rosenbluthRatio * std::exp(equilibriumTerm) *
                                     std::exp(-system.beta * energyDifference.potentialEnergy());
  if (random.uniform() >= acceptanceProbability)
  {
    return std::nullopt;
  }

  system.mc_moves_statistics.addConstructed(move);
  system.mc_moves_statistics.addAccepted(move);

  deleteReactionSideFractionalMolecules(system, reaction, reaction.fractionalSideIsReactants);
  std::vector<std::vector<std::size_t>>& targetIds =
      reactantsToProducts ? reaction.productFractionalMoleculeIds : reaction.reactantFractionalMoleculeIds;
  insertSerialSideFractionalMolecules(system, reaction, growData->molecules, newStoichiometry, reaction.currentLambda,
                                      targetIds);
  reaction.fractionalSideIsReactants = !reaction.fractionalSideIsReactants;
  system.syncReactionLambdaBin(reaction);
  acceptChargedEwaldMove(system);

  return energyDifference;
}

[[nodiscard]] std::optional<RunningEnergy> serialWholeMoleculeReactionMove(RandomNumber& random, System& system,
                                                                           Reaction& reaction, Move::Types move,
                                                                           bool useCBMC) noexcept
{
  const bool forward = reaction.fractionalSideIsReactants;
  const std::vector<std::size_t>& removeStoichiometry =
      forward ? reaction.reactantStoichiometry : reaction.productStoichiometry;
  const std::vector<std::size_t>& insertStoichiometry =
      forward ? reaction.productStoichiometry : reaction.reactantStoichiometry;
  const std::vector<std::size_t>& newFractionalStoichiometry =
      forward ? reaction.productStoichiometry : reaction.reactantStoichiometry;

  if (totalStoichiometry(removeStoichiometry) == 0 || totalStoichiometry(insertStoichiometry) == 0)
  {
    return std::nullopt;
  }

  std::vector<std::pair<std::size_t, std::size_t>> selectedMolecules;
  if (!selectRandomIntegerMolecules(random, system, removeStoichiometry, selectedMolecules))
  {
    return std::nullopt;
  }

  std::vector<std::pair<std::size_t, std::size_t>> excludeMolecules = selectedMolecules;
  appendAllReactionFractionalMoleculeExclusions(system, excludeMolecules);
  std::optional<MoleculeGroupGrowData> growData =
      growMoleculeGroupInsertion(random, system, insertStoichiometry, excludeMolecules);
  if (!growData)
  {
    return std::nullopt;
  }

  std::optional<MoleculeGroupRetraceData> retraceData =
      retraceMoleculeGroupDeletion(random, system, selectedMolecules);
  if (!retraceData)
  {
    return std::nullopt;
  }

  std::vector<std::pair<std::size_t, std::size_t>> excludeForFractionalGrow = selectedMolecules;
  appendAllReactionFractionalMoleculeExclusions(system, excludeForFractionalGrow);

  std::optional<MoleculeGroupGrowData> fractionalGrowData = growMoleculeGroupInsertion(
      random, system, newFractionalStoichiometry, excludeForFractionalGrow, reaction.currentLambda);
  if (!fractionalGrowData)
  {
    return std::nullopt;
  }

  std::vector<Atom> newIntegerAtoms;
  for (const ChainGrowData& data : growData->molecules)
  {
    newIntegerAtoms.insert(newIntegerAtoms.end(), data.atom.begin(), data.atom.end());
  }

  std::vector<Atom> newFractionalAtoms;
  std::uint32_t trialMoleculeId = 1'000'000;
  for (const ChainGrowData& data : fractionalGrowData->molecules)
  {
    std::vector<Atom> scaledAtoms(data.atom.begin(), data.atom.end());
    applySerialFractionalScaling(scaledAtoms, reaction.currentLambda);
    for (Atom& atom : scaledAtoms)
    {
      atom.moleculeId = trialMoleculeId;
    }
    ++trialMoleculeId;
    newFractionalAtoms.insert(newFractionalAtoms.end(), scaledAtoms.begin(), scaledAtoms.end());
  }

  std::vector<Atom> newAtoms = newIntegerAtoms;
  newAtoms.insert(newAtoms.end(), newFractionalAtoms.begin(), newFractionalAtoms.end());

  std::vector<Atom> oldFractionalAtoms = collectActiveSerialFractionalAtoms(system, reaction);
  std::vector<Atom> oldIntegerAtoms;
  for (const auto& [componentId, moleculeId] : selectedMolecules)
  {
    std::span<Atom> molecule = system.spanOfMolecule(componentId, moleculeId);
    oldIntegerAtoms.insert(oldIntegerAtoms.end(), molecule.begin(), molecule.end());
  }

  std::vector<Atom> oldAtoms = oldFractionalAtoms;
  oldAtoms.insert(oldAtoms.end(), oldIntegerAtoms.begin(), oldIntegerAtoms.end());

  const bool splitEwald = system.forceField.useCharge;
  std::optional<RunningEnergy> swapDifference =
      computeGroupSwapEnergyDifference(system, newAtoms, oldAtoms, true, !splitEwald);
  if (!swapDifference)
  {
    return std::nullopt;
  }

  const RunningEnergy internalDifference =
      internalEnergyFromGrowMolecules(system, growData->molecules) -
      internalEnergyFromSelectedMolecules(system, selectedMolecules) +
      internalEnergyFromScaledGrowMolecules(system, fractionalGrowData->molecules, reaction.currentLambda) -
      internalEnergyFromActiveSerialFractional(system, reaction);

  const ReactionMoveKind moveKind =
      forward ? ReactionMoveKind::ForwardInsert : ReactionMoveKind::BackwardDelete;
  const double equilibriumTerm = computeReactionEquilibriumLogTerm(system, reaction, moveKind);

  RunningEnergy energyDifference{};
  if (useCBMC && !splitEwald && system.framework.has_value())
  {
    std::optional<RunningEnergy> removeFractional =
        computeGroupSwapEnergyDifference(system, {}, oldFractionalAtoms, true, true);
    if (!removeFractional)
    {
      return std::nullopt;
    }
    energyDifference = (growData->energies - retraceData->energies) + fractionalGrowData->energies +
                       removeFractional.value() - internalEnergyFromActiveSerialFractional(system, reaction);
  }
  else
  {
    energyDifference = swapDifference.value() + internalDifference;
  }
  RunningEnergy energyFourierInteger{};
  if (splitEwald)
  {
    const std::vector<std::pair<std::complex<double>, std::complex<double>>> totalEikBeforeInteger =
        system.totalEik;
    energyFourierInteger =
        computeChargedGroupEwaldDifference(system, newIntegerAtoms, oldIntegerAtoms);
    system.totalEik = totalEikBeforeInteger;

    energyDifference += computeChargedWholeMoleculeEwaldDifference(system, newAtoms, oldAtoms, newIntegerAtoms,
                                                                   oldIntegerAtoms, newFractionalAtoms,
                                                                   oldFractionalAtoms);
  }

  double rosenbluthRatio = 1.0;
  if (useCBMC)
  {
    rosenbluthRatio = growData->RosenbluthWeight / retraceData->RosenbluthWeight;
    rosenbluthRatio *= fractionalGrowData->RosenbluthWeight;
  }

  double acceptanceProbability = 0.0;
  if (useCBMC)
  {
    if (splitEwald)
    {
      const double correctionFactorEwald = std::exp(-system.beta * energyFourierInteger.potentialEnergy());
      acceptanceProbability = correctionFactorEwald * rosenbluthRatio * std::exp(equilibriumTerm);
    }
    else
    {
      acceptanceProbability = rosenbluthRatio * std::exp(equilibriumTerm);
    }
  }
  else
  {
    acceptanceProbability =
        rosenbluthRatio * std::exp(equilibriumTerm) * std::exp(-system.beta * energyDifference.potentialEnergy());
  }

  if (random.uniform() >= acceptanceProbability)
  {
    return std::nullopt;
  }

  system.mc_moves_statistics.addConstructed(move);
  system.mc_moves_statistics.addAccepted(move);

  if (forward)
  {
    acceptReactionForwardInsert(system, reaction, reaction.currentLambda, selectedMolecules, growData->molecules);
  }
  else
  {
    acceptReactionBackwardDelete(system, reaction, reaction.currentLambda, selectedMolecules, growData->molecules);
  }

  deleteReactionSideFractionalMolecules(system, reaction, reaction.fractionalSideIsReactants);
  std::vector<std::vector<std::size_t>>& targetIds =
      forward ? reaction.productFractionalMoleculeIds : reaction.reactantFractionalMoleculeIds;
  insertSerialSideFractionalMolecules(system, reaction, fractionalGrowData->molecules, newFractionalStoichiometry,
                                      reaction.currentLambda, targetIds);
  reaction.fractionalSideIsReactants = !reaction.fractionalSideIsReactants;
  system.syncReactionLambdaBin(reaction);
  acceptChargedEwaldMove(system);

  return energyDifference;
}

}  // namespace MC_Moves::ReactionCommon
