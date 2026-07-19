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
import cbmc_interactions;
import randomnumbers;
import system;
import forcefield;
import units;
import potential_coulomb_real_space;
import running_energy;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import property_lambda_probability_histogram;
import interactions_external_field;
import reaction;
import mc_moves_move_types;

namespace MC_Moves::ReactionCommon
{

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

void acceptChargedEwaldMove(System& system) noexcept
{
  if (!system.forceField.usesEwaldFourier())
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
  const std::vector<std::vector<Atom>> newMolecules = splitAtomsByMoleculeId(newAtoms);
  const std::vector<std::vector<Atom>> oldMolecules = splitAtomsByMoleculeId(oldAtoms);

  RunningEnergy ewaldCombined = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
      system.simulationBox, newAtoms, oldAtoms, system.netCharge);

  // The combined difference treats each atom set as one big molecule and therefore adds spurious
  // cross-molecule exclusion (erf) pairs, to both the energy and the per-group dU/dlambda. Undo
  // the all-pairs exclusion terms and re-add the true intra-molecular ones.
  auto accumulateExclusion = [&system](std::span<const Atom> atoms, double sign, RunningEnergy& accumulator)
  {
    const double alpha = system.forceField.EwaldAlpha;
    for (std::size_t i = 0; i != atoms.size(); ++i)
    {
      for (std::size_t j = i + 1; j != atoms.size(); ++j)
      {
        double3 dr = atoms[i].position - atoms[j].position;
        dr = system.simulationBox.applyPeriodicBoundaryConditions(dr);
        const double r = std::sqrt(double3::dot(dr, dr));
        const double scalingTotal = atoms[i].scalingCoulomb * atoms[j].scalingCoulomb;
        const Potentials::EwaldExclusionFactors exclusion =
            Potentials::ewaldExclusionFactors(alpha, scalingTotal, r);
        const double prefactor = sign * Units::CoulombicConversionFactor * atoms[i].charge * atoms[j].charge;
        accumulator.ewald_exclusion += scalingTotal * prefactor * exclusion.potential;
        accumulator.addDudlambdaEwald(atoms[i].groupId, atoms[j].groupId, atoms[i].scalingCoulomb,
                                      atoms[j].scalingCoulomb, prefactor * exclusion.dUdlambda);
      }
    }
  };
  // remove the all-pairs terms the combined call added (+pairs(old), -pairs(new))
  accumulateExclusion(oldAtoms, -1.0, ewaldCombined);
  accumulateExclusion(newAtoms, +1.0, ewaldCombined);
  // add the intra-molecular exclusion terms per molecule
  for (const std::vector<Atom>& oldMolecule : oldMolecules)
  {
    accumulateExclusion(oldMolecule, +1.0, ewaldCombined);
  }
  for (const std::vector<Atom>& newMolecule : newMolecules)
  {
    accumulateExclusion(newMolecule, -1.0, ewaldCombined);
  }

  // rebuild trialEik by replacing the molecules one at a time, keeping storedEik untouched
  const std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> storedEikSnapshot = system.storedEik;
  std::vector<std::pair<std::complex<double>, std::array<std::complex<double>, 4>>> workingStoredEik = storedEikSnapshot;
  for (const std::vector<Atom>& oldMolecule : oldMolecules)
  {
    (void)Interactions::energyDifferenceEwaldFourier(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                     workingStoredEik, system.trialEik, system.forceField,
                                                     system.simulationBox, std::span<const Atom>{}, oldMolecule);
    workingStoredEik = system.trialEik;
  }
  for (const std::vector<Atom>& newMolecule : newMolecules)
  {
    (void)Interactions::energyDifferenceEwaldFourier(system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
                                                     workingStoredEik, system.trialEik, system.forceField,
                                                     system.simulationBox, newMolecule, std::span<const Atom>{});
    workingStoredEik = system.trialEik;
  }
  system.storedEik = storedEikSnapshot;

  return ewaldCombined;
}

[[nodiscard]] std::optional<RunningEnergy> computeGroupSwapEnergyDifference(
    System& system, std::span<const Atom> newAtoms, std::span<const Atom> oldAtoms, bool includeTailCorrections = true,
    bool includeEwaldCorrections = true, std::span<const Atom> excludeFromBackground = {}) noexcept;

void appendAllReactionFractionalMoleculeExclusions(const System& system,
                                                   std::vector<std::pair<std::size_t, std::size_t>>& exclude) noexcept
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

void applyLinearReactionScaling(std::span<Atom> atoms, bool isReactant, double lambda,
                                std::uint8_t dUdlambdaGroupId) noexcept
{
  // staged schedule from scaling.ixx: VDW switches on for lambda in [0, 0.5], Coulomb in [0.5, 1].
  // Reactants are coupled with (1 - lambda) so that their electrostatics vanish first.
  const double effectiveLambda = isReactant ? (1.0 - lambda) : lambda;
  for (Atom& atom : atoms)
  {
    atom.setScalingToFractional(effectiveLambda, dUdlambdaGroupId);
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
    for (std::size_t moleculeId = 0; moleculeId < system.numberOfIntegerMoleculesPerComponent[componentId];
         ++moleculeId)
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
    std::span<const std::pair<std::size_t, std::size_t>> excludeMolecules, double scaling, bool isFractional,
    std::uint8_t dUdlambdaGroupId, bool useCBMC) noexcept
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

  // Determine cutoff distances based on whether dual cutoff is used (only used for the CBMC growth).
  const double cutOffFrameworkVDW =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffFrameworkVDW;
  const double cutOffMoleculeVDW =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffMoleculeVDW;
  const double cutOffCoulomb =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffCoulomb;

  for (std::size_t componentId = 0; componentId < stoichiometry.size(); ++componentId)
  {
    for (std::size_t n = 0; n < stoichiometry[componentId]; ++n)
    {
      Component& component = system.components[componentId];
      const std::size_t selectedMolecule = system.numberOfMolecules() + result.molecules.size();

      std::optional<ChainGrowData> growData;
      if (useCBMC)
      {
        const CBMC::GrowContext growContext{system.hasExternalField, system.forceField, system.simulationBox,
                                            system.interpolationGrids, system.externalFieldInterpolationGrid,
                                            system.framework, system.spanOfFrameworkAtoms(), background, system.beta,
                                            cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb};
        growData = CBMC::growMoleculeSwapInsertion(random, growContext, component, componentId, selectedMolecule,
                                                   scaling, dUdlambdaGroupId, isFractional);

        if (growData && system.forceField.useDualCutOff)
        {
          // Dual cut-off scheme: correct the grown molecule from the inner cut-off to the full
          // cut-offs, using the same background (previously grown group members included) as the growth.
          std::optional<RunningEnergy> correctionNew =
              CBMC::computeDualCutOffCorrection(growContext, component, growData->atoms);
          if (!correctionNew.has_value())
          {
            return std::nullopt;
          }

          growData->energies += correctionNew.value();
          growData->RosenbluthWeight *= std::exp(-system.beta * correctionNew->potentialEnergy());
        }
      }
      else
      {
        // Non-CBMC path (ReactionCFCMC / ReactionConventionalCFCMC), same as SwapCFCMC:
        // rigid reference geometry or an ideal-gas Boltzmann conformation, then random COM/orientation.
        auto [molecule, atoms] = system.equilibratedIdealGasMoleculeRandomInBox(random, componentId);
        for (Atom& atom : atoms)
        {
          atom.moleculeId = static_cast<std::uint32_t>(selectedMolecule);
          atom.componentId = static_cast<std::uint8_t>(componentId);
          if (isFractional)
          {
            atom.setScalingToFractional(scaling, dUdlambdaGroupId);
          }
          else
          {
            atom.setScalingToInteger();
          }
        }
        growData.emplace(molecule, std::move(atoms), RunningEnergy{}, 1.0, 0.0);
      }

      if (!growData)
      {
        return std::nullopt;
      }

      result.RosenbluthWeight *= growData->RosenbluthWeight;
      result.energies += growData->energies;

      for (const Atom& atom : growData->atoms)
      {
        background.push_back(atom);
      }
      result.molecules.push_back(std::move(*growData));
    }
  }

  return result;
}

[[nodiscard]] std::optional<MoleculeGroupRetraceData> retraceMoleculeGroupDeletion(
    RandomNumber& random, System& system, std::span<const std::pair<std::size_t, std::size_t>> selectedMolecules,
    std::span<const std::pair<std::size_t, std::size_t>> excludeMolecules, bool useCBMC) noexcept
{
  MoleculeGroupRetraceData result;

  std::vector<std::size_t> groupGlobalMoleculeIds;
  groupGlobalMoleculeIds.reserve(selectedMolecules.size());
  for (const auto& [componentId, moleculeId] : selectedMolecules)
  {
    groupGlobalMoleculeIds.push_back(system.moleculeIndexOfComponent(componentId, moleculeId));
  }

  std::unordered_set<std::size_t> extraExcludedGlobalMoleculeIds;
  extraExcludedGlobalMoleculeIds.reserve(excludeMolecules.size());
  for (const auto& [componentId, moleculeId] : excludeMolecules)
  {
    extraExcludedGlobalMoleculeIds.insert(system.moleculeIndexOfComponent(componentId, moleculeId));
  }

  // Determine cutoff distances based on whether dual cutoff is used (only used for the CBMC retrace).
  const double cutOffFrameworkVDW =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffFrameworkVDW;
  const double cutOffMoleculeVDW =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffMoleculeVDW;
  const double cutOffCoulomb =
      system.forceField.useDualCutOff ? system.forceField.dualCutOff : system.forceField.cutOffCoulomb;

  // Detailed balance: the reverse move grows the group sequentially in the order of 'selectedMolecules'
  // (growMoleculeGroupInsertion), where molecule k only sees the group members 0..k-1. The retrace must
  // reproduce these nested environments: molecule k is retraced with the members k+1..m-1 excluded from
  // the background. This is equivalent to retracing in the opposite order of growth and removing each
  // molecule from the background after its weight has been computed. The intra-group interactions are
  // thereby contained once in the Rosenbluth weights and retrace energies.
  for (std::size_t k = 0; k < selectedMolecules.size(); ++k)
  {
    const auto& [componentId, moleculeId] = selectedMolecules[k];

    std::unordered_set<std::size_t> excludedGlobalMoleculeIds = extraExcludedGlobalMoleculeIds;
    excludedGlobalMoleculeIds.insert(groupGlobalMoleculeIds.begin() + static_cast<std::ptrdiff_t>(k) + 1,
                                     groupGlobalMoleculeIds.end());

    std::vector<Atom> background;
    background.reserve(system.spanOfMoleculeAtoms().size());
    for (const Atom& atom : system.spanOfMoleculeAtoms())
    {
      if (!excludedGlobalMoleculeIds.contains(static_cast<std::size_t>(atom.moleculeId)))
      {
        background.push_back(atom);
      }
    }

    const Component& component = system.components[componentId];
    std::span<Atom> moleculeAtoms = system.spanOfMolecule(componentId, moleculeId);

    ChainRetraceData retraceData{RunningEnergy{}, 1.0, 0.0};
    if (useCBMC)
    {
      const CBMC::GrowContext retraceContext{system.hasExternalField, system.forceField, system.simulationBox,
                                             system.interpolationGrids, system.externalFieldInterpolationGrid,
                                             system.framework, system.spanOfFrameworkAtoms(), background, system.beta,
                                             cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb};
      try
      {
        retraceData = CBMC::retraceMoleculeSwapDeletion(random, retraceContext, component, moleculeAtoms);
      }
      catch (const std::runtime_error&)
      {
        return std::nullopt;
      }

      if (system.forceField.useDualCutOff)
      {
        // Dual cut-off scheme: correct the retraced molecule from the inner cut-off to the full
        // cut-offs, using the same nested background as the retrace.
        std::vector<Atom> moleculeCopy(moleculeAtoms.begin(), moleculeAtoms.end());
        std::optional<RunningEnergy> correctionOld =
            CBMC::computeDualCutOffCorrection(retraceContext, component, moleculeCopy);
        if (!correctionOld.has_value())
        {
          return std::nullopt;
        }

        retraceData.energies += correctionOld.value();
        retraceData.RosenbluthWeight *= std::exp(-system.beta * correctionOld->potentialEnergy());
      }
    }

    result.RosenbluthWeight *= retraceData.RosenbluthWeight;
    result.energies += retraceData.energies;
    result.molecules.push_back(std::move(retraceData));
  }

  // A vanishing retrace weight means the current configuration itself is in (near-)overlap: the
  // acceptance ratio would divide by zero and force-accept a move whose retrace energy does not
  // describe the actual state (overlapping trials are discarded), corrupting the running energies.
  // Reject instead, mirroring the minimumRosenbluthFactor guard of the insertion paths.
  if (!(result.RosenbluthWeight >= system.forceField.minimumRosenbluthFactor))
  {
    return std::nullopt;
  }

  return result;
}

[[nodiscard]] double idealGasRosenbluthWeightProduct(const System& system,
                                                     std::span<const std::size_t> stoichiometry) noexcept
{
  // reference ideal-gas Rosenbluth weights: with W/W_IG ratios the full ideal-gas partition functions
  // q_i can be used in the acceptance rules of reaction moves (Rosch and Maginn, eqs. 18-24)
  double weight = 1.0;
  for (std::size_t componentId = 0; componentId < stoichiometry.size(); ++componentId)
  {
    for (std::size_t k = 0; k < stoichiometry[componentId]; ++k)
    {
      weight *= system.components[componentId].idealGasRosenbluthWeight.value_or(1.0);
    }
  }
  return weight;
}

void setReactionFractionalScaling(System& system, Reaction& reaction, double lambda) noexcept
{
  for (std::size_t componentId = 0; componentId < reaction.reactantFractionalMoleculeIds.size(); ++componentId)
  {
    for (const std::size_t moleculeId : reaction.reactantFractionalMoleculeIds[componentId])
    {
      applyLinearReactionScaling(system.spanOfMolecule(componentId, moleculeId), true, lambda,
                                 reaction.dUdlambdaGroup(true));
    }
    for (const std::size_t moleculeId : reaction.productFractionalMoleculeIds[componentId])
    {
      applyLinearReactionScaling(system.spanOfMolecule(componentId, moleculeId), false, lambda,
                                 reaction.dUdlambdaGroup(false));
    }
  }
}

[[nodiscard]] double computeReactionEquilibriumLogTerm(const System& system, const Reaction& reaction,
                                                       ReactionMoveKind moveKind) noexcept
{
  // Rosch and Maginn eq. 24 / eq. 27: prod_i [N_i! / (N_i + nu_i delta)!] q_i^(nu_i delta), evaluated
  // with the whole-molecule counts N_i of the current (old) configuration and the full ideal-gas
  // molecular partition function q_i (proportional to the volume).
  if (moveKind == ReactionMoveKind::LambdaChange)
  {
    return 0.0;
  }

  const bool forward = (moveKind == ReactionMoveKind::ForwardInsert);
  const std::vector<std::size_t>& removeStoichiometry =
      forward ? reaction.reactantStoichiometry : reaction.productStoichiometry;
  const std::vector<std::size_t>& insertStoichiometry =
      forward ? reaction.productStoichiometry : reaction.reactantStoichiometry;

  double term = 0.0;
  const double logVolume = std::log(system.simulationBox.volume);

  for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
  {
    const Component& component = system.components[componentId];
    const std::size_t N = system.numberOfIntegerMoleculesPerComponent[componentId];

    // side losing whole molecules: N! / (N - nu)! divided by (q V)^nu
    for (std::size_t j = 0; j < removeStoichiometry[componentId]; ++j)
    {
      term += std::log(static_cast<double>(N - j));
    }
    term -= static_cast<double>(removeStoichiometry[componentId]) * (component.lnPartitionFunction + logVolume);

    // side gaining whole molecules: (q V)^nu times N! / (N + nu)!
    const std::size_t NAfterRemoval = N - removeStoichiometry[componentId];
    for (std::size_t j = 0; j < insertStoichiometry[componentId]; ++j)
    {
      term -= std::log(static_cast<double>(NAfterRemoval + j + 1));
    }
    term += static_cast<double>(insertStoichiometry[componentId]) * (component.lnPartitionFunction + logVolume);
  }

  return term;
}

[[nodiscard]] std::optional<RunningEnergy> computeGroupSwapEnergyDifference(
    System& system, std::span<const Atom> newAtoms, std::span<const Atom> oldAtoms, bool includeTailCorrections,
    bool includeEwaldCorrections, std::span<const Atom> excludeFromBackground) noexcept
{
  RunningEnergy energyDifference;

  std::unordered_set<std::size_t> oldGlobalMoleculeIds;
  oldGlobalMoleculeIds.reserve(oldAtoms.size() + excludeFromBackground.size());
  for (const Atom& atom : oldAtoms)
  {
    oldGlobalMoleculeIds.insert(static_cast<std::size_t>(atom.moleculeId));
  }
  for (const Atom& atom : excludeFromBackground)
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
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
        system.simulationBox, newAtoms, oldAtoms, system.netCharge);
    energyDifference += ewaldDifference;
  }

  if (includeTailCorrections)
  {
    energyDifference += Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                                                background, newAtoms, oldAtoms) +
                        Interactions::computeFrameworkMoleculeTailEnergyDifference(
                            system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), newAtoms, oldAtoms);
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

void deleteSelectedMolecules(System& system,
                             std::span<const std::pair<std::size_t, std::size_t>> selectedMolecules) noexcept
{
  std::vector<std::pair<std::size_t, std::size_t>> sorted(selectedMolecules.begin(), selectedMolecules.end());
  std::ranges::sort(sorted,
                    [](const auto& a, const auto& b)
                    {
                      if (a.first != b.first)
                      {
                        return a.first < b.first;
                      }
                      return a.second > b.second;
                    });

  for (const auto& [componentId, moleculeId] : sorted)
  {
    // swap the doomed molecule to the last slot of its component so that deleteMolecule pops the tail
    const std::size_t lastIntegerIndex = system.numberOfMoleculesPerComponent[componentId] - 1;
    if (moleculeId != lastIntegerIndex)
    {
      std::span<Atom> moleculeA = system.spanOfMolecule(componentId, moleculeId);
      std::span<Atom> moleculeB = system.spanOfMolecule(componentId, lastIntegerIndex);
      std::swap_ranges(moleculeA.begin(), moleculeA.end(), moleculeB.begin());
      std::swap(system.moleculeData[system.moleculeIndexOfComponent(componentId, moleculeId)],
                system.moleculeData[system.moleculeIndexOfComponent(componentId, lastIntegerIndex)]);
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
      std::vector<Atom> acceptedAtoms(data.atoms.begin(), data.atoms.end());
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

void applySerialFractionalScaling(std::span<Atom> atoms, double lambda, std::uint8_t dUdlambdaGroupId) noexcept
{
  // staged schedule from scaling.ixx: VDW switches on for lambda in [0, 0.5], Coulomb in [0.5, 1]
  for (Atom& atom : atoms)
  {
    atom.setScalingToFractional(lambda, dUdlambdaGroupId);
  }
}

void setSerialReactionFractionalScaling(System& system, Reaction& reaction, double lambda) noexcept
{
  std::vector<std::vector<std::size_t>>& activeIds = reaction.fractionalSideIsReactants
                                                         ? reaction.reactantFractionalMoleculeIds
                                                         : reaction.productFractionalMoleculeIds;
  for (std::size_t componentId = 0; componentId < activeIds.size(); ++componentId)
  {
    for (const std::size_t moleculeId : activeIds[componentId])
    {
      // serial mode: both sides share one dU/dlambda group (only one side is fractional at a time)
      applySerialFractionalScaling(system.spanOfMolecule(componentId, moleculeId), lambda,
                                   reaction.lambda.dUdlambdaGroupId);
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

  std::ranges::sort(toDelete,
                    [](const auto& a, const auto& b)
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

void insertSerialSideFractionalMolecules(System& system, Reaction& reaction, std::span<const ChainGrowData> growData,
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
      const std::size_t moleculeIndex = system.serialReactionFractionalMoleculeIndex(reaction.id, componentId, n);
      system.insertSerialReactionFractionalMolecule(componentId, moleculeIndex, data.molecule, data.atoms, lambda,
                                                    reaction.lambda.dUdlambdaGroupId);
      targetIds[componentId].push_back(moleculeIndex);
    }
  }
  system.syncReactionFractionalMoleculeIndices();
}

[[nodiscard]] std::vector<Atom> collectActiveSerialFractionalAtoms(System& system, Reaction& reaction) noexcept
{
  std::vector<Atom> atoms;
  std::vector<std::vector<std::size_t>>& activeIds = reaction.fractionalSideIsReactants
                                                         ? reaction.reactantFractionalMoleculeIds
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

[[nodiscard]] PropertyLambdaProbabilityHistogram& activeSerialLambdaHistogram(Reaction& reaction) noexcept
{
  return reaction.fractionalSideIsReactants ? reaction.lambda : reaction.lambdaProductSide;
}

[[nodiscard]] PropertyLambdaProbabilityHistogram& inactiveSerialLambdaHistogram(Reaction& reaction) noexcept
{
  return reaction.fractionalSideIsReactants ? reaction.lambdaProductSide : reaction.lambda;
}

[[nodiscard]] double serialBiasFactor(const PropertyLambdaProbabilityHistogram& histogram, double lambda) noexcept
{
  const std::size_t bin =
      std::min(static_cast<std::size_t>(static_cast<double>(histogram.numberOfSamplePoints) * lambda),
               histogram.numberOfSamplePoints > 0 ? histogram.numberOfSamplePoints - 1 : 0uz);
  return histogram.biasFactor[bin];
}

// sum of intramolecular energies, one molecule per distinct moleculeId in 'atoms'
[[nodiscard]] RunningEnergy internalEnergyFromAtomGroups(System& system, std::span<const Atom> atoms) noexcept
{
  RunningEnergy internal;
  for (const std::vector<Atom>& molecule : splitAtomsByMoleculeId(atoms))
  {
    const std::size_t componentId = static_cast<std::size_t>(molecule.front().componentId);
    internal += system.components[componentId].intraMolecularPotentials.computeInternalEnergies(molecule);
  }
  return internal;
}

// Monolithic serial Rx/CFC move (Rosch and Maginn): one function containing the three sub-moves.
// The sub-move is chosen randomly based on the lambda switch point, or forced via 'forcedKind'
// (used by the tests to exercise a specific sub-move deterministically).
[[nodiscard]] std::optional<RunningEnergy> serialReactionMove(RandomNumber& random, System& system, Reaction& reaction,
                                                              Move::Types move, bool useCBMC,
                                                              std::optional<SerialMoveKind> forcedKind) noexcept
{
  // the active fractional side must have molecules to operate on
  const std::vector<std::size_t>& activeStoichiometry =
      reaction.fractionalSideIsReactants ? reaction.reactantStoichiometry : reaction.productStoichiometry;
  if (totalStoichiometry(activeStoichiometry) == 0)
  {
    return std::nullopt;
  }

  // sub-move selection (paper section 3.3): below the switch point mix lambda-changes with
  // fractional reactions, above it with whole-molecule reactions
  SerialMoveKind kind = SerialMoveKind::LambdaChange;
  if (forcedKind.has_value())
  {
    kind = forcedKind.value();
  }
  else if (reaction.lambdaSwitchPoint <= 1.0)
  {
    if (reaction.currentLambda < reaction.lambdaSwitchPoint)
    {
      kind = random.uniform() < 0.5 ? SerialMoveKind::LambdaChange : SerialMoveKind::FractionalReaction;
    }
    else
    {
      kind = random.uniform() < 0.5 ? SerialMoveKind::LambdaChange : SerialMoveKind::WholeMoleculeReaction;
    }
  }

  //================================================================================================================
  // Sub-move 1: lambda change of the active fractional molecules (paper eq. S26)
  //================================================================================================================
  if (kind == SerialMoveKind::LambdaChange)
  {
    PropertyLambdaProbabilityHistogram& lambdaHistogram = activeSerialLambdaHistogram(reaction);
    const double lambdaOld = reaction.currentLambda;
    const std::size_t oldBin =
        std::min(static_cast<std::size_t>(static_cast<double>(lambdaHistogram.numberOfSamplePoints) * lambdaOld),
                 lambdaHistogram.numberOfSamplePoints > 0 ? lambdaHistogram.numberOfSamplePoints - 1 : 0uz);
    const double biasOld = lambdaHistogram.biasFactor[oldBin];

    const double maximumChange =
        reaction.fractionalSideIsReactants ? reaction.maximumLambdaChange : reaction.maximumLambdaChangeProducts;
    const double lambdaNew = lambdaOld + (2.0 * random.uniform() - 1.0) * maximumChange;
    if (lambdaNew < 0.0 || lambdaNew > 1.0)
    {
      return std::nullopt;
    }

    const std::size_t newBin =
        std::min(static_cast<std::size_t>(static_cast<double>(lambdaHistogram.numberOfSamplePoints) * lambdaNew),
                 lambdaHistogram.numberOfSamplePoints > 0 ? lambdaHistogram.numberOfSamplePoints - 1 : 0uz);
    const double biasNew = lambdaHistogram.biasFactor[newBin];

    // energy difference of rescaling the active fractional molecules to the new lambda
    std::optional<RunningEnergy> scalingDifference = RunningEnergy{};
    const std::vector<Atom> oldAtoms = collectActiveSerialFractionalAtoms(system, reaction);
    if (!oldAtoms.empty())
    {
      std::vector<Atom> newAtoms = oldAtoms;
      applySerialFractionalScaling(newAtoms, lambdaNew, reaction.lambda.dUdlambdaGroupId);

      const bool useChargedEwald = system.forceField.useCharge;
      scalingDifference = computeGroupSwapEnergyDifference(system, newAtoms, oldAtoms, true, !useChargedEwald);
      if (!scalingDifference)
      {
        return std::nullopt;
      }
      if (useChargedEwald)
      {
        scalingDifference = scalingDifference.value() + computeChargedGroupEwaldDifference(system, newAtoms, oldAtoms);
      }
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

  //================================================================================================================
  // Sub-move 2: reaction for fractional molecules (paper eq. S29 / eq. 7): remove the fractional
  // molecules of the current side and insert fractional molecules of the other side at the same
  // coupling lambda; delta flips while lambda, the number of whole molecules, and all other positions
  // stay the same. The new molecules are grown with the Rosenbluth scheme (allowed by the paper for
  // the internal configuration), so detailed balance requires the ratio W_new / W_old where W_old
  // comes from retracing the removed fractional molecules. The V^(+/-nu) q^(+/-nu) factors are part
  // of the acceptance rule (equilibriumTerm below).
  //================================================================================================================
  if (kind == SerialMoveKind::FractionalReaction)
  {
    const bool reactantsToProducts = reaction.fractionalSideIsReactants;
    const std::vector<std::size_t>& newStoichiometry =
        reactantsToProducts ? reaction.productStoichiometry : reaction.reactantStoichiometry;

    if (totalStoichiometry(newStoichiometry) == 0)
    {
      return std::nullopt;
    }

    // the fractional molecules of the current side that will be removed
    std::vector<std::pair<std::size_t, std::size_t>> removedMolecules;
    const std::vector<std::vector<std::size_t>>& activeIds =
        reactantsToProducts ? reaction.reactantFractionalMoleculeIds : reaction.productFractionalMoleculeIds;
    for (std::size_t componentId = 0; componentId < activeIds.size(); ++componentId)
    {
      for (const std::size_t moleculeId : activeIds[componentId])
      {
        removedMolecules.emplace_back(componentId, moleculeId);
      }
    }

    std::vector<Atom> oldAtoms = collectActiveSerialFractionalAtoms(system, reaction);

    // grow the new fractional molecules at coupling lambda (staged scaling applied via Atom::setScaling)
    std::optional<MoleculeGroupGrowData> growData =
        growMoleculeGroupInsertion(random, system, newStoichiometry, removedMolecules, reaction.currentLambda, true,
                                   reaction.lambda.dUdlambdaGroupId, useCBMC);
    if (!growData)
    {
      return std::nullopt;
    }

    // retrace the removed fractional molecules at their current coupling
    std::optional<MoleculeGroupRetraceData> retraceData =
        retraceMoleculeGroupDeletion(random, system, removedMolecules, {}, useCBMC);
    if (!retraceData)
    {
      return std::nullopt;
    }

    std::vector<Atom> newAtoms;
    for (const ChainGrowData& data : growData->molecules)
    {
      newAtoms.insert(newAtoms.end(), data.atoms.begin(), data.atoms.end());
    }

    const bool useChargedEwald = system.forceField.useCharge;
    RunningEnergy tailDifference{};
    RunningEnergy fourierDifference{};
    RunningEnergy energyDifference{};
    if (useCBMC)
    {
      tailDifference = computeGroupSwapTailEnergyDifference(system, newAtoms, oldAtoms);
      if (useChargedEwald)
      {
        fourierDifference = computeChargedGroupEwaldDifference(system, newAtoms, oldAtoms);
      }
      energyDifference =
          (growData->energies - retraceData->energies) + fourierDifference + tailDifference;
    }
    else
    {
      std::optional<RunningEnergy> completeDifference =
          computeGroupSwapEnergyDifference(system, newAtoms, oldAtoms, true, !useChargedEwald);
      if (!completeDifference)
      {
        return std::nullopt;
      }
      energyDifference = completeDifference.value() +
                         (internalEnergyFromAtomGroups(system, newAtoms) -
                          internalEnergyFromAtomGroups(system, oldAtoms));
      if (useChargedEwald)
      {
        energyDifference += computeChargedGroupEwaldDifference(system, newAtoms, oldAtoms);
      }
    }

    // V^(+/-nu) q^(+/-nu) equilibrium factors: the disappearing side loses, the appearing side gains
    double equilibriumTerm = 0.0;
    {
      const double logVolume = std::log(system.simulationBox.volume);
      const std::vector<std::size_t>& lostStoichiometry =
          reactantsToProducts ? reaction.reactantStoichiometry : reaction.productStoichiometry;
      const std::vector<std::size_t>& gainedStoichiometry =
          reactantsToProducts ? reaction.productStoichiometry : reaction.reactantStoichiometry;
      for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
      {
        const Component& component = system.components[componentId];
        equilibriumTerm -=
            static_cast<double>(lostStoichiometry[componentId]) * (component.lnPartitionFunction + logVolume);
        equilibriumTerm +=
            static_cast<double>(gainedStoichiometry[componentId]) * (component.lnPartitionFunction + logVolume);
      }
    }

    // biased sampling: exp[W(lambda, delta_new) - W(lambda, delta_old)] (paper section 3.3)
    const double biasOld = serialBiasFactor(activeSerialLambdaHistogram(reaction), reaction.currentLambda);
    const double biasNew = serialBiasFactor(inactiveSerialLambdaHistogram(reaction), reaction.currentLambda);

    double acceptanceProbability;
    if (useCBMC)
    {
      const double correctionFactor =
          std::exp(-system.beta * (fourierDifference.potentialEnergy() + tailDifference.potentialEnergy()));
      // grow and retrace referenced to the ideal-gas Rosenbluth weights, so that the full ideal-gas
      // partition functions q_i can be used in the equilibrium term
      const std::vector<std::size_t>& oldStoichiometry =
          reactantsToProducts ? reaction.reactantStoichiometry : reaction.productStoichiometry;
      const double idealGasNew = idealGasRosenbluthWeightProduct(system, newStoichiometry);
      const double idealGasOld = idealGasRosenbluthWeightProduct(system, oldStoichiometry);
      acceptanceProbability =
          ((growData->RosenbluthWeight / idealGasNew) / (retraceData->RosenbluthWeight / idealGasOld)) *
          correctionFactor * std::exp(equilibriumTerm + biasNew - biasOld);
    }
    else
    {
      acceptanceProbability =
          std::exp(-system.beta * energyDifference.potentialEnergy() + equilibriumTerm + biasNew - biasOld);
    }
    if (random.uniform() >= acceptanceProbability)
    {
      return std::nullopt;
    }

    system.mc_moves_statistics.addConstructed(move);
    system.mc_moves_statistics.addAccepted(move);

    deleteReactionSideFractionalMolecules(system, reaction, reaction.fractionalSideIsReactants);
    // flip delta before inserting so that syncReactionFractionalMoleculeIndices assigns the ids of the new side
    reaction.fractionalSideIsReactants = !reaction.fractionalSideIsReactants;
    std::vector<std::vector<std::size_t>>& targetIds =
        reactantsToProducts ? reaction.productFractionalMoleculeIds : reaction.reactantFractionalMoleculeIds;
    insertSerialSideFractionalMolecules(system, reaction, growData->molecules, newStoichiometry, reaction.currentLambda,
                                        targetIds);
    system.syncReactionLambdaBin(reaction);
    acceptChargedEwaldMove(system);

    return energyDifference;
  }

  //================================================================================================================
  // Sub-move 3: reaction for whole molecules (paper eq. S36 / eq. 9): all positions and lambda stay
  // the same. The fractional molecules of the current side are transformed in place into whole
  // molecules, and randomly selected whole molecules of the other side are transformed in place into
  // the new fractional molecules at the same lambda; delta flips. The acceptance rule contains only
  // the factorial terms and exp(-beta dU); ideal-gas partition functions do not appear because the
  // total number of (whole + fractional) molecules of each component is unchanged.
  //================================================================================================================
  const bool fractionalsAreReactants = reaction.fractionalSideIsReactants;
  const std::vector<std::size_t>& promoteStoichiometry =
      fractionalsAreReactants ? reaction.reactantStoichiometry : reaction.productStoichiometry;
  const std::vector<std::size_t>& demoteStoichiometry =
      fractionalsAreReactants ? reaction.productStoichiometry : reaction.reactantStoichiometry;

  if (totalStoichiometry(promoteStoichiometry) == 0 || totalStoichiometry(demoteStoichiometry) == 0)
  {
    return std::nullopt;
  }

  // automatically rejected when there are not enough whole molecules to turn into fractional molecules
  std::vector<std::pair<std::size_t, std::size_t>> selectedMolecules;
  if (!selectRandomIntegerMolecules(random, system, demoteStoichiometry, selectedMolecules))
  {
    return std::nullopt;
  }

  std::vector<Atom> oldFractionalAtoms = collectActiveSerialFractionalAtoms(system, reaction);
  std::vector<Atom> oldIntegerAtoms;
  for (const auto& [componentId, moleculeId] : selectedMolecules)
  {
    std::span<Atom> molecule = system.spanOfMolecule(componentId, moleculeId);
    oldIntegerAtoms.insert(oldIntegerAtoms.end(), molecule.begin(), molecule.end());
  }

  // new state: identical positions, only the coupling changes
  std::vector<Atom> newIntegerAtoms = oldFractionalAtoms;
  for (Atom& atom : newIntegerAtoms)
  {
    atom.setScalingToInteger();
  }
  std::vector<Atom> newFractionalAtoms = oldIntegerAtoms;
  applySerialFractionalScaling(newFractionalAtoms, reaction.currentLambda, reaction.lambda.dUdlambdaGroupId);

  std::vector<Atom> oldAtoms = oldFractionalAtoms;
  oldAtoms.insert(oldAtoms.end(), oldIntegerAtoms.begin(), oldIntegerAtoms.end());
  std::vector<Atom> newAtoms = newIntegerAtoms;
  newAtoms.insert(newAtoms.end(), newFractionalAtoms.begin(), newFractionalAtoms.end());

  const bool useChargedEwald = system.forceField.useCharge;
  std::optional<RunningEnergy> swapDifference =
      computeGroupSwapEnergyDifference(system, newAtoms, oldAtoms, true, !useChargedEwald);
  if (!swapDifference)
  {
    return std::nullopt;
  }

  const RunningEnergy internalDifference =
      internalEnergyFromAtomGroups(system, newAtoms) - internalEnergyFromAtomGroups(system, oldAtoms);
  RunningEnergy energyDifference = swapDifference.value() + internalDifference;

  if (useChargedEwald)
  {
    energyDifference += computeChargedGroupEwaldDifference(system, newAtoms, oldAtoms);
  }

  // factorial terms of eq. S36: the demoted side loses nu whole molecules, the promoted side gains nu
  double factorialTerm = 0.0;
  for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
  {
    const std::size_t N = system.numberOfIntegerMoleculesPerComponent[componentId];
    for (std::size_t k = 0; k < demoteStoichiometry[componentId]; ++k)
    {
      factorialTerm += std::log(static_cast<double>(N - k));
    }
    const std::size_t NAfterDemote = N - demoteStoichiometry[componentId];
    for (std::size_t k = 0; k < promoteStoichiometry[componentId]; ++k)
    {
      factorialTerm -= std::log(static_cast<double>(NAfterDemote + k + 1));
    }
  }

  // biased sampling: exp[W(lambda, delta_new) - W(lambda, delta_old)] (paper section 3.3)
  const double biasOld = serialBiasFactor(activeSerialLambdaHistogram(reaction), reaction.currentLambda);
  const double biasNew = serialBiasFactor(inactiveSerialLambdaHistogram(reaction), reaction.currentLambda);

  const double acceptanceProbability =
      std::exp(factorialTerm + biasNew - biasOld) * std::exp(-system.beta * energyDifference.potentialEnergy());
  if (random.uniform() >= acceptanceProbability)
  {
    return std::nullopt;
  }

  system.mc_moves_statistics.addConstructed(move);
  system.mc_moves_statistics.addAccepted(move);

  // capture the promoted (fractional -> whole) molecules before any bookkeeping changes the indices
  struct TransformedMolecule
  {
    std::size_t componentId;
    Molecule molecule;
    std::vector<Atom> atoms;
  };
  std::vector<TransformedMolecule> promotedMolecules;
  const std::vector<std::vector<std::size_t>>& activeIds =
      fractionalsAreReactants ? reaction.reactantFractionalMoleculeIds : reaction.productFractionalMoleculeIds;
  for (std::size_t componentId = 0; componentId < activeIds.size(); ++componentId)
  {
    for (const std::size_t moleculeId : activeIds[componentId])
    {
      std::span<Atom> molecule = system.spanOfMolecule(componentId, moleculeId);
      std::vector<Atom> atoms(molecule.begin(), molecule.end());
      for (Atom& atom : atoms)
      {
        atom.setScalingToInteger();
      }
      promotedMolecules.push_back({componentId, *system.indexForMolecule(componentId, moleculeId), std::move(atoms)});
    }
  }

  // capture the demoted (whole -> fractional) molecules, ordered by component to match the slot layout
  std::vector<ChainGrowData> demotedMolecules;
  for (const auto& [componentId, moleculeId] : selectedMolecules)
  {
    std::span<Atom> molecule = system.spanOfMolecule(componentId, moleculeId);
    demotedMolecules.emplace_back(*system.indexForMolecule(componentId, moleculeId),
                                  std::vector<Atom>(molecule.begin(), molecule.end()), RunningEnergy{}, 1.0, 0.0);
  }

  // remove the selected whole molecules and re-insert the promoted molecules as whole molecules
  deleteSelectedMolecules(system, selectedMolecules);
  for (TransformedMolecule& promotedMolecule : promotedMolecules)
  {
    system.insertMolecule(promotedMolecule.componentId, promotedMolecule.molecule, std::move(promotedMolecule.atoms));
  }

  // replace the fractional molecules of the old side by the demoted molecules on the new side
  deleteReactionSideFractionalMolecules(system, reaction, fractionalsAreReactants);
  // flip delta before inserting so that syncReactionFractionalMoleculeIndices assigns the ids of the new side
  reaction.fractionalSideIsReactants = !reaction.fractionalSideIsReactants;
  std::vector<std::vector<std::size_t>>& targetIds =
      fractionalsAreReactants ? reaction.productFractionalMoleculeIds : reaction.reactantFractionalMoleculeIds;
  insertSerialSideFractionalMolecules(system, reaction, demotedMolecules, demoteStoichiometry, reaction.currentLambda,
                                      targetIds);
  system.syncReactionLambdaBin(reaction);
  acceptChargedEwaldMove(system);

  return energyDifference;
}

// Boundary-crossing chemical reaction of parallel Rx/CFC (Rosch and Maginn, eq. 27). Following the
// paper, the transformation happens in place: the (nearly fully coupled) fractional molecules of the
// promoted side are turned into whole molecules, randomly selected whole molecules of the other side
// are turned into the new fractional molecules at the wrapped lambda, the nearly decoupled fractional
// molecules of the demoted side are deleted, and new nearly decoupled fractional molecules of the
// promoted side are grown. All couplings follow the staged schedule of scaling.ixx via
// Atom::setScaling: products are coupled with lambda, reactants with (1 - lambda).
[[nodiscard]] static std::optional<RunningEnergy> parallelReactionBoundaryMove(RandomNumber& random, System& system,
                                                                               Reaction& reaction, Move::Types move,
                                                                               double lambdaNew, bool forward,
                                                                               bool useCBMC) noexcept
{
  const std::vector<std::size_t>& demoteStoichiometry =
      forward ? reaction.reactantStoichiometry : reaction.productStoichiometry;
  const std::vector<std::size_t>& growStoichiometry =
      forward ? reaction.productStoichiometry : reaction.reactantStoichiometry;
  const std::vector<std::vector<std::size_t>>& promotedIds =
      forward ? reaction.productFractionalMoleculeIds : reaction.reactantFractionalMoleculeIds;
  const std::vector<std::vector<std::size_t>>& ghostIds =
      forward ? reaction.reactantFractionalMoleculeIds : reaction.productFractionalMoleculeIds;

  if (totalStoichiometry(demoteStoichiometry) == 0 || totalStoichiometry(growStoichiometry) == 0)
  {
    return std::nullopt;
  }

  // automatically rejected when there are not enough whole molecules to demote
  std::vector<std::pair<std::size_t, std::size_t>> selectedMolecules;
  if (!selectRandomIntegerMolecules(random, system, demoteStoichiometry, selectedMolecules))
  {
    return std::nullopt;
  }

  std::vector<std::pair<std::size_t, std::size_t>> ghostMolecules;
  std::vector<std::pair<std::size_t, std::size_t>> promotedMoleculeIds;
  for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
  {
    for (const std::size_t moleculeId : ghostIds[componentId])
    {
      ghostMolecules.emplace_back(componentId, moleculeId);
    }
    for (const std::size_t moleculeId : promotedIds[componentId])
    {
      promotedMoleculeIds.emplace_back(componentId, moleculeId);
    }
  }

  auto collectAtoms = [&system](const std::vector<std::pair<std::size_t, std::size_t>>& molecules)
  {
    std::vector<Atom> atoms;
    for (const auto& [componentId, moleculeId] : molecules)
    {
      std::span<Atom> molecule = system.spanOfMolecule(componentId, moleculeId);
      atoms.insert(atoms.end(), molecule.begin(), molecule.end());
    }
    return atoms;
  };
  const std::vector<Atom> promotedOldAtoms = collectAtoms(promotedMoleculeIds);
  const std::vector<Atom> ghostOldAtoms = collectAtoms(ghostMolecules);
  const std::vector<Atom> demotedOldAtoms = collectAtoms(selectedMolecules);

  // retrace of the deleted ghost fractional molecules against the unmodified old configuration; this
  // provides the Rosenbluth weight of the reverse move
  std::optional<MoleculeGroupRetraceData> retraceData =
      retraceMoleculeGroupDeletion(random, system, ghostMolecules, {}, useCBMC);
  if (!retraceData)
  {
    return std::nullopt;
  }

  // in-place transformations: promoted fractional molecules become whole; the selected whole molecules
  // become the new fractional molecules of the demoted side at the wrapped lambda
  std::vector<Atom> promotedNewAtoms = promotedOldAtoms;
  for (Atom& atom : promotedNewAtoms)
  {
    atom.setScalingToInteger();
  }
  // forward demotes whole reactant molecules into reactant fractionals, backward demotes products
  std::vector<Atom> demotedNewAtoms = demotedOldAtoms;
  applyLinearReactionScaling(demotedNewAtoms, forward, lambdaNew, reaction.dUdlambdaGroup(forward));

  std::vector<Atom> rescaleNewAtoms = promotedNewAtoms;
  rescaleNewAtoms.insert(rescaleNewAtoms.end(), demotedNewAtoms.begin(), demotedNewAtoms.end());
  std::vector<Atom> rescaleOldAtoms = promotedOldAtoms;
  rescaleOldAtoms.insert(rescaleOldAtoms.end(), demotedOldAtoms.begin(), demotedOldAtoms.end());

  // energy difference of the in-place transformations; the deleted ghost molecules are excluded from
  // the background because their removal is accounted for by the retrace Rosenbluth weight
  std::optional<RunningEnergy> rescaleDifference =
      computeGroupSwapEnergyDifference(system, rescaleNewAtoms, rescaleOldAtoms, false, false, ghostOldAtoms);
  if (!rescaleDifference)
  {
    return std::nullopt;
  }

  auto writeAtomsToSpans =
      [&system](const std::vector<std::pair<std::size_t, std::size_t>>& molecules, std::span<const Atom> source)
  {
    std::size_t offset = 0;
    for (const auto& [componentId, moleculeId] : molecules)
    {
      std::span<Atom> molecule = system.spanOfMolecule(componentId, moleculeId);
      std::copy_n(source.begin() + static_cast<std::ptrdiff_t>(offset), molecule.size(), molecule.begin());
      offset += molecule.size();
    }
  };
  auto restoreOldScalings = [&]()
  {
    writeAtomsToSpans(promotedMoleculeIds, promotedOldAtoms);
    writeAtomsToSpans(selectedMolecules, demotedOldAtoms);
  };

  // apply the new scalings in memory so that the CBMC growth of the new ghost fractional molecules
  // sees the final configuration
  writeAtomsToSpans(promotedMoleculeIds, promotedNewAtoms);
  writeAtomsToSpans(selectedMolecules, demotedNewAtoms);

  const double growScaling = forward ? lambdaNew : (1.0 - lambdaNew);
  // forward grows product-side ghost fractionals, backward grows reactant-side ones
  std::optional<MoleculeGroupGrowData> growData =
      growMoleculeGroupInsertion(random, system, growStoichiometry, ghostMolecules, growScaling, true,
                                 reaction.dUdlambdaGroup(!forward), useCBMC);
  if (!growData)
  {
    restoreOldScalings();
    return std::nullopt;
  }

  std::vector<Atom> grownAtoms;
  for (const ChainGrowData& data : growData->molecules)
  {
    grownAtoms.insert(grownAtoms.end(), data.atoms.begin(), data.atoms.end());
  }

  std::vector<Atom> fullNewAtoms = rescaleNewAtoms;
  fullNewAtoms.insert(fullNewAtoms.end(), grownAtoms.begin(), grownAtoms.end());
  std::vector<Atom> fullOldAtoms = rescaleOldAtoms;
  fullOldAtoms.insert(fullOldAtoms.end(), ghostOldAtoms.begin(), ghostOldAtoms.end());

  RunningEnergy tailDifference{};
  RunningEnergy fourierDifference{};
  RunningEnergy energyDifference{};
  if (useCBMC)
  {
    tailDifference = computeGroupSwapTailEnergyDifference(system, fullNewAtoms, fullOldAtoms);
    if (system.forceField.useCharge)
    {
      fourierDifference = computeChargedGroupEwaldDifference(system, fullNewAtoms, fullOldAtoms);
    }
    energyDifference =
        rescaleDifference.value() + (growData->energies - retraceData->energies) + tailDifference + fourierDifference;
  }
  else
  {
    std::optional<RunningEnergy> completeDifference = computeGroupSwapEnergyDifference(
        system, fullNewAtoms, fullOldAtoms, true, !system.forceField.useCharge);
    if (!completeDifference)
    {
      restoreOldScalings();
      return std::nullopt;
    }
    energyDifference = completeDifference.value() +
                       (internalEnergyFromAtomGroups(system, fullNewAtoms) -
                        internalEnergyFromAtomGroups(system, fullOldAtoms));
    if (system.forceField.useCharge)
    {
      energyDifference += computeChargedGroupEwaldDifference(system, fullNewAtoms, fullOldAtoms);
    }
  }

  const ReactionMoveKind moveKind = forward ? ReactionMoveKind::ForwardInsert : ReactionMoveKind::BackwardDelete;
  const double equilibriumTerm = computeReactionEquilibriumLogTerm(system, reaction, moveKind);

  // one-dimensional weight function eta(lambda) of parallel Rx/CFC
  const double biasOld = serialBiasFactor(reaction.lambda, reaction.currentLambda);
  const double biasNew = serialBiasFactor(reaction.lambda, lambdaNew);

  double acceptanceProbability;
  if (useCBMC)
  {
    // grow and retrace referenced to the ideal-gas Rosenbluth weights, so that the full ideal-gas
    // partition functions q_i can be used in the acceptance rule
    const double idealGasGrow = idealGasRosenbluthWeightProduct(system, growStoichiometry);
    const double idealGasGhost = idealGasRosenbluthWeightProduct(system, demoteStoichiometry);
    const double rosenbluthRatio =
        (growData->RosenbluthWeight / idealGasGrow) / (retraceData->RosenbluthWeight / idealGasGhost);
    acceptanceProbability =
        rosenbluthRatio * std::exp(equilibriumTerm + biasNew - biasOld) *
        std::exp(-system.beta * (rescaleDifference->potentialEnergy() + tailDifference.potentialEnergy() +
                                 fourierDifference.potentialEnergy()));
  }
  else
  {
    acceptanceProbability =
        std::exp(equilibriumTerm + biasNew - biasOld - system.beta * energyDifference.potentialEnergy());
  }

  if (random.uniform() >= acceptanceProbability)
  {
    restoreOldScalings();
    return std::nullopt;
  }

  system.mc_moves_statistics.addConstructed(move);
  system.mc_moves_statistics.addAccepted(move);

  // capture the promoted (fractional -> whole) molecules before any bookkeeping changes the indices
  struct TransformedMolecule
  {
    std::size_t componentId;
    Molecule molecule;
    std::vector<Atom> atoms;
  };
  std::vector<TransformedMolecule> promotedMolecules;
  for (const auto& [componentId, moleculeId] : promotedMoleculeIds)
  {
    std::span<Atom> molecule = system.spanOfMolecule(componentId, moleculeId);
    promotedMolecules.push_back({componentId, *system.indexForMolecule(componentId, moleculeId),
                                 std::vector<Atom>(molecule.begin(), molecule.end())});
  }

  // capture the demoted (whole -> fractional) molecules, ordered by component to match the slot layout
  std::vector<TransformedMolecule> demotedMolecules;
  for (const auto& [componentId, moleculeId] : selectedMolecules)
  {
    std::span<Atom> molecule = system.spanOfMolecule(componentId, moleculeId);
    demotedMolecules.push_back({componentId, *system.indexForMolecule(componentId, moleculeId),
                                std::vector<Atom>(molecule.begin(), molecule.end())});
  }

  // remove the demoted whole molecules from the integer region
  deleteSelectedMolecules(system, selectedMolecules);

  // remove the fractional molecules of this reaction on both sides (descending ids per component)
  {
    std::vector<std::pair<std::size_t, std::size_t>> toDelete;
    for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
    {
      for (const std::size_t moleculeId : reaction.reactantFractionalMoleculeIds[componentId])
      {
        toDelete.emplace_back(componentId, moleculeId);
      }
      for (const std::size_t moleculeId : reaction.productFractionalMoleculeIds[componentId])
      {
        toDelete.emplace_back(componentId, moleculeId);
      }
    }
    std::ranges::sort(toDelete,
                      [](const auto& a, const auto& b)
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
    for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
    {
      reaction.reactantFractionalMoleculeIds[componentId].clear();
      reaction.productFractionalMoleculeIds[componentId].clear();
    }
  }

  // insert the promoted molecules as whole molecules
  for (TransformedMolecule& promoted : promotedMolecules)
  {
    system.insertMolecule(promoted.componentId, promoted.molecule, std::move(promoted.atoms));
  }

  // refill the fractional slots in ascending slot order: for a forward reaction the reactant slots
  // receive the demoted whole molecules and the product slots the freshly grown molecules (backward:
  // the mirror image)
  std::size_t demotedIndex = 0;
  std::size_t grownIndex = 0;
  for (std::size_t componentId = 0; componentId < system.components.size(); ++componentId)
  {
    for (std::size_t k = 0; k < reaction.reactantStoichiometry[componentId]; ++k)
    {
      const std::size_t slotIndex = system.parallelReactionFractionalMoleculeIndex(reaction.id, componentId, false, k);
      if (forward)
      {
        TransformedMolecule& demoted = demotedMolecules[demotedIndex++];
        system.insertReactionFractionalMolecule(componentId, slotIndex, demoted.molecule, std::move(demoted.atoms),
                                                true, lambdaNew, reaction.dUdlambdaGroup(true));
      }
      else
      {
        const ChainGrowData& grown = growData->molecules[grownIndex++];
        system.insertReactionFractionalMolecule(componentId, slotIndex, grown.molecule, grown.atoms, true, lambdaNew,
                                                reaction.dUdlambdaGroup(true));
      }
    }
    for (std::size_t k = 0; k < reaction.productStoichiometry[componentId]; ++k)
    {
      const std::size_t slotIndex = system.parallelReactionFractionalMoleculeIndex(reaction.id, componentId, true, k);
      if (forward)
      {
        const ChainGrowData& grown = growData->molecules[grownIndex++];
        system.insertReactionFractionalMolecule(componentId, slotIndex, grown.molecule, grown.atoms, false, lambdaNew,
                                                reaction.dUdlambdaGroup(false));
      }
      else
      {
        TransformedMolecule& demoted = demotedMolecules[demotedIndex++];
        system.insertReactionFractionalMolecule(componentId, slotIndex, demoted.molecule, std::move(demoted.atoms),
                                                false, lambdaNew, reaction.dUdlambdaGroup(false));
      }
    }
  }

  system.syncReactionFractionalMoleculeIndices();
  reaction.currentLambda = lambdaNew;
  system.syncReactionLambdaBin(reaction);
  acceptChargedEwaldMove(system);

  return energyDifference;
}

[[nodiscard]] std::optional<RunningEnergy> parallelReactionMove(RandomNumber& random, System& system,
                                                                Move::Types move, bool useCBMC) noexcept
{
  if (system.reactions.list.empty())
  {
    return std::nullopt;
  }

  system.mc_moves_statistics.addTrial(move);

  std::vector<Reaction*> parallelReactions;
  parallelReactions.reserve(system.reactions.list.size());
  for (Reaction& candidate : system.reactions.list)
  {
    if (candidate.reactionMove == move &&
        candidate.reactantFractionalMoleculeIds.size() == system.components.size())
    {
      parallelReactions.push_back(&candidate);
    }
  }
  if (parallelReactions.empty())
  {
    return std::nullopt;
  }

  Reaction& reaction = *parallelReactions[random.uniform_integer(0, parallelReactions.size() - 1)];

  const double lambdaOld = reaction.currentLambda;
  double lambdaNew = lambdaOld + (2.0 * random.uniform() - 1.0) * reaction.maximumLambdaChange;

  // boundary crossing: attempt a chemical reaction with the wrapped lambda (Rosch and Maginn)
  if (lambdaNew > 1.0)
  {
    return parallelReactionBoundaryMove(random, system, reaction, move, lambdaNew - 1.0, true, useCBMC);
  }
  if (lambdaNew < 0.0)
  {
    return parallelReactionBoundaryMove(random, system, reaction, move, lambdaNew + 1.0, false, useCBMC);
  }

  // lambda change within [0, 1] (Rosch and Maginn, eq. 26)
  const double biasOld = serialBiasFactor(reaction.lambda, lambdaOld);
  const double biasNew = serialBiasFactor(reaction.lambda, lambdaNew);

  // energy difference of rescaling all reaction fractional molecules to the new lambda:
  // reactants are coupled with (1 - lambda), products with lambda
  std::optional<RunningEnergy> scalingDifference = RunningEnergy{};
  {
    std::vector<Atom> oldAtoms;
    std::vector<Atom> newAtoms;
    for (std::size_t componentId = 0; componentId < reaction.reactantFractionalMoleculeIds.size(); ++componentId)
    {
      for (const std::size_t moleculeId : reaction.reactantFractionalMoleculeIds[componentId])
      {
        std::span<Atom> molecule = system.spanOfMolecule(componentId, moleculeId);
        oldAtoms.insert(oldAtoms.end(), molecule.begin(), molecule.end());
        std::vector<Atom> trialMolecule(molecule.begin(), molecule.end());
        applyLinearReactionScaling(trialMolecule, true, lambdaNew, reaction.dUdlambdaGroup(true));
        newAtoms.insert(newAtoms.end(), trialMolecule.begin(), trialMolecule.end());
      }
      for (const std::size_t moleculeId : reaction.productFractionalMoleculeIds[componentId])
      {
        std::span<Atom> molecule = system.spanOfMolecule(componentId, moleculeId);
        oldAtoms.insert(oldAtoms.end(), molecule.begin(), molecule.end());
        std::vector<Atom> trialMolecule(molecule.begin(), molecule.end());
        applyLinearReactionScaling(trialMolecule, false, lambdaNew, reaction.dUdlambdaGroup(false));
        newAtoms.insert(newAtoms.end(), trialMolecule.begin(), trialMolecule.end());
      }
    }

    if (!oldAtoms.empty())
    {
      scalingDifference = computeGroupSwapEnergyDifference(system, newAtoms, oldAtoms, true, true);
      if (!scalingDifference)
      {
        return std::nullopt;
      }
    }
  }

  const double acceptanceProbability =
      std::exp(-system.beta * scalingDifference->potentialEnergy() + biasNew - biasOld);
  if (random.uniform() >= acceptanceProbability)
  {
    return std::nullopt;
  }

  system.mc_moves_statistics.addConstructed(move);
  system.mc_moves_statistics.addAccepted(move);

  Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);
  setReactionFractionalScaling(system, reaction, lambdaNew);
  reaction.currentLambda = lambdaNew;
  system.syncReactionLambdaBin(reaction);

  return scalingDifference;
}

}  // namespace MC_Moves::ReactionCommon
