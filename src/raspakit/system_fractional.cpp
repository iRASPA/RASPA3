module;

module system;

import std;

import randomnumbers;
import atom;
import component;
import reaction;
import reactions;
import mc_moves_move_types;
import mc_moves_probabilities;
import property_lambda_probability_histogram;
import running_energy;
import forcefield;
import simulationbox;
import framework;
import cbmc;
import cbmc_chain_data;
import interpolation_energy_grid;
import interactions_external_field;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_polarization;
import scaling;

// System fractional molecules: CFCMC slots, lambda histograms, and reaction bookkeeping.

[[nodiscard]] std::optional<RunningEnergy> conventionalReactionInitializationEnergy(
    System& system, std::size_t componentId, std::span<const Atom> atoms) noexcept
{
  std::optional<RunningEnergy> external = Interactions::computeExternalFieldEnergyDifference(
      system.hasExternalField, system.forceField, system.simulationBox, system.externalFieldInterpolationGrid, atoms,
      {});
  std::optional<RunningEnergy> framework = Interactions::computeFrameworkMoleculeEnergyDifference(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
      system.spanOfFrameworkAtoms(), atoms, {});
  std::optional<RunningEnergy> intermolecular = Interactions::computeInterMolecularEnergyDifference(
      system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), atoms, {});
  if (!external || !framework || !intermolecular)
  {
    return std::nullopt;
  }
  return external.value() + framework.value() + intermolecular.value() +
         system.components[componentId].intraMolecularPotentials.computeInternalEnergies(atoms);
}

// Whether the component pins the given lambda coordinate at a fixed bin (fixed-lambda
// thermodynamic integration, "SimulationType": "ThermodynamicIntegration"). The corresponding
// fractional molecules then exist without any lambda-changing move being enabled.
static bool pinsFixedLambda(const Component &component, Component::FixedLambdaCoordinate coordinate)
{
  return component.fixedLambdaBin.has_value() && component.fixedLambdaCoordinate == coordinate;
}

void System::determineFractionalComponents()
{
  for (std::size_t i = 0; i < components.size(); ++i)
  {
    numberOfFractionalMoleculesPerComponent[i] = 0;
    numberOfGCFractionalMoleculesPerComponent_CFCMC[i] = 0;
    numberOfPairGCFractionalMoleculesPerComponent_CFCMC[i] = 0;
    numberOfPairSwapFractionalMoleculesPerComponent_CFCMC[i] = 0;
    numberOfPairSwapCBFractionalMoleculesPerComponent_CFCMC[i] = 0;
    numberOfGroupSwapFractionalMoleculesPerComponent_CFCMC[i] = 0;
    numberOfGroupSwapCBFractionalMoleculesPerComponent_CFCMC[i] = 0;
    numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC[i] = 0;
    numberOfGibbsFractionalMoleculesPerComponent_CFCMC[i] = 0;
    numberOfParallelReactionFractionalMoleculesPerComponent_CFCMC[i] = 0;
    numberOfSerialReactionFractionalMoleculesPerComponent_CFCMC[i] = 0;

    // fixed-lambda thermodynamic integration pins a fractional molecule in the GC slot without
    // any lambda-changing move being enabled
    const bool needsGCSlot = components[i].mc_moves_probabilities.getProbability(Move::Types::SwapCFCMC) > 0.0 ||
                             components[i].mc_moves_probabilities.getProbability(Move::Types::WidomCFCMC) > 0.0 ||
                             pinsFixedLambda(components[i], Component::FixedLambdaCoordinate::GC);
    if (needsGCSlot)
    {
      numberOfGCFractionalMoleculesPerComponent_CFCMC[i] = 1;
      components[i].hasFractionalMolecule = true;
    }

    if (components[i].mc_moves_probabilities.getProbability(Move::Types::SwapCBCFCMC) > 0.0 ||
        components[i].mc_moves_probabilities.getProbability(Move::Types::WidomCBCFCMC) > 0.0)
    {
      numberOfPairGCFractionalMoleculesPerComponent_CFCMC[i] = 1;
      components[i].hasFractionalMolecule = true;
    }

    if (components[i].mc_moves_probabilities.getProbability(Move::Types::GibbsSwapCFCMC) > 0.0 ||
        components[i].mc_moves_probabilities.getProbability(Move::Types::GibbsSwapCBCFCMC) > 0.0)
    {
      numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC[i] = 1;
      components[i].hasFractionalMolecule = true;
    }

    if (components[i].mc_moves_probabilities.getProbability(Move::Types::GibbsConventionalCFCMC) > 0.0 ||
        components[i].mc_moves_probabilities.getProbability(Move::Types::GibbsConventionalCBCFCMC) > 0.0)
    {
      numberOfGibbsFractionalMoleculesPerComponent_CFCMC[i] = 1;
      components[i].hasFractionalMolecule = true;
    }
  }

  // ion-pair CFCMC moves need a fractional molecule on both components of the pair
  // (done in a separate pass because a component allocates the slot of its pair partner)
  for (std::size_t i = 0; i < components.size(); ++i)
  {
    if (!components[i].pairComponentId.has_value() || components[i].pairComponentId.value() >= components.size())
    {
      continue;
    }
    const std::size_t partner = components[i].pairComponentId.value();

    if (components[i].mc_moves_probabilities.getProbability(Move::Types::PairSwapCFCMC) > 0.0 ||
        pinsFixedLambda(components[i], Component::FixedLambdaCoordinate::PairSwap))
    {
      numberOfPairSwapFractionalMoleculesPerComponent_CFCMC[i] = 1;
      numberOfPairSwapFractionalMoleculesPerComponent_CFCMC[partner] = 1;
      components[i].hasFractionalMolecule = true;
      components[partner].hasFractionalMolecule = true;
    }

    if (components[i].mc_moves_probabilities.getProbability(Move::Types::PairSwapCBCFCMC) > 0.0)
    {
      numberOfPairSwapCBFractionalMoleculesPerComponent_CFCMC[i] = 1;
      numberOfPairSwapCBFractionalMoleculesPerComponent_CFCMC[partner] = 1;
      components[i].hasFractionalMolecule = true;
      components[partner].hasFractionalMolecule = true;
    }
  }

  // group CFCMC moves need one fractional molecule for the central component plus one per satellite
  // occurrence (a component occurring twice among the satellites, e.g. Cl of CaCl2, gets two slots)
  for (std::size_t i = 0; i < components.size(); ++i)
  {
    if (components[i].groupComponentIds.empty())
    {
      continue;
    }

    if (components[i].mc_moves_probabilities.getProbability(Move::Types::GroupSwapCFCMC) > 0.0 ||
        pinsFixedLambda(components[i], Component::FixedLambdaCoordinate::GroupSwap))
    {
      numberOfGroupSwapFractionalMoleculesPerComponent_CFCMC[i] += 1;
      components[i].hasFractionalMolecule = true;
      for (std::size_t satelliteComponentId : components[i].groupComponentIds)
      {
        if (satelliteComponentId < components.size())
        {
          numberOfGroupSwapFractionalMoleculesPerComponent_CFCMC[satelliteComponentId] += 1;
          components[satelliteComponentId].hasFractionalMolecule = true;
        }
      }
    }

    if (components[i].mc_moves_probabilities.getProbability(Move::Types::GroupSwapCBCFCMC) > 0.0)
    {
      numberOfGroupSwapCBFractionalMoleculesPerComponent_CFCMC[i] += 1;
      components[i].hasFractionalMolecule = true;
      for (std::size_t satelliteComponentId : components[i].groupComponentIds)
      {
        if (satelliteComponentId < components.size())
        {
          numberOfGroupSwapCBFractionalMoleculesPerComponent_CFCMC[satelliteComponentId] += 1;
          components[satelliteComponentId].hasFractionalMolecule = true;
        }
      }
    }
  }

  numberOfReactionFractionalMoleculesPerComponent_CFCMC.assign(
      reactions.list.size(), std::vector<std::size_t>(components.size(), 0));
}

std::size_t System::indexOfGCFractionalMoleculesPerComponent_CFCMC(
    [[maybe_unused]] std::size_t selectedComponent) const noexcept
{
  return 0;
}

std::size_t System::indexOfPairGCFractionalMoleculesPerComponent_CFCMC(std::size_t selectedComponent) const noexcept
{
  return numberOfGCFractionalMoleculesPerComponent_CFCMC[selectedComponent];
}

std::size_t System::indexOfPairSwapFractionalMoleculesPerComponent_CFCMC(std::size_t selectedComponent) const noexcept
{
  return numberOfGCFractionalMoleculesPerComponent_CFCMC[selectedComponent] +
         numberOfPairGCFractionalMoleculesPerComponent_CFCMC[selectedComponent];
}

std::size_t System::indexOfPairSwapCBFractionalMoleculesPerComponent_CFCMC(std::size_t selectedComponent) const noexcept
{
  return indexOfPairSwapFractionalMoleculesPerComponent_CFCMC(selectedComponent) +
         numberOfPairSwapFractionalMoleculesPerComponent_CFCMC[selectedComponent];
}

std::size_t System::indexOfGroupSwapFractionalMoleculesPerComponent_CFCMC(std::size_t selectedComponent) const noexcept
{
  return indexOfPairSwapCBFractionalMoleculesPerComponent_CFCMC(selectedComponent) +
         numberOfPairSwapCBFractionalMoleculesPerComponent_CFCMC[selectedComponent];
}

std::size_t System::indexOfGroupSwapCBFractionalMoleculesPerComponent_CFCMC(
    std::size_t selectedComponent) const noexcept
{
  return indexOfGroupSwapFractionalMoleculesPerComponent_CFCMC(selectedComponent) +
         numberOfGroupSwapFractionalMoleculesPerComponent_CFCMC[selectedComponent];
}

std::size_t System::indexOfGibbsSwapFractionalMoleculesPerComponent_CFCMC(std::size_t selectedComponent) const noexcept
{
  return indexOfGroupSwapCBFractionalMoleculesPerComponent_CFCMC(selectedComponent) +
         numberOfGroupSwapCBFractionalMoleculesPerComponent_CFCMC[selectedComponent];
}

std::size_t System::indexOfParallelReactionFractionalMoleculesPerComponent_CFCMC(
    std::size_t selectedComponent) const noexcept
{
  return indexOfGibbsSwapFractionalMoleculesPerComponent_CFCMC(selectedComponent) +
         numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC[selectedComponent];
}

std::size_t System::indexOfSerialReactionFractionalMoleculesPerComponent_CFCMC(
    std::size_t selectedComponent) const noexcept
{
  return indexOfParallelReactionFractionalMoleculesPerComponent_CFCMC(selectedComponent) +
         numberOfParallelReactionFractionalMoleculesPerComponent_CFCMC[selectedComponent];
}

std::size_t System::indexOfGibbsConventionalFractionalMoleculesPerComponent_CFCMC(
    std::size_t selectedComponent) const noexcept
{
  return indexOfSerialReactionFractionalMoleculesPerComponent_CFCMC(selectedComponent) +
         numberOfSerialReactionFractionalMoleculesPerComponent_CFCMC[selectedComponent];
}

std::size_t System::indexOfFractionalMoleculeForMove(Move::Types move, std::size_t selectedComponent,
                                                     std::size_t subIndex) const noexcept
{
  switch (move)
  {
    case Move::Types::SwapCFCMC:
    case Move::Types::WidomCFCMC:
      return indexOfGCFractionalMoleculesPerComponent_CFCMC(selectedComponent) + subIndex;
    case Move::Types::SwapCBCFCMC:
    case Move::Types::WidomCBCFCMC:
      return indexOfPairGCFractionalMoleculesPerComponent_CFCMC(selectedComponent) + subIndex;
    case Move::Types::PairSwapCFCMC:
      return indexOfPairSwapFractionalMoleculesPerComponent_CFCMC(selectedComponent) + subIndex;
    case Move::Types::PairSwapCBCFCMC:
      return indexOfPairSwapCBFractionalMoleculesPerComponent_CFCMC(selectedComponent) + subIndex;
    case Move::Types::GroupSwapCFCMC:
      return indexOfGroupSwapFractionalMoleculesPerComponent_CFCMC(selectedComponent) + subIndex;
    case Move::Types::GroupSwapCBCFCMC:
      return indexOfGroupSwapCBFractionalMoleculesPerComponent_CFCMC(selectedComponent) + subIndex;
    case Move::Types::GibbsSwapCFCMC:
    case Move::Types::GibbsSwapCBCFCMC:
      return indexOfGibbsSwapFractionalMoleculesPerComponent_CFCMC(selectedComponent) + subIndex;
    case Move::Types::GibbsConventionalCFCMC:
    case Move::Types::GibbsConventionalCBCFCMC:
      return indexOfGibbsConventionalFractionalMoleculesPerComponent_CFCMC(selectedComponent) + subIndex;
    default:
      return indexOfGCFractionalMoleculesPerComponent_CFCMC(selectedComponent) + subIndex;
  }
}

std::size_t System::parallelReactionFractionalMoleculeIndex(std::size_t reactionId, std::size_t componentId,
                                                            bool isProduct, std::size_t localIndex) const noexcept
{
  std::size_t index = indexOfParallelReactionFractionalMoleculesPerComponent_CFCMC(componentId);
  for (std::size_t r = 0; r < reactionId; ++r)
  {
    if (!reactions.list[r].isParallelRxCFC())
    {
      continue;
    }
    index += reactions.list[r].reactantStoichiometry[componentId] +
             reactions.list[r].productStoichiometry[componentId];
  }

  const Reaction& reaction = reactions.list[reactionId];
  if (isProduct)
  {
    index += reaction.reactantStoichiometry[componentId] + localIndex;
  }
  else
  {
    index += localIndex;
  }
  return index;
}

std::size_t System::serialReactionFractionalMoleculeIndex(std::size_t reactionId, std::size_t componentId,
                                                          std::size_t localIndex) const noexcept
{
  // the number of physically present serial-reaction fractional molecules depends on which side of
  // each reaction is currently fractional, so the offset must use the active-side stoichiometry
  std::size_t index = indexOfSerialReactionFractionalMoleculesPerComponent_CFCMC(componentId);
  for (std::size_t r = 0; r < reactionId; ++r)
  {
    const Reaction& reaction = reactions.list[r];
    if (!reaction.isSerialRxCFC() || numberOfReactionFractionalMoleculesPerComponent_CFCMC[r][componentId] == 0)
    {
      continue;
    }
    index += reaction.fractionalSideIsReactants ? reaction.reactantStoichiometry[componentId]
                                                : reaction.productStoichiometry[componentId];
  }
  return index + localIndex;
}

void System::precomputeReactionFractionalLayout() noexcept
{
  numberOfParallelReactionFractionalMoleculesPerComponent_CFCMC.assign(components.size(), 0);
  numberOfSerialReactionFractionalMoleculesPerComponent_CFCMC.assign(components.size(), 0);
  numberOfReactionFractionalMoleculesPerComponent_CFCMC.assign(
      reactions.list.size(), std::vector<std::size_t>(components.size(), 0));

  const bool useSerial = usesSerialReactionCFCMC();
  const bool useParallel = usesParallelReactionCFCMC();

  if (!useParallel && !useSerial)
  {
    return;
  }

  for (std::size_t reactionId = 0; reactionId < reactions.list.size(); ++reactionId)
  {
    const Reaction& reaction = reactions.list[reactionId];
    for (std::size_t componentId = 0; componentId < components.size(); ++componentId)
    {
      if (reaction.isSerialRxCFC())
      {
        if (!useSerial)
        {
          continue;
        }
        const std::size_t count = std::max(reaction.reactantStoichiometry[componentId],
                                           reaction.productStoichiometry[componentId]);
        numberOfSerialReactionFractionalMoleculesPerComponent_CFCMC[componentId] += count;
        numberOfReactionFractionalMoleculesPerComponent_CFCMC[reactionId][componentId] = count;
      }
      else if (reaction.isParallelRxCFC())
      {
        if (!useParallel)
        {
          continue;
        }
        const std::size_t count = reaction.reactantStoichiometry[componentId] +
                                  reaction.productStoichiometry[componentId];
        numberOfParallelReactionFractionalMoleculesPerComponent_CFCMC[componentId] += count;
        numberOfReactionFractionalMoleculesPerComponent_CFCMC[reactionId][componentId] = count;
      }
    }
  }
}

void System::syncReactionFractionalMoleculeIndices() noexcept
{
  const bool useSerial = usesSerialReactionCFCMC();
  const bool useParallel = usesParallelReactionCFCMC();
  for (std::size_t reactionId = 0; reactionId < reactions.list.size(); ++reactionId)
  {
    Reaction& reaction = reactions.list[reactionId];
    reaction.reactantFractionalMoleculeIds.assign(components.size(), {});
    reaction.productFractionalMoleculeIds.assign(components.size(), {});

    if (reaction.isSerialRxCFC() && useSerial)
    {
      std::vector<std::vector<std::size_t>>& activeIds =
          reaction.fractionalSideIsReactants ? reaction.reactantFractionalMoleculeIds
                                           : reaction.productFractionalMoleculeIds;
      const std::vector<std::size_t>& stoichiometry =
          reaction.fractionalSideIsReactants ? reaction.reactantStoichiometry : reaction.productStoichiometry;
      for (std::size_t componentId = 0; componentId < components.size(); ++componentId)
      {
        for (std::size_t k = 0; k < stoichiometry[componentId]; ++k)
        {
          activeIds[componentId].push_back(serialReactionFractionalMoleculeIndex(reactionId, componentId, k));
        }
      }
      continue;
    }

    if (!reaction.isParallelRxCFC() || !useParallel)
    {
      continue;
    }

    for (std::size_t componentId = 0; componentId < components.size(); ++componentId)
    {
      for (std::size_t k = 0; k < reaction.reactantStoichiometry[componentId]; ++k)
      {
        reaction.reactantFractionalMoleculeIds[componentId].push_back(
            parallelReactionFractionalMoleculeIndex(reactionId, componentId, false, k));
      }
      for (std::size_t k = 0; k < reaction.productStoichiometry[componentId]; ++k)
      {
        reaction.productFractionalMoleculeIds[componentId].push_back(
            parallelReactionFractionalMoleculeIndex(reactionId, componentId, true, k));
      }
    }
  }
}

bool System::usesReactionConventionalCFCMC() const noexcept
{
  return usesParallelReactionCFCMC() || usesSerialReactionCFCMC();
}

bool System::usesParallelReactionCFCMC() const noexcept
{
  return mc_moves_probabilities.getProbability(Move::Types::ReactionConventionalCFCMC) > 0.0 ||
         mc_moves_probabilities.getProbability(Move::Types::ReactionConventionalCBCFCMC) > 0.0;
}

bool System::usesSerialReactionCFCMC() const noexcept
{
  return mc_moves_probabilities.getProbability(Move::Types::ReactionCFCMC) > 0.0 ||
         mc_moves_probabilities.getProbability(Move::Types::ReactionCBCFCMC) > 0.0;
}

bool System::usesGibbsConventionalCFCMC() const noexcept
{
  for (const Component& component : components)
  {
    if (component.mc_moves_probabilities.getProbability(Move::Types::GibbsConventionalCFCMC) > 0.0 ||
        component.mc_moves_probabilities.getProbability(Move::Types::GibbsConventionalCBCFCMC) > 0.0)
    {
      return true;
    }
  }
  return false;
}

void System::initializeGibbsConventionalFractionalMolecules() noexcept
{
  if (!usesGibbsConventionalCFCMC())
  {
    return;
  }

  for (std::size_t componentId = 0; componentId < components.size(); ++componentId)
  {
    if (numberOfGibbsFractionalMoleculesPerComponent_CFCMC[componentId] == 0)
    {
      continue;
    }

    const std::size_t fractionalMoleculeIndex =
        indexOfGibbsConventionalFractionalMoleculesPerComponent_CFCMC(componentId);
    std::span<Atom> fractionalMolecule = spanOfMolecule(componentId, fractionalMoleculeIndex);
    const std::uint8_t groupId = components[componentId].lambdaGibbs.dUdlambdaGroupId;
    for (Atom& atom : fractionalMolecule)
    {
      atom.setScalingToFractional(0.5, groupId);
    }

    PropertyLambdaProbabilityHistogram& lambdaGibbs = components[componentId].lambdaGibbs;
    if (lambdaGibbs.numberOfSamplePoints > 0)
    {
      lambdaGibbs.setCurrentBin(lambdaGibbs.numberOfSamplePoints / 2);
    }
  }
}

// Fixed-lambda thermodynamic integration: pin the fractional molecule(s) of each component with a
// 'fixedLambdaBin' at lambda = bin * delta. Depending on the pinned coordinate this is the single
// GC fractional molecule, both molecules of an ion-pair, or all members of a CFCMC group. The bin
// never changes during the simulation, so all dU/dlambda samples accumulate in this single bin.
void System::initializeFixedLambdaFractionalMolecules() noexcept
{
  auto rescale = [this](std::size_t componentId, std::size_t slotIndex, double lambda, std::uint8_t groupId)
  {
    for (Atom& atom : spanOfMolecule(componentId, slotIndex))
    {
      atom.setScalingToFractional(lambda, groupId);
    }
  };

  for (std::size_t componentId = 0; componentId < components.size(); ++componentId)
  {
    Component& component = components[componentId];
    if (!component.fixedLambdaBin.has_value())
    {
      continue;
    }

    PropertyLambdaProbabilityHistogram& histogram = component.fixedLambdaHistogram();
    histogram.setCurrentBin(component.fixedLambdaBin.value());
    const double lambda = histogram.lambdaValue();

    switch (component.fixedLambdaCoordinate)
    {
      case Component::FixedLambdaCoordinate::GC:
      {
        rescale(componentId, indexOfGCFractionalMoleculesPerComponent_CFCMC(componentId), lambda,
                histogram.dUdlambdaGroupId);
        break;
      }
      case Component::FixedLambdaCoordinate::PairSwap:
      {
        // both fractional molecules of the ion-pair are coupled to the driving component's histogram
        const std::size_t partner = component.pairComponentId.value();
        rescale(componentId, indexOfPairSwapFractionalMoleculesPerComponent_CFCMC(componentId), lambda,
                histogram.dUdlambdaGroupId);
        rescale(partner, indexOfPairSwapFractionalMoleculesPerComponent_CFCMC(partner), lambda,
                histogram.dUdlambdaGroupId);
        break;
      }
      case Component::FixedLambdaCoordinate::GroupSwap:
      {
        // all fractional molecules of the group (the central component and every satellite
        // occurrence) are coupled to the central component's histogram
        for (std::size_t member = 0; member < components.size(); ++member)
        {
          if (groupSwapLambdaDriver(member, Move::Types::GroupSwapCFCMC) != componentId)
          {
            continue;
          }
          const std::size_t begin = indexOfGroupSwapFractionalMoleculesPerComponent_CFCMC(member);
          for (std::size_t k = 0; k < numberOfGroupSwapFractionalMoleculesPerComponent_CFCMC[member]; ++k)
          {
            rescale(member, begin + k, lambda, histogram.dUdlambdaGroupId);
          }
        }
        break;
      }
    }
  }
}

// Serial Gibbs CFCMC: only the system that owns the active fractional molecule contributes dUdlambda.
// Fractional molecules in inactive systems are switched off (lambda=0) and must have groupId=0,
// otherwise the soft-core derivative still yields a spurious dUdlambda contribution.
void System::initializeGibbsSwapFractionalMoleculeGroupIds() noexcept
{
  for (std::size_t componentId = 0; componentId < components.size(); ++componentId)
  {
    if (numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC[componentId] == 0)
    {
      continue;
    }

    const std::size_t fractionalMoleculeIndex = indexOfGibbsSwapFractionalMoleculesPerComponent_CFCMC(componentId);
    std::span<Atom> fractionalMolecule = spanOfMolecule(componentId, fractionalMoleculeIndex);
    const std::uint8_t groupId =
        containsTheFractionalMolecule ? components[componentId].lambdaGC.dUdlambdaGroupId : std::uint8_t{0};
    for (Atom& atom : fractionalMolecule)
    {
      atom.groupId = groupId;
    }
  }
}

bool System::hasReactionFractionalMolecules() const noexcept
{
  for (const std::vector<std::size_t>& perReaction : numberOfReactionFractionalMoleculesPerComponent_CFCMC)
  {
    for (std::size_t count : perReaction)
    {
      if (count > 0)
      {
        return true;
      }
    }
  }
  return false;
}

PropertyLambdaProbabilityHistogram& System::activeReactionLambdaHistogram(Reaction& reaction) noexcept
{
  if (!reaction.isSerialRxCFC())
  {
    return reaction.lambda;
  }
  return reaction.fractionalSideIsReactants ? reaction.lambda : reaction.lambdaProductSide;
}

const PropertyLambdaProbabilityHistogram& System::activeReactionLambdaHistogram(const Reaction& reaction) const noexcept
{
  if (!reaction.isSerialRxCFC())
  {
    return reaction.lambda;
  }
  return reaction.fractionalSideIsReactants ? reaction.lambda : reaction.lambdaProductSide;
}

// dU/dlambda of a reaction lambda coordinate. Serial Rx/CFC couples the active side directly at
// lambda, so the plain per-group derivative applies. Parallel Rx/CFC couples products at lambda
// and reactants at (1 - lambda); the chain rule d(1 - lambda)/dlambda = -1 makes the reactant
// contribution enter with a minus sign, evaluated at the reactant coupling (1 - lambda).
double System::reactionDUdlambda(const Reaction& reaction) const noexcept
{
  const double lambda = reaction.currentLambda;
  if (reaction.isParallelRxCFC())
  {
    return currentDUdlambda(lambda, reaction.dUdlambdaGroup(false)) -
           currentDUdlambda(1.0 - lambda, reaction.dUdlambdaGroup(true));
  }
  const PropertyLambdaProbabilityHistogram& histogram = activeReactionLambdaHistogram(reaction);
  return currentDUdlambda(lambda, histogram.dUdlambdaGroupId);
}

double System::currentDUdlambda(double lambda, std::size_t groupId) const noexcept
{
  double value = runningEnergies.dudlambda(lambda, groupId);
  if (forceField.computePolarization && groupId != 0)
  {
    value += Scaling::scalingCoulombDerivative(lambda) *
             Interactions::computePolarizationDUdlambda(forceField, simulationBox, spanOfMoleculeAtoms(),
                                                        spanOfMoleculeElectricField(), groupId);
  }
  return value;
}

void System::initializeReactionLambdaHistograms(std::size_t numberOfBlocks, std::size_t numberOfLambdaBins)
{
  for (Reaction& reaction : reactions.list)
  {
    reaction.lambda = PropertyLambdaProbabilityHistogram(numberOfBlocks, numberOfLambdaBins);
    reaction.lambdaProductSide = PropertyLambdaProbabilityHistogram(numberOfBlocks, numberOfLambdaBins);
    reaction.lambda.computeDUdlambda = true;
    reaction.lambdaProductSide.computeDUdlambda = true;
  }
}

void System::syncReactionLambdaBin(Reaction& reaction) noexcept
{
  auto setBinFromLambda = [](PropertyLambdaProbabilityHistogram& histogram, double lambda)
  {
    const std::size_t bin =
        std::min(static_cast<std::size_t>(static_cast<double>(histogram.numberOfSamplePoints) * lambda),
                 histogram.numberOfSamplePoints > 0 ? histogram.numberOfSamplePoints - 1 : 0uz);
    histogram.setCurrentBin(bin);
  };

  if (!reaction.isSerialRxCFC())
  {
    setBinFromLambda(reaction.lambda, reaction.currentLambda);
    return;
  }

  setBinFromLambda(reaction.lambda, reaction.fractionalSideIsReactants ? reaction.currentLambda : 0.0);
  setBinFromLambda(reaction.lambdaProductSide, reaction.fractionalSideIsReactants ? 0.0 : reaction.currentLambda);
}

void System::syncReactionLambdaBins() noexcept
{
  for (Reaction& reaction : reactions.list)
  {
    syncReactionLambdaBin(reaction);
  }
}

void System::reactionLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase phase) noexcept
{
  if (reactions.list.empty() || !usesReactionConventionalCFCMC())
  {
    return;
  }

  const bool hasFractionals = hasReactionFractionalMolecules();

  for (Reaction& reaction : reactions.list)
  {
    if (phase == PropertyLambdaProbabilityHistogram::WangLandauPhase::Sample)
    {
      syncReactionLambdaBin(reaction);
      activeReactionLambdaHistogram(reaction).WangLandauIteration(phase, hasFractionals);
      continue;
    }

    reaction.lambda.WangLandauIteration(phase, hasFractionals);
    if (reaction.isSerialRxCFC())
    {
      reaction.lambdaProductSide.WangLandauIteration(phase, hasFractionals);
    }
  }
}

void System::reactionLambdaSampleOccupancy() noexcept
{
  if (reactions.list.empty() || !usesReactionConventionalCFCMC())
  {
    return;
  }

  const bool hasFractionals = hasReactionFractionalMolecules();
  for (Reaction& reaction : reactions.list)
  {
    if (reaction.isSerialRxCFC())
    {
      // sample both sides so that each histogram's occupancy measures the fraction of samples for
      // which that side held the fractional molecules
      reaction.lambda.sampleOccupancy(hasFractionals && reaction.fractionalSideIsReactants);
      reaction.lambdaProductSide.sampleOccupancy(hasFractionals && !reaction.fractionalSideIsReactants);
      continue;
    }
    activeReactionLambdaHistogram(reaction).sampleOccupancy(hasFractionals);
  }
}

void System::reactionLambdaClearBookkeeping() noexcept
{
  for (Reaction& reaction : reactions.list)
  {
    reaction.lambda.clear();
    if (reaction.isSerialRxCFC())
    {
      reaction.lambdaProductSide.clear();
    }
  }
}

void System::reactionLambdaFinalize() noexcept
{
  if (reactions.list.empty() || !usesReactionConventionalCFCMC())
  {
    return;
  }

  const bool hasFractionals = hasReactionFractionalMolecules();
  for (Reaction& reaction : reactions.list)
  {
    reaction.lambda.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Finalize, hasFractionals);
    if (reaction.isSerialRxCFC())
    {
      reaction.lambdaProductSide.WangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase::Finalize,
                                                     hasFractionals);
    }
  }
}

double System::reactionLambdaMinBias() const noexcept
{
  double minBias = std::numeric_limits<double>::max();
  for (const Reaction& reaction : reactions.list)
  {
    const PropertyLambdaProbabilityHistogram& histogram = activeReactionLambdaHistogram(reaction);
    if (!histogram.biasFactor.empty())
    {
      const double currentMinBias =
          *std::min_element(histogram.biasFactor.cbegin(), histogram.biasFactor.cend());
      minBias = currentMinBias < minBias ? currentMinBias : minBias;
    }
  }
  return minBias;
}

void System::reactionLambdaNormalize(double minBias) noexcept
{
  for (Reaction& reaction : reactions.list)
  {
    reaction.lambda.normalize(minBias);
    if (reaction.isSerialRxCFC())
    {
      reaction.lambdaProductSide.normalize(minBias);
    }
  }
}

// A component "drives" a pair-swap lambda histogram if it is the lower-index component of an
// ion-pair and the corresponding CFCMC pair move is enabled. The partner's histogram stays unused;
// both fractional molecules of the pair are coupled to the driving component's histogram.
bool System::componentDrivesPairSwapLambda(std::size_t componentId, Move::Types move) const noexcept
{
  const Component& component = components[componentId];
  if (!component.pairComponentId.has_value())
  {
    return false;
  }
  const std::size_t partner = component.pairComponentId.value();
  if (partner >= components.size() || componentId >= partner)
  {
    return false;
  }
  // a fixed lambda-bin (thermodynamic integration at constant lambda) drives the conventional
  // pair-swap histogram without any move being enabled
  return component.mc_moves_probabilities.getProbability(move) > 0.0 ||
         (move == Move::Types::PairSwapCFCMC &&
          pinsFixedLambda(component, Component::FixedLambdaCoordinate::PairSwap));
}

// A component "drives" a group-swap lambda histogram if it has a group definition (GroupComponents)
// and the corresponding group CFCMC move is enabled. All fractional molecules of the group (central
// and satellites) are coupled to the driving component's histogram.
bool System::componentDrivesGroupSwapLambda(std::size_t componentId, Move::Types move) const noexcept
{
  const Component& component = components[componentId];
  if (component.groupComponentIds.empty())
  {
    return false;
  }
  // a fixed lambda-bin (thermodynamic integration at constant lambda) drives the conventional
  // group-swap histogram without any move being enabled
  return component.mc_moves_probabilities.getProbability(move) > 0.0 ||
         (move == Move::Types::GroupSwapCFCMC &&
          pinsFixedLambda(component, Component::FixedLambdaCoordinate::GroupSwap));
}

// The driving component of the group fractional slots held by 'componentId': either the component
// itself (when it drives a group CFCMC move) or the first component listing it among its satellites.
// A component can participate in at most one CFCMC group per move flavor.
std::optional<std::size_t> System::groupSwapLambdaDriver(std::size_t componentId, Move::Types move) const noexcept
{
  if (componentDrivesGroupSwapLambda(componentId, move))
  {
    return componentId;
  }
  for (std::size_t i = 0; i < components.size(); ++i)
  {
    if (!componentDrivesGroupSwapLambda(i, move))
    {
      continue;
    }
    if (std::find(components[i].groupComponentIds.begin(), components[i].groupComponentIds.end(), componentId) !=
        components[i].groupComponentIds.end())
    {
      return i;
    }
  }
  return std::nullopt;
}

// Assigns sequential 1-based thermodynamic-integration group ids to all lambda histograms with
// computeDUdlambda enabled. Atoms of the corresponding fractional molecules carry this id in
// Atom::groupId, and RunningEnergy accumulates dU/dlambda separately per group, so up to
// maximumNumberOfDUDlambdaGroups lambda coordinates can be followed simultaneously.
void System::assignDUdlambdaGroups()
{
  std::uint8_t nextGroupId = 1;
  auto assign = [&nextGroupId](PropertyLambdaProbabilityHistogram& histogram)
  {
    histogram.dUdlambdaGroupId = 0;
    if (!histogram.computeDUdlambda) return;
    if (nextGroupId > maximumNumberOfDUDlambdaGroups)
    {
      throw std::runtime_error(
          std::format("Too many thermodynamic integrations: at most {} lambda coordinates can be "
                      "followed simultaneously (Atom::groupId slots)\n",
                      maximumNumberOfDUDlambdaGroups));
    }
    histogram.dUdlambdaGroupId = nextGroupId;
    ++nextGroupId;
  };

  for (std::size_t i = 0; i < components.size(); ++i)
  {
    // only histograms that actually drive a fractional molecule consume a group slot
    const bool hasGCStyleSlot = numberOfGCFractionalMoleculesPerComponent_CFCMC[i] > 0 ||
                                numberOfPairGCFractionalMoleculesPerComponent_CFCMC[i] > 0 ||
                                numberOfGibbsSwapFractionalMoleculesPerComponent_CFCMC[i] > 0;
    if (hasGCStyleSlot)
    {
      assign(components[i].lambdaGC);
    }
    else
    {
      components[i].lambdaGC.dUdlambdaGroupId = 0;
    }

    if (numberOfGibbsFractionalMoleculesPerComponent_CFCMC[i] > 0)
    {
      assign(components[i].lambdaGibbs);
    }
    else
    {
      components[i].lambdaGibbs.dUdlambdaGroupId = 0;
    }

    // the pair-swap histograms live on the driving (lower-index) component of the pair
    if (componentDrivesPairSwapLambda(i, Move::Types::PairSwapCFCMC))
    {
      assign(components[i].lambdaPairSwap);
    }
    else
    {
      components[i].lambdaPairSwap.dUdlambdaGroupId = 0;
    }
    if (componentDrivesPairSwapLambda(i, Move::Types::PairSwapCBCFCMC))
    {
      assign(components[i].lambdaPairSwapCB);
    }
    else
    {
      components[i].lambdaPairSwapCB.dUdlambdaGroupId = 0;
    }

    // the group-swap histograms live on the driving component of the group
    if (componentDrivesGroupSwapLambda(i, Move::Types::GroupSwapCFCMC))
    {
      assign(components[i].lambdaGroupSwap);
    }
    else
    {
      components[i].lambdaGroupSwap.dUdlambdaGroupId = 0;
    }
    if (componentDrivesGroupSwapLambda(i, Move::Types::GroupSwapCBCFCMC))
    {
      assign(components[i].lambdaGroupSwapCB);
    }
    else
    {
      components[i].lambdaGroupSwapCB.dUdlambdaGroupId = 0;
    }
  }

  // Reaction lambda coordinates. Serial Rx/CFC has one fractional side at a time, coupled directly
  // at lambda, so one group shared between both side histograms suffices. Parallel Rx/CFC couples
  // reactants at (1 - lambda) and products at lambda simultaneously; the two opposite chain-rule
  // factors require separate groups (combined at sampling time in reactionDUdlambda()).
  for (Reaction& reaction : reactions.list)
  {
    reaction.lambda.dUdlambdaGroupId = 0;
    reaction.lambdaProductSide.dUdlambdaGroupId = 0;
    if (!reaction.lambda.computeDUdlambda)
    {
      continue;
    }
    if (reaction.isSerialRxCFC())
    {
      assign(reaction.lambda);
      reaction.lambdaProductSide.dUdlambdaGroupId = reaction.lambda.dUdlambdaGroupId;
    }
    else if (reaction.isParallelRxCFC())
    {
      assign(reaction.lambda);
      assign(reaction.lambdaProductSide);
    }
  }
}

// The 1-based dU/dlambda group id carried by the fractional molecule in the given slot (0 when the
// slot is not tracked). The pair-swap slots follow the histogram of the driving (lower-index)
// component of the pair; the GC, pair-GC, and Gibbs-swap slots follow the component's lambdaGC
// histogram; the Gibbs-conventional slot follows lambdaGibbs.
std::uint8_t System::fractionalSlotDUdlambdaGroupId(std::size_t componentId, std::size_t slotIndex) const noexcept
{
  std::size_t driver = componentId;
  if (components[componentId].pairComponentId.has_value() &&
      components[componentId].pairComponentId.value() < components.size())
  {
    driver = std::min(componentId, components[componentId].pairComponentId.value());
  }

  const std::size_t pairSwapBegin = indexOfPairSwapFractionalMoleculesPerComponent_CFCMC(componentId);
  const std::size_t pairSwapEnd = pairSwapBegin + numberOfPairSwapFractionalMoleculesPerComponent_CFCMC[componentId];
  const std::size_t pairSwapCBEnd = pairSwapEnd + numberOfPairSwapCBFractionalMoleculesPerComponent_CFCMC[componentId];

  if (slotIndex >= pairSwapBegin && slotIndex < pairSwapEnd)
  {
    return components[driver].lambdaPairSwap.dUdlambdaGroupId;
  }
  if (slotIndex >= pairSwapEnd && slotIndex < pairSwapCBEnd)
  {
    return components[driver].lambdaPairSwapCB.dUdlambdaGroupId;
  }

  const std::size_t groupSwapBegin = indexOfGroupSwapFractionalMoleculesPerComponent_CFCMC(componentId);
  const std::size_t groupSwapEnd = groupSwapBegin + numberOfGroupSwapFractionalMoleculesPerComponent_CFCMC[componentId];
  const std::size_t groupSwapCBEnd =
      groupSwapEnd + numberOfGroupSwapCBFractionalMoleculesPerComponent_CFCMC[componentId];

  if (slotIndex >= groupSwapBegin && slotIndex < groupSwapEnd)
  {
    const std::optional<std::size_t> groupDriver = groupSwapLambdaDriver(componentId, Move::Types::GroupSwapCFCMC);
    return groupDriver.has_value() ? components[groupDriver.value()].lambdaGroupSwap.dUdlambdaGroupId
                                   : std::uint8_t{0};
  }
  if (slotIndex >= groupSwapEnd && slotIndex < groupSwapCBEnd)
  {
    const std::optional<std::size_t> groupDriver = groupSwapLambdaDriver(componentId, Move::Types::GroupSwapCBCFCMC);
    return groupDriver.has_value() ? components[groupDriver.value()].lambdaGroupSwapCB.dUdlambdaGroupId
                                   : std::uint8_t{0};
  }

  const std::size_t gibbsConventionalBegin =
      indexOfGibbsConventionalFractionalMoleculesPerComponent_CFCMC(componentId);
  const std::size_t gibbsConventionalEnd =
      gibbsConventionalBegin + numberOfGibbsFractionalMoleculesPerComponent_CFCMC[componentId];
  if (slotIndex >= gibbsConventionalBegin && slotIndex < gibbsConventionalEnd)
  {
    return components[componentId].lambdaGibbs.dUdlambdaGroupId;
  }

  return components[componentId].lambdaGC.dUdlambdaGroupId;
}

// The pairSwapLambda* bookkeeping functions cover both the ion-pair and the group CFCMC lambda
// histograms (they share the same lifecycle: Wang-Landau, occupancy, normalization, biasing files).
void System::pairSwapLambdaWangLandauIteration(PropertyLambdaProbabilityHistogram::WangLandauPhase phase) noexcept
{
  for (std::size_t i = 0; i < components.size(); ++i)
  {
    if (componentDrivesPairSwapLambda(i, Move::Types::PairSwapCFCMC))
    {
      components[i].lambdaPairSwap.WangLandauIteration(phase, containsTheFractionalMolecule);
    }
    if (componentDrivesPairSwapLambda(i, Move::Types::PairSwapCBCFCMC))
    {
      components[i].lambdaPairSwapCB.WangLandauIteration(phase, containsTheFractionalMolecule);
    }
    if (componentDrivesGroupSwapLambda(i, Move::Types::GroupSwapCFCMC))
    {
      components[i].lambdaGroupSwap.WangLandauIteration(phase, containsTheFractionalMolecule);
    }
    if (componentDrivesGroupSwapLambda(i, Move::Types::GroupSwapCBCFCMC))
    {
      components[i].lambdaGroupSwapCB.WangLandauIteration(phase, containsTheFractionalMolecule);
    }
  }
}

void System::pairSwapLambdaSampleOccupancy() noexcept
{
  for (std::size_t i = 0; i < components.size(); ++i)
  {
    if (componentDrivesPairSwapLambda(i, Move::Types::PairSwapCFCMC))
    {
      components[i].lambdaPairSwap.sampleOccupancy(containsTheFractionalMolecule);
    }
    if (componentDrivesPairSwapLambda(i, Move::Types::PairSwapCBCFCMC))
    {
      components[i].lambdaPairSwapCB.sampleOccupancy(containsTheFractionalMolecule);
    }
    if (componentDrivesGroupSwapLambda(i, Move::Types::GroupSwapCFCMC))
    {
      components[i].lambdaGroupSwap.sampleOccupancy(containsTheFractionalMolecule);
    }
    if (componentDrivesGroupSwapLambda(i, Move::Types::GroupSwapCBCFCMC))
    {
      components[i].lambdaGroupSwapCB.sampleOccupancy(containsTheFractionalMolecule);
    }
  }
}

void System::pairSwapLambdaClearBookkeeping() noexcept
{
  for (std::size_t i = 0; i < components.size(); ++i)
  {
    if (componentDrivesPairSwapLambda(i, Move::Types::PairSwapCFCMC))
    {
      components[i].lambdaPairSwap.clear();
    }
    if (componentDrivesPairSwapLambda(i, Move::Types::PairSwapCBCFCMC))
    {
      components[i].lambdaPairSwapCB.clear();
    }
    if (componentDrivesGroupSwapLambda(i, Move::Types::GroupSwapCFCMC))
    {
      components[i].lambdaGroupSwap.clear();
    }
    if (componentDrivesGroupSwapLambda(i, Move::Types::GroupSwapCBCFCMC))
    {
      components[i].lambdaGroupSwapCB.clear();
    }
  }
}

double System::pairSwapLambdaMinBias() const noexcept
{
  double minBias = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i < components.size(); ++i)
  {
    if (componentDrivesPairSwapLambda(i, Move::Types::PairSwapCFCMC))
    {
      minBias = std::min(minBias, *std::min_element(components[i].lambdaPairSwap.biasFactor.cbegin(),
                                                    components[i].lambdaPairSwap.biasFactor.cend()));
    }
    if (componentDrivesPairSwapLambda(i, Move::Types::PairSwapCBCFCMC))
    {
      minBias = std::min(minBias, *std::min_element(components[i].lambdaPairSwapCB.biasFactor.cbegin(),
                                                    components[i].lambdaPairSwapCB.biasFactor.cend()));
    }
    if (componentDrivesGroupSwapLambda(i, Move::Types::GroupSwapCFCMC))
    {
      minBias = std::min(minBias, *std::min_element(components[i].lambdaGroupSwap.biasFactor.cbegin(),
                                                    components[i].lambdaGroupSwap.biasFactor.cend()));
    }
    if (componentDrivesGroupSwapLambda(i, Move::Types::GroupSwapCBCFCMC))
    {
      minBias = std::min(minBias, *std::min_element(components[i].lambdaGroupSwapCB.biasFactor.cbegin(),
                                                    components[i].lambdaGroupSwapCB.biasFactor.cend()));
    }
  }
  return minBias;
}

void System::pairSwapLambdaNormalize(double minBias) noexcept
{
  for (std::size_t i = 0; i < components.size(); ++i)
  {
    if (componentDrivesPairSwapLambda(i, Move::Types::PairSwapCFCMC))
    {
      components[i].lambdaPairSwap.normalize(minBias);
    }
    if (componentDrivesPairSwapLambda(i, Move::Types::PairSwapCBCFCMC))
    {
      components[i].lambdaPairSwapCB.normalize(minBias);
    }
    if (componentDrivesGroupSwapLambda(i, Move::Types::GroupSwapCFCMC))
    {
      components[i].lambdaGroupSwap.normalize(minBias);
    }
    if (componentDrivesGroupSwapLambda(i, Move::Types::GroupSwapCBCFCMC))
    {
      components[i].lambdaGroupSwapCB.normalize(minBias);
    }
  }
}

void System::pairSwapLambdaWriteBiasingFiles(std::size_t systemId)
{
  for (std::size_t i = 0; i < components.size(); ++i)
  {
    if (componentDrivesPairSwapLambda(i, Move::Types::PairSwapCFCMC))
    {
      components[i].lambdaPairSwap.writeBiasingFile(
          std::format("bias_factors/lambda_pair_bias_{}.s{}.json", components[i].name, systemId));
    }
    if (componentDrivesPairSwapLambda(i, Move::Types::PairSwapCBCFCMC))
    {
      components[i].lambdaPairSwapCB.writeBiasingFile(
          std::format("bias_factors/lambda_pair_cb_bias_{}.s{}.json", components[i].name, systemId));
    }
    if (componentDrivesGroupSwapLambda(i, Move::Types::GroupSwapCFCMC))
    {
      components[i].lambdaGroupSwap.writeBiasingFile(
          std::format("bias_factors/lambda_group_bias_{}.s{}.json", components[i].name, systemId));
    }
    if (componentDrivesGroupSwapLambda(i, Move::Types::GroupSwapCBCFCMC))
    {
      components[i].lambdaGroupSwapCB.writeBiasingFile(
          std::format("bias_factors/lambda_group_cb_bias_{}.s{}.json", components[i].name, systemId));
    }
  }
}

void System::reactionLambdaSampleProductionHistograms(std::size_t blockIndex, double weight) noexcept
{
  if (reactions.list.empty() || !usesReactionConventionalCFCMC())
  {
    return;
  }

  const bool hasFractionals = hasReactionFractionalMolecules();
  const double totalDensity = static_cast<double>(numberOfIntegerMolecules()) / simulationBox.volume;

  for (Reaction& reaction : reactions.list)
  {
    PropertyLambdaProbabilityHistogram& histogram = activeReactionLambdaHistogram(reaction);
    syncReactionLambdaBin(reaction);
    const double dudlambda = reactionDUdlambda(reaction);
    histogram.sampleHistogram(blockIndex, totalDensity, dudlambda, hasFractionals, weight);
  }
}

void System::createReactionFractionalMolecules()
{
  if (reactions.list.empty())
  {
    return;
  }

  if (!components.empty())
  {
    initializeReactionLambdaHistograms(components.front().lambdaGC.numberOfBlocks,
                                     components.front().lambdaGC.numberOfSamplePoints);
  }

  // the histograms above are fresh objects: (re-)assign the thermodynamic-integration group ids
  // before any fractional molecule is inserted (insertion tags the atoms with these ids)
  assignDUdlambdaGroups();

  const bool useParallel = usesParallelReactionCFCMC();
  const bool useSerial = usesSerialReactionCFCMC();

  precomputeReactionFractionalLayout();

  if (useParallel)
  {
    createParallelReactionFractionalMolecules();
  }
  if (useSerial)
  {
    createSerialReactionFractionalMolecules();
  }

  syncReactionFractionalMoleculeIndices();
  syncReactionLambdaBins();
}

void System::incrementReactionFractionalMoleculeIds(std::size_t componentId) noexcept
{
  for (Reaction& otherReaction : reactions.list)
  {
    if (otherReaction.reactantFractionalMoleculeIds.size() != components.size())
    {
      continue;
    }
    for (std::size_t& moleculeId : otherReaction.reactantFractionalMoleculeIds[componentId])
    {
      ++moleculeId;
    }
    for (std::size_t& moleculeId : otherReaction.productFractionalMoleculeIds[componentId])
    {
      ++moleculeId;
    }
  }
}

void System::createParallelReactionFractionalMolecules()
{
  if (reactions.list.empty())
  {
    return;
  }

  if (!usesParallelReactionCFCMC())
  {
    return;
  }

  RandomNumber random(1200);

  for (std::size_t reactionId = 0; reactionId < reactions.list.size(); ++reactionId)
  {
    Reaction& reaction = reactions.list[reactionId];
    if (!reaction.isParallelRxCFC())
    {
      continue;
    }

    reaction.reactantFractionalMoleculeIds.assign(components.size(), {});
    reaction.productFractionalMoleculeIds.assign(components.size(), {});
    const bool useCBMC = reaction.reactionMove == Move::Types::ReactionConventionalCBCFCMC;

    for (std::size_t componentId = 0; componentId < components.size(); ++componentId)
    {
      // grow at the actual effective coupling (reactants: 1 - lambda, products: lambda) so that the
      // overlap check is meaningful for fractional molecules that start (nearly) fully coupled
      const double reactantScaling = 1.0 - reaction.currentLambda;
      for (std::size_t k = 0; k < reaction.reactantStoichiometry[componentId]; ++k)
      {
        std::optional<ChainGrowData> growData = std::nullopt;
        do
        {
          if (useCBMC)
          {
            growData = CBMC::growMoleculeSwapInsertion(
                random,
                CBMC::GrowContext{hasExternalField, forceField, simulationBox, interpolationGrids,
                                  externalFieldInterpolationGrid, framework, spanOfFrameworkAtoms(),
                                  spanOfMoleculeAtoms(), beta, forceField.cutOffFrameworkVDW,
                                  forceField.cutOffMoleculeVDW, forceField.cutOffCoulomb},
                components[componentId], componentId, numberOfMolecules(), reactantScaling,
                reaction.dUdlambdaGroup(true), true);
          }
          else
          {
            auto [molecule, atoms] = equilibratedIdealGasMoleculeRandomInBox(random, componentId);
            for (Atom& atom : atoms)
            {
              atom.moleculeId = static_cast<std::uint32_t>(numberOfMolecules());
              atom.componentId = static_cast<std::uint8_t>(componentId);
              atom.setScalingToFractional(reactantScaling, reaction.dUdlambdaGroup(true));
            }
            if (std::optional<RunningEnergy> energy =
                    conventionalReactionInitializationEnergy(*this, componentId, atoms))
            {
              growData.emplace(molecule, std::move(atoms), energy.value(), 1.0, 0.0);
            }
          }
        } while (!growData || growData->energies.potentialEnergy() > forceField.energyOverlapCriteria ||
                 insideBlockedPockets(components[componentId], growData->atoms));

        const std::size_t moleculeIndex =
            parallelReactionFractionalMoleculeIndex(reactionId, componentId, false, k);
        insertReactionFractionalMolecule(componentId, moleculeIndex, growData->molecule, growData->atoms, true,
                                         reaction.currentLambda, reaction.dUdlambdaGroup(true));
      }

      const double productScaling = reaction.currentLambda;
      for (std::size_t k = 0; k < reaction.productStoichiometry[componentId]; ++k)
      {
        std::optional<ChainGrowData> growData = std::nullopt;
        do
        {
          if (useCBMC)
          {
            growData = CBMC::growMoleculeSwapInsertion(
                random,
                CBMC::GrowContext{hasExternalField, forceField, simulationBox, interpolationGrids,
                                  externalFieldInterpolationGrid, framework, spanOfFrameworkAtoms(),
                                  spanOfMoleculeAtoms(), beta, forceField.cutOffFrameworkVDW,
                                  forceField.cutOffMoleculeVDW, forceField.cutOffCoulomb},
                components[componentId], componentId, numberOfMolecules(), productScaling,
                reaction.dUdlambdaGroup(false), true);
          }
          else
          {
            auto [molecule, atoms] = equilibratedIdealGasMoleculeRandomInBox(random, componentId);
            for (Atom& atom : atoms)
            {
              atom.moleculeId = static_cast<std::uint32_t>(numberOfMolecules());
              atom.componentId = static_cast<std::uint8_t>(componentId);
              atom.setScalingToFractional(productScaling, reaction.dUdlambdaGroup(false));
            }
            if (std::optional<RunningEnergy> energy =
                    conventionalReactionInitializationEnergy(*this, componentId, atoms))
            {
              growData.emplace(molecule, std::move(atoms), energy.value(), 1.0, 0.0);
            }
          }
        } while (!growData || growData->energies.potentialEnergy() > forceField.energyOverlapCriteria ||
                 insideBlockedPockets(components[componentId], growData->atoms));

        const std::size_t moleculeIndex = parallelReactionFractionalMoleculeIndex(reactionId, componentId, true, k);
        insertReactionFractionalMolecule(componentId, moleculeIndex, growData->molecule, growData->atoms, false,
                                         reaction.currentLambda, reaction.dUdlambdaGroup(false));
      }
    }
  }
}

void System::createSerialReactionFractionalMolecules()
{
  if (reactions.list.empty())
  {
    return;
  }

  if (!usesSerialReactionCFCMC())
  {
    return;
  }

  RandomNumber random(1200);

  for (std::size_t reactionId = 0; reactionId < reactions.list.size(); ++reactionId)
  {
    Reaction& reaction = reactions.list[reactionId];
    if (!reaction.isSerialRxCFC())
    {
      continue;
    }

    reaction.fractionalSideIsReactants = true;
    reaction.currentLambda = 0.0;
    reaction.reactantFractionalMoleculeIds.assign(components.size(), {});
    reaction.productFractionalMoleculeIds.assign(components.size(), {});
    const bool useCBMC = reaction.reactionMove == Move::Types::ReactionCBCFCMC;

    for (std::size_t componentId = 0; componentId < components.size(); ++componentId)
    {
      for (std::size_t k = 0; k < reaction.reactantStoichiometry[componentId]; ++k)
      {
        std::optional<ChainGrowData> growData = std::nullopt;
        do
        {
          // grow at full coupling so that the overlap check below is meaningful: at the actual initial
          // coupling (lambda = 0) the molecule is invisible and could be placed on top of another
          // molecule, which produces an astronomically high energy as soon as lambda is raised
          if (useCBMC)
          {
            growData = CBMC::growMoleculeSwapInsertion(
                random,
                CBMC::GrowContext{hasExternalField, forceField, simulationBox, interpolationGrids,
                                  externalFieldInterpolationGrid, framework, spanOfFrameworkAtoms(),
                                  spanOfMoleculeAtoms(), beta, forceField.cutOffFrameworkVDW,
                                  forceField.cutOffMoleculeVDW, forceField.cutOffCoulomb},
                components[componentId], componentId, numberOfMolecules(), 1.0, reaction.lambda.dUdlambdaGroupId,
                true);
          }
          else
          {
            auto [molecule, atoms] = equilibratedIdealGasMoleculeRandomInBox(random, componentId);
            for (Atom& atom : atoms)
            {
              atom.moleculeId = static_cast<std::uint32_t>(numberOfMolecules());
              atom.componentId = static_cast<std::uint8_t>(componentId);
              atom.setScalingToFractional(1.0, reaction.lambda.dUdlambdaGroupId);
            }
            if (std::optional<RunningEnergy> energy =
                    conventionalReactionInitializationEnergy(*this, componentId, atoms))
            {
              growData.emplace(molecule, std::move(atoms), energy.value(), 1.0, 0.0);
            }
          }
        } while (!growData || growData->energies.potentialEnergy() > forceField.energyOverlapCriteria ||
                 insideBlockedPockets(components[componentId], growData->atoms));

        const std::size_t moleculeIndex = serialReactionFractionalMoleculeIndex(reactionId, componentId, k);
        insertSerialReactionFractionalMolecule(componentId, moleculeIndex, growData->molecule, growData->atoms, 0.0,
                                               reaction.lambda.dUdlambdaGroupId);
      }
    }
  }
}
