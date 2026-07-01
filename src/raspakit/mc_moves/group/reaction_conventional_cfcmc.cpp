module;

module mc_moves_reaction_conventional_cfcmc;

import std;

import component;
import atom;
import molecule;
import double3;
import simulationbox;
import cbmc;
import cbmc_chain_data;
import randomnumbers;
import forcefield;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import property_lambda_probability_histogram;
import interactions_external_field;
import reaction;
import mc_moves_move_types;
import running_energy;
import system;
import mc_moves_reaction_common;

std::optional<RunningEnergy> MC_Moves::reactionMove_ConventionalCFCMC(RandomNumber& random, System& system)
{
  const Move::Types move = Move::Types::ReactionConventionalCFCMC;

  if (system.reactions.list.empty())
  {
    return std::nullopt;
  }

  system.mc_moves_statistics.addTrial(move);

  std::vector<Reaction*> parallelReactions;
  parallelReactions.reserve(system.reactions.list.size());
  for (Reaction& reaction : system.reactions.list)
  {
    if (!reaction.serialRxCFC && reaction.reactantFractionalMoleculeIds.size() == system.components.size())
    {
      parallelReactions.push_back(&reaction);
    }
  }
  if (parallelReactions.empty())
  {
    return std::nullopt;
  }

  Reaction& reaction = *parallelReactions[static_cast<std::size_t>(
      random.uniform_integer(0, static_cast<int>(parallelReactions.size()) - 1))];

  PropertyLambdaProbabilityHistogram& lambdaHistogram = reaction.lambda;
  const double lambdaOld = reaction.currentLambda;
  const std::size_t oldBin = std::min(
      static_cast<std::size_t>(lambdaHistogram.numberOfSamplePoints * lambdaOld),
      lambdaHistogram.numberOfSamplePoints > 0 ? lambdaHistogram.numberOfSamplePoints - 1 : 0);
  const double biasOld = lambdaHistogram.biasFactor[oldBin];

  const double vNew = (2.0 * random.uniform() - 1.0) * reaction.maximumLambdaChange;
  double lambdaNew = lambdaOld + vNew;

  ReactionCommon::ReactionMoveKind moveKind = ReactionCommon::ReactionMoveKind::LambdaChange;
  if (lambdaNew > 1.0)
  {
    moveKind = ReactionCommon::ReactionMoveKind::ForwardInsert;
    lambdaNew -= 1.0;
  }
  else if (lambdaNew < 0.0)
  {
    moveKind = ReactionCommon::ReactionMoveKind::BackwardDelete;
    lambdaNew += 1.0;
  }

  const std::size_t newBin = std::min(
      static_cast<std::size_t>(lambdaHistogram.numberOfSamplePoints * lambdaNew),
      lambdaHistogram.numberOfSamplePoints > 0 ? lambdaHistogram.numberOfSamplePoints - 1 : 0);
  const double biasNew = lambdaHistogram.biasFactor[newBin];

  RunningEnergy energyDifference;
  RunningEnergy energyFourierDifference{};
  double rosenbluthRatio = 1.0;
  std::vector<std::pair<std::size_t, std::size_t>> selectedMolecules;
  std::optional<ReactionCommon::MoleculeGroupGrowData> growData;
  std::vector<Atom> newAtoms;
  std::vector<Atom> oldAtoms;

  if (moveKind == ReactionCommon::ReactionMoveKind::LambdaChange)
  {
    std::optional<RunningEnergy> scalingDifference =
        ReactionCommon::computeReactionFractionalScalingEnergyDifference(system, reaction, lambdaOld, lambdaNew);
    if (!scalingDifference)
    {
      return std::nullopt;
    }
    energyDifference = scalingDifference.value();
  }
  else
  {
    // Boundary crossing: the fractional molecules are rescaled from lambdaOld to the wrapped lambdaNew
    // (forward: reactants come back on, products go nearly off; backward: the reverse) and simultaneously
    // a group of integer molecules is swapped.
    const bool forward = (moveKind == ReactionCommon::ReactionMoveKind::ForwardInsert);
    const std::vector<std::size_t>& deleteStoichiometry =
        forward ? reaction.reactantStoichiometry : reaction.productStoichiometry;
    const std::vector<std::size_t>& insertStoichiometry =
        forward ? reaction.productStoichiometry : reaction.reactantStoichiometry;

    if (!ReactionCommon::selectRandomIntegerMolecules(random, system, deleteStoichiometry, selectedMolecules))
    {
      return std::nullopt;
    }

    // non-Ewald part of rescaling the fractionals to lambdaNew; the Ewald Fourier part is computed
    // below in a single combined difference together with the integer-molecule swap
    std::vector<Atom> fractionalOldAtoms =
        ReactionCommon::collectReactionFractionalAtomsAtLambda(system, reaction, lambdaOld);
    std::optional<RunningEnergy> scalingDifference = ReactionCommon::computeReactionFractionalScalingEnergyDifference(
        system, reaction, lambdaOld, lambdaNew, false);
    if (!scalingDifference)
    {
      return std::nullopt;
    }
    energyDifference = scalingDifference.value();

    // apply the new scaling in memory so that the CBMC grow and retrace see the correct (final) fractional
    // background; otherwise the energy bookkeeping misses the cross-terms between the swapped molecules and
    // the fractional molecules, and grown molecules can be placed on top of fractional molecules
    ReactionCommon::setReactionFractionalScaling(system, reaction, lambdaNew);

    growData = ReactionCommon::growMoleculeGroupInsertion(random, system, insertStoichiometry, selectedMolecules);
    if (!growData)
    {
      ReactionCommon::setReactionFractionalScaling(system, reaction, lambdaOld);
      return std::nullopt;
    }

    std::optional<ReactionCommon::MoleculeGroupRetraceData> retraceData =
        ReactionCommon::retraceMoleculeGroupDeletion(random, system, selectedMolecules);
    if (!retraceData)
    {
      ReactionCommon::setReactionFractionalScaling(system, reaction, lambdaOld);
      return std::nullopt;
    }

    for (const ChainGrowData& data : growData->molecules)
    {
      newAtoms.insert(newAtoms.end(), data.atom.begin(), data.atom.end());
    }
    for (const auto& [componentId, moleculeId] : selectedMolecules)
    {
      std::span<Atom> molecule = system.spanOfMolecule(componentId, moleculeId);
      oldAtoms.insert(oldAtoms.end(), molecule.begin(), molecule.end());
    }

    // single Ewald Fourier difference for the combined change: fractional rescaling + molecule swap
    std::vector<Atom> ewaldNewAtoms = newAtoms;
    std::vector<Atom> fractionalNewAtoms =
        ReactionCommon::collectReactionFractionalAtomsAtLambda(system, reaction, lambdaNew);
    ewaldNewAtoms.insert(ewaldNewAtoms.end(), fractionalNewAtoms.begin(), fractionalNewAtoms.end());
    std::vector<Atom> ewaldOldAtoms = oldAtoms;
    ewaldOldAtoms.insert(ewaldOldAtoms.end(), fractionalOldAtoms.begin(), fractionalOldAtoms.end());

    energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
        system.simulationBox, ewaldNewAtoms, ewaldOldAtoms);
    energyDifference += (growData->energies - retraceData->energies) + energyFourierDifference;
  }

  RunningEnergy tailEnergyDifference{};
  if (moveKind != ReactionCommon::ReactionMoveKind::LambdaChange)
  {
    // scaling-aware tail difference of the molecule swap (the tail change of the fractional rescaling is
    // already contained in the scaling energy difference above)
    tailEnergyDifference = ReactionCommon::computeGroupSwapTailEnergyDifference(system, newAtoms, oldAtoms);
    energyDifference += tailEnergyDifference;
  }

  const double equilibriumTerm = moveKind == ReactionCommon::ReactionMoveKind::LambdaChange
                                     ? 0.0
                                     : ReactionCommon::computeReactionEquilibriumLogTerm(system, reaction, moveKind);

  double acceptanceProbability = 0.0;
  if (moveKind == ReactionCommon::ReactionMoveKind::LambdaChange)
  {
    acceptanceProbability = std::exp(-system.beta * energyDifference.potentialEnergy() + biasNew - biasOld);
  }
  else
  {
    acceptanceProbability =
        rosenbluthRatio * std::exp(equilibriumTerm + biasNew - biasOld) *
        std::exp(-system.beta * energyDifference.potentialEnergy());
  }

  if (random.uniform() >= acceptanceProbability)
  {
    if (moveKind != ReactionCommon::ReactionMoveKind::LambdaChange)
    {
      ReactionCommon::setReactionFractionalScaling(system, reaction, lambdaOld);
    }
    return std::nullopt;
  }

  system.mc_moves_statistics.addConstructed(move);
  system.mc_moves_statistics.addAccepted(move);

  Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);

  if (moveKind == ReactionCommon::ReactionMoveKind::ForwardInsert)
  {
    ReactionCommon::acceptReactionForwardInsert(system, reaction, lambdaNew, selectedMolecules, growData->molecules);
  }
  else if (moveKind == ReactionCommon::ReactionMoveKind::BackwardDelete)
  {
    ReactionCommon::acceptReactionBackwardDelete(system, reaction, lambdaNew, selectedMolecules, growData->molecules);
  }
  ReactionCommon::setReactionFractionalScaling(system, reaction, lambdaNew);

  reaction.currentLambda = lambdaNew;
  lambdaHistogram.setCurrentBin(newBin);

  return energyDifference;
}
