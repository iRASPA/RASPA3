module;

module mc_moves_reaction;

import std;

import atom;
import randomnumbers;
import system;
import running_energy;
import reaction;
import interactions_intermolecular;
import interactions_ewald;
import cbmc_chain_data;
import mc_moves_move_types;
import mc_moves_reaction_common;

std::optional<RunningEnergy> MC_Moves::reactionMove_CBMC(RandomNumber& random, System& system)
{
  const Move::Types move = Move::Types::ReactionCBMC;

  if (system.reactions.list.empty())
  {
    return std::nullopt;
  }

  system.mc_moves_statistics.addTrial(move);

  const std::size_t reactionIndex =
      random.uniform_integer(0, static_cast<int>(system.reactions.list.size()) - 1);
  Reaction& reaction = system.reactions.list[reactionIndex];

  const bool forward = random.uniform() < 0.5;
  const std::vector<std::size_t>& removeStoichiometry =
      forward ? reaction.reactantStoichiometry : reaction.productStoichiometry;
  const std::vector<std::size_t>& insertStoichiometry =
      forward ? reaction.productStoichiometry : reaction.reactantStoichiometry;

  if (ReactionCommon::totalStoichiometry(removeStoichiometry) == 0 ||
      ReactionCommon::totalStoichiometry(insertStoichiometry) == 0)
  {
    return std::nullopt;
  }

  std::vector<std::pair<std::size_t, std::size_t>> selectedMolecules;
  if (!ReactionCommon::selectRandomIntegerMolecules(random, system, removeStoichiometry, selectedMolecules))
  {
    return std::nullopt;
  }

  std::vector<std::pair<std::size_t, std::size_t>> excludeMolecules = selectedMolecules;
  ReactionCommon::appendAllReactionFractionalMoleculeExclusions(system, excludeMolecules);
  std::optional<ReactionCommon::MoleculeGroupGrowData> growData =
      ReactionCommon::growMoleculeGroupInsertion(random, system, insertStoichiometry, excludeMolecules);
  if (!growData)
  {
    return std::nullopt;
  }

  std::optional<ReactionCommon::MoleculeGroupRetraceData> retraceData =
      ReactionCommon::retraceMoleculeGroupDeletion(random, system, selectedMolecules);
  if (!retraceData)
  {
    return std::nullopt;
  }

  std::vector<Atom> newAtoms;
  for (const ChainGrowData& data : growData->molecules)
  {
    newAtoms.insert(newAtoms.end(), data.atom.begin(), data.atom.end());
  }

  std::vector<Atom> oldAtoms;
  for (const auto& [componentId, moleculeId] : selectedMolecules)
  {
    std::span<Atom> molecule = system.spanOfMolecule(componentId, moleculeId);
    oldAtoms.insert(oldAtoms.end(), molecule.begin(), molecule.end());
  }

  RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
      system.simulationBox, newAtoms, oldAtoms, system.netCharge);

  RunningEnergy tailEnergyDifference = Interactions::computeInterMolecularTailEnergyDifferenceReaction(
      system.forceField, system.simulationBox, system.totalNumberOfPseudoAtoms, reaction.reactantStoichiometry,
      reaction.productStoichiometry, system.components, forward);

  const ReactionCommon::ReactionMoveKind moveKind =
      forward ? ReactionCommon::ReactionMoveKind::ForwardInsert
              : ReactionCommon::ReactionMoveKind::BackwardDelete;
  const double equilibriumTerm = ReactionCommon::computeReactionEquilibriumLogTerm(system, reaction, moveKind);

  const double correctionFactorEwald = std::exp(-system.beta * energyFourierDifference.potentialEnergy());
  const double rosenbluthNew =
      growData->RosenbluthWeight * std::exp(-system.beta * tailEnergyDifference.potentialEnergy());
  const double rosenbluthOld = retraceData->RosenbluthWeight;

  const double acceptanceProbability = correctionFactorEwald * (rosenbluthNew / rosenbluthOld) * std::exp(equilibriumTerm);

  if (random.uniform() >= acceptanceProbability)
  {
    return std::nullopt;
  }

  system.mc_moves_statistics.addConstructed(move);
  system.mc_moves_statistics.addAccepted(move);

  Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);

  ReactionCommon::deleteSelectedMolecules(system, selectedMolecules);
  ReactionCommon::insertGrownMolecules(system, growData->molecules, insertStoichiometry);

  return (growData->energies - retraceData->energies) + energyFourierDifference + tailEnergyDifference;
}
