module;

module mc_moves_reaction;

import std;

import atom;
import randomnumbers;
import system;
import running_energy;
import reaction;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import cbmc_results;
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

  std::vector<Reaction*> matchingReactions;
  matchingReactions.reserve(system.reactions.list.size());
  for (Reaction& candidate : system.reactions.list)
  {
    if (candidate.reactionMove == move)
    {
      matchingReactions.push_back(&candidate);
    }
  }
  if (matchingReactions.empty())
  {
    return std::nullopt;
  }
  Reaction& reaction = *matchingReactions[random.uniform_integer(0uz, matchingReactions.size() - 1uz)];

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

  // The new molecules are grown in the system without the molecules they replace. Everything else, including the
  // fractional molecules of CFCMC reactions (coupled with their lambda), is part of the Hamiltonian and stays in
  // the background of the growth and of the retrace; the reverse move grows in the same environment.
  std::optional<ReactionCommon::MoleculeGroupGrowData> growData =
      ReactionCommon::growMoleculeGroupInsertion(random, system, insertStoichiometry, selectedMolecules);
  if (!growData)
  {
    return std::nullopt;
  }

  std::optional<ReactionCommon::MoleculeGroupRetraceData> retraceData =
      ReactionCommon::retraceMoleculeGroupDeletion(random, system, selectedMolecules, {});
  if (!retraceData)
  {
    return std::nullopt;
  }

  std::vector<Atom> newAtoms;
  for (const CBMC::GrowResult& data : growData->molecules)
  {
    newAtoms.insert(newAtoms.end(), data.atoms.begin(), data.atoms.end());
  }

  std::vector<Atom> oldAtoms;
  for (const auto& [componentId, moleculeId] : selectedMolecules)
  {
    std::span<Atom> molecule = system.spanOfMolecule(componentId, moleculeId);
    oldAtoms.insert(oldAtoms.end(), molecule.begin(), molecule.end());
  }

  RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
      system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.trialEik, system.forceField,
      system.simulationBox, system.components, newAtoms, oldAtoms, system.netCharge);

  // Tail corrections of the molecule-molecule and of the framework-molecule interactions: both change with the
  // pseudo-atom composition of the adsorbate phase. Pairwise with the scalings of the present molecules, so that
  // the coupling of coexisting fractional molecules (and their dU/dlambda) is honoured.
  RunningEnergy tailEnergyDifference =
      ReactionCommon::computeGroupSwapTailEnergyDifference(system, newAtoms, oldAtoms);

  const ReactionCommon::ReactionMoveKind moveKind =
      forward ? ReactionCommon::ReactionMoveKind::ForwardInsert
              : ReactionCommon::ReactionMoveKind::BackwardDelete;
  const double equilibriumTerm = ReactionCommon::computeReactionEquilibriumLogTerm(system, reaction, moveKind);

  // Ewald energy difference in log space (its Boltzmann factor overflows or underflows on its own for
  // charged molecules)
  const double logCorrectionFactorEwald = -system.beta * energyFourierDifference.potentialEnergy();
  // grow and retrace referenced to the ideal-gas Rosenbluth weights, so that the full ideal-gas
  // partition functions q_i can be used in the equilibrium term (Rosch and Maginn, eqs. 18-24)
  const double idealGasInsert = ReactionCommon::idealGasRosenbluthWeightProduct(system, insertStoichiometry);
  const double idealGasRemove = ReactionCommon::idealGasRosenbluthWeightProduct(system, removeStoichiometry);
  // Rosenbluth weights through their exact logarithms (raw weights of long chains underflow to zero).
  const double acceptanceProbability =
      std::exp(logCorrectionFactorEwald + growData->logRosenbluthWeight - std::log(idealGasInsert) -
               system.beta * tailEnergyDifference.potentialEnergy() - retraceData->logRosenbluthWeight +
               std::log(idealGasRemove) + equilibriumTerm);

  if (random.uniform() >= acceptanceProbability)
  {
    return std::nullopt;
  }

  system.mc_moves_statistics.addConstructed(move);
  system.mc_moves_statistics.addAccepted(move);

  Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.trialEik);

  ReactionCommon::deleteSelectedMolecules(system, selectedMolecules);
  ReactionCommon::insertGrownMolecules(system, growData->molecules, insertStoichiometry);

  return (growData->energies - retraceData->energies) + energyFourierDifference + tailEnergyDifference;
}
