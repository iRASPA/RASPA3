module;

module dctst_moves;

import std;

import system;
import component;
import randomnumbers;
import running_energy;
import mc_moves_move_types;
import mc_moves_translation;
import mc_moves_random_translation;
import mc_moves_rotation;
import mc_moves_random_rotation;
import mc_moves_reinsertion;
import mc_moves_reinsertion_multi_trial_cbmc;
import mc_moves_partial_reinsertion;
import mc_moves_pivot;
import mc_moves_pivot_cbmc;
import mc_moves_crankshaft;
import mc_moves_reptation;
import mc_moves_concerted_rotation;
import mc_moves_bead_displacement;
import mc_moves_bead_flip;

bool DcTSTMoves::supported(Move::Types move)
{
  switch (move)
  {
    case Move::Types::Translation:
    case Move::Types::RandomTranslation:
    case Move::Types::Rotation:
    case Move::Types::RandomRotation:
    case Move::Types::ReinsertionCBMC:
    case Move::Types::ReinsertionMultiTrialCBMC:
    case Move::Types::PartialReinsertionCBMC:
    case Move::Types::Pivot:
    case Move::Types::PivotCBMC:
    case Move::Types::Crankshaft:
    case Move::Types::Reptation:
    case Move::Types::ConcertedRotation:
    case Move::Types::BeadDisplacement:
    case Move::Types::BeadFlip:
      return true;
    default:
      return false;
  }
}

std::optional<RunningEnergy> DcTSTMoves::perform(RandomNumber& random, System& system, std::size_t selectedComponent,
                                                 std::size_t selectedMolecule, Move::Types move)
{
  switch (move)
  {
    case Move::Types::Translation:
      return MC_Moves::translationMove(random, system, selectedComponent, selectedMolecule);
    case Move::Types::RandomTranslation:
      return MC_Moves::randomTranslationMove(random, system, selectedComponent, selectedMolecule);
    case Move::Types::Rotation:
      return MC_Moves::rotationMove(random, system, selectedComponent, selectedMolecule);
    case Move::Types::RandomRotation:
      return MC_Moves::randomRotationMove(random, system, selectedComponent, selectedMolecule);
    case Move::Types::ReinsertionCBMC:
      return MC_Moves::reinsertionMove(random, system, selectedComponent, selectedMolecule);
    case Move::Types::ReinsertionMultiTrialCBMC:
      return MC_Moves::reinsertionMoveMultiTrialCBMC(random, system, selectedComponent, selectedMolecule);
    case Move::Types::PartialReinsertionCBMC:
      return MC_Moves::partialReinsertionMove(random, system, selectedComponent, selectedMolecule);
    case Move::Types::Pivot:
      return MC_Moves::pivotMove(random, system, selectedComponent, selectedMolecule);
    case Move::Types::PivotCBMC:
      return MC_Moves::pivotCBMCMove(random, system, selectedComponent, selectedMolecule);
    case Move::Types::Crankshaft:
      return MC_Moves::crankshaftMove(random, system, selectedComponent, selectedMolecule);
    case Move::Types::Reptation:
      return MC_Moves::reptationMove(random, system, selectedComponent, selectedMolecule);
    case Move::Types::ConcertedRotation:
      return MC_Moves::concertedRotationMove(random, system, selectedComponent, selectedMolecule);
    case Move::Types::BeadDisplacement:
      return MC_Moves::beadDisplacementMove(random, system, selectedComponent, selectedMolecule);
    case Move::Types::BeadFlip:
      return MC_Moves::beadFlipMove(random, system, selectedComponent, selectedMolecule);
    default:
      break;
  }
  throw std::logic_error(std::format("dcTST: move '{}' is not a single-molecule move of the constrained ensemble\n",
                                     Move::moveNames[std::to_underlying(move)]));
}

void DcTSTMoves::validate(const System& system, std::size_t taggedComponent, std::string_view driverName,
                          double taggedReinsertionProbability, bool occupancyConstraint)
{
  for (const Component& component : system.components)
  {
    for (std::size_t index = 0; index < std::to_underlying(Move::Types::Count); ++index)
    {
      const Move::Types move = static_cast<Move::Types>(index);
      if (!(component.mc_moves_probabilities.getProbability(move) > 0.0)) continue;
      if (!supported(move))
      {
        throw std::runtime_error(std::format(
            "dcTST: '{}' is a fixed-loading (NVT) method with a constraint on one molecule and performs the "
            "single-molecule moves only (translation, rotation, reinsertion, partial reinsertion, pivot, crankshaft, "
            "reptation, concerted rotation, bead displacement and bead flip); remove move '{}' of component '{}'\n",
            driverName, Move::moveNames[index], component.name));
      }
      // the occupancy constraint undoes an accepted move that changes the occupancy of cage A or B by putting
      // the atoms back; a reptation also relabels its atoms, which that undo does not reverse
      if (occupancyConstraint && move == Move::Types::Reptation)
      {
        throw std::runtime_error(std::format(
            "dcTST: '{}' with 'CageOccupancies' does not support move '{}' (component '{}')\n", driverName,
            Move::moveNames[index], component.name));
      }
    }
  }
  // the system moves (volume, Gibbs volume, hybrid MC, reactions, cross-link moves) as well: run in the pilot stage
  // they rescaled the box behind the coordinate grid, moved all molecules at once or changed the molecule counts
  for (std::size_t index = 0; index < std::to_underlying(Move::Types::Count); ++index)
  {
    const Move::Types move = static_cast<Move::Types>(index);
    if (!(system.mc_moves_probabilities.getProbability(move) > 0.0) || supported(move)) continue;
    throw std::runtime_error(
        std::format("dcTST: '{}' is a fixed-loading (NVT) method with a constraint on one molecule and performs the "
                    "single-molecule moves only; remove the system move '{}'\n",
                    driverName, Move::moveNames[index]));
  }
  if (taggedComponent < system.components.size())
  {
    const Component& tagged = system.components[taggedComponent];
    if (tagged.atoms.size() > 1uz && !(taggedReinsertionProbability > 0.0))
    {
      throw std::runtime_error(std::format(
          "dcTST: '{}': the reinsertion probability of the tagged molecule ('ReinsertionInPlaneProbability', "
          "'ReinsertionTaggedProbability') must be positive for the polyatomic tagged component '{}': the tagged "
          "molecule is otherwise only translated, and its orientation and conformation never change\n",
          driverName, tagged.name));
    }
  }
}
