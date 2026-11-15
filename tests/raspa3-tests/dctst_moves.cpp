#include <gtest/gtest.h>

import std;

import forcefield;
import simulationbox;
import component;
import system;
import mc_moves_move_types;
import mc_moves_probabilities;
import dctst_moves;

// The moves of the molecules around the tagged one in the dcTST drivers: every single-molecule move is performed
// under its own name, and every other move is refused by the input reader and the constructors alike.

namespace
{

System makeSystem(const std::vector<std::pair<Move::Types, double>>& moves)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component co2 = Component::makeCO2(forceField, 0, true);
  for (std::size_t index = 0; index < std::to_underlying(Move::Types::Count); ++index)
  {
    co2.mc_moves_probabilities.setProbability(static_cast<Move::Types>(index), 0.0);
  }
  for (const auto& [move, probability] : moves) co2.mc_moves_probabilities.setProbability(move, probability);
  return System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {co2}, {}, {4}, 5);
}

std::string refusal(const System& system, double reinsertion = 0.5, bool occupancies = false)
{
  try
  {
    DcTSTMoves::validate(system, 0, "DynamicallyCorrectedTST", reinsertion, occupancies);
  }
  catch (const std::runtime_error& error)
  {
    return error.what();
  }
  return {};
}

}  // namespace

TEST(DCTST_MOVES, the_conformational_moves_are_performed_and_everything_else_is_refused)
{
  for (const Move::Types move :
       {Move::Types::Translation, Move::Types::RandomTranslation, Move::Types::Rotation, Move::Types::RandomRotation,
        Move::Types::ReinsertionCBMC, Move::Types::ReinsertionMultiTrialCBMC, Move::Types::PartialReinsertionCBMC,
        Move::Types::Pivot, Move::Types::PivotCBMC, Move::Types::Crankshaft, Move::Types::Reptation,
        Move::Types::ConcertedRotation, Move::Types::BeadDisplacement, Move::Types::BeadFlip})
  {
    EXPECT_TRUE(DcTSTMoves::supported(move)) << Move::moveNames[std::to_underlying(move)];
    EXPECT_TRUE(refusal(makeSystem({{Move::Types::Translation, 1.0}, {move, 1.0}})).empty())
        << Move::moveNames[std::to_underlying(move)];
  }

  // particle-number, identity, fractional-molecule, test-particle and connectivity moves
  for (const Move::Types move :
       {Move::Types::Swap, Move::Types::SwapCBMC, Move::Types::SwapMultiTrialCBMC, Move::Types::SwapCFCMC,
        Move::Types::SwapCBCFCMC, Move::Types::PairSwap, Move::Types::PairSwapCBMC, Move::Types::PairSwapCFCMC,
        Move::Types::GroupSwap, Move::Types::GroupSwapCFCMC, Move::Types::IdentityChangeCBMC, Move::Types::Widom,
        Move::Types::WidomCFCMC, Move::Types::DoubleBridging, Move::Types::IntramolecularDoubleRebridging,
        Move::Types::TetheredProtonHop, Move::Types::GibbsSwapCBMC})
  {
    EXPECT_FALSE(DcTSTMoves::supported(move)) << Move::moveNames[std::to_underlying(move)];
    const std::string message = refusal(makeSystem({{Move::Types::Translation, 1.0}, {move, 1.0}}));
    EXPECT_NE(message.find(Move::moveNames[std::to_underlying(move)]), std::string::npos)
        << Move::moveNames[std::to_underlying(move)] << ": '" << message << "'";
  }
}

TEST(DCTST_MOVES, the_occupancy_constraint_refuses_the_move_its_undo_does_not_cover)
{
  const System system = makeSystem({{Move::Types::Translation, 1.0}, {Move::Types::Reptation, 1.0}});
  EXPECT_TRUE(refusal(system, 0.5, false).empty());
  EXPECT_NE(refusal(system, 0.5, true).find("CageOccupancies"), std::string::npos);
}

TEST(DCTST_MOVES, a_polyatomic_tagged_molecule_needs_its_reinsertion)
{
  const System system = makeSystem({{Move::Types::Translation, 1.0}});
  EXPECT_TRUE(refusal(system, 0.3).empty());
  EXPECT_NE(refusal(system, 0.0).find("reinsertion probability"), std::string::npos);
}
