#include <gtest/gtest.h>

import std;

import atom;
import forcefield;
import component;
import system;
import simulationbox;
import running_energy;
import randomnumbers;
import move_statistics;
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_gibbs_swap_cbmc;
import mc_moves_gibbs_swap_multi_trial_cbmc;

// A Gibbs transfer out of an empty box is an attempt: the transition matrix records it as a rejected transfer,
// and the move statistics have to count it too (as the grand-canonical deletion counts an attempt on an empty
// component), or the two disagree on the number of attempts.

namespace
{

System makeBox(std::size_t molecules, Move::Types move)
{
  const ForceField forceField({{"CH4", false, 16.04246, 0.0, 0.0, 6, false}}, {{158.5, 3.72}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  MCMoveProbabilities probabilities;
  probabilities.setProbability(move, 1.0);
  const Component methane(forceField, "methane", 190.564, 4599200, 0.01142,
                          {Atom({0, 0, 0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21, probabilities,
                          std::nullopt, false);
  return System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {methane}, {}, {molecules}, 5);
}

}  // namespace

TEST(GIBBS_EMPTY_BOX_STATISTICS, a_transfer_out_of_an_empty_box_is_counted_as_a_trial)
{
  for (const Move::Types move : {Move::Types::GibbsSwapCBMC, Move::Types::GibbsSwapMultiTrialCBMC})
  {
    System systemA = makeBox(3, move);
    System systemB = makeBox(0, move);  // nothing to transfer out of B
    RandomNumber random(7);
    for (std::size_t attempt = 0; attempt < 5; ++attempt)
    {
      const auto result = (move == Move::Types::GibbsSwapCBMC)
                              ? MC_Moves::GibbsSwapMove_CBMC(random, systemA, systemB, 0)
                              : MC_Moves::GibbsSwapMoveMultiTrialCBMC(random, systemA, systemB, 0);
      EXPECT_FALSE(result.has_value());
    }
    EXPECT_DOUBLE_EQ(std::get<MoveStatistics<double>>(systemA.components[0].mc_moves_statistics[move]).counts, 5.0)
        << Move::moveNames[std::to_underlying(move)];
    EXPECT_DOUBLE_EQ(std::get<MoveStatistics<double>>(systemB.components[0].mc_moves_statistics[move]).counts, 5.0)
        << Move::moveNames[std::to_underlying(move)];
  }
}
