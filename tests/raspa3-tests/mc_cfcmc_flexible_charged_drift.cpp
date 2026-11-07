#include <gtest/gtest.h>

#include "../test_support.hpp"
#include "molecule_fixtures.hpp"

import std;

import atom;
import forcefield;
import component;
import system;
import monte_carlo;
import simulationbox;
import running_energy;
import mc_moves_move_types;
import mc_moves_probabilities;

// CFCMC moves rescale the charges of fractional molecules in place. For a flexible molecule with charged
// non-excluded intramolecular pairs (the 1-4 pair of butane here) the Ewald Fourier sum counts such a pair with
// the scaled charges, and the intramolecular Coulomb term carries the compensating, equally scaled, kernel: the
// two change by opposite amounts when the scaling changes, and a move that updates only the Fourier part drifts
// from a full recomputation in both terms.

namespace
{

ForceField makeChargedAlkaneForceField()
{
  return ForceField({{"CH3", false, 15.04, 0.25, 0.0, 6, false}, {"CH2", false, 14.03, -0.25, 0.0, 6, false}},
                    {{98.0, 3.75}, {46.0, 3.95}}, ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true,
                    false, true);
}

Component makeChargedButane(const ForceField& forceField, const MCMoveProbabilities& probabilities)
{
  TemporaryFile file("butane.json", molecule_fixtures::kButaneJson);
  Component butane(Component::Type::Adsorbate, 0, forceField, "butane", file.stemPath().string(), 5, 21, probabilities,
                   std::nullopt, false);
  butane.lnPartitionFunction = 1.0;
  return butane;
}

void runAndCheckDrift(const MCMoveProbabilities& probabilities)
{
  const ForceField forceField = makeChargedAlkaneForceField();
  std::vector<System> systems{System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {},
                                     {makeChargedButane(forceField, probabilities)}, {}, {6}, 5)};
  ASSERT_TRUE(systems[0].forceField.usesEwaldFourier());
  ASSERT_NE(systems[0].computeTotalEnergies().intraCoul, 0.0);  // the charged 1-4 pairs are exercised

  MonteCarlo mc({20, 0, 3, 3, 1000, 10000, 5000, 5000}, systems, 17uz, 5, false);
  mc.run();

  for (System& s : mc.systems)
  {
    const RunningEnergy recomputed = s.computeTotalEnergies();
    const RunningEnergy drift = s.runningEnergies - recomputed;
    EXPECT_NEAR(drift.intraCoul, 0.0, 1e-6);
    EXPECT_NEAR(drift.ewald_fourier, 0.0, 1e-6);
    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
  }
}

}  // namespace

TEST(MC_CFCMC_FLEXIBLE_CHARGED_DRIFT, swap_cfcmc)
{
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::SwapCFCMC, 1.0);
  runAndCheckDrift(probabilities);
}

TEST(MC_CFCMC_FLEXIBLE_CHARGED_DRIFT, swap_cbcfcmc)
{
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::SwapCBCFCMC, 1.0);
  runAndCheckDrift(probabilities);
}

TEST(MC_CFCMC_FLEXIBLE_CHARGED_DRIFT, pair_swap_cfcmc_on_two_charged_butanes)
{
  const ForceField forceField = makeChargedAlkaneForceField();
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::PairSwapCFCMC, 1.0);
  Component a = makeChargedButane(forceField, probabilities);
  Component b = makeChargedButane(forceField, probabilities);
  b.name = "butane-b";
  a.pairComponentId = 1;
  a.maximumPairDistance = 6.0;
  std::vector<System> systems{
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {a, b}, {}, {4, 4}, 5)};
  MonteCarlo mc({20, 0, 3, 3, 1000, 10000, 5000, 5000}, systems, 19uz, 5, false);
  mc.run();
  for (System& s : mc.systems)
  {
    const RunningEnergy drift = s.runningEnergies - s.computeTotalEnergies();
    EXPECT_NEAR(drift.intraCoul, 0.0, 1e-6);
    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
  }
}

TEST(MC_CFCMC_FLEXIBLE_CHARGED_DRIFT, group_swap_cfcmc)
{
  const ForceField forceField = makeChargedAlkaneForceField();
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::GroupSwapCFCMC, 1.0);
  Component a = makeChargedButane(forceField, probabilities);
  Component b = makeChargedButane(forceField, MCMoveProbabilities{});
  b.name = "butane-b";
  a.groupComponentIds = {1};
  a.maximumGroupDistance = 6.0;
  std::vector<System> systems{
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {a, b}, {}, {4, 4}, 5)};
  MonteCarlo mc({20, 0, 3, 3, 1000, 10000, 5000, 5000}, systems, 23uz, 5, false);
  mc.run();
  for (System& s : mc.systems)
  {
    const RunningEnergy drift = s.runningEnergies - s.computeTotalEnergies();
    EXPECT_NEAR(drift.intraCoul, 0.0, 1e-6);
    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
  }
}

TEST(MC_CFCMC_FLEXIBLE_CHARGED_DRIFT, gibbs_swap_cfcmc_between_two_boxes)
{
  const ForceField forceField = makeChargedAlkaneForceField();
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::GibbsSwapCFCMC, 1.0);
  const Component butane = makeChargedButane(forceField, probabilities);
  std::vector<System> systems{
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {butane}, {}, {5}, 5),
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {butane}, {}, {5}, 5)};
  MonteCarlo mc({20, 0, 3, 3, 1000, 10000, 5000, 5000}, systems, 29uz, 5, false);
  mc.run();
  for (System& s : mc.systems)
  {
    const RunningEnergy drift = s.runningEnergies - s.computeTotalEnergies();
    EXPECT_NEAR(drift.intraCoul, 0.0, 1e-6);
    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
  }
}

namespace
{
void runGibbsDrift(Move::Types move, std::size_t seed)
{
  const ForceField forceField = makeChargedAlkaneForceField();
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(move, 1.0);
  const Component butane = makeChargedButane(forceField, probabilities);
  std::vector<System> systems{
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {butane}, {}, {5}, 5),
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {butane}, {}, {5}, 5)};
  MonteCarlo mc({20, 0, 3, 3, 1000, 10000, 5000, 5000}, systems, seed, 5, false);
  mc.run();
  for (System& s : mc.systems)
  {
    const RunningEnergy drift = s.runningEnergies - s.computeTotalEnergies();
    EXPECT_NEAR(drift.intraCoul, 0.0, 1e-6);
    EXPECT_NEAR(drift.potentialEnergy(), 0.0, 1e-6);
  }
}
}  // namespace

TEST(MC_CFCMC_FLEXIBLE_CHARGED_DRIFT, gibbs_swap_cbcfcmc_between_two_boxes)
{
  runGibbsDrift(Move::Types::GibbsSwapCBCFCMC, 31);
}

TEST(MC_CFCMC_FLEXIBLE_CHARGED_DRIFT, gibbs_conventional_cfcmc_between_two_boxes)
{
  runGibbsDrift(Move::Types::GibbsConventionalCFCMC, 37);
}
