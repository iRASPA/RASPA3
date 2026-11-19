#include <gtest/gtest.h>

import std;

import randomnumbers;
import forcefield;
import component;
import system;
import simulationbox;
import running_energy;
import move_statistics;
import mc_moves_statistics;
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_identity_change;
import mc_moves_identity_switch;
import mc_moves_gibbs_identity_change;

// The dispatcher books the overall count and the total time of a component move on the component the move was
// drawn for. The identity moves involve a second component (the partner of the switch, or the old component
// after the coin flip of the identity change), but their trials, constructions and acceptances belong to the
// same component as the overall count, or the printed statistics of the two components disagree and the sum
// over the components counts a switch twice.

namespace
{
System makeTwoBeadSystem(std::size_t beadsA, std::size_t beadsB)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, false);
  const std::size_t typeCH4 = forceField.findPseudoAtom("CH4").value();

  Component beadA = Component::makeIon(forceField, 0, "beadA", typeCH4, 0.0);
  beadA.identitySwitches = {1};
  beadA.identityChanges = {1};
  beadA.molFraction = 0.5;
  beadA.idealGasRosenbluthWeight = 1.0;

  Component beadB = Component::makeIon(forceField, 1, "beadB", typeCH4, 0.0);
  beadB.molFraction = 0.5;
  beadB.idealGasRosenbluthWeight = 1.0;

  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {beadA, beadB}, {},
                {beadsA, beadsB}, 5);
  system.runningEnergies = system.computeTotalEnergies();
  return system;
}
}  // namespace

TEST(MC_IDENTITY_MOVE_STATISTICS, identity_switch_counts_on_the_component_it_was_drawn_for)
{
  System system = makeTwoBeadSystem(30, 10);
  RandomNumber random(7);
  const std::size_t numberOfTrials = 200;
  for (std::size_t i = 0; i != numberOfTrials; ++i)
  {
    const std::optional<RunningEnergy> energy = MC_Moves::identitySwitchMove(random, system, 0);
    if (energy.has_value()) system.runningEnergies += energy.value();
  }
  const MoveStatistics<double>& switchesA =
      std::get<MoveStatistics<double>>(system.components[0].mc_moves_statistics[Move::Types::IdentitySwitchCBMC]);
  const MoveStatistics<double>& switchesB =
      std::get<MoveStatistics<double>>(system.components[1].mc_moves_statistics[Move::Types::IdentitySwitchCBMC]);
  EXPECT_EQ(switchesA.totalCounts, static_cast<double>(numberOfTrials));
  EXPECT_EQ(switchesB.totalCounts, 0.0);
}

TEST(MC_IDENTITY_MOVE_STATISTICS, identity_change_counts_on_the_component_it_was_drawn_for)
{
  System system = makeTwoBeadSystem(30, 30);
  RandomNumber random(7);
  const std::size_t numberOfTrials = 200;
  for (std::size_t i = 0; i != numberOfTrials; ++i)
  {
    const std::optional<RunningEnergy> energy = MC_Moves::identityChangeMove(random, system, 0);
    if (energy.has_value()) system.runningEnergies += energy.value();
  }
  // the coin flip makes the partner's molecule the one that changes half of the time; the trial is still one
  // of component 0
  const MoveStatistics<double>& changesA =
      std::get<MoveStatistics<double>>(system.components[0].mc_moves_statistics[Move::Types::IdentityChangeCBMC]);
  const MoveStatistics<double>& changesB =
      std::get<MoveStatistics<double>>(system.components[1].mc_moves_statistics[Move::Types::IdentityChangeCBMC]);
  EXPECT_EQ(changesA.totalCounts, static_cast<double>(numberOfTrials));
  EXPECT_EQ(changesB.totalCounts, 0.0);
  EXPECT_GT(changesA.totalAccepted, 0.0);
}

// An attempt that finds no molecule to change is an attempt, and is counted before it is turned down, as the swap
// moves count theirs. The trials were counted after the test for an empty component, so that the acceptance ratio
// was taken over the attempts that found a molecule only.
TEST(MC_IDENTITY_MOVE_STATISTICS, an_attempt_on_an_empty_component_is_counted)
{
  const std::size_t numberOfTrials = 50;
  const auto trials = [](const System& system, std::size_t component, Move::Types move)
  { return std::get<MoveStatistics<double>>(system.components[component].mc_moves_statistics[move]).totalCounts; };
  {
    System system = makeTwoBeadSystem(0, 10);
    RandomNumber random(7);
    for (std::size_t i = 0; i != numberOfTrials; ++i)
    {
      EXPECT_FALSE(MC_Moves::identitySwitchMove(random, system, 0).has_value());
    }
    EXPECT_EQ(trials(system, 0, Move::Types::IdentitySwitchCBMC), static_cast<double>(numberOfTrials));
  }
  {
    System system = makeTwoBeadSystem(0, 0);
    RandomNumber random(7);
    for (std::size_t i = 0; i != numberOfTrials; ++i)
    {
      EXPECT_FALSE(MC_Moves::identityChangeMove(random, system, 0).has_value());
    }
    EXPECT_EQ(trials(system, 0, Move::Types::IdentityChangeCBMC), static_cast<double>(numberOfTrials));
  }
  {
    // the Gibbs identity change counts a trial on the component of each box that would change
    System boxI = makeTwoBeadSystem(0, 0);
    System boxII = makeTwoBeadSystem(0, 0);
    boxI.components[0].gibbsIdentityChanges = {1};
    boxII.components[0].gibbsIdentityChanges = {1};
    RandomNumber random(7);
    for (std::size_t i = 0; i != numberOfTrials; ++i)
    {
      EXPECT_FALSE(MC_Moves::GibbsIdentityChangeMove_CBMC(random, boxI, boxII, 0).has_value());
    }
    double total = 0.0;
    for (const System* box : {&boxI, &boxII})
    {
      for (std::size_t component : {0uz, 1uz}) total += trials(*box, component, Move::Types::GibbsIdentityChangeCBMC);
    }
    EXPECT_EQ(total, 2.0 * static_cast<double>(numberOfTrials));
  }
}
