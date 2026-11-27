#include <gtest/gtest.h>

import std;

import double3;
import randomnumbers;
import atom;
import forcefield;
import component;
import system;
import simulationbox;
import running_energy;
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_deletion_cbmc;
import mc_moves_insertion_cbmc;
import cbmc;

// The Ewald, tail and polarization energy differences of the swap-type moves enter the acceptance rule in
// log space. Their separate Boltzmann factor is not representable when |beta dU| exceeds ~709 (it becomes
// 0 or +inf and its logarithm -inf or +inf), while the complete exponent of the acceptance rule is of a
// normal size when the other factors (the Rosenbluth weight, the fugacity prefactor) compensate. This is
// arranged here with a single Lennard-Jones atom whose (absurdly large) epsilon gives a tail correction of
// beta dU ~ -900 on deletion, compensated by a tiny fugacity in the prefactor; with an intermediate
// Boltzmann factor the acceptance probability collapses to zero instead of exp(-200).
TEST(MC_ACCEPTANCE_LOG_SPACE, deletion_with_an_enormous_tail_correction_keeps_the_acceptance_finite)
{
  const ForceField forceField({{"X", false, 16.0, 0.0, 0.0, 6, false}}, {{6.0e8, 3.0}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, false, true, false);
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::SwapCBMC, 1.0);
  const Component x(forceField, "X", 190.0, 4.5e6, 0.01, {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false)},
                    {}, {}, 5, 21, probabilities, std::nullopt, false);

  System system(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1e4, 1.0, {}, {x}, {}, {1}, 5);
  system.runningEnergies = system.computeTotalEnergies();

  // the tail energy of the lone atom (its self term) is removed by the deletion: dU_tail = -U_tail > 0
  const double betaTail = system.beta * system.runningEnergies.tail;
  ASSERT_LT(betaTail, -750.0);   // the Boltzmann factor exp(-beta dU_tail) alone underflows to zero
  ASSERT_GT(betaTail, -1300.0);  // but the complete exponent below is representable

  // prefactor N / (beta f V) of the deletion with N = 1: choose the fugacity such that its logarithm is +700
  const double volume = system.simulationBox.volume;
  system.pressure = std::exp(-700.0) / (system.beta * volume);
  const Component& component = system.components[0];
  const double fugacity = component.molFraction * component.fugacityCoefficient.value_or(1.0) * system.pressure;
  const double logPrefactor = -std::log(system.beta * fugacity * volume);
  ASSERT_NEAR(logPrefactor, 700.0, 1e-2);  // the fugacity coefficient of the equation of state is close to one

  // the retrace in the otherwise empty box sees zero energy at every trial position: log W = 0
  {
    RandomNumber retraceRandom(5);
    const CBMC::RetraceResult retrace = CBMC::retraceMolecule(retraceRandom, system.makeGrowContext(),
                                                              system.components[0], system.spanOfMolecule(0, 0));
    ASSERT_EQ(retrace.logRosenbluthWeight, 0.0);
  }
  const double expectedLogPacc = betaTail + logPrefactor;
  ASSERT_LT(expectedLogPacc, -100.0);
  ASSERT_GT(expectedLogPacc, -300.0);

  RandomNumber random(5);
  const auto [energy, acceptance] = MC_Moves::deletionMoveCBMC(random, system, 0, 0);
  const double Pacc = acceptance.x;
  ASSERT_GT(Pacc, 0.0);
  EXPECT_TRUE(std::isfinite(Pacc));
  EXPECT_NEAR(std::log(Pacc), expectedLogPacc, 1e-6);
  EXPECT_FALSE(energy.has_value());  // exp(-200) is never accepted
}

// A reservoir of no fugacity (an external pressure of zero, a desorption into vacuum, or a mole fraction of zero): an
// insertion from it is never accepted, a deletion into it always. The logarithm of the prefactor beta f V / (N + 1)
// was log(0) and the prefactor N / (beta f V) of the deletion a division by zero, infinities the fast-math flags of
// the build assume never occur. (NaN read off the printed value: the fast-math build folds the comparisons.)
TEST(MC_ACCEPTANCE_LOG_SPACE, swap_moves_with_a_reservoir_of_no_fugacity)
{
  const ForceField forceField({{"X", false, 16.0, 0.0, 0.0, 6, false}}, {{100.0, 3.0}},
                              ForceField::MixingRule::Lorentz_Berthelot, 9.0, 9.0, 9.0, false, false, false);
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::SwapCBMC, 1.0);
  const Component x(forceField, "X", 190.0, 4.5e6, 0.01, {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false)}, {},
                    {}, 5, 21, probabilities, std::nullopt, false);
  System system(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {x}, {}, {3}, 5);
  system.pressure = 0.0;
  system.runningEnergies = system.computeTotalEnergies();

  RandomNumber random(7);
  for (std::size_t attempt = 0; attempt < 20; ++attempt)
  {
    const auto [energy, acceptance] = MC_Moves::insertionMoveCBMC(random, system, 0);
    EXPECT_FALSE(energy.has_value());
    EXPECT_EQ(std::to_string(acceptance.z), std::to_string(0.0));
  }
  const auto [energy, acceptance] = MC_Moves::deletionMoveCBMC(random, system, 0, 0);
  EXPECT_TRUE(energy.has_value());
  EXPECT_EQ(std::to_string(acceptance.x).find("nan"), std::string::npos);
  EXPECT_EQ(std::to_string(acceptance.x).find("inf"), std::string::npos);
  EXPECT_GE(acceptance.x, 1.0);
}
