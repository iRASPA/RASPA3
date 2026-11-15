#include <gtest/gtest.h>

import std;

import int3;
import forcefield;
import simulationbox;
import atom;
import framework;
import component;
import system;
import units;
import mc_moves_move_types;
import mc_moves_probabilities;
import mc_moves_parallel_tempering_swap;
import property_energy_histogram;
import isotherm_bet;
import dctst_wham;
import reweighted_histogram;

// Regression tests for the WHAM solvers, the reweighted-histogram guards and the auto-BET setup helpers.

namespace
{

System makeCO2InFAU(double temperature, std::size_t numberOfMolecules = 0)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  const Framework framework = Framework::makeFAU(forceField, int3(1, 1, 1));
  Component co2 = Component::makeCO2(forceField, 0, true);
  co2.mc_moves_probabilities.setProbability(Move::Types::Translation, 1.0);
  co2.mc_moves_probabilities.setProbability(Move::Types::Swap, 1.0);
  return System(forceField, std::nullopt, false, temperature, 1e5, 1.0, {framework}, {co2}, {}, {numberOfMolecules},
                5);
}

System makeCO2InBox(double temperature)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component co2 = Component::makeCO2(forceField, 0, true);
  co2.mc_moves_probabilities.setProbability(Move::Types::Translation, 1.0);
  co2.mc_moves_probabilities.setProbability(Move::Types::Swap, 1.0);
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, temperature, 1e5, 1.0, {}, {co2}, {}, {0}, 5);
  system.tmmc.maxMacrostate = 20;
  return system;
}

ReweightedHistogramParameters quietParameters()
{
  ReweightedHistogramParameters parameters;
  parameters.numberOfInitializationCycles = 0;
  parameters.numberOfProductionCycles = 10;
  parameters.writeBinaryRestartEvery = 0;
  return parameters;
}

}  // namespace

// State 0 has no samples and a forbidden bias over every bin (a dcTST pilot without cycles in occupancy mode):
// the gauge was fixed on it, moved every other g by its ~1e30 and the iteration alternated without converging.
TEST(WHAM_AND_BET_SETUP_FIXES, dctst_wham_gauges_on_the_first_sampled_state)
{
  const std::size_t numberOfBins = 10;
  std::vector<std::vector<double>> counts(3, std::vector<double>(numberOfBins, 0.0));
  std::vector<std::vector<double>> bias(3, std::vector<double>(numberOfBins, 0.0));
  for (std::size_t b = 0; b < numberOfBins; ++b)
  {
    bias[0][b] = 1e30;
    bias[2][b] = 0.1 * static_cast<double>(b);
    counts[1][b] = 100.0;
    counts[2][b] = 100.0 * std::exp(-bias[2][b]);
  }

  const DcTSTWham::Result result = DcTSTWham::solve(counts, bias, 1e-10, 10000);
  EXPECT_TRUE(result.converged);
  EXPECT_TRUE(result.connected);
  EXPECT_LT(result.iterations, 10000uz);
  EXPECT_NEAR(result.freeEnergies[1], 0.0, 1e-12);
}

// Two states whose biases leave no bin in common: their free-energy offset is undetermined. The iteration used to
// stop on its first step with a zero residual and report convergence.
TEST(WHAM_AND_BET_SETUP_FIXES, dctst_wham_reports_disconnected_states)
{
  const std::size_t numberOfBins = 10;
  std::vector<std::vector<double>> counts(2, std::vector<double>(numberOfBins, 0.0));
  std::vector<std::vector<double>> bias(2, std::vector<double>(numberOfBins, 0.0));
  for (std::size_t b = 0; b < numberOfBins; ++b)
  {
    const bool left = b < numberOfBins / 2;
    bias[0][b] = left ? 0.0 : 2000.0;
    bias[1][b] = left ? 2000.0 : 0.0;
    counts[0][b] = left ? 50.0 : 0.0;
    counts[1][b] = left ? 0.0 : 50.0;
  }

  const DcTSTWham::Result disconnected = DcTSTWham::solve(counts, bias, 1e-10, 10000);
  EXPECT_FALSE(disconnected.connected);
  EXPECT_FALSE(disconnected.converged);

  // a moderate bias overlaps, and the same layout is then one connected group
  for (std::size_t b = 0; b < numberOfBins; ++b)
  {
    const bool left = b < numberOfBins / 2;
    bias[0][b] = left ? 0.0 : 2.0;
    bias[1][b] = left ? 2.0 : 0.0;
    counts[0][b] = left ? 50.0 : 50.0 * std::exp(-2.0);
    counts[1][b] = left ? 50.0 * std::exp(-2.0) : 50.0;
  }
  const DcTSTWham::Result connected = DcTSTWham::solve(counts, bias, 1e-10, 10000);
  EXPECT_TRUE(connected.connected);
  EXPECT_TRUE(connected.converged);
}

// A flexible molecule samples beta f / Z_intra(beta): across temperatures Z_intra does not cancel and the swap is
// refused, at one temperature (a pressure ladder) it does and the swap goes ahead.
TEST(WHAM_AND_BET_SETUP_FIXES, flexible_molecules_swap_only_at_equal_temperature)
{
  System cold = makeCO2InFAU(300.0);
  System coldOther = makeCO2InFAU(300.0);
  System hot = makeCO2InFAU(350.0);
  for (System* system : {&cold, &coldOther, &hot})
  {
    system->components[0].rigid = false;
  }
  EXPECT_TRUE(MC_Moves::parallelTemperingUnsupportedReason(cold, false).has_value());
  EXPECT_TRUE(MC_Moves::ParallelTemperingLogAcceptance(cold, coldOther).has_value());
  EXPECT_FALSE(MC_Moves::ParallelTemperingLogAcceptance(cold, hot).has_value());
}

TEST(WHAM_AND_BET_SETUP_FIXES, reweighted_histogram_refuses_what_it_can_not_reweight)
{
  // a flexible adsorbate across temperatures
  {
    System system = makeCO2InBox(300.0);
    system.components[0].rigid = false;
    EXPECT_THROW(ReweightedHistogram(system, {300.0, 350.0}, {1e4}, quietParameters()), std::runtime_error);
    // at one temperature it is fine
    EXPECT_NO_THROW(ReweightedHistogram(system, {300.0}, {1e4, 1e5}, quietParameters()));
  }
  // volume moves
  {
    System system = makeCO2InBox(300.0);
    system.mc_moves_probabilities.setProbability(Move::Types::VolumeChange, 0.1);
    EXPECT_THROW(ReweightedHistogram(system, {300.0}, {1e4, 1e5}, quietParameters()), std::runtime_error);
  }
  // no insertion/deletion move on a pressure ladder
  {
    System system = makeCO2InBox(300.0);
    system.components[0].mc_moves_probabilities.setProbability(Move::Types::Swap, 0.0);
    EXPECT_THROW(ReweightedHistogram(system, {300.0}, {1e4, 1e5}, quietParameters()), std::runtime_error);
  }
}

// The energy-histogram bins were laid out from the lower limit over |upper - lower|: a reversed range binned the
// interval above the lower limit instead of the one between the limits.
TEST(WHAM_AND_BET_SETUP_FIXES, energy_histogram_refuses_a_reversed_range)
{
  EXPECT_THROW(PropertyEnergyHistogram(5, 128, {-1000.0, -5000.0}, 1, std::nullopt), std::invalid_argument);
  EXPECT_THROW(PropertyEnergyHistogram(5, 128, {-1000.0, -1000.0}, 1, std::nullopt), std::invalid_argument);
  EXPECT_NO_THROW(PropertyEnergyHistogram(5, 128, {-5000.0, -1000.0}, 1, std::nullopt));
}

// A 'FugacityCoefficient' given in the input is kept by the pressure pin of the WHAM and TMMC replicas.
TEST(WHAM_AND_BET_SETUP_FIXES, pressure_pin_keeps_a_given_fugacity_coefficient)
{
  System system = makeCO2InFAU(300.0);
  Component& co2 = system.components[0];
  co2.criticalTemperature = 304.1282;
  co2.criticalPressure = 7377300.0;
  co2.acentricFactor = 0.22394;

  co2.fugacityCoefficient = 1.0;
  co2.fugacityCoefficientGiven = true;
  pinSystemPengRobinsonPressure(system, 5e6);
  ASSERT_TRUE(system.components[0].fugacityCoefficient.has_value());
  EXPECT_EQ(system.components[0].fugacityCoefficient.value(), 1.0);

  system.components[0].fugacityCoefficientGiven = false;
  pinSystemPengRobinsonPressure(system, 5e6);
  ASSERT_TRUE(system.components[0].fugacityCoefficient.has_value());
  EXPECT_LT(system.components[0].fugacityCoefficient.value(), 0.9);
}

// The Henry scout runs on an empty host at the temperature it is asked at.
TEST(WHAM_AND_BET_SETUP_FIXES, scout_helpers_pin_the_temperature_and_empty_the_host)
{
  System system = makeCO2InFAU(300.0, 4);
  ASSERT_TRUE(system.components[0].swappable);
  ASSERT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 4uz);

  pinSystemTemperature(system, 77.0);
  EXPECT_EQ(system.temperature, 77.0);
  EXPECT_DOUBLE_EQ(system.beta, 1.0 / (Units::KB * 77.0));
  EXPECT_EQ(system.forceField.temperature, 77.0);

  removeSwappableIntegerMolecules(system);
  EXPECT_EQ(system.numberOfIntegerMoleculesPerComponent[0], 0uz);
  EXPECT_EQ(system.numberOfMoleculesPerComponent[0], 0uz);
  EXPECT_TRUE(system.spanOfMoleculeAtoms().empty());
}
