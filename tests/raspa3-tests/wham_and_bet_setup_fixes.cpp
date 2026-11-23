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
import property_lambda_probability_histogram;
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

// A flexible molecule samples beta f / Z_intra(beta): at one temperature (a pressure ladder) Z_intra cancels and
// the swap goes ahead; across temperatures the ratio is needed, and a pair of replicas that was not given it
// ('Preparation::prepareIdealGasReferences') is refused. Under solute tempering it is refused up front.
TEST(WHAM_AND_BET_SETUP_FIXES, flexible_molecules_swap_only_at_equal_temperature)
{
  System cold = makeCO2InFAU(300.0);
  System coldOther = makeCO2InFAU(300.0);
  System hot = makeCO2InFAU(350.0);
  for (System* system : {&cold, &coldOther, &hot})
  {
    system->components[0].rigid = false;
  }
  EXPECT_FALSE(MC_Moves::parallelTemperingUnsupportedReason(cold, false).has_value());
  EXPECT_TRUE(MC_Moves::parallelTemperingUnsupportedReason(cold, true).has_value());
  EXPECT_TRUE(MC_Moves::ParallelTemperingLogAcceptance(cold, coldOther).has_value());
  EXPECT_FALSE(MC_Moves::ParallelTemperingLogAcceptance(cold, hot).has_value());
}

TEST(WHAM_AND_BET_SETUP_FIXES, reweighted_histogram_refuses_what_it_can_not_reweight)
{
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

// The density of states returned with the free energies is the one those free energies give,
// Omega_b = M_b / sum_k n_k exp(g_k - u_k(b)). It was computed from the free energies of the iteration before (and
// before that iteration's gauge shift), and with no iteration at all it was never computed.
TEST(WHAM_AND_BET_SETUP_FIXES, dctst_wham_density_belongs_to_the_returned_free_energies)
{
  const std::size_t numberOfBins = 12;
  std::vector<std::vector<double>> counts(3, std::vector<double>(numberOfBins, 0.0));
  std::vector<std::vector<double>> bias(3, std::vector<double>(numberOfBins, 0.0));
  for (std::size_t b = 0; b < numberOfBins; ++b)
  {
    const double x = static_cast<double>(b);
    bias[1][b] = 0.3 * x;
    bias[2][b] = 0.05 * (x - 8.0) * (x - 8.0);
    counts[0][b] = 50.0 + 3.0 * x;
    counts[1][b] = 200.0 * std::exp(-bias[1][b]);
    counts[2][b] = 90.0 * std::exp(-bias[2][b]) + 1.0;
  }
  auto expectedLogDensity = [&](const DcTSTWham::Result& result, std::size_t b)
  {
    double binCount = 0.0;
    double denominator = 0.0;
    for (std::size_t k = 0; k < 3; ++k)
    {
      binCount += counts[k][b];
      const double stateCount = std::accumulate(counts[k].begin(), counts[k].end(), 0.0);
      denominator += stateCount * std::exp(result.freeEnergies[k] - bias[k][b]);
    }
    return std::log(binCount) - std::log(denominator);
  };

  // stopped far from convergence, so that one iteration of lag is visible
  for (const std::size_t iterations : {0uz, 1uz, 3uz})
  {
    const DcTSTWham::Result result = DcTSTWham::solve(counts, bias, 1e-14, iterations);
    EXPECT_EQ(result.iterations, iterations);
    EXPECT_FALSE(result.converged);
    for (std::size_t b = 0; b < numberOfBins; ++b)
    {
      EXPECT_NEAR(result.logDensity[b], expectedLogDensity(result, b), 1e-10) << iterations << " iterations, bin " << b;
    }
  }
}

// The probe of the BET scouts is the adsorbate, which fixed extra-framework cations may precede. The scouts took
// component 0: they asked the cation for its BET probe properties and counted the cations as the loading.
TEST(WHAM_AND_BET_SETUP_FIXES, bet_scout_probes_the_adsorbate_behind_a_cation)
{
  const ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Component cation = Component::makeCO2(forceField, 0, true);
  cation.name = "cation";
  cation.type = Component::Type::Cation;
  Component probe = Component::makeCO2(forceField, 1, true);
  probe.mc_moves_probabilities.setProbability(Move::Types::Translation, 1.0);
  probe.mc_moves_probabilities.setProbability(Move::Types::Swap, 1.0);
  probe.crossSection = 16.2;
  probe.liquidVolume = 57.7;
  probe.saturationPressure = 101325.0;
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {cation, probe}, {}, {3, 0},
                5);
  ASSERT_EQ(system.tmmcComponent(), 1uz);

  NitrogenBETScoutPoint point;
  ASSERT_NO_THROW(point = scoutNitrogenBETOccupancy(system, 2.0e6, 1000));
  // about 13 molecules of an ideal gas at this pressure, and not the three cations
  EXPECT_GT(point.meanOccupancy, 6.0);
  EXPECT_LT(point.meanOccupancy, 40.0);
  EXPECT_GT(point.peakOccupancy, 6uz);
}

// A bin that holds samples but is forbidden (a reduced potential of 1e30) in every sampled state: its largest log
// weight, about -1e30, is far above the 'log of zero', so the bin got the density ln M + 1e30 and entered every
// free energy as the difference of two numbers of that size.
TEST(WHAM_AND_BET_SETUP_FIXES, dctst_wham_leaves_out_a_bin_forbidden_in_every_state)
{
  const std::vector<std::vector<double>> counts{{100.0, 50.0, 7.0}, {40.0, 110.0, 3.0}};
  const std::vector<std::vector<double>> reducedPotential{{0.0, 0.0, 1e30}, {0.5, -0.5, 1e30}};
  const DcTSTWham::Result result = DcTSTWham::solve(counts, reducedPotential, 1e-10, 10000);

  EXPECT_LE(result.logDensity[2], DcTSTWham::logZeroThreshold);
  EXPECT_GT(result.logDensity[0], -50.0);
  EXPECT_LT(result.logDensity[0], 50.0);
  EXPECT_LT(std::abs(result.freeEnergies[1]), 5.0);
  EXPECT_TRUE(result.converged);

  // the same two states without the forbidden bin
  const DcTSTWham::Result reference =
      DcTSTWham::solve({{100.0, 50.0}, {40.0, 110.0}}, {{0.0, 0.0}, {0.5, -0.5}}, 1e-10, 10000);
  EXPECT_NEAR(result.freeEnergies[1], reference.freeEnergies[1], 1e-8);
}

// Wang-Landau only lowers the bias: after an equilibration the raw values are hundreds below zero and the weight
// exp(-bias) of a production sample overflows. Each histogram is shifted by its own minimum.
TEST(WHAM_AND_BET_SETUP_FIXES, lambda_bias_is_normalized_to_its_own_minimum)
{
  PropertyLambdaProbabilityHistogram histogram(5, 3);
  histogram.biasFactor = {-800.0, -750.0, -790.0};
  histogram.normalizeToMinimum();
  EXPECT_DOUBLE_EQ(histogram.biasFactor[0], 0.0);
  EXPECT_DOUBLE_EQ(histogram.biasFactor[1], 50.0);
  EXPECT_DOUBLE_EQ(histogram.biasFactor[2], 10.0);
  for (std::size_t bin = 0; bin < 3uz; ++bin)
  {
    histogram.currentBin = bin;
    EXPECT_LE(histogram.weight(), 1.0);
    EXPECT_GE(histogram.weight(), 0.0);
  }

  System system = makeCO2InBox(300.0);
  system.components[0].lambdaGC.biasFactor.assign(system.components[0].lambdaGC.biasFactor.size(), -900.0);
  system.components[0].lambdaGC.biasFactor.front() = -1000.0;
  system.components[0].lambdaPairSwap.biasFactor.assign(system.components[0].lambdaPairSwap.biasFactor.size(), -5.0);
  system.normalizeLambdaBiasFactors();
  EXPECT_DOUBLE_EQ(system.components[0].lambdaGC.biasFactor.front(), 0.0);
  EXPECT_DOUBLE_EQ(system.components[0].lambdaGC.biasFactor.back(), 100.0);
  // by its own minimum, and not by the -1000 of the other histogram (a weight of exp(-995) is zero)
  EXPECT_DOUBLE_EQ(system.components[0].lambdaPairSwap.biasFactor.front(), 0.0);
}

// The bimodal window of a first-order transition is a hundredth wide in ln f for a few hundred molecules; 200 scan
// points over ten units are 0.05 apart. No two of them were bimodal, and no coexistence was found.
TEST(WHAM_AND_BET_SETUP_FIXES, equal_weight_point_is_found_inside_a_window_narrower_than_the_scan)
{
  const double transition = 1.23456;
  auto evaluate = [&](double x) -> MacrostateScanPoint
  {
    if (std::abs(x - transition) < 0.004) return {500.0 * (x - transition), 205.0 + 49000.0 * (x - transition)};
    return {std::nullopt, x < transition ? 10.0 + 0.1 * x : 400.0 + 0.1 * x};
  };
  const std::optional<double> located = locateEqualWeightPoint(-5.0, 5.0, 200, evaluate);
  ASSERT_TRUE(located.has_value());
  EXPECT_NEAR(located.value(), transition, 1e-9);

  // two bimodal scan points on either side are still the bracket
  auto wide = [&](double x) -> MacrostateScanPoint
  {
    if (std::abs(x - transition) < 0.5) return {500.0 * (x - transition), 205.0 + 390.0 * (x - transition)};
    return {std::nullopt, x < transition ? 10.0 : 400.0};
  };
  const std::optional<double> locatedWide = locateEqualWeightPoint(-5.0, 5.0, 200, wide);
  ASSERT_TRUE(locatedWide.has_value());
  EXPECT_NEAR(locatedWide.value(), transition, 1e-9);

  // a loading that does not rise anywhere has no transition, and neither has an empty range
  auto flat = [](double) -> MacrostateScanPoint { return {std::nullopt, 50.0}; };
  EXPECT_FALSE(locateEqualWeightPoint(-5.0, 5.0, 200, flat).has_value());
  EXPECT_FALSE(locateEqualWeightPoint(1.0, 1.0, 200, evaluate).has_value());
}

// The swap divides the activity of each replica by its own ideal-gas Rosenbluth weight, as its insertions do: a
// weight was refused outright, and a pressure ladder of ReweightedHistogram with one ran without an accepted swap.
// A flexible molecule is still refused across temperatures, with or without a weight.
//
// The weight enters where the insertions divide by it and nowhere else: the swaps that grow the molecule do, and
// the conventional swap, which draws the conformation, samples beta f whatever weight was given. So the replicas
// here insert with the growing swap, and the same weights on replicas with the conventional one change nothing.
TEST(WHAM_AND_BET_SETUP_FIXES, an_ideal_gas_rosenbluth_weight_enters_the_swap_through_the_activity)
{
  System cold = makeCO2InFAU(300.0, 2);
  System coldOther = makeCO2InFAU(300.0, 5);
  System hot = makeCO2InFAU(350.0, 5);
  auto insertWith = [&](Move::Types wanted, Move::Types other)
  {
    for (System* system : {&cold, &coldOther, &hot})
    {
      system->components[0].mc_moves_probabilities.setProbability(wanted, 1.0);
      system->components[0].mc_moves_probabilities.setProbability(other, 0.0);
    }
  };
  insertWith(Move::Types::SwapCBMC, Move::Types::Swap);
  const double reference = MC_Moves::ParallelTemperingLogAcceptance(cold, coldOther).value();
  cold.components[0].idealGasRosenbluthWeight = 0.5;
  coldOther.components[0].idealGasRosenbluthWeight = 0.25;
  hot.components[0].idealGasRosenbluthWeight = 0.25;
  // (N_B - N_A) ln(a_A / a_B) with a = beta f / W_IG
  const std::optional<double> weighted = MC_Moves::ParallelTemperingLogAcceptance(cold, coldOther);
  ASSERT_TRUE(weighted.has_value());
  EXPECT_NEAR(weighted.value() - reference, 3.0 * std::log(0.25 / 0.5), 1e-9);
  EXPECT_TRUE(MC_Moves::ParallelTemperingLogAcceptance(cold, hot).has_value());

  // with the conventional swap the weights are not what the insertions divide by, and the exchange leaves them out
  insertWith(Move::Types::Swap, Move::Types::SwapCBMC);
  const std::optional<double> conventional = MC_Moves::ParallelTemperingLogAcceptance(cold, coldOther);
  ASSERT_TRUE(conventional.has_value());
  EXPECT_NEAR(conventional.value(), reference, 1e-9);

  for (System* system : {&cold, &coldOther, &hot}) system->components[0].rigid = false;
  EXPECT_TRUE(MC_Moves::ParallelTemperingLogAcceptance(cold, coldOther).has_value());
  EXPECT_FALSE(MC_Moves::ParallelTemperingLogAcceptance(cold, hot).has_value());
}

// A state that takes no part stays at zero: the gauge shift moved it along with the others, by the shift of the
// last iteration (visible in a solve that is stopped before it converged).
TEST(WHAM_AND_BET_SETUP_FIXES, dctst_wham_does_not_gauge_a_state_that_takes_no_part)
{
  const DcTSTWham::Result result = DcTSTWham::solve({{0.0, 0.0}, {100.0, 50.0}, {40.0, 110.0}},
                                                    {{1e30, 1e30}, {0.0, 0.0}, {0.5, -0.5}}, 1e-10, 1);
  EXPECT_FALSE(result.converged);
  EXPECT_EQ(result.freeEnergies[0], 0.0);
}

// A fractional molecule pinned at a lambda is never decoupled, and the samples are taken at lambda = 0 only.
TEST(WHAM_AND_BET_SETUP_FIXES, reweighted_histogram_refuses_a_fixed_lambda)
{
  System system = makeCO2InBox(300.0);
  system.components[0].fixedLambdaBin = 5uz;
  EXPECT_THROW(ReweightedHistogram(system, {300.0}, {1e4, 1e5}, quietParameters()), std::runtime_error);
}
