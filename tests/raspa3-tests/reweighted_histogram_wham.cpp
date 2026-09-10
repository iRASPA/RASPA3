#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import atom;
import averages;
import component;
import forcefield;
import isotherm_bet;
import property_loading;
import randomnumbers;
import reweighted_histogram;
import simulationbox;
import system;
import units;

namespace
{

class ScopedCurrentPath
{
 public:
  explicit ScopedCurrentPath(const std::filesystem::path& path) : original_(std::filesystem::current_path())
  {
    std::filesystem::current_path(path);
  }

  ~ScopedCurrentPath() { std::filesystem::current_path(original_); }

 private:
  std::filesystem::path original_;
};

System makeIdealParticleSystem()
{
  const ForceField forceField({{"X", false, 1.0, 0.0, 0.0, 1, false}}, {{0.0, 1.0}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  const Component particle(forceField, "particle", 1.0, 1.0, 0.0,
                           {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21);
  System system(forceField, SimulationBox(20.0, 20.0, 20.0), false, 300.0, 1e4, 1.0, {}, {particle}, {}, {0}, 5);
  system.components[0].fugacityCoefficient = 1.0;
  return system;
}

ReweightedHistogram makeSingleStateDriver()
{
  System system = makeIdealParticleSystem();
  system.tmmc.maxMacrostate = 50;
  ReweightedHistogramParameters parameters;
  parameters.numberOfPreInitializationCycles = 0;
  parameters.numberOfInitializationCycles = 0;
  parameters.numberOfEquilibrationCycles = 0;
  parameters.numberOfProductionCycles = 10;
  parameters.numberOfBlocks = 5;
  parameters.parallelTemperingSwapEvery = 0;
  parameters.writeBinaryRestartEvery = 0;
  parameters.reweightingTemperatures = {300.0};
  parameters.reweightingPressureRange = std::make_pair(1e4, 1e4);
  parameters.reweightingNumberOfPressures = 1;
  return ReweightedHistogram(std::move(system), {300.0}, {1e4}, parameters);
}

void addBalancedSamples(ReweightedHistogram& simulation)
{
  simulation.reweightingSamples.assign(1, {});
  for (std::size_t block = 0; block < 5uz; ++block)
  {
    for (std::size_t sample = 0; sample < 20uz; ++sample)
    {
      simulation.reweightingSamples.front().push_back(
          ReweightedHistogram::Sample{.energy = 0.0, .numberOfMolecules = 0, .block = static_cast<std::uint32_t>(block)});
      simulation.reweightingSamples.front().push_back(
          ReweightedHistogram::Sample{.energy = 0.0, .numberOfMolecules = 1, .block = static_cast<std::uint32_t>(block)});
    }
  }
}

}  // namespace

TEST(REWEIGHTED_HISTOGRAM_WHAM, recovers_sampled_mean_and_reports_convergence)
{
  TemporaryDirectory workspace;
  ScopedCurrentPath currentPath(workspace.path());

  ReweightedHistogram simulation = makeSingleStateDriver();
  addBalancedSamples(simulation);
  simulation.computeBET = true;
  simulation.setup();
  simulation.performReweightingAnalysis();

  EXPECT_TRUE(simulation.whamConverged);
  EXPECT_LE(simulation.whamResidual, 1.0e-8);
  EXPECT_GT(simulation.whamIterations, 0uz);
  EXPECT_EQ(simulation.whamUnconvergedBlocks, 0uz);
  EXPECT_TRUE(simulation.whamUnconvergedBlockDetails.empty());
  ASSERT_TRUE(simulation.outputJson["output"]["reweighting"]["unconvergedBlockDetails"].is_array());
  EXPECT_TRUE(simulation.outputJson["output"]["reweighting"]["unconvergedBlockDetails"].empty());
  ASSERT_TRUE(simulation.outputJson["output"]["reweighting"]["states"].is_array());
  ASSERT_FALSE(simulation.outputJson["output"]["reweighting"]["states"].empty());
  EXPECT_NEAR(simulation.outputJson["output"]["reweighting"]["states"][0]["reweightedLoading"].get<double>(), 0.5,
              1.0e-8);
  EXPECT_TRUE(std::filesystem::exists("wham/reweighted_isotherm_300.reweighted_histogram.txt"));
  EXPECT_TRUE(std::filesystem::exists("wham/density_of_states.reweighted_histogram.txt"));
  EXPECT_TRUE(std::filesystem::exists("wham/reweighted_free_energies.reweighted_histogram.txt"));
}

TEST(REWEIGHTED_HISTOGRAM_WHAM, empty_samples_skip_analysis)
{
  TemporaryDirectory workspace;
  ScopedCurrentPath currentPath(workspace.path());

  ReweightedHistogram simulation = makeSingleStateDriver();
  simulation.reweightingSamples.assign(1, {});
  simulation.setup();
  simulation.performReweightingAnalysis();

  EXPECT_TRUE(simulation.whamConverged);
  EXPECT_EQ(simulation.whamIterations, 0uz);
  EXPECT_TRUE(simulation.reweightedIsotherms.empty());
}

TEST(REWEIGHTED_HISTOGRAM_WHAM, molecules_per_cell_is_box_total_without_framework)
{
  TemporaryDirectory workspace;
  ScopedCurrentPath currentPath(workspace.path());

  ReweightedHistogram simulation = makeSingleStateDriver();
  simulation.reweightingSamples.assign(1, {});
  for (std::size_t block = 0; block < 5uz; ++block)
  {
    for (std::size_t sample = 0; sample < 10uz; ++sample)
    {
      simulation.reweightingSamples.front().push_back(
          ReweightedHistogram::Sample{.energy = 0.0, .numberOfMolecules = 4, .block = static_cast<std::uint32_t>(block)});
    }
  }
  simulation.setup();
  simulation.performReweightingAnalysis();

  ASSERT_FALSE(simulation.reweightedIsotherms.empty());
  ASSERT_FALSE(simulation.reweightedIsotherms.front().points.empty());
  EXPECT_NEAR(simulation.reweightedIsotherms.front().points.front().moleculesPerCell, 4.0, 1.0e-12);
}

TEST(REWEIGHTED_HISTOGRAM_WHAM, inner_steps_equal_filling_ceiling)
{
  ReweightedHistogram simulation = makeSingleStateDriver();
  EXPECT_EQ(simulation.numberOfStepsPerCycle, 50uz);
  EXPECT_EQ(simulation.numberOfStepsPerCycle, simulation.systems.front().tmmc.maxMacrostate);
}

TEST(REWEIGHTED_HISTOGRAM_WHAM, lambda0_block_error_needs_three_blocks)
{
  const double mean = 2.0;
  const std::vector<double> twoBlocks{1.5, 2.5};
  EXPECT_DOUBLE_EQ(blockErrorEstimate(twoBlocks, mean), 0.0);

  const std::vector<double> fiveBlocks{1.6, 1.8, 2.0, 2.2, 2.4};
  EXPECT_GT(blockErrorEstimate(fiveBlocks, mean), 0.0);
}

TEST(REWEIGHTED_HISTOGRAM_WHAM, weight_ignores_lambda_bias_without_fractional_molecule)
{
  System system = makeIdealParticleSystem();
  ASSERT_FALSE(system.components[0].hasFractionalMolecule);
  EXPECT_FALSE(system.lambdaWangLandauIsActive(0));

  system.components[0].lambdaGC.biasFactor.assign(system.components[0].lambdaGC.numberOfSamplePoints, -1.0e4);
  EXPECT_DOUBLE_EQ(system.weight(), 1.0);
  EXPECT_TRUE(std::isfinite(system.weight()));

  system.components[0].hasFractionalMolecule = true;
  EXPECT_TRUE(system.lambdaWangLandauIsActive(0));
  system.components[0].lambdaGC.biasFactor[0] = std::log(2.0);
  EXPECT_NEAR(system.weight(), 0.5, 1.0e-12);
}

TEST(REWEIGHTED_HISTOGRAM_WHAM, empty_loading_average_is_finite)
{
  const SimulationBox box(10.0, 10.0, 10.0);
  const LoadingData empty(1, {0uz}, box);
  EXPECT_TRUE(std::isfinite(empty.inverseNumberDensities[0]));
  EXPECT_DOUBLE_EQ(empty.inverseNumberDensities[0], 0.0);

  PropertyLoading loadings(5, 1);
  for (std::size_t block = 0; block < 5uz; ++block)
  {
    loadings.addSample(block, empty, 1.0);
    loadings.addSample(block, LoadingData(1, {8uz}, box), 1.0);
  }

  const auto [mean, error] = loadings.average();
  EXPECT_TRUE(std::isfinite(mean.numberOfMolecules[0]));
  EXPECT_TRUE(std::isfinite(mean.inverseNumberDensities[0]));
  EXPECT_NEAR(mean.numberOfMolecules[0], 4.0, 1.0e-12);
}

TEST(REWEIGHTED_HISTOGRAM_WHAM, equal_occupancy_keeps_log_spacing_for_henry_isotherm)
{
  const std::vector<double> pressures{1.0e-2, 1.0e-1, 1.0, 1.0e1, 1.0e2};
  std::vector<double> occupancies(pressures.size());
  for (std::size_t i = 0; i < pressures.size(); ++i)
  {
    occupancies[i] = std::log(pressures[i]);
  }

  const std::vector<double> placed = placePressuresByEqualOccupancy(pressures, occupancies);
  ASSERT_EQ(placed.size(), pressures.size());
  EXPECT_DOUBLE_EQ(placed.front(), pressures.front());
  EXPECT_DOUBLE_EQ(placed.back(), pressures.back());
  for (std::size_t i = 0; i < pressures.size(); ++i)
  {
    EXPECT_NEAR(placed[i], pressures[i], 1.0e-6 * pressures[i]);
  }
}

TEST(REWEIGHTED_HISTOGRAM_WHAM, equal_occupancy_clusters_rungs_in_type_I_step)
{
  const std::vector<double> pressures{1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0, 1.0e1, 1.0e2, 1.0e5};
  const std::vector<double> occupancies{0.0, 0.0, 0.1, 9.9, 10.0, 10.0, 10.0, 10.0};

  const std::vector<double> placed = placePressuresByEqualOccupancy(pressures, occupancies);
  ASSERT_EQ(placed.size(), pressures.size());
  EXPECT_DOUBLE_EQ(placed.front(), pressures.front());
  EXPECT_DOUBLE_EQ(placed.back(), pressures.back());

  std::size_t inJump = 0uz;
  for (std::size_t i = 1; i + 1 < placed.size(); ++i)
  {
    if (placed[i] >= 1.0e-3 && placed[i] <= 1.0)
    {
      ++inJump;
    }
  }
  EXPECT_GE(inJump, 4uz);
}

TEST(REWEIGHTED_HISTOGRAM_WHAM, equal_occupancy_keeps_flat_ladder)
{
  const std::vector<double> pressures{1.0, 10.0, 100.0};
  const std::vector<double> occupancies{5.0, 5.0, 5.0};
  const std::vector<double> placed = placePressuresByEqualOccupancy(pressures, occupancies);
  ASSERT_EQ(placed, pressures);
}

TEST(REWEIGHTED_HISTOGRAM_WHAM, langmuir_equal_occupancy_keeps_endpoints)
{
  const double henry = 600.0;
  const double nSat = 23.0;
  const double lo = 1.0e-4;
  const double hi = 1.0e5;
  const std::vector<double> placed = placePressuresByLangmuirEqualOccupancy(lo, hi, 8uz, henry, nSat);
  ASSERT_EQ(placed.size(), 8uz);
  EXPECT_DOUBLE_EQ(placed.front(), lo);
  EXPECT_DOUBLE_EQ(placed.back(), hi);
  for (std::size_t i = 1; i < placed.size(); ++i)
  {
    EXPECT_GT(placed[i], placed[i - 1]);
  }
}

TEST(REWEIGHTED_HISTOGRAM_WHAM, langmuir_equal_occupancy_clusters_rungs_at_half_filling)
{
  const double henry = 600.0;
  const double nSat = 23.0;
  const double halfFilling = nSat / henry;
  const std::vector<double> placed =
      placePressuresByLangmuirEqualOccupancy(1.0e-4, 1.0e5, 16uz, henry, nSat);

  std::size_t nearStep = 0uz;
  for (std::size_t i = 1; i + 1 < placed.size(); ++i)
  {
    if (placed[i] >= 0.1 * halfFilling && placed[i] <= 10.0 * halfFilling)
    {
      ++nearStep;
    }
  }
  EXPECT_GE(nearStep, 3uz);

  const double maxLogGap =
      nitrogenBETMaxLogSpacingMultiplier * std::log(placed.back() / placed.front()) / 15.0;
  EXPECT_GT(placed[placed.size() - 2], std::exp(std::log(placed.back()) - maxLogGap * 1.000001));
}

TEST(REWEIGHTED_HISTOGRAM_WHAM, langmuir_log_gap_at_most_twice_equal_spacing)
{
  const std::vector<double> placed =
      placePressuresByLangmuirEqualOccupancy(1.0e-4, 1.0e5, 16uz, 600.0, 23.0);
  ASSERT_EQ(placed.size(), 16uz);

  const double maxLogGap =
      nitrogenBETMaxLogSpacingMultiplier * std::log(placed.back() / placed.front()) / 15.0;
  double largest = 0.0;
  for (std::size_t i = 1; i < placed.size(); ++i)
  {
    largest = std::max(largest, std::log(placed[i] / placed[i - 1]));
    EXPECT_LE(std::log(placed[i] / placed[i - 1]), maxLogGap * (1.0 + 1.0e-12));
  }
  EXPECT_GT(largest, 0.5 * maxLogGap);
  EXPECT_GT(placed[placed.size() - 2], 1.0e2);
}

TEST(REWEIGHTED_HISTOGRAM_WHAM, langmuir_equal_occupancy_falls_back_without_henry)
{
  const std::vector<double> placed = placePressuresByLangmuirEqualOccupancy(1.0e-2, 1.0e2, 5uz, 0.0, 20.0);
  ASSERT_EQ(placed.size(), 5uz);
  EXPECT_DOUBLE_EQ(placed.front(), 1.0e-2);
  EXPECT_DOUBLE_EQ(placed.back(), 1.0e2);
  EXPECT_NEAR(placed[2], 1.0, 1.0e-12);
}

TEST(REWEIGHTED_HISTOGRAM_WHAM, preisotherm_keeps_langmuir_when_probes_match)
{
  const double nSat = 23.0;
  const double henry = 600.0;
  const double b = henry / nSat;
  std::vector<SimulatedIsothermPoint> probes;
  for (double theta : nitrogenBETLangmuirProbeCoverages)
  {
    const double P = langmuirPressureAtCoverage(theta, b);
    probes.push_back(SimulatedIsothermPoint{P, nSat * theta});
  }

  const NitrogenBETPreIsothermFit fit = fitNitrogenBETPreIsotherm(henry, nSat, 90.0, probes);
  EXPECT_EQ(fit.model, NitrogenBETIsothermModel::Langmuir);
  EXPECT_FALSE(fit.langmuirRejected);
  EXPECT_TRUE(fit.hillPlotFlat);
  EXPECT_NEAR(fit.exponent, 1.0, 1.0e-8);
}

TEST(REWEIGHTED_HISTOGRAM_WHAM, preisotherm_selects_sips_when_m_is_not_one)
{
  const double nSat = 23.4;
  const double b = 8.5;
  const double m = 0.63;
  const double henry = 600.0;
  std::vector<SimulatedIsothermPoint> probes;
  for (double thetaL : nitrogenBETLangmuirProbeCoverages)
  {
    const double P = langmuirPressureAtCoverage(thetaL, henry / 22.6);
    const double x = std::pow(b * P, m);
    probes.push_back(SimulatedIsothermPoint{P, nSat * x / (1.0 + x)});
  }

  const NitrogenBETPreIsothermFit fit = fitNitrogenBETPreIsotherm(henry, 22.6, 90.0, probes);
  EXPECT_EQ(fit.model, NitrogenBETIsothermModel::LangmuirFreundlich);
  EXPECT_TRUE(fit.langmuirRejected);
  EXPECT_TRUE(fit.hillPlotFlat);
  EXPECT_NEAR(fit.exponent, m, 0.08);
}

TEST(REWEIGHTED_HISTOGRAM_WHAM, preisotherm_selects_toth_when_hill_plot_drifts)
{
  const double n1 = 16.0;
  const double n2 = 7.0;
  const double nSat = n1 + n2;
  const double henry = 600.0;
  const double b1 = henry / n1;
  const double b2 = 1.0e-4;
  std::vector<SimulatedIsothermPoint> probes;
  for (double thetaL : nitrogenBETLangmuirProbeCoverages)
  {
    const double P = langmuirPressureAtCoverage(thetaL, henry / nSat);
    const double n = n1 * (b1 * P) / (1.0 + b1 * P) + n2 * (b2 * P) / (1.0 + b2 * P);
    probes.push_back(SimulatedIsothermPoint{P, n});
  }

  const NitrogenBETPreIsothermFit fit = fitNitrogenBETPreIsotherm(henry, nSat, 90.0, probes);
  EXPECT_EQ(fit.model, NitrogenBETIsothermModel::Toth);
  EXPECT_FALSE(fit.hillPlotFlat);
  EXPECT_LT(fit.exponent, 0.85);
}

TEST(REWEIGHTED_HISTOGRAM_WHAM, sips_fisher_shifts_interior_rungs_up)
{
  NitrogenBETPreIsothermFit langmuir;
  langmuir.model = NitrogenBETIsothermModel::Langmuir;
  langmuir.saturationLoadingPerCell = 22.6;
  langmuir.affinity = 600.0 / 22.6;
  langmuir.exponent = 1.0;

  NitrogenBETPreIsothermFit sips;
  sips.model = NitrogenBETIsothermModel::LangmuirFreundlich;
  sips.saturationLoadingPerCell = 23.4;
  sips.affinity = 8.5;
  sips.exponent = 0.63;

  const std::vector<double> langmuirPlaced =
      placePressuresByPreIsothermFisher(1.0e-4, 1.0e5, 16uz, langmuir);
  const std::vector<double> sipsPlaced = placePressuresByPreIsothermFisher(1.0e-4, 1.0e5, 16uz, sips);
  ASSERT_EQ(langmuirPlaced.size(), 16uz);
  ASSERT_EQ(sipsPlaced.size(), 16uz);

  double langmuirLogSum = 0.0;
  double sipsLogSum = 0.0;
  for (std::size_t i = 1; i + 1 < langmuirPlaced.size(); ++i)
  {
    langmuirLogSum += std::log(langmuirPlaced[i]);
    sipsLogSum += std::log(sipsPlaced[i]);
  }
  EXPECT_GT(sipsLogSum, langmuirLogSum);

  std::size_t sipsAboveOne = 0uz;
  std::size_t langmuirAboveOne = 0uz;
  for (std::size_t i = 1; i + 1 < sipsPlaced.size(); ++i)
  {
    if (sipsPlaced[i] >= 1.0)
    {
      ++sipsAboveOne;
    }
    if (langmuirPlaced[i] >= 1.0)
    {
      ++langmuirAboveOne;
    }
  }
  EXPECT_GE(sipsAboveOne, langmuirAboveOne);
}
