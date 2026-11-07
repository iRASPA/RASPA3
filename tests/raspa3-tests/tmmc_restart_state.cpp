#include <gtest/gtest.h>

import std;

import archive;
import atom;
import forcefield;
import component;
import system;
import simulationbox;
import parallel_tmmc;
import reweighted_histogram;
import isotherm_bet;
import mc_moves_move_types;
import mc_moves_probabilities;

// The BET auto-setup state of the transition-matrix and reweighted-histogram drivers (the auto flags, the
// scouted pressure plan and the filling ceiling) is produced by scouts that are not rerun on a binary restart,
// so it has to come back from the restart file.

namespace
{

System makeSystem()
{
  const ForceField forceField({{"CH4", false, 16.04246, 0.0, 0.0, 6, false}}, {{158.5, 3.72}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Swap, 1.0);
  const Component methane(forceField, "methane", 190.564, 4599200, 0.01142,
                          {Atom({0, 0, 0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21, probabilities,
                          std::nullopt, false);
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {methane}, {}, {5}, 5);
  system.tmmc.minMacrostate = 0;
  system.tmmc.maxMacrostate = 20;
  return system;
}

NitrogenBETPressurePlan makePlan()
{
  NitrogenBETPressurePlan plan;
  plan.henryCoefficientPerCell = 2.5e-3;
  plan.packingCapacity = 48.0;
  plan.lowestPressure = 0.75;
  plan.highestPressure = 101325.0;
  plan.samplingPressures = {0.75, 10.0, 1000.0, 101325.0};
  plan.whamReweightingNumberOfPressures = 321;
  plan.tmmcReweightingNumberOfPressures = 123;
  plan.probe = {16.2, 57.7, 101325.0};
  return plan;
}

NitrogenBETFillingCeiling makeCeiling()
{
  NitrogenBETFillingCeiling ceiling;
  ceiling.maxMacrostate = 44;
  ceiling.meanOccupancy = 39.5;
  ceiling.standardDeviation = 1.25;
  ceiling.peakOccupancy = 43;
  ceiling.cycles = 777;
  ceiling.converged = true;
  ceiling.gurvichCap = 48;
  ceiling.saturationPressure = 101325.0;
  return ceiling;
}

template <typename Driver>
void roundTrip(const Driver& original, Driver& restored, std::string_view name)
{
  const std::filesystem::path path = std::filesystem::temp_directory_path() / std::format("raspa3_{}.bin", name);
  {
    std::ofstream stream(path, std::ios::binary);
    Archive<std::ofstream> archive(stream);
    archive << original;
  }
  {
    std::ifstream stream(path, std::ios::binary);
    Archive<std::ifstream> archive(stream);
    archive >> restored;
  }
  std::filesystem::remove(path);
}

void expectPlan(const std::optional<NitrogenBETPressurePlan>& plan)
{
  ASSERT_TRUE(plan.has_value());
  const NitrogenBETPressurePlan expected = makePlan();
  EXPECT_DOUBLE_EQ(plan->henryCoefficientPerCell, expected.henryCoefficientPerCell);
  EXPECT_DOUBLE_EQ(plan->packingCapacity, expected.packingCapacity);
  EXPECT_DOUBLE_EQ(plan->lowestPressure, expected.lowestPressure);
  EXPECT_DOUBLE_EQ(plan->highestPressure, expected.highestPressure);
  EXPECT_EQ(plan->samplingPressures, expected.samplingPressures);
  EXPECT_EQ(plan->whamReweightingNumberOfPressures, expected.whamReweightingNumberOfPressures);
  EXPECT_EQ(plan->tmmcReweightingNumberOfPressures, expected.tmmcReweightingNumberOfPressures);
  EXPECT_DOUBLE_EQ(plan->probe.crossSection, expected.probe.crossSection);
  EXPECT_DOUBLE_EQ(plan->probe.liquidVolume, expected.probe.liquidVolume);
  EXPECT_DOUBLE_EQ(plan->probe.saturationPressure, expected.probe.saturationPressure);
}

void expectCeiling(const std::optional<NitrogenBETFillingCeiling>& ceiling)
{
  ASSERT_TRUE(ceiling.has_value());
  const NitrogenBETFillingCeiling expected = makeCeiling();
  EXPECT_EQ(ceiling->maxMacrostate, expected.maxMacrostate);
  EXPECT_DOUBLE_EQ(ceiling->meanOccupancy, expected.meanOccupancy);
  EXPECT_DOUBLE_EQ(ceiling->standardDeviation, expected.standardDeviation);
  EXPECT_EQ(ceiling->peakOccupancy, expected.peakOccupancy);
  EXPECT_EQ(ceiling->cycles, expected.cycles);
  EXPECT_EQ(ceiling->converged, expected.converged);
  EXPECT_EQ(ceiling->gurvichCap, expected.gurvichCap);
  EXPECT_DOUBLE_EQ(ceiling->saturationPressure, expected.saturationPressure);
}

}  // namespace

TEST(TMMC_RESTART_STATE, parallel_tmmc_restores_the_bet_auto_setup_state)
{
  ParallelTMMCParameters parameters;
  parameters.numberOfWindows = 2;
  parameters.randomSeed = 3;

  ParallelTMMC original(makeSystem(), {300.0}, parameters);
  original.computeBET = true;
  original.autoReweightingPressureRange = true;
  original.autoMacroStateMaximum = true;
  original.nitrogenBETPressurePlan = makePlan();
  original.nitrogenBETFillingCeiling = makeCeiling();

  ParallelTMMC restored(makeSystem(), {300.0}, parameters);
  ASSERT_FALSE(restored.computeBET);
  ASSERT_FALSE(restored.nitrogenBETPressurePlan.has_value());
  roundTrip(original, restored, "parallel_tmmc_restart_state");

  EXPECT_TRUE(restored.computeBET);
  EXPECT_TRUE(restored.autoReweightingPressureRange);
  EXPECT_TRUE(restored.autoMacroStateMaximum);
  expectPlan(restored.nitrogenBETPressurePlan);
  expectCeiling(restored.nitrogenBETFillingCeiling);
  EXPECT_EQ(restored.systems.size(), original.systems.size());
}

TEST(TMMC_RESTART_STATE, reweighted_histogram_restores_the_bet_auto_setup_state)
{
  ReweightedHistogramParameters parameters;
  parameters.randomSeed = 5;

  ReweightedHistogram original(makeSystem(), {300.0}, {1e4, 1e5}, parameters);
  original.computeBET = true;
  original.autoExternalPressures = true;
  original.autoReweightingPressureRange = true;
  original.autoMacroStateMaximum = true;
  original.nitrogenBETPressurePlan = makePlan();
  original.nitrogenBETFillingCeiling = makeCeiling();

  ReweightedHistogram restored(makeSystem(), {300.0}, {1e4, 1e5}, parameters);
  ASSERT_FALSE(restored.autoExternalPressures);
  ASSERT_FALSE(restored.nitrogenBETFillingCeiling.has_value());
  roundTrip(original, restored, "reweighted_histogram_restart_state");

  EXPECT_TRUE(restored.computeBET);
  EXPECT_TRUE(restored.autoExternalPressures);
  EXPECT_TRUE(restored.autoReweightingPressureRange);
  EXPECT_TRUE(restored.autoMacroStateMaximum);
  expectPlan(restored.nitrogenBETPressurePlan);
  expectCeiling(restored.nitrogenBETFillingCeiling);
  EXPECT_EQ(restored.systems.size(), original.systems.size());
}
