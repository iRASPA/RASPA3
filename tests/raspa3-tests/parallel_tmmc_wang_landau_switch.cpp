#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import atom;
import forcefield;
import component;
import system;
import simulationbox;
import parallel_tmmc;
import mc_moves_move_types;
import mc_moves_probabilities;

// Production of the parallel transition-matrix driver starts with the statistics emptied and Wang-Landau
// still driving the walk; after a quarter of production the bias is rebuilt as -ln Pi from the production
// histogram and Wang-Landau retires. A short production (fewer than four cycles) put that switch at cycle
// zero, on the just-emptied histogram: nothing could be rebuilt, Wang-Landau was retired anyway, and its
// un-annealed bias stayed frozen. The switch has to wait for at least one cycle of statistics.

namespace
{

class ScopedCurrentPath
{
 public:
  explicit ScopedCurrentPath(const std::filesystem::path& path) : original_(std::filesystem::current_path())
  {
    std::filesystem::create_directories(path);
    std::filesystem::current_path(path);
  }

  ~ScopedCurrentPath()
  {
    std::error_code ignored;
    const std::filesystem::path used = std::filesystem::current_path();
    std::filesystem::current_path(original_);
    std::filesystem::remove_all(used, ignored);
  }

 private:
  std::filesystem::path original_;
};

System makeSwappingSystem()
{
  const ForceField forceField({{"CH4", false, 16.04246, 0.0, 0.0, 6, false}}, {{158.5, 3.72}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  MCMoveProbabilities probabilities;
  probabilities.setProbability(Move::Types::Translation, 1.0);
  probabilities.setProbability(Move::Types::Swap, 1.0);
  const Component methane(forceField, "methane", 190.564, 4599200, 0.01142,
                          {Atom({0, 0, 0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21, probabilities,
                          std::nullopt, false);
  System system(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e5, 1.0, {}, {methane}, {}, {4}, 5);
  system.tmmc.minMacrostate = 0;
  system.tmmc.maxMacrostate = 12;
  return system;
}

}  // namespace

TEST(PARALLEL_TMMC_WANG_LANDAU_SWITCH, a_short_production_still_hands_the_bias_to_the_transition_matrix)
{
  ScopedCurrentPath scratch(std::filesystem::temp_directory_path() /
                            std::format("raspa3_test_{}_parallel_tmmc_switch", test_support::uniqueSuffix()));

  ParallelTMMCParameters parameters;
  parameters.numberOfPreInitializationCycles = 0;
  parameters.numberOfInitializationCycles = 2;
  parameters.numberOfEquilibrationCycles = 2;
  parameters.numberOfProductionCycles = 2;
  parameters.numberOfWindows = 2;
  parameters.numberOfBlocks = 2;
  parameters.randomSeed = 29;

  ParallelTMMC simulation(makeSwappingSystem(), {300.0}, parameters);
  for (const System& walker : simulation.systems) ASSERT_TRUE(walker.tmmc.useWangLandau);

  simulation.run();

  for (std::size_t walkerId = 0; walkerId < simulation.systems.size(); ++walkerId)
  {
    const System& walker = simulation.systems[walkerId];
    EXPECT_FALSE(walker.tmmc.useWangLandau) << "walker " << walkerId << " never retired Wang-Landau";
    // the bias was rebuilt from a histogram with visits in it: it is -ln Pi, not the frozen Wang-Landau bias
    ASSERT_EQ(walker.tmmc.bias.size(), walker.tmmc.lnpi.size());
    for (std::size_t i = 0; i < walker.tmmc.bias.size(); ++i)
    {
      EXPECT_DOUBLE_EQ(walker.tmmc.bias[i], -walker.tmmc.lnpi[i]) << "walker " << walkerId << ", state " << i;
    }
  }
}
