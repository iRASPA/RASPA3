#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import randomnumbers;
import atom;
import json;
import forcefield;
import component;
import system;
import simulationbox;
import monte_carlo_transition_matrix;
import simulation_schedule;
import mc_moves_move_types;
import mc_moves_probabilities;
import property_number_of_molecules_histogram;

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
}  // namespace

TEST(TMMC_PROPERTY_OUTPUT, sampled_properties_and_averages_are_written)
{
  // The transition-matrix driver samples the properties of its systems every production cycle; the files of those
  // properties and the averages (text and JSON) must be written as in the Monte Carlo driver.
  TemporaryDirectory workspace;
  ScopedCurrentPath currentPath(workspace.path());

  const ForceField forceField({{"CH4", false, 16.04246, 0.0, 0.0, 6, false}}, {{158.5, 3.72}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  MCMoveProbabilities moves;
  moves.setProbability(Move::Types::Translation, 1.0);
  moves.setProbability(Move::Types::Swap, 1.0);
  const Component methane(forceField, "methane", 190.564, 4599200, 0.01142,
                          {Atom({0, 0, 0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21, moves, std::nullopt,
                          false);
  std::vector<System> systems{
      System(forceField, SimulationBox(30.0, 30.0, 30.0), false, 300.0, 1e4, 1.0, {}, {methane}, {}, {2}, 5)};
  System& system = systems[0];
  system.tmmc.doTMMC = true;
  system.tmmc.useBias = true;
  system.tmmc.useTMBias = true;
  system.tmmc.rejectOutOfBound = true;
  system.tmmc.minMacrostate = 0;
  system.tmmc.maxMacrostate = 10;
  system.tmmc.initialize();
  system.averageNumberOfMoleculesHistogram = PropertyNumberOfMoleculesHistogram(5, 1, {0, 11}, 1, 5);

  RandomNumber random(7);
  const SimulationSchedule schedule{10, 0, 0, 0, 5, 1000, 1000, 1000};
  MonteCarloTransitionMatrix tmmc(schedule, systems, random, 5);
  tmmc.run();

  EXPECT_TRUE(std::filesystem::exists("number_of_molecules_histogram/number_of_molecules_histogram_histogram.s0.txt"));

  ASSERT_FALSE(tmmc.outputJsonFileNames.empty());
  std::ifstream json(tmmc.outputJsonFileNames[0]);
  const nlohmann::json output = nlohmann::json::parse(json);
  ASSERT_TRUE(output.contains("properties"));
  EXPECT_TRUE(output["properties"].contains("averageEnergies"));
  EXPECT_TRUE(output["properties"].contains("averagePressure"));
}
