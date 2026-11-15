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

// The statistics file written at the end of production held the ln Pi of the last rebuild, at a multiple of
// 'TMMCUpdateEvery' steps (a million by default): all zeros for a shorter run, next to the final collection matrix.
// And the collection matrix carried the transitions of initialization and equilibration into production.
TEST(TMMC_PROPERTY_OUTPUT, statistics_file_holds_the_final_lnpi_of_production_only)
{
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

  RandomNumber random(7);
  const SimulationSchedule schedule{50, 0, 40, 40, 1000, 0, 1000, 1000};
  MonteCarloTransitionMatrix tmmc(schedule, systems, random, 5);
  tmmc.run();

  std::ifstream file("tmmc/tmmc_statistics.txt");
  ASSERT_TRUE(file.good());
  std::size_t performed = 0uz;
  std::size_t visits = 0uz;
  bool anyLnPi = false;
  for (std::string line; std::getline(file, line);)
  {
    if (line.starts_with("# performed:"))
    {
      std::istringstream stream(line.substr(12));
      stream >> performed;
      continue;
    }
    if (line.empty() || line[0] == '#' || line[0] == 'N') continue;
    std::istringstream stream(line);
    std::size_t N, histogram;
    double down, stay, up, bias, lnpi, forward, reverse;
    stream >> N >> down >> stay >> up >> bias >> lnpi >> forward >> reverse >> histogram;
    ASSERT_FALSE(stream.fail()) << line;
    visits += histogram;
    anyLnPi = anyLnPi || lnpi != 0.0;
  }
  // the ln Pi of the collection matrix as it stands at the end, not the zeros of before the first rebuild
  EXPECT_TRUE(anyLnPi);
  // every visit in the file is a production step: the earlier stages were dropped at the start of production
  EXPECT_GT(performed, 0uz);
  EXPECT_EQ(visits, performed);
}
