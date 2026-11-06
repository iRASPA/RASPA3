#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import atom;
import component;
import forcefield;
import parallel_tmmc;
import simulationbox;
import system;

// The window layout of the parallel transition-matrix driver when the box cannot be packed as full as the
// macrostate range asks. The range then ends at the packing, the windows above it are dropped, and what is
// left has to be one consistent layout: every walker inside the range, every window with room to move in,
// and the same windows at every temperature.

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

// Particles ten Ångström across in a box of twenty-five: thirty or so can be placed, fewer than any of the
// macrostate ranges asked for here runs to.
System makeCrowdedSystem(std::size_t largestMacrostate)
{
  const ForceField forceField({{"X", false, 1.0, 0.0, 0.0, 1, false}}, {{100.0, 10.0}},
                              ForceField::MixingRule::Lorentz_Berthelot, 12.0, 12.0, 12.0, true, false, false);
  const Component particle(forceField, "particle", 1.0, 1.0, 0.0,
                           {Atom({0.0, 0.0, 0.0}, 0.0, 1.0, 0, 0, 0, false, false)}, {}, {}, 5, 21);
  System system(forceField, SimulationBox(25.0, 25.0, 25.0), false, 300.0, 1e4, 1.0, {}, {particle}, {}, {0}, 5);
  system.components[0].fugacityCoefficient = 1.0;
  system.tmmc.minMacrostate = 0;
  system.tmmc.maxMacrostate = largestMacrostate;
  return system;
}

}  // namespace


TEST(PARALLEL_TMMC_WINDOWS, a_range_cut_short_by_the_packing_leaves_one_consistent_layout)
{
  ScopedCurrentPath scratch(std::filesystem::temp_directory_path() /
                            std::format("raspa3_test_{}_parallel_tmmc_windows", test_support::uniqueSuffix()));

  // The packing falls in the first window of the widest range, and two or three windows up in the others,
  // where the windows below the last have upper bounds of their own to be brought down.
  for (auto [largestMacrostate, seed] : std::vector<std::pair<std::size_t, std::size_t>>{
           {800uz, 11uz}, {800uz, 12uz}, {160uz, 13uz}, {160uz, 14uz}, {120uz, 15uz}, {120uz, 16uz}, {120uz, 17uz}})
  {
    ParallelTMMCParameters parameters;
    parameters.numberOfWindows = 8;
    parameters.randomSeed = seed;
    ParallelTMMC simulation(makeCrowdedSystem(largestMacrostate), {250.0, 300.0, 350.0}, parameters);
    ASSERT_EQ(simulation.maxMacrostate, largestMacrostate);
    ASSERT_EQ(simulation.numberOfWindows, 8uz);
    const std::vector<std::size_t> laidOut = simulation.windowBoundaries;

    simulation.setup();

    // The range was cut short, and what is left of the partition is what the packing leaves of it.
    EXPECT_LT(simulation.maxMacrostate, largestMacrostate / 2uz) << "seed " << seed;
    EXPECT_EQ(simulation.windowBoundaries, tmmcWindowBoundariesWithinPacking(laidOut, simulation.maxMacrostate))
        << "seed " << seed;
    EXPECT_GT(simulation.maxMacrostate, simulation.minMacrostate) << "seed " << seed;
    ASSERT_LT(simulation.numberOfWindows, 8uz) << "seed " << seed;
    ASSERT_GE(simulation.numberOfWindows, 1uz) << "seed " << seed;

    // The partition is of the range that is left, each part of it with something in it.
    ASSERT_EQ(simulation.windowBoundaries.size(), simulation.numberOfWindows + 1uz) << "seed " << seed;
    EXPECT_EQ(simulation.windowBoundaries.front(), simulation.minMacrostate) << "seed " << seed;
    EXPECT_EQ(simulation.windowBoundaries.back(), simulation.maxMacrostate) << "seed " << seed;
    for (std::size_t window = 0; window < simulation.numberOfWindows; ++window)
    {
      EXPECT_LT(simulation.windowBoundaries[window], simulation.windowBoundaries[window + 1uz])
          << "seed " << seed << ", window " << window;
    }

    // One walker to a temperature and a window, and no others.
    ASSERT_EQ(simulation.numberOfWalkers, 3uz * simulation.numberOfWindows) << "seed " << seed;
    ASSERT_EQ(simulation.systems.size(), simulation.numberOfWalkers) << "seed " << seed;

    for (std::size_t temperature = 0; temperature < 3uz; ++temperature)
    {
      for (std::size_t window = 0; window < simulation.numberOfWindows; ++window)
      {
        const System& walker = simulation.systems[simulation.walkerIndex(temperature, window)];
        const std::string where = std::format("seed {}, temperature {}, window {}", seed, temperature, window);

        // Inside the range, with more than the one state, and holding its segment of the partition.
        EXPECT_GE(walker.tmmc.minMacrostate, simulation.minMacrostate) << where;
        EXPECT_LE(walker.tmmc.maxMacrostate, simulation.maxMacrostate) << where;
        EXPECT_LT(walker.tmmc.minMacrostate, walker.tmmc.maxMacrostate) << where;
        EXPECT_LE(walker.tmmc.minMacrostate, simulation.windowBoundaries[window]) << where;
        EXPECT_GE(walker.tmmc.maxMacrostate, simulation.windowBoundaries[window + 1uz]) << where;

        // The walker is where its window is.
        const std::size_t molecules = walker.numberOfIntegerMoleculesPerComponent[0];
        EXPECT_GE(molecules, walker.tmmc.minMacrostate) << where;
        EXPECT_LE(molecules, walker.tmmc.maxMacrostate) << where;

        // The same window at every temperature.
        const System& first = simulation.systems[simulation.walkerIndex(0uz, window)];
        EXPECT_EQ(walker.tmmc.minMacrostate, first.tmmc.minMacrostate) << where;
        EXPECT_EQ(walker.tmmc.maxMacrostate, first.tmmc.maxMacrostate) << where;
      }
    }
  }
}


// What the packing leaves of a partition. A window is kept for the part of the range that is its own, the
// segment between its boundary and the next, and not for the overlap it has with the window below: one whose
// segment starts at the packing or above it has nothing of its own left, and closed at the packing it would
// be a window of one macrostate, or of a few that the window below holds already.
TEST(PARALLEL_TMMC_WINDOWS, a_window_with_nothing_of_its_own_under_the_packing_is_dropped)
{
  const std::vector<std::size_t> boundaries{0, 100, 200, 300, 400};
  using Boundaries = std::vector<std::size_t>;

  // well inside the second window
  EXPECT_EQ(tmmcWindowBoundariesWithinPacking(boundaries, 150uz), (Boundaries{0, 100, 150}));
  // one macrostate into it
  EXPECT_EQ(tmmcWindowBoundariesWithinPacking(boundaries, 101uz), (Boundaries{0, 100, 101}));
  // on its boundary: the second window would hold that one macrostate and no other of its own
  EXPECT_EQ(tmmcWindowBoundariesWithinPacking(boundaries, 100uz), (Boundaries{0, 100}));
  // under its boundary, where only its overlap with the first reaches
  EXPECT_EQ(tmmcWindowBoundariesWithinPacking(boundaries, 80uz), (Boundaries{0, 80}));
  // in the last window but one, and in the first
  EXPECT_EQ(tmmcWindowBoundariesWithinPacking(boundaries, 399uz), (Boundaries{0, 100, 200, 300, 399}));
  EXPECT_EQ(tmmcWindowBoundariesWithinPacking(boundaries, 1uz), (Boundaries{0, 1}));

  // what is left is a partition still: it starts where it did, ends at the packing, and rises all the way
  for (std::size_t packed = 1uz; packed < 400uz; ++packed)
  {
    const Boundaries kept = tmmcWindowBoundariesWithinPacking(boundaries, packed);
    ASSERT_GE(kept.size(), 2uz);
    EXPECT_EQ(kept.front(), 0uz);
    EXPECT_EQ(kept.back(), packed);
    for (std::size_t index = 0; index + 1uz < kept.size(); ++index) EXPECT_LT(kept[index], kept[index + 1uz]);
  }
}
