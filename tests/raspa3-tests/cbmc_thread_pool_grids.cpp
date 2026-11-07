#include <gtest/gtest.h>

import std;

import int3;
import uint3;
import double3;
import atom;
import forcefield;
import framework;
import component;
import system;
import running_energy;
import threadpool;
import cbmc;
import cbmc_external_energy;

// Under the thread pool the CBMC framework energy of a context without a framework cell list (a box too small
// for one, say) went to the threaded pair kernel, which does not read the interpolation grids; the serial path,
// the full evaluation and the moves read them. A grid-reading atom then got another framework energy under the
// thread pool than without it, and the running energies of the CBMC moves drifted from a full recomputation.

namespace
{

using Pool = ThreadPool::ThreadPool<ThreadPool::details::default_function_type, std::jthread>;

class ScopedThreadingType
{
 public:
  explicit ScopedThreadingType(ThreadPool::ThreadingType type) : previous_(Pool::instance().getThreadingType())
  {
    Pool::instance().init(4, type);
  }
  ~ScopedThreadingType() { Pool::instance().init(4, previous_); }

 private:
  ThreadPool::ThreadingType previous_;
};

}  // namespace

TEST(CBMC_THREAD_POOL_GRIDS, a_grid_reading_atom_gets_the_grid_energy_under_the_thread_pool)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  forceField.settings.gridPseudoAtomIndices = {2};  // methane
  forceField.settings.numberOfVDWGridPoints = uint3(16, 16, 16);
  forceField.settings.numberOfCoulombGridPoints = uint3(16, 16, 16);
  forceField.settings.numberOfGridTestPoints = 10;
  const Framework framework = Framework::makeFAU(forceField, int3(1, 1, 1));
  const Component methane = Component(forceField, "methane", 190.564, 4599200, 0.01142,
                                      {Atom({0, 0, 0}, 0.0, 1.0, 0, 2, 0, false, false)}, {}, {}, 5, 21);
  System system(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {framework}, {methane}, {}, {0}, 5);
  std::ostringstream log;
  system.createFrameworkInterpolationGrids(log);  // the drivers build the grids explicitly
  ASSERT_TRUE(system.interpolationGrids[2].has_value());

  CBMC::GrowContext context = system.makeGrowContext();
  context.frameworkCellList = nullptr;  // the path without a framework cell list
  const Atom probe(double3(3.1, 7.4, 11.9), 0.0, 1.0, 0, 2, 0, false, false);

  std::optional<RunningEnergy> serial;
  std::optional<RunningEnergy> threaded;
  {
    ScopedThreadingType threading(ThreadPool::ThreadingType::Serial);
    serial = CBMC::computeFrameworkMoleculeEnergy(context, std::span<const Atom>(&probe, 1));
  }
  {
    ScopedThreadingType threading(ThreadPool::ThreadingType::ThreadPool);
    threaded = CBMC::computeFrameworkMoleculeEnergy(context, std::span<const Atom>(&probe, 1));
  }
  ASSERT_TRUE(serial.has_value());
  ASSERT_TRUE(threaded.has_value());
  EXPECT_NE(serial->frameworkMoleculeVDW, 0.0);
  EXPECT_DOUBLE_EQ(threaded->frameworkMoleculeVDW, serial->frameworkMoleculeVDW);
  EXPECT_DOUBLE_EQ(threaded->potentialEnergy(), serial->potentialEnergy());
}
