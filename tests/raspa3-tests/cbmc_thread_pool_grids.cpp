#include <gtest/gtest.h>

#include "../test_support.hpp"

import std;

import int3;
import uint3;
import double3;
import atom;
import forcefield;
import vdwparameters;
import framework;
import interpolation_energy_grid;
import component;
import system;
import running_energy;
import threadpool;
import cbmc;
import cbmc_external_energy;
import interactions_framework_molecule;

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

// The Coulomb grid tabulates the real-space Ewald energy of a unit charge, and the full evaluation and the moves read
// no grid under another charge method (the grids are built all the same). The CBMC framework energy must follow the
// same rule, with and without a framework cell list, or the growth sees another framework energy than the rest of
// the code and the running energies of its moves drift from a recomputation.
TEST(CBMC_THREAD_POOL_GRIDS, no_grid_is_read_under_a_charge_method_other_than_ewald)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  const std::size_t sodium = *forceField.findPseudoAtom("Na+");
  forceField.settings.gridPseudoAtomIndices = {sodium};
  forceField.settings.numberOfVDWGridPoints = uint3(16, 16, 16);
  forceField.settings.numberOfCoulombGridPoints = uint3(16, 16, 16);
  forceField.settings.numberOfGridTestPoints = 10;
  const Framework framework = Framework::makeFAU(forceField, int3(1, 1, 1));
  const Component ion =
      Component(forceField, "sodium", 190.564, 4599200, 0.01142,
                {Atom({0, 0, 0}, 1.0, 1.0, 0, static_cast<std::uint16_t>(sodium), 0, false, false)}, {}, {}, 5, 21);
  System system(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {framework}, {ion}, {}, {0}, 5);
  std::ostringstream log;
  system.createFrameworkInterpolationGrids(log);
  ASSERT_TRUE(system.interpolationGrids[sodium].has_value());
  ASSERT_TRUE(system.interpolationGrids.back().has_value());

  const Atom probe(double3(3.1, 7.4, 11.9), 1.0, 1.0, 0, static_cast<std::uint16_t>(sodium), 0, false, false);
  const std::span<const Atom> probeSpan(&probe, 1);
  for (const ForceField::ChargeMethod method : {ForceField::ChargeMethod::Ewald, ForceField::ChargeMethod::Wolf})
  {
    SCOPED_TRACE(method == ForceField::ChargeMethod::Ewald ? "Ewald" : "Wolf");
    system.forceField.chargeMethod = method;
    const RunningEnergy expected =
        Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, system.interpolationGrids,
                                                     system.framework, system.spanOfFrameworkAtoms(), probeSpan);
    EXPECT_NE(expected.frameworkMoleculeCharge, 0.0);
    for (const bool withCellList : {true, false})
    {
      SCOPED_TRACE(withCellList ? "with a framework cell list" : "without a framework cell list");
      CBMC::GrowContext context = system.makeGrowContext();
      if (!withCellList) context.frameworkCellList = nullptr;
      const std::optional<RunningEnergy> energy = CBMC::computeFrameworkMoleculeEnergy(context, probeSpan);
      ASSERT_TRUE(energy.has_value());
      EXPECT_NEAR(energy->frameworkMoleculeVDW, expected.frameworkMoleculeVDW,
                  1e-6 * std::abs(expected.frameworkMoleculeVDW) + 1e-8);
      EXPECT_NEAR(energy->frameworkMoleculeCharge, expected.frameworkMoleculeCharge,
                  1e-6 * std::abs(expected.frameworkMoleculeCharge) + 1e-8);
    }
  }
}

// A framework interpolation grid holds the framework as it was built. A flexible framework moves (hybrid MC, MD), and
// the gradient evaluates its pairs then; the energies and the energy differences read the frozen grid instead, so
// the energy of a moved framework was that of the original one. With mobile framework atoms no routine reads a grid.
TEST(CBMC_THREAD_POOL_GRIDS, no_grid_is_read_for_a_framework_with_mobile_atoms)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  const std::size_t methaneType = *forceField.findPseudoAtom("CH4");
  forceField.settings.gridPseudoAtomIndices = {methaneType};
  forceField.settings.numberOfVDWGridPoints = uint3(16, 16, 16);
  forceField.settings.numberOfCoulombGridPoints = uint3(16, 16, 16);
  forceField.settings.numberOfGridTestPoints = 10;
  Framework framework = Framework::makeFAU(forceField, int3(1, 1, 1));
  TemporaryFile definition("raspa3-fau-flexible-grid.json", R"({"Type": "Flexible"})");
  framework.readFrameworkDefinition(forceField, definition.path().string());
  ASSERT_TRUE(framework.hasMobileAtoms());
  const Component methane = Component(
      forceField, "methane", 190.564, 4599200, 0.01142,
      {Atom({0, 0, 0}, 0.0, 1.0, 0, static_cast<std::uint16_t>(methaneType), 0, false, false)}, {}, {}, 5, 21);
  System system(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {framework}, {methane}, {}, {0}, 5);
  std::ostringstream log;
  system.createFrameworkInterpolationGrids(log);
  ASSERT_TRUE(system.interpolationGrids[methaneType].has_value());
  const std::vector<std::optional<InterpolationEnergyGrid>> noGrids(system.interpolationGrids.size(), std::nullopt);

  const Atom probe(double3(3.1, 7.4, 11.9), 0.0, 1.0, 0, static_cast<std::uint16_t>(methaneType), 0, false, false);
  const std::span<const Atom> probeSpan(&probe, 1);
  auto withGrids = [&]
  {
    return Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox,
                                                        system.interpolationGrids, system.framework,
                                                        system.spanOfFrameworkAtoms(), probeSpan);
  };
  auto byPairs = [&]
  {
    return Interactions::computeFrameworkMoleculeEnergy(system.forceField, system.simulationBox, noGrids,
                                                        system.framework, system.spanOfFrameworkAtoms(), probeSpan);
  };
  const double before = byPairs().frameworkMoleculeVDW;

  // the framework moves: every atom by 0.3 A along x
  for (Atom& atom : system.spanOfFrameworkAtoms()) atom.position.x += 0.3;
  const double moved = byPairs().frameworkMoleculeVDW;
  ASSERT_GT(std::abs(moved - before), 1e-3 * std::abs(before));
  EXPECT_NEAR(withGrids().frameworkMoleculeVDW, moved, 1e-10 * std::abs(moved));

  // the energy differences and the CBMC growth see the moved framework as well
  const std::optional<RunningEnergy> difference = Interactions::computeFrameworkMoleculeEnergyDifference(
      system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
      system.spanOfFrameworkAtoms(), probeSpan, {});
  ASSERT_TRUE(difference.has_value());
  EXPECT_NEAR(difference->frameworkMoleculeVDW, moved, 1e-10 * std::abs(moved));
  CBMC::GrowContext context = system.makeGrowContext();
  const std::optional<RunningEnergy> grown = CBMC::computeFrameworkMoleculeEnergy(context, probeSpan);
  ASSERT_TRUE(grown.has_value());
  EXPECT_NEAR(grown->frameworkMoleculeVDW, moved, 1e-10 * std::abs(moved));
}

// The VDW interpolation grid tabulates the 12-6 Lennard-Jones form: a framework pair of another form was gridded as
// the Lennard-Jones potential of its first two parameters, while the pairwise routes evaluated the real potential.
TEST(CBMC_THREAD_POOL_GRIDS, a_non_lennard_jones_framework_pair_is_not_gridded)
{
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  forceField.settings.gridPseudoAtomIndices = {2};  // methane
  forceField.settings.numberOfVDWGridPoints = uint3(8, 8, 8);
  forceField.settings.numberOfCoulombGridPoints = uint3(8, 8, 8);
  forceField.settings.numberOfGridTestPoints = 0;
  const std::size_t oxygen = forceField.findPseudoAtom("O").value();
  const VDWParameters buckingham(VDWParameters::Type::BuckingHam, {1.0e5, 3.5, 1.0e3});
  forceField(2, oxygen) = buckingham;
  forceField(oxygen, 2) = buckingham;
  const Framework framework = Framework::makeFAU(forceField, int3(1, 1, 1));
  const Component methane = Component(forceField, "methane", 190.564, 4599200, 0.01142,
                                      {Atom({0, 0, 0}, 0.0, 1.0, 0, 2, 0, false, false)}, {}, {}, 5, 21);
  System system(forceField, std::nullopt, false, 300.0, 1e5, 1.0, {framework}, {methane}, {}, {0}, 5);
  std::ostringstream log;
  EXPECT_THROW(system.createFrameworkInterpolationGrids(log), std::runtime_error);
}
