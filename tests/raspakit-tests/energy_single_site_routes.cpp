#include <gtest/gtest.h>

import std;

import uint3;
import double3;
import units;
import unit_cell;
import crystal;
import pair_interactions;
import blocking_spheres;
import opencl;
import energy_backend;
import energy_isosurface;
import energy_shared_energy_backend;
import energy_shared_isosurface;
import energy_shared_probe_energy_grid;
import energy_shared_blocking_mask;
import energy_void_fraction;
import energy_opencl_surface_area;
import energy_opencl_void_fraction;

// The routes that take one probe atom and build their own energies: Widom insertion on the processor, and the
// surface area and the void fraction on the device. None of them goes through the backend, which is where a run's
// blocked pockets are written onto a field, so each has to be handed the pockets and has to use them; and the
// Widom route has to be told the temperature it is averaging at.
//
// They are tested from here and not from the structural library's own tests because they read the unit system,
// which the engine defines.

namespace
{

constexpr double edge = 20.0;
constexpr double sigma = 3.0;
constexpr double epsilon = 100.0;

PairInteractions probeAndAtom()
{
  PairInteractions interactions;
  interactions.numberOfTypes = 2;
  interactions.names = {"X", "O"};
  interactions.charges = {0.0, 0.0};
  interactions.parameters = {PairParameters{sigma, epsilon, 0.0}, PairParameters{sigma, epsilon, 0.0},
                             PairParameters{sigma, epsilon, 0.0}, PairParameters{sigma, epsilon, 0.0}};
  interactions.cutOffVDW = 9.0;
  interactions.cutOffCoulomb = 9.0;
  return interactions;
}

// One atom at the middle of a cubic cell wide enough that it meets none of its own copies.
Crystal oneAtom(const std::string &name)
{
  Crystal framework;
  framework.name = name;
  framework.unitCell = UnitCell(edge, edge, edge);
  framework.mass = 16.0;
  const double3 fractional(0.5, 0.5, 0.5);
  framework.atoms.push_back(CrystalAtom{fractional * edge, 0.0, 1});
  framework.fractionalPositions.push_back(fractional);
  return framework;
}

// A sphere five Ångström from the atom, inside the shell where the atom attracts, and steep enough that a point
// inside it is worth nothing.
const BlockingSphere pocket{double3(0.75, 0.5, 0.5), 1.5};
constexpr double steep = 1.0e6;
constexpr double ceiling = 1.0e9;

double sphereShareOfCell()
{
  return 4.0 / 3.0 * std::numbers::pi * pocket.radius * pocket.radius * pocket.radius / (edge * edge * edge);
}

// The field the processor backend makes for the probe, with the pocket written onto it or not.
ProbeEnergyGrid processorField(const Crystal &framework, uint3 gridSize, bool blocked)
{
  ProbeEnergyGrid grid = cpuEnergyBackend().probeEnergyGrid(probeAndAtom(), framework, "X", gridSize);
  if (blocked)
  {
    const std::vector<BlockingSphere> spheres{pocket};
    applyBlockingSpheresToEnergy(grid.energy, grid.gridSize, grid.unitCell, spheres, steep, ceiling);
  }
  return grid;
}

}  // namespace

// A table for the routes that read sizes: every size scaled, the probe's with the atoms', and nothing else.
TEST(energy_single_site_routes, a_table_of_scaled_sizes_scales_the_sizes_and_nothing_else)
{
  PairInteractions interactions = probeAndAtom();
  interactions.parameters[1] = PairParameters{3.4, 80.0, 0.25};
  interactions.parameters[2] = PairParameters{3.4, 80.0, 0.25};
  const double factor = std::pow(2.0, 1.0 / 6.0);
  const PairInteractions scaled = interactions.withSizesScaledBy(factor);

  ASSERT_EQ(scaled.parameters.size(), interactions.parameters.size());
  for (std::size_t pair = 0; pair < interactions.parameters.size(); ++pair)
  {
    EXPECT_DOUBLE_EQ(scaled.parameters[pair].sizeParameter, factor * interactions.parameters[pair].sizeParameter);
    EXPECT_DOUBLE_EQ(scaled.parameters[pair].strengthParameter, interactions.parameters[pair].strengthParameter);
    EXPECT_DOUBLE_EQ(scaled.parameters[pair].shift, interactions.parameters[pair].shift);
  }
  EXPECT_EQ(scaled.names, interactions.names);
  EXPECT_EQ(scaled.charges, interactions.charges);
  EXPECT_EQ(scaled.findType("X"), interactions.findType("X"));
  EXPECT_DOUBLE_EQ(scaled.cutOffVDW, interactions.cutOffVDW);
}

// The average of exp(-U/kT) is of a temperature. The atom attracts, so the colder the probe the more the shell
// round the atom is worth, and the average rises: it was the same at every temperature asked for, being taken at
// room temperature whatever was asked.
TEST(energy_single_site_routes, widom_insertion_is_at_the_temperature_it_is_given)
{
  Crystal framework = oneAtom("widom_temperature");
  EnergyVoidFraction warm;
  warm.run(probeAndAtom(), framework, "X", 80, 500, 300.0);
  EnergyVoidFraction cold;
  cold.run(probeAndAtom(), framework, "X", 80, 500, 60.0);
  std::filesystem::remove(framework.name + ".energy.vf.cpu.txt");

  EXPECT_EQ(warm.numberOfInsertions, 40000uz);

  // against the same average summed over the processor's field
  const ProbeEnergyGrid grid = processorField(framework, uint3{40, 40, 40}, false);
  auto summed = [&](double temperature)
  {
    double sum = 0.0;
    for (float energy : grid.energy) sum += std::exp(-static_cast<double>(energy) / (Units::KB * temperature));
    return sum / static_cast<double>(grid.energy.size());
  };
  EXPECT_NEAR(warm.boltzmannAverage, summed(300.0), 0.02);
  EXPECT_NEAR(cold.boltzmannAverage, summed(60.0), 0.03);
  EXPECT_GT(summed(60.0), summed(300.0) + 0.08)
      << "the two temperatures have to differ for this to have asked anything";
  EXPECT_GT(cold.boltzmannAverage, warm.boltzmannAverage + 0.05);

  EXPECT_THROW(warm.run(probeAndAtom(), framework, "X", 1, 1, 0.0), std::runtime_error);
}

// A point in a pocket the probe cannot reach is worth nothing. A sphere at the corner of the cell, out of the atom's
// reach where every point is otherwise worth one, takes its share of the cell out of the average, and that share
// is a sixth of it here, far more than the scatter of two runs.
TEST(energy_single_site_routes, widom_insertion_leaves_out_a_blocked_pocket)
{
  Crystal framework = oneAtom("widom_blocked");
  const BlockingSphere corner{double3(0.0, 0.0, 0.0), 7.0};
  const std::vector<BlockingSphere> spheres{corner};
  const double share = 4.0 / 3.0 * std::numbers::pi * corner.radius * corner.radius * corner.radius / (edge * edge * edge);

  EnergyVoidFraction open;
  open.run(probeAndAtom(), framework, "X", 80, 500, 300.0);
  EnergyVoidFraction blocked;
  blocked.run(probeAndAtom(), framework, "X", 80, 500, 300.0, spheres, steep, ceiling);
  std::filesystem::remove(framework.name + ".energy.vf.cpu.txt");

  EXPECT_GT(share, 0.15);
  EXPECT_NEAR(open.boltzmannAverage - blocked.boltzmannAverage, share, 0.02);
}

// The device's own surface area, with the pocket written onto its field, is the surface the processor's field
// has with the pocket written onto it: the sheets round the atom and the wall of the pocket as well.
TEST(energy_single_site_routes, the_device_surface_area_walls_off_a_blocked_pocket)
{
  OpenCL::initialize();
  if (!OpenCL::clContext.has_value()) GTEST_SKIP() << "no OpenCL device";

  Crystal framework = oneAtom("device_surface_blocked");
  const uint3 gridSize{64, 64, 64};
  const double isoValue = -2.0;  // inside the shell where the atom attracts, which the pocket is in
  const std::vector<BlockingSphere> spheres{pocket};

  EnergyOpenCLSurfaceArea open;
  open.run(probeAndAtom(), framework, isoValue, "X", gridSize);
  EnergyOpenCLSurfaceArea blocked;
  blocked.run(probeAndAtom(), framework, isoValue, "X", gridSize, spheres, steep, ceiling);
  std::filesystem::remove(framework.name + ".energy.sa.gpu.txt");

  auto areaOnProcessor = [&](bool withThePocket)
  {
    const ProbeEnergyGrid grid = processorField(framework, gridSize, withThePocket);
    return EnergyIsosurface::areaOfIsosurface(framework, grid.energy, gridSize, isoValue).area;
  };
  const double openOnProcessor = areaOnProcessor(false);
  const double blockedOnProcessor = areaOnProcessor(true);

  EXPECT_NEAR(open.area, openOnProcessor, 0.01 * openOnProcessor);
  EXPECT_NEAR(blocked.area, blockedOnProcessor, 0.01 * blockedOnProcessor);

  // the wall of the pocket is surface the open field does not have
  const double wall = 4.0 * std::numbers::pi * pocket.radius * pocket.radius;
  EXPECT_NEAR(blocked.area - open.area, wall, 0.25 * wall);
}

TEST(energy_single_site_routes, the_device_void_fraction_leaves_out_a_blocked_pocket)
{
  OpenCL::initialize();
  if (!OpenCL::clContext.has_value()) GTEST_SKIP() << "no OpenCL device";

  Crystal framework = oneAtom("device_void_blocked");
  const uint3 gridSize{40, 40, 40};
  const double temperature = 300.0;
  const std::vector<BlockingSphere> spheres{pocket};

  EnergyOpenCLVoidFraction open;
  open.run(probeAndAtom(), framework, "X", gridSize, temperature);
  EnergyOpenCLVoidFraction blocked;
  blocked.run(probeAndAtom(), framework, "X", gridSize, temperature, spheres, steep, ceiling);
  std::filesystem::remove(framework.name + ".energy.vf.gpu.txt");

  auto summed = [&](bool withThePocket)
  {
    const ProbeEnergyGrid grid = processorField(framework, gridSize, withThePocket);
    double sum = 0.0;
    for (float energy : grid.energy) sum += std::exp(-static_cast<double>(energy) / (Units::KB * temperature));
    return sum / static_cast<double>(grid.energy.size());
  };
  EXPECT_NEAR(open.voidFraction, summed(false), 1.0e-3 * summed(false));
  EXPECT_NEAR(blocked.voidFraction, summed(true), 1.0e-3 * summed(true));
  EXPECT_GT(open.voidFraction - blocked.voidFraction, 0.5 * sphereShareOfCell());
}


// A caller that hands over the pockets and says nothing of what a point inside one is worth gets the usual
// figures, the ramp the backend writes with and the ceiling every field is held at. The arguments are nought
// when they are left out, and nought taken at its word is a pocket whose every point is worth exactly
// nothing: open void, weighing one, where the pocket was to weigh nothing.
TEST(energy_single_site_routes, pockets_handed_over_without_an_energy_are_blocked_all_the_same)
{
  const double usualRamp = blockedEnergyPerAngstromInKelvin * Units::KelvinToEnergy;
  const double usualCeiling = probeEnergyCeilingInKelvin * Units::KelvinToEnergy;

  // By insertion, with the sphere at the corner where every point is otherwise worth one.
  {
    Crystal framework = oneAtom("widom_blocked_by_default");
    const BlockingSphere corner{double3(0.0, 0.0, 0.0), 7.0};
    const std::vector<BlockingSphere> spheres{corner};
    const double share =
        4.0 / 3.0 * std::numbers::pi * corner.radius * corner.radius * corner.radius / (edge * edge * edge);

    EnergyVoidFraction open;
    open.run(probeAndAtom(), framework, "X", 80, 500, 300.0);
    EnergyVoidFraction blocked;
    blocked.run(probeAndAtom(), framework, "X", 80, 500, 300.0, spheres);
    std::filesystem::remove(framework.name + ".energy.vf.cpu.txt");

    // The usual ramp is a thousand kelvin to the Ångström, so at this temperature a point keeps a weight of
    // exp(-depth / skin) with a skin of three tenths of an Ångström, and what the sphere takes out of the
    // average is its volume less the integral of that over it.
    const double skin = 300.0 / blockedEnergyPerAngstromInKelvin;
    const double radius = corner.radius;
    const double kept = 4.0 * std::numbers::pi *
                        (skin * radius * radius - 2.0 * skin * skin * radius +
                         2.0 * skin * skin * skin * (1.0 - std::exp(-radius / skin))) /
                        (edge * edge * edge);
    ASSERT_GT(share - kept, 0.15);
    EXPECT_NEAR(open.boltzmannAverage - blocked.boltzmannAverage, share - kept, 0.02);
  }

  OpenCL::initialize();
  if (!OpenCL::clContext.has_value()) GTEST_SKIP() << "no OpenCL device";

  const std::vector<BlockingSphere> spheres{pocket};
  {
    Crystal framework = oneAtom("device_void_blocked_by_default");
    const uint3 gridSize{40, 40, 40};

    EnergyOpenCLVoidFraction open;
    open.run(probeAndAtom(), framework, "X", gridSize, 300.0);
    EnergyOpenCLVoidFraction byDefault;
    byDefault.run(probeAndAtom(), framework, "X", gridSize, 300.0, spheres);
    EnergyOpenCLVoidFraction asUsual;
    asUsual.run(probeAndAtom(), framework, "X", gridSize, 300.0, spheres, usualRamp, usualCeiling);
    std::filesystem::remove(framework.name + ".energy.vf.gpu.txt");

    // most of the pocket's share of the cell, the rest being the skin of it that still has some weight
    EXPECT_DOUBLE_EQ(byDefault.voidFraction, asUsual.voidFraction);
    EXPECT_LT(byDefault.voidFraction, open.voidFraction - 0.5 * sphereShareOfCell());
  }
  {
    Crystal framework = oneAtom("device_surface_blocked_by_default");
    const uint3 gridSize{64, 64, 64};
    const double isoValue = -2.0;

    EnergyOpenCLSurfaceArea open;
    open.run(probeAndAtom(), framework, isoValue, "X", gridSize);
    EnergyOpenCLSurfaceArea byDefault;
    byDefault.run(probeAndAtom(), framework, isoValue, "X", gridSize, spheres);
    EnergyOpenCLSurfaceArea asUsual;
    asUsual.run(probeAndAtom(), framework, isoValue, "X", gridSize, spheres, usualRamp, usualCeiling);
    std::filesystem::remove(framework.name + ".energy.sa.gpu.txt");

    const double wall = 4.0 * std::numbers::pi * pocket.radius * pocket.radius;
    EXPECT_DOUBLE_EQ(byDefault.area, asUsual.area);
    EXPECT_NEAR(byDefault.area - open.area, wall, 0.25 * wall);
  }
}
