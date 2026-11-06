#include <gtest/gtest.h>

import std;

import uint3;
import double3;
import unit_cell;
import crystal;
import units;
import pair_interactions;
import blocking_spheres;
import energy_shared_blocking_mask;
import energy_shared_probe_energy_grid;
import energy_void_fraction;
import energy_opencl_void_fraction;
import energy_opencl_surface_area;
import opencl;

// The single-site device routes for the void fraction and the surface area. The device works over a range
// padded out to whole work-groups, so the grid is followed by points that are not on it; those must weigh
// nothing in the average. And a grid point can lie on an atom's centre, where the pair's separation is zero;
// an atom with no strength of its own has to count for nothing there rather than poison the whole field.

namespace
{

constexpr double edge = 10.0;
constexpr double sigma = 3.0;
constexpr double epsilon = 100.0;
constexpr double cutOff = 4.0;
constexpr double temperature = 298.0;

// What every builder of an energy field holds an overlap at, in the units the fields are in, and the
// temperature at which a point held there weighs one over e. Nothing is measured at a temperature like
// that; it is where the height of the ceiling can be read off an average. The unit system is the engine's
// and is set up when the program starts, so the first is asked for when it is wanted and not before.
double fieldCeiling() { return probeEnergyCeilingInKelvin * Units::KelvinToEnergy; }
constexpr double ceilingTemperature = probeEnergyCeilingInKelvin;

// Three types: the probe, an atom with strength, and one without.
PairInteractions probeAndTwoAtoms()
{
  PairInteractions interactions;
  interactions.numberOfTypes = 3;
  interactions.names = {"He", "O", "H"};
  interactions.charges = {0.0, 0.0, 0.0};
  interactions.parameters = {PairParameters{sigma, epsilon, 0.0}, PairParameters{sigma, epsilon, 0.0},
                             PairParameters{2.5, 0.0, 0.0},       PairParameters{sigma, epsilon, 0.0},
                             PairParameters{sigma, epsilon, 0.0}, PairParameters{2.5, 0.0, 0.0},
                             PairParameters{2.5, 0.0, 0.0},       PairParameters{2.5, 0.0, 0.0},
                             PairParameters{2.5, 0.0, 0.0}};
  interactions.cutOffVDW = cutOff;
  interactions.cutOffCoulomb = cutOff;
  return interactions;
}

// One atom with strength in the middle of the cell. The origin, which is the first grid point and the
// position every padding point is given, is a cutoff and more away from it.
Crystal atomInTheMiddle(const std::string &name)
{
  Crystal framework;
  framework.name = name;
  framework.unitCell = UnitCell(edge, edge, edge);
  framework.mass = 16.0;
  framework.atoms.push_back(CrystalAtom{double3(0.5 * edge, 0.5 * edge, 0.5 * edge), 0.0, 1});
  framework.fractionalPositions.push_back(double3(0.5, 0.5, 0.5));
  return framework;
}

// The same, with a strengthless atom on the origin, exactly on the first grid point.
Crystal withStrengthlessAtomOnTheOrigin(const std::string &name)
{
  Crystal framework = atomInTheMiddle(name);
  framework.atoms.push_back(CrystalAtom{double3(0.0, 0.0, 0.0), 0.0, 2});
  framework.fractionalPositions.push_back(double3(0.0, 0.0, 0.0));
  return framework;
}

// The average of exp(-U/kT) over the grid, summed in double on the host from the one atom with strength.
double expectedVoidFraction(const uint3 &gridSize)
{
  const double beta = 1.0 / (Units::KB * temperature);
  double sum = 0.0;
  for (std::size_t k = 0; k < gridSize.z; ++k)
  {
    for (std::size_t j = 0; j < gridSize.y; ++j)
    {
      for (std::size_t i = 0; i < gridSize.x; ++i)
      {
        double3 ds(static_cast<double>(i) / static_cast<double>(gridSize.x) - 0.5,
                   static_cast<double>(j) / static_cast<double>(gridSize.y) - 0.5,
                   static_cast<double>(k) / static_cast<double>(gridSize.z) - 0.5);
        ds.x -= std::rint(ds.x);
        ds.y -= std::rint(ds.y);
        ds.z -= std::rint(ds.z);
        double rr = edge * edge * double3::dot(ds, ds);
        double energy = 0.0;
        if (rr < cutOff * cutOff)
        {
          double ratio3 = std::pow(sigma * sigma / rr, 3.0);
          energy = std::min(4.0 * epsilon * ratio3 * (ratio3 - 1.0), fieldCeiling());
        }
        sum += std::exp(-beta * energy);
      }
    }
  }
  return sum / (static_cast<double>(gridSize.x) * static_cast<double>(gridSize.y) * static_cast<double>(gridSize.z));
}

bool noDevice()
{
  OpenCL::initialize();
  return !OpenCL::clContext.has_value();
}

// More atoms than the device is handed in one batch, which is 4096: two with strength, both at the middle of
// the cell, and strengthless ones to fill the room between them. With `together` the two with strength come
// first and share the first batch; without it one opens the first batch and the other is alone in the second.
Crystal twoAtomsInTheMiddle(const std::string &name, bool together)
{
  Crystal framework;
  framework.name = name;
  framework.unitCell = UnitCell(edge, edge, edge);
  framework.mass = 16.0;
  auto add = [&](double3 fractional, std::size_t type)
  {
    framework.atoms.push_back(CrystalAtom{edge * fractional, 0.0, type});
    framework.fractionalPositions.push_back(fractional);
  };
  add(double3(0.5, 0.5, 0.5), 1);
  if (together) add(double3(0.5, 0.5, 0.5), 1);
  for (std::size_t filler = 0; filler < 4095; ++filler) add(double3(0.0, 0.0, 0.0), 2);
  if (!together) add(double3(0.5, 0.5, 0.5), 1);
  return framework;
}

// A cell with nothing in it.
Crystal emptyCell(const std::string &name)
{
  Crystal framework;
  framework.name = name;
  framework.unitCell = UnitCell(edge, edge, edge);
  framework.mass = 0.0;
  return framework;
}

// The last line of a report that is not a comment and not about curvature, and whether the report is there.
std::optional<std::string> firstResultOf(const std::string &fileName)
{
  std::ifstream report(fileName);
  if (!report.is_open()) return std::nullopt;
  std::string line;
  while (std::getline(report, line))
  {
    if (!line.empty() && line[0] != '#') return line;
  }
  return std::string{};
}

}  // namespace


// An odd number of grid points is not a multiple of any even work-group size, so the range is padded.
TEST(energy_legacy_void_fraction, padding_points_weigh_nothing_in_the_average)
{
  if (noDevice()) GTEST_SKIP() << "no OpenCL device";
  const uint3 gridSize{9, 9, 9};
  PairInteractions interactions = probeAndTwoAtoms();
  Crystal framework = atomInTheMiddle("legacy-void-fraction-padding");

  EnergyOpenCLVoidFraction vf;
  vf.run(interactions, framework, "He", gridSize, temperature);
  std::filesystem::remove(framework.name + ".energy.vf.gpu.txt");

  const double expected = expectedVoidFraction(gridSize);
  EXPECT_NEAR(vf.voidFraction, expected, 1.0e-4 * expected);
}


TEST(energy_legacy_void_fraction, a_strengthless_atom_on_a_grid_point_counts_for_nothing)
{
  if (noDevice()) GTEST_SKIP() << "no OpenCL device";
  const uint3 gridSize{8, 8, 8};
  PairInteractions interactions = probeAndTwoAtoms();
  Crystal framework = withStrengthlessAtomOnTheOrigin("legacy-void-fraction-origin");

  EnergyOpenCLVoidFraction vf;
  vf.run(interactions, framework, "He", gridSize, temperature);
  std::filesystem::remove(framework.name + ".energy.vf.gpu.txt");

  const double expected = expectedVoidFraction(gridSize);
  ASSERT_TRUE(std::isfinite(vf.voidFraction));
  EXPECT_NEAR(vf.voidFraction, expected, 1.0e-4 * expected);
}


TEST(energy_legacy_void_fraction, a_strengthless_atom_on_a_grid_point_leaves_the_surface_alone)
{
  if (noDevice()) GTEST_SKIP() << "no OpenCL device";
  const uint3 gridSize{32, 32, 32};
  // On the repulsive side of the well, so that the surface is one sphere about the atom with strength and
  // every point beyond the cutoff, the origin among them, lies below it.
  const double isoValue = 10.0;
  PairInteractions interactions = probeAndTwoAtoms();
  Crystal alone = atomInTheMiddle("legacy-surface-area-alone");
  Crystal withStrengthless = withStrengthlessAtomOnTheOrigin("legacy-surface-area-origin");

  EnergyOpenCLSurfaceArea without;
  without.run(interactions, alone, isoValue, "He", gridSize);
  std::filesystem::remove(alone.name + ".energy.sa.gpu.txt");

  EnergyOpenCLSurfaceArea with;
  with.run(interactions, withStrengthless, isoValue, "He", gridSize);
  std::filesystem::remove(withStrengthless.name + ".energy.sa.gpu.txt");

  ASSERT_TRUE(std::isfinite(without.area));
  ASSERT_GT(without.area, 0.0);
  EXPECT_TRUE(std::isfinite(with.area));
  EXPECT_DOUBLE_EQ(with.area, without.area);
}


// A device may report a work-group size that is not a power of two, and the range has to be padded to a multiple of
// it whatever it is, or the device refuses the work. The padding points are not on the grid, so the area is the same.
TEST(energy_legacy_void_fraction, the_range_is_padded_to_a_work_group_size_that_is_no_power_of_two)
{
  if (noDevice()) GTEST_SKIP() << "no OpenCL device";

  PairInteractions interactions = probeAndTwoAtoms();
  Crystal framework = atomInTheMiddle("energy_legacy_workgroup_test");
  const uint3 gridSize(16, 16, 16);

  EnergyOpenCLSurfaceArea asReported;
  asReported.run(interactions, framework, 1000.0, "He", gridSize);
  ASSERT_GT(asReported.area, 0.0);

  for (std::size_t workGroupSize : {3uz, 7uz, 12uz})
  {
    ASSERT_LE(workGroupSize, asReported.energyGridWorkGroupSize) << "the device has to allow this group";

    EnergyOpenCLSurfaceArea route;
    route.energyGridWorkGroupSize = workGroupSize;
    EXPECT_NO_THROW(route.run(interactions, framework, 1000.0, "He", gridSize)) << "group of " << workGroupSize;
    EXPECT_DOUBLE_EQ(route.area, asReported.area) << "group of " << workGroupSize;
  }

  std::filesystem::remove("energy_legacy_workgroup_test.energy.sa.gpu.txt");
}


// The same for the void fraction, which has more to get right. Two kernels run over the one range, each with a
// group size of its own, so the range has to be a multiple of both; and the second adds each group up by folding
// it onto itself, which has to reach every point of a group of any size. A group of twelve on a grid of twelve
// cubed is the case that says nothing when it goes wrong: twelve divides the range, the device takes the work, and
// a fold that halves twelve to six to three to one leaves every third partial sum out of the average.
TEST(energy_legacy_void_fraction, the_average_is_the_same_whatever_the_two_work_group_sizes)
{
  if (noDevice()) GTEST_SKIP() << "no OpenCL device";

  PairInteractions interactions = probeAndTwoAtoms();
  Crystal framework = atomInTheMiddle("energy_legacy_vf_workgroup_test");
  const uint3 gridSize(12, 12, 12);

  EnergyOpenCLVoidFraction asReported;
  asReported.run(interactions, framework, "He", gridSize, temperature);
  ASSERT_GT(asReported.voidFraction, 0.0);

  const std::vector<std::pair<std::size_t, std::size_t>> groupSizes{{12, 12}, {3, 7}, {7, 3}, {5, 12}, {16, 6}, {1, 1}};
  for (auto [ofTheGrid, ofTheSum] : groupSizes)
  {
    ASSERT_LE(ofTheGrid, asReported.energyGridWorkGroupSize) << "the device has to allow this group";
    ASSERT_LE(ofTheSum, asReported.energyVoidFractionWorkGroupSize) << "the device has to allow this group";

    EnergyOpenCLVoidFraction route;
    route.energyGridWorkGroupSize = ofTheGrid;
    route.energyVoidFractionWorkGroupSize = ofTheSum;
    EXPECT_NO_THROW(route.run(interactions, framework, "He", gridSize, temperature))
        << "groups of " << ofTheGrid << " and " << ofTheSum;

    // the partial sums are single precision and are added up in another order
    EXPECT_NEAR(route.voidFraction, asReported.voidFraction, 1.0e-5 * asReported.voidFraction)
        << "groups of " << ofTheGrid << " and " << ofTheSum;
  }

  std::filesystem::remove("energy_legacy_vf_workgroup_test.energy.vf.gpu.txt");
}


// The field the device builds for these two routes is held at the ceiling every other builder holds its
// field at, and not at a figure of its own a thousand times lower. Eight grid points, one of them on the
// atom and the rest beyond its cutoff: at the temperature of the ceiling the one weighs one over e and the
// rest weigh one.
TEST(energy_legacy_void_fraction, an_overlap_is_held_at_the_ceiling_every_field_is_held_at)
{
  if (noDevice()) GTEST_SKIP() << "no OpenCL device";
  const uint3 gridSize{2, 2, 2};
  PairInteractions interactions = probeAndTwoAtoms();
  Crystal framework = atomInTheMiddle("legacy-void-fraction-ceiling");

  EnergyOpenCLVoidFraction vf;
  vf.run(interactions, framework, "He", gridSize, ceilingTemperature);
  std::filesystem::remove(framework.name + ".energy.vf.gpu.txt");

  EXPECT_NEAR(vf.voidFraction, (7.0 + std::exp(-1.0)) / 8.0, 1.0e-5);
}


// The atoms are handed to the device 4096 at a time and the field is their sum, so it cannot matter which
// of them share a batch. Two atoms on the one grid point are an overlap whichever way they are dealt out, and
// an overlap is at the ceiling: not at twice it for the two having come in different batches, which is what
// holding each batch down by itself made of them.
TEST(energy_legacy_void_fraction, the_field_is_the_same_however_the_atoms_fall_into_batches)
{
  if (noDevice()) GTEST_SKIP() << "no OpenCL device";
  const uint3 gridSize{2, 2, 2};
  PairInteractions interactions = probeAndTwoAtoms();

  for (bool together : {true, false})
  {
    Crystal framework = twoAtomsInTheMiddle("legacy-void-fraction-batches", together);
    ASSERT_GT(framework.atoms.size(), 4096uz);

    EnergyOpenCLVoidFraction vf;
    vf.run(interactions, framework, "He", gridSize, ceilingTemperature);
    std::filesystem::remove(framework.name + ".energy.vf.gpu.txt");

    EXPECT_NEAR(vf.voidFraction, (7.0 + std::exp(-1.0)) / 8.0, 1.0e-5)
        << (together ? "the two in one batch" : "the two in two batches");
  }
}


// The same of the surface area's field, which is built by a kernel of its own. A level a little under the
// ceiling is crossed, close in to the atom, and a level over it is crossed nowhere, however many atoms are
// on the one spot and however they were dealt out.
TEST(energy_legacy_void_fraction, the_surface_area_field_reaches_the_ceiling_and_stops_there)
{
  if (noDevice()) GTEST_SKIP() << "no OpenCL device";
  const uint3 gridSize{64, 64, 64};
  PairInteractions interactions = probeAndTwoAtoms();

  // Where one atom's energy comes to half the ceiling: a sphere under an Ångström across.
  const double radius = sigma * std::pow(4.0 * epsilon / (0.5 * fieldCeiling()), 1.0 / 12.0);
  ASSERT_GT(radius, 4.0 * edge / 64.0);
  ASSERT_LT(radius, 1.0);
  const double sphere = 4.0 * std::numbers::pi * radius * radius;

  Crystal alone = atomInTheMiddle("legacy-surface-area-ceiling");
  EnergyOpenCLSurfaceArea under;
  under.run(interactions, alone, 0.5 * fieldCeiling(), "He", gridSize);
  EXPECT_NEAR(under.area, sphere, 0.2 * sphere);

  EnergyOpenCLSurfaceArea over;
  over.run(interactions, alone, 1.5 * fieldCeiling(), "He", gridSize);
  EXPECT_EQ(over.area, 0.0);
  std::filesystem::remove(alone.name + ".energy.sa.gpu.txt");

  for (bool together : {true, false})
  {
    Crystal framework = twoAtomsInTheMiddle("legacy-surface-area-batches", together);
    EnergyOpenCLSurfaceArea route;
    route.run(interactions, framework, 1.5 * fieldCeiling(), "He", gridSize);
    std::filesystem::remove(framework.name + ".energy.sa.gpu.txt");
    EXPECT_EQ(route.area, 0.0) << (together ? "the two in one batch" : "the two in two batches");
  }
}


// A cell with no atoms in it is all void: the average of exp(-U/kT) is one, on the device as by insertion on
// the processor, and the report is written as for any other cell. The device route used to do nothing at
// all, which left the average at the nought it starts at and no report behind it.
TEST(energy_legacy_void_fraction, a_cell_with_no_atoms_is_all_void)
{
  PairInteractions interactions = probeAndTwoAtoms();
  Crystal framework = emptyCell("legacy-void-fraction-empty");

  EnergyVoidFraction widom;
  widom.run(interactions, framework, "He", 10, 100, temperature);
  EXPECT_DOUBLE_EQ(widom.boltzmannAverage, 1.0);
  std::filesystem::remove(framework.name + ".energy.vf.cpu.txt");

  if (noDevice()) GTEST_SKIP() << "no OpenCL device";

  const std::string report = framework.name + ".energy.vf.gpu.txt";
  std::filesystem::remove(report);

  EnergyOpenCLVoidFraction vf;
  vf.run(interactions, framework, "He", uint3{9, 9, 9}, temperature);
  EXPECT_DOUBLE_EQ(vf.voidFraction, 1.0);
  const std::optional<std::string> written = firstResultOf(report);
  ASSERT_TRUE(written.has_value()) << "no report was written";
  EXPECT_DOUBLE_EQ(std::stod(written.value()), 1.0);

  // A pocket handed over with it is still left out, at the usual thousand kelvin for every Ångström a point
  // is inside it: all of the sphere but a skin a few tenths of an Ångström thick, where a probe this warm
  // still has some weight.
  const double radius = 3.0;
  const std::size_t points = 40;
  const std::vector<BlockingSphere> spheres{BlockingSphere{double3(0.5, 0.5, 0.5), radius}};
  double sum = 0.0;
  for (std::size_t k = 0; k < points; ++k)
  {
    for (std::size_t j = 0; j < points; ++j)
    {
      for (std::size_t i = 0; i < points; ++i)
      {
        const double3 fromCentre(static_cast<double>(i) / static_cast<double>(points) - 0.5,
                                 static_cast<double>(j) / static_cast<double>(points) - 0.5,
                                 static_cast<double>(k) / static_cast<double>(points) - 0.5);
        const double depth = radius - edge * fromCentre.length();
        sum += depth > 0.0 ? std::exp(-depth * blockedEnergyPerAngstromInKelvin / temperature) : 1.0;
      }
    }
  }
  const double expected = sum / static_cast<double>(points * points * points);
  const double share = 4.0 / 3.0 * std::numbers::pi * radius * radius * radius / (edge * edge * edge);
  ASSERT_LT(expected, 1.0 - 0.6 * share);

  EnergyOpenCLVoidFraction blocked;
  blocked.run(interactions, framework, "He", uint3{40, 40, 40}, temperature, spheres);
  EXPECT_NEAR(blocked.voidFraction, expected, 1.0e-4);
  std::filesystem::remove(report);
}


// And it has no surface: the area is nought, and the report says so rather than not being written, with
// nought for the area to a gram of a cell that weighs nothing.
TEST(energy_legacy_void_fraction, a_cell_with_no_atoms_has_no_surface)
{
  if (noDevice()) GTEST_SKIP() << "no OpenCL device";
  PairInteractions interactions = probeAndTwoAtoms();
  Crystal framework = emptyCell("legacy-surface-area-empty");
  const std::string report = framework.name + ".energy.sa.gpu.txt";
  std::filesystem::remove(report);

  EnergyOpenCLSurfaceArea sa;
  sa.run(interactions, framework, 10.0, "He", uint3{16, 16, 16});
  EXPECT_EQ(sa.area, 0.0);

  std::ifstream written(report);
  ASSERT_TRUE(written.is_open()) << "no report was written";
  std::vector<std::string> results;
  for (std::string line; std::getline(written, line);)
  {
    if (!line.empty() && line[0] != '#') results.push_back(line);
  }
  ASSERT_GE(results.size(), 3uz);
  EXPECT_EQ(results[0], "0 [A^2]");
  EXPECT_EQ(results[1], "0 [m^2/g]");
  EXPECT_EQ(results[2], "0 [m^2/cm^3]");
  written.close();
  std::filesystem::remove(report);
}
