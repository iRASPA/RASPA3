#include <gtest/gtest.h>

import std;

import double3;
import unit_cell;
import opencl;
import sampled_structure;
import mc_pore_size_distribution;
import mc_opencl_pore_size_distribution;

// The sampled pore-size distribution, on the processor and on the device. A single sphere in a large box has a void
// fraction known in closed form, and every point of the void is at least covered by the sphere about itself, so what
// the curve must add up to is known however few centres are drawn.

namespace
{

constexpr double edge = 20.0;
constexpr double atomRadius = 2.0;

SampledStructure oneSphere(const std::string &name)
{
  SampledStructure structure;
  structure.name = name;
  structure.unitCell = UnitCell(edge, edge, edge);
  structure.positions.push_back(double3(10.0, 10.0, 10.0));
  structure.radii.push_back(atomRadius);
  return structure;
}

double voidFraction()
{
  return 1.0 - 4.0 / 3.0 * std::numbers::pi * atomRadius * atomRadius * atomRadius / (edge * edge * edge);
}

// Column 1 is the diameter and column 2 the distribution, from the report the run writes.
std::vector<std::pair<double, double>> readCurve(const std::string &file)
{
  std::vector<std::pair<double, double>> curve;
  std::ifstream stream(file);
  std::string line;
  while (std::getline(stream, line))
  {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream row(line);
    double diameter = 0.0;
    double distribution = 0.0;
    row >> diameter >> distribution;
    curve.push_back({diameter, distribution});
  }
  std::filesystem::remove(file);
  return curve;
}

template <typename Route>
void checkRoute(Route &route, const SampledStructure &structure, const std::string &file, std::size_t bins,
                std::size_t iterations)
{
  const double standardError = std::sqrt(voidFraction() * (1.0 - voidFraction()) / static_cast<double>(iterations));
  const double range = 2.5;  // short of most of the room in the box, so nearly every point is past the last bin

  route.run(structure, iterations, std::size_t{1}, range);

  // Every point of the void is at or above every diameter short of its own, past the end of the range or not, so at
  // the first bin the cumulative curve is the void fraction, and it does not depend on what the one centre covered.
  const double atFirst = route.histogram_cummulative[0] / static_cast<double>(iterations);
  EXPECT_NEAR(atFirst, voidFraction(), std::max(5.0 * standardError, 1.0 / static_cast<double>(iterations)));

  // What is past the range counts at every bin and is drawn in none.
  double drawn = 0.0;
  for (double count : route.histogram) drawn += count;
  EXPECT_DOUBLE_EQ(route.histogram_cummulative[bins - 1], route.histogram[bins - 1] + (route.histogram_cummulative[0] - drawn));
  EXPECT_LT(drawn, route.histogram_cummulative[0]);

  // A density over the diameter: the height is the share per Å, and the width of a bin is twice that of the radii.
  const double width = 2.0 * range / static_cast<double>(bins);
  std::vector<std::pair<double, double>> curve = readCurve(file);
  ASSERT_EQ(curve.size(), bins);
  double integral = 0.0;
  for (std::size_t bin = 0; bin < bins; ++bin)
  {
    EXPECT_NEAR(curve[bin].first, width * (static_cast<double>(bin) + 0.5), 1.0e-9);
    EXPECT_NEAR(curve[bin].second, route.histogram[bin] / static_cast<double>(iterations) / width, 1.0e-9);
    integral += curve[bin].second * width;
  }
  EXPECT_NEAR(integral, drawn / static_cast<double>(iterations), 1.0e-9);

  // A second run is a run of its own and does not add to the first.
  route.run(structure, iterations, std::size_t{1}, range);
  EXPECT_NEAR(route.histogram_cummulative[0] / static_cast<double>(iterations), voidFraction(),
              std::max(5.0 * standardError, 1.0 / static_cast<double>(iterations)));
  std::filesystem::remove(file);
}

}  // namespace


TEST(mc_pore_size_distribution, the_curve_adds_up_on_the_processor)
{
  MC_PoreSizeDistribution route(25);
  checkRoute(route, oneSphere("mc_psd_test_cpu"), "mc_psd_test_cpu.mc.psd.cpu.txt", 25, 4000);
}

TEST(mc_pore_size_distribution, the_curve_adds_up_on_the_device)
{
  OpenCL::initialize();
  if (!OpenCL::clContext.has_value()) GTEST_SKIP() << "no OpenCL device";
  MC_OpenCL_PoreSizeDistribution route(25);
  checkRoute(route, oneSphere("mc_psd_test_gpu"), "mc_psd_test_gpu.mc.psd.gpu.txt", 25, 500);
}

// With no centres drawn at all a point still has the sphere about itself.
TEST(mc_pore_size_distribution, a_point_is_covered_by_the_sphere_about_itself)
{
  MC_PoreSizeDistribution route(25);
  constexpr std::size_t iterations = 4000;
  route.run(oneSphere("mc_psd_test_none"), iterations, std::size_t{0}, 12.0);
  std::filesystem::remove("mc_psd_test_none.mc.psd.cpu.txt");
  const double standardError = std::sqrt(voidFraction() * (1.0 - voidFraction()) / static_cast<double>(iterations));
  EXPECT_NEAR(route.histogram_cummulative[0] / static_cast<double>(iterations), voidFraction(),
              std::max(5.0 * standardError, 1.0 / static_cast<double>(iterations)));
}

// No atoms, so a point has as much room as a double can say, which is past any range and must not be turned into an
// index. No bins and no range have nothing to draw and must not try.
TEST(mc_pore_size_distribution, an_empty_box_and_an_empty_range_are_handled)
{
  SampledStructure empty;
  empty.name = "mc_psd_test_empty";
  empty.unitCell = UnitCell(edge, edge, edge);

  MC_PoreSizeDistribution route(5);
  route.run(empty, 100, std::size_t{1}, 10.0);
  for (double value : route.histogram_cummulative) EXPECT_DOUBLE_EQ(value, 100.0);
  for (double value : route.histogram) EXPECT_DOUBLE_EQ(value, 0.0);

  MC_PoreSizeDistribution noBins(0);
  noBins.run(oneSphere("mc_psd_test_empty"), 10, std::size_t{1}, 10.0);
  EXPECT_TRUE(noBins.histogram.empty());

  MC_PoreSizeDistribution noRange(5);
  noRange.run(oneSphere("mc_psd_test_empty"), 10, std::size_t{1}, 0.0);
  for (double value : noRange.histogram) EXPECT_DOUBLE_EQ(value, 0.0);

  std::filesystem::remove("mc_psd_test_empty.mc.psd.cpu.txt");
}

// The same on the device, and that the device asks for no more centres than it was told to.
TEST(mc_pore_size_distribution, no_centres_and_an_empty_box_on_the_device)
{
  OpenCL::initialize();
  if (!OpenCL::clContext.has_value()) GTEST_SKIP() << "no OpenCL device";

  MC_OpenCL_PoreSizeDistribution route(5);
  route.run(oneSphere("mc_psd_test_device"), 200, std::size_t{0}, 12.0);
  EXPECT_NEAR(route.histogram_cummulative[0] / 200.0, voidFraction(), 0.05);

  SampledStructure empty;
  empty.name = "mc_psd_test_device";
  empty.unitCell = UnitCell(edge, edge, edge);
  route.run(empty, 20, std::size_t{1}, 10.0);
  for (double value : route.histogram_cummulative) EXPECT_DOUBLE_EQ(value, 20.0);
  for (double value : route.histogram) EXPECT_DOUBLE_EQ(value, 0.0);

  std::filesystem::remove("mc_psd_test_device.mc.psd.gpu.txt");
}

namespace
{

// One small atom to a cell whose a and b are at sixty degrees: columns of atoms on a triangular net, ten apart, in
// layers ten apart. The point farthest from all of them is over the middle of a triangle and half way between two
// layers, 10/sqrt(3) from three columns and five from each layer, so no sphere that misses the atoms has a radius
// above sqrt(100/3 + 25) less the atom's own, and no point of the void is in a larger one.
constexpr double skewedAtomRadius = 1.0;
const double largestFreeRadius = std::sqrt(100.0 / 3.0 + 25.0) - skewedAtomRadius;

SampledStructure oneAtomInASkewedCell(const std::string &name)
{
  SampledStructure structure;
  structure.name = name;
  structure.unitCell =
      UnitCell(10.0, 10.0, 10.0, 0.5 * std::numbers::pi, 0.5 * std::numbers::pi, std::numbers::pi / 3.0);
  structure.positions.push_back(double3(0.0, 0.0, 0.0));
  structure.radii.push_back(skewedAtomRadius);
  return structure;
}

// In a skewed cell the copy of an atom nearest a point in fractions of the cell is not always the nearest in
// Ångström, and a sphere measured against that copy comes out wider than the room there is. So nothing is drawn
// above the largest radius the cell has room for, and nothing past the range, which is wider still.
template <typename Route>
void checkNothingIsWiderThanTheRoom(Route &route, const SampledStructure &structure, const std::string &file,
                                    std::size_t iterations)
{
  constexpr std::size_t bins = 45;
  constexpr double range = 9.0;
  route.run(structure, iterations, std::size_t{64}, range);
  std::filesystem::remove(file);

  const std::size_t firstBinPastTheRoom = static_cast<std::size_t>(std::ceil(largestFreeRadius / (range / bins))) + 1;
  ASSERT_LT(firstBinPastTheRoom, bins);
  double drawn = 0.0;
  double pastTheRoom = 0.0;
  for (std::size_t bin = 0; bin < bins; ++bin)
  {
    drawn += route.histogram[bin];
    if (bin >= firstBinPastTheRoom) pastTheRoom += route.histogram[bin];
  }
  EXPECT_EQ(pastTheRoom, 0.0) << "points were put in spheres wider than any the cell has room for";
  EXPECT_DOUBLE_EQ(drawn, route.histogram_cummulative[0]) << "and none can be past a range that is wider still";
  EXPECT_GT(drawn, 0.9 * static_cast<double>(iterations));
}

}  // namespace

TEST(mc_pore_size_distribution, no_sphere_is_wider_than_the_room_in_a_skewed_cell_on_the_processor)
{
  MC_PoreSizeDistribution route(45);
  checkNothingIsWiderThanTheRoom(route, oneAtomInASkewedCell("mc_psd_skewed_cpu"), "mc_psd_skewed_cpu.mc.psd.cpu.txt",
                                 2000);
}

TEST(mc_pore_size_distribution, no_sphere_is_wider_than_the_room_in_a_skewed_cell_on_the_device)
{
  OpenCL::initialize();
  if (!OpenCL::clContext.has_value()) GTEST_SKIP() << "no OpenCL device";
  MC_OpenCL_PoreSizeDistribution route(45);
  checkNothingIsWiderThanTheRoom(route, oneAtomInASkewedCell("mc_psd_skewed_gpu"), "mc_psd_skewed_gpu.mc.psd.gpu.txt",
                                 500);
}
