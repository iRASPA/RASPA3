#include <gtest/gtest.h>

import std;

import double3;
import double3x3;
import unit_cell;
import randomnumbers;
import opencl;
import sampled_structure;
import sampling_backend;
import mc_opencl_backend;

// The device against the structure itself. `SampledStructure` walks the periodic bins out from a point for as far as
// anything could still matter, so it finds the right image of every atom in any cell and is what the kernels have to
// agree with, to the precision of the single-precision floats they work in.

namespace
{

constexpr double agreement = 2.0e-3;

SampledStructure structureOf(const UnitCell &cell, const std::vector<double3> &fractional, double radius)
{
  SampledStructure structure;
  structure.name = "sampling_backend_opencl_test";
  structure.unitCell = cell;
  for (const double3 &position : fractional)
  {
    structure.positions.push_back(cell.cell * position);
    structure.radii.push_back(radius);
  }
  return structure;
}

std::vector<SampledWay> devicesWays(const SampledStructure &structure, const std::vector<double3> &origins,
                                    const std::vector<double3> &displacements)
{
  std::vector<SampledWay> ways(origins.size());
  samplingBackendOpenCL().straightWays(structure, origins, displacements, ways);
  return ways;
}

}  // namespace


// A long segment in a line the atom's lattice does not lie along. The atom's images fall on the line near one end and
// elsewhere near the other, and the image nearest the middle of the segment is one of those off the line: judged by
// that one image alone the segment clears the atom by a long way, when it runs straight through another.
TEST(sampling_backend_opencl, a_long_segment_is_blocked_by_an_image_far_from_its_middle)
{
  OpenCL::initialize();
  if (!samplingOpenCLAvailable()) GTEST_SKIP() << "no OpenCL device";

  const UnitCell box(10.0, 10.0, 30.0);
  const double3 onTheLine(-7.0, -3.5, 5.0);
  SampledStructure structure = structureOf(box, {box.inverseCell * onTheLine}, 1.0);

  const std::vector<double3> origins{double3(-8.0, -4.0, 5.0)};
  const std::vector<double3> displacements{double3(16.0, 8.0, 0.0)};

  const SegmentBottleneck expected = structure.segmentBottleneck(origins[0], displacements[0]);
  ASSERT_NEAR(expected.radius, -1.0, 1.0e-9) << "the segment runs through the centre of the atom";

  std::vector<SampledWay> ways = devicesWays(structure, origins, displacements);
  EXPECT_NEAR(ways[0].radius, -1.0, agreement);
  EXPECT_NEAR((ways[0].position - expected.position).length(), 0.0, agreement);
}

// Any segment the roadmap may ask about, in cells square and skewed.
TEST(sampling_backend_opencl, straight_ways_and_clearances_agree_with_the_structure)
{
  OpenCL::initialize();
  if (!samplingOpenCLAvailable()) GTEST_SKIP() << "no OpenCL device";

  const std::vector<std::pair<const char *, UnitCell>> cells = {
      {"cube", UnitCell(10.0, 10.0, 10.0)},
      {"long", UnitCell(8.0, 11.0, 24.0)},
      {"skew", UnitCell(double3x3(double3(10.0, 0.0, 0.0), double3(4.0, 9.0, 0.0), double3(3.0, 3.5, 11.0)))},
      {"very skew", UnitCell(double3x3(double3(9.0, 0.0, 0.0), double3(8.0, 5.0, 0.0), double3(7.5, 4.0, 4.0)))}};

  RandomNumber random{std::optional<std::size_t>(17)};
  for (const auto &[name, cell] : cells)
  {
    std::vector<double3> fractional;
    for (std::size_t i = 0; i < 24; ++i) fractional.push_back(double3(random.uniform(), random.uniform(), random.uniform()));
    SampledStructure structure = structureOf(cell, fractional, 1.2);

    const double longest = std::max({cell.lengthA, cell.lengthB, cell.lengthC});
    std::vector<double3> origins;
    std::vector<double3> displacements;
    for (std::size_t i = 0; i < 4000; ++i)
    {
      origins.push_back(cell.cell * double3(random.uniform(), random.uniform(), random.uniform()));
      displacements.push_back(longest * random.uniform() * random.randomVectorOnUnitSphere());
    }

    std::vector<SampledWay> ways = devicesWays(structure, origins, displacements);
    std::size_t waysApart = 0;
    double worstWay = 0.0;
    for (std::size_t i = 0; i < origins.size(); ++i)
    {
      const double difference = std::abs(ways[i].radius - structure.segmentBottleneck(origins[i], displacements[i]).radius);
      worstWay = std::max(worstWay, difference);
      if (difference > agreement) ++waysApart;
    }
    EXPECT_EQ(waysApart, 0uz) << name << ": worst by " << worstWay;

    std::vector<double> clearances(origins.size());
    samplingBackendOpenCL().clearances(structure, origins, clearances);
    std::size_t clearancesApart = 0;
    double worstClearance = 0.0;
    for (std::size_t i = 0; i < origins.size(); ++i)
    {
      const double difference = std::abs(clearances[i] - structure.clearance(origins[i]));
      worstClearance = std::max(worstClearance, difference);
      if (difference > agreement) ++clearancesApart;
    }
    EXPECT_EQ(clearancesApart, 0uz) << name << ": worst by " << worstClearance;
  }
}
