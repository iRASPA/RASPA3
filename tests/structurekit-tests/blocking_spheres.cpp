#include <gtest/gtest.h>

import std;

import double3;
import unit_cell;
import pore_accessibility;
import voronoi_blocking_spheres;

// The sampled blocking spheres, and the allowances made for the sampling.
//
// A sphere may reach no nearer the channel than the nearest accessible sample, less how far that sample can
// stand off the channel's edge, which is the spacing of the accessible samples. The samples are thrown over the
// whole cell, so that spacing is set by the cell's volume and the number thrown, and not by how much of the
// cell turned out to be accessible.

namespace
{

struct Lattice
{
  UnitCell box{20.0, 20.0, 20.0};
  std::vector<double3> fractionalPositions;
  std::vector<double> radii;

  void add(const double3& cartesian, double radius)
  {
    this->fractionalPositions.push_back(double3::fract(this->box.inverseCell * cartesian));
    this->radii.push_back(radius);
  }
};

// A sealed cage of `count` spheres on a shell about `centre`, placed on the Fibonacci spiral.
std::vector<double3> fibonacciShell(const double3& centre, double shellRadius, std::size_t count)
{
  std::vector<double3> positions;
  const double golden = std::numbers::pi * (1.0 + std::sqrt(5.0));
  for (std::size_t i = 0; i < count; ++i)
  {
    double fraction = (static_cast<double>(i) + 0.5) / static_cast<double>(count);
    double cosinePolar = 1.0 - 2.0 * fraction;
    double sinePolar = std::sqrt(std::max(0.0, 1.0 - cosinePolar * cosinePolar));
    double azimuth = golden * (static_cast<double>(i) + 0.5);
    positions.push_back(centre + double3(shellRadius * sinePolar * std::cos(azimuth),
                                         shellRadius * sinePolar * std::sin(azimuth), shellRadius * cosinePolar));
  }
  return positions;
}

const double3 cageCentre(10.0, 10.0, 10.0);
constexpr double atomRadius = 1.8;
constexpr double probe = 0.5;

// A cage in the middle of the cell, with the void round it as the channel. With `filledBeyond` set, the
// corners of the cell farther than that from the cage's centre are packed solid, which takes away a good
// part of the channel but none of it anywhere near the cage.
Lattice cageInACell(std::optional<double> filledBeyond)
{
  Lattice lattice;
  for (const double3& position : fibonacciShell(cageCentre, 6.0, 128)) lattice.add(position, atomRadius);
  if (filledBeyond.has_value())
  {
    for (int i = 0; i < 10; ++i)
      for (int j = 0; j < 10; ++j)
        for (int k = 0; k < 10; ++k)
        {
          double3 position(2.0 * i, 2.0 * j, 2.0 * k);
          double3 delta = lattice.box.applyPeriodicBoundaryConditions(position - cageCentre);
          if (delta.length() >= filledBeyond.value()) lattice.add(position, atomRadius);
        }
  }
  return lattice;
}

std::vector<BlockingSphere> spheresOf(const Lattice& lattice, std::size_t samples)
{
  PoreAccessibility accessibility =
      PoreAccessibility::create(lattice.box, lattice.fractionalPositions, lattice.radii, probe);
  return computeBlockingSpheres(accessibility, samples);
}

}  // namespace


TEST(blocking_spheres, the_channel_margin_is_the_spacing_of_the_samples_and_not_of_the_accessible_ones)
{
  // Few enough samples that the margin, and not the far wall of the cavity, is what bounds the sphere.
  const std::size_t samples = 300;

  std::vector<BlockingSphere> open = spheresOf(cageInACell(std::nullopt), samples);
  std::vector<BlockingSphere> cramped = spheresOf(cageInACell(12.0), samples);

  ASSERT_FALSE(open.empty());
  ASSERT_EQ(cramped.size(), open.size());

  // The same samples land in the cavity and the same ones land just outside it, so the spheres are the same.
  // Taking the spacing from the count of accessible samples alone instead would have widened the margin by the
  // cube root of the part of the cell the corners took away, and held the sphere short by that.
  for (std::size_t i = 0; i < open.size(); ++i)
  {
    EXPECT_NEAR(cramped[i].radius, open[i].radius, 1.0e-6) << "sphere " << i;
    EXPECT_NEAR((cramped[i].centerFractional - open[i].centerFractional).length(), 0.0, 1.0e-6) << "sphere " << i;
  }

  // And the sphere is bounded by the channel, not by the free ball at its centre or by its own far wall: it
  // reaches past the cavity's room, and stops short of the cage's outer surface.
  const double room = 6.0 - atomRadius - probe;
  EXPECT_GT(open.front().radius, room + 1.0);
  EXPECT_LT(open.front().radius, 6.0 + atomRadius + probe);
}
