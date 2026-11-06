#include <gtest/gtest.h>

import std;

import double3;
import unit_cell;
import crystal;
import pair_interactions;
import voronoi_blocking_spheres;
import apollonius_pore_size_distribution;

// Whether the default pore-size curve was measured on more than the framework's atoms: whether a pocket
// sealed to the probe got a ball. The spheres the curve was measured on do not say so by their number, and
// not by being a different list from the atoms either. An atom wholly inside another is left off that list
// whether or not anything was filled, having no surface of its own, and a ball that swallows an atom takes
// one off as it adds itself. What says it is a sphere on the list that is not an atom.
//
// The lattice is simple cubic with one wide atom to a cell. At a radius of 3.3 Å on an edge of 6 Å the
// faces of the cell are open to the probe and nothing is sealed; at 4 Å they are shut and the room at the
// body centre is a pocket.

namespace
{

constexpr double edge = 6.0;
constexpr double openRadius = 3.3;
constexpr double sealedRadius = 4.0;
constexpr double smallRadius = 0.4;
constexpr double probeRadius = 0.5;

// A wide atom, a small one, and the probe.
PairInteractions interactionsFor(double wideRadius)
{
  const std::array<double, 3> radii{wideRadius, smallRadius, probeRadius};
  PairInteractions interactions;
  interactions.numberOfTypes = 3;
  interactions.names = {"O", "H", "He"};
  interactions.charges = {0.0, 0.0, 0.0};
  for (std::size_t i = 0; i < 3; ++i)
  {
    for (std::size_t j = 0; j < 3; ++j)
    {
      interactions.parameters.push_back(PairParameters{radii[i] + radii[j], 100.0, 0.0});
    }
  }
  return interactions;
}

// The wide atom at `fractional`, and with `withBuried` a small atom well inside it.
Crystal lattice(const double3& fractional, bool withBuried)
{
  Crystal framework;
  framework.name = "pockets_filled_test";
  framework.unitCell = UnitCell(edge, edge, edge);
  framework.mass = 16.0;
  framework.atoms.push_back(CrystalAtom{framework.unitCell.cell * fractional, 0.0, 0});
  framework.fractionalPositions.push_back(fractional);
  if (withBuried)
  {
    const double3 inside = fractional + double3(0.1, 0.05, 0.0);
    framework.atoms.push_back(CrystalAtom{framework.unitCell.cell * inside, 0.0, 1});
    framework.fractionalPositions.push_back(inside);
  }
  return framework;
}

struct ReportCleanup
{
  ~ReportCleanup()
  {
    for (const char* kind : {"psd", "psd-peaks"})
    {
      std::error_code ignored;
      std::filesystem::remove(std::format("pockets_filled_test.apollonius.{}.txt", kind), ignored);
    }
  }
};

ApolloniusPoreSizeDistribution sweep(double wideRadius, const double3& fractional, bool withBuried)
{
  ApolloniusPoreSizeDistribution driver;
  driver.run(interactionsFor(wideRadius), lattice(fractional, withBuried), "He", 8.0, 40, 1, 0.0, true);
  return driver;
}

}  // namespace

TEST(apollonius_pockets_filled, nothing_sealed_is_nothing_filled)
{
  ReportCleanup cleanup;
  const ApolloniusPoreSizeDistribution driver = sweep(openRadius, double3(0.0, 0.0, 0.0), false);
  EXPECT_EQ(driver.sweptSpheres.size(), 1uz);
  EXPECT_FALSE(driver.pocketsFilled);
}

TEST(apollonius_pockets_filled, an_atom_given_outside_the_cell_is_still_only_an_atom)
{
  // The positions go to the diagram as they are given, and come back on the list as they were given.
  ReportCleanup cleanup;
  const ApolloniusPoreSizeDistribution driver = sweep(openRadius, double3(1.25, -0.5, 2.75), false);
  EXPECT_EQ(driver.sweptSpheres.size(), 1uz);
  EXPECT_FALSE(driver.pocketsFilled);
}

TEST(apollonius_pockets_filled, an_atom_buried_in_another_is_left_off_the_list_and_fills_nothing)
{
  ReportCleanup cleanup;
  const ApolloniusPoreSizeDistribution driver = sweep(openRadius, double3(0.0, 0.0, 0.0), true);

  // Two atoms, one sphere: the small one has no surface of its own.
  ASSERT_EQ(driver.sweptSpheres.size(), 1uz);
  EXPECT_DOUBLE_EQ(driver.sweptSpheres[0].radius, openRadius);
  EXPECT_FALSE(driver.pocketsFilled);
}

TEST(apollonius_pockets_filled, a_sealed_pocket_is_filled_whether_or_not_an_atom_is_buried)
{
  ReportCleanup cleanup;
  for (bool withBuried : {false, true})
  {
    const ApolloniusPoreSizeDistribution driver = sweep(sealedRadius, double3(0.0, 0.0, 0.0), withBuried);

    // The atom and the ball over the pocket at the body centre.
    ASSERT_EQ(driver.sweptSpheres.size(), 2uz) << "buried atom: " << withBuried;
    EXPECT_TRUE(driver.pocketsFilled) << "buried atom: " << withBuried;

    std::size_t atBodyCentre = 0;
    for (const BlockingSphere& sphere : driver.sweptSpheres)
    {
      const double3 fromCentre = sphere.centerFractional - double3(0.5, 0.5, 0.5);
      const double3 wrapped(fromCentre.x - std::round(fromCentre.x), fromCentre.y - std::round(fromCentre.y),
                            fromCentre.z - std::round(fromCentre.z));
      if (wrapped.length() < 1.0e-3) ++atBodyCentre;
    }
    EXPECT_EQ(atBodyCentre, 1uz) << "buried atom: " << withBuried;
  }
}
