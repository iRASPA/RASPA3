#include <gtest/gtest.h>

import std;

import units;
import connectivity_table;
import bend_potential;

// The bonded terms enumerated from the connectivity of a molecule, and the printing of a bend.

namespace
{
ConnectivityTable ring(std::size_t n)
{
  ConnectivityTable table(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    const std::size_t j = (i + 1) % n;
    table[i, j] = true;
    table[j, i] = true;
  }
  return table;
}
}  // namespace

// A three-membered ring has three bends and no torsion: a torsion needs four different atoms, and on the ring the
// neighbour that would extend a bend is the far end of the bend itself.
TEST(connectivity_table, three_membered_ring_has_no_torsions)
{
  const ConnectivityTable triangle = ring(3);
  EXPECT_EQ(triangle.findAllBonds().size(), 3uz);
  EXPECT_EQ(triangle.findAllBends().size(), 3uz);
  EXPECT_TRUE(triangle.findAllTorsions().empty());
}

TEST(connectivity_table, ring_torsions_have_four_different_atoms)
{
  for (const std::size_t n : {4uz, 5uz, 6uz})
  {
    const ConnectivityTable cycle = ring(n);
    const std::vector<std::array<std::size_t, 4>> torsions = cycle.findAllTorsions();
    EXPECT_EQ(torsions.size(), n) << "ring of " << n;
    for (const std::array<std::size_t, 4>& torsion : torsions)
    {
      std::array<std::size_t, 4> sorted = torsion;
      std::ranges::sort(sorted);
      EXPECT_EQ(std::ranges::unique(sorted).begin(), sorted.end()) << "ring of " << n;
    }
  }
  // a chain of five: two torsions
  ConnectivityTable chain(5);
  for (std::size_t i = 0; i + 1 < 5; ++i)
  {
    chain[i, i + 1] = true;
    chain[i + 1, i] = true;
  }
  EXPECT_EQ(chain.findAllTorsions().size(), 2uz);
}

// The printed equilibrium angle of every bend with one is in degrees.
TEST(bend_potential, print_reports_the_equilibrium_angle_in_degrees)
{
  const BendPotential harmonic({0, 1, 2}, BendType::Harmonic, {62500.0, 112.0});
  const BendPotential quartic({0, 1, 2}, BendType::Quartic, {62500.0, 112.0, 100.0, 10.0});
  const BendPotential cffQuartic({0, 1, 2}, BendType::CFF_Quartic, {62500.0, 112.0, 100.0, 10.0});
  for (const BendPotential& bend : {harmonic, quartic, cffQuartic})
  {
    const std::string text = bend.print();
    EXPECT_NE(text.find("p_1=112 [degrees]"), std::string::npos) << text;
  }
}
