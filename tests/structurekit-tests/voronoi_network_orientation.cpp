#include <gtest/gtest.h>

import std;

import double3;
import double3x3;
import unit_cell;
import randomnumbers;
import voronoi_network;
import voronoi_pore_diameters;

// The network of a crystal is the crystal's: it does not depend on which way the cell is turned, nor on the
// order the atoms are listed in. The one place it could is where the vertices of the diagram are gathered into
// nodes, vertices nearer each other than a tolerance being one node, and it used to in two ways.
//
// Nearer was asked of each Cartesian component of the difference, which is a cube round the node, reaching
// 1.7 times as far along its diagonals as along the axes: two vertices 0.03 Å apart were one node in a cell
// turned one way and two in the same cell turned another.
//
// And each vertex went to the first node already made that was near enough. Of three vertices in a row, each
// near the next and the ends not near each other, the middle one then goes with whichever end came first, and
// the node stands where its first vertex stood; the order is that of the atoms in the input and of the corners
// in each cell. The same atoms listed in another order gave other nodes in nearly half of the structures here,
// and a free sphere diameter that differed by up to 0.007 Å.
//
// Both take vertices a few hundredths of an Ångström apart, which is what a nearly symmetric structure has:
// here a simple cubic lattice, eight cells of it, with each atom moved by up to a few hundredths, so that the
// one vertex where eight cells meet comes apart into a cluster of them about as wide as the tolerance.

namespace
{

constexpr double edge = 8.0;

// The same turn of the axes for every structure: about z, then about x, then about y, by angles of no
// particular kind.
double3 turned(const double3& v)
{
  const double aboutZ = 0.7;
  const double aboutX = 0.4;
  const double aboutY = 1.1;
  const double3 first(std::cos(aboutZ) * v.x - std::sin(aboutZ) * v.y, std::sin(aboutZ) * v.x + std::cos(aboutZ) * v.y,
                      v.z);
  const double3 second(first.x, std::cos(aboutX) * first.y - std::sin(aboutX) * first.z,
                       std::sin(aboutX) * first.y + std::cos(aboutX) * first.z);
  return double3(std::cos(aboutY) * second.x + std::sin(aboutY) * second.z, second.y,
                 -std::sin(aboutY) * second.x + std::cos(aboutY) * second.z);
}

std::vector<double3> nearlyCubic(double amplitude, std::size_t seed)
{
  RandomNumber random{std::optional<std::size_t>(seed)};
  std::vector<double3> fractionalPositions;
  for (int x = 0; x < 2; ++x)
  {
    for (int y = 0; y < 2; ++y)
    {
      for (int z = 0; z < 2; ++z)
      {
        fractionalPositions.push_back(double3((x + 0.25) / 2.0 + amplitude * (random.uniform() - 0.5) / edge,
                                              (y + 0.25) / 2.0 + amplitude * (random.uniform() - 0.5) / edge,
                                              (z + 0.25) / 2.0 + amplitude * (random.uniform() - 0.5) / edge));
      }
    }
  }
  return fractionalPositions;
}

// What two networks of one crystal have to have in common.
void expectTheSameNetwork(const VoronoiNetwork& one, const VoronoiNetwork& other, const std::string& which)
{
  EXPECT_EQ(one.nodes.size(), other.nodes.size()) << which;
  EXPECT_EQ(one.edges.size(), other.edges.size()) << which;

  const PoreDiameters mine = PoreDiameters::compute(one);
  const PoreDiameters theirs = PoreDiameters::compute(other);
  EXPECT_NEAR(mine.includedSphereDiameter, theirs.includedSphereDiameter, 1.0e-9) << which;
  EXPECT_NEAR(mine.freeSphereDiameter, theirs.freeSphereDiameter, 1.0e-9) << which;
  EXPECT_NEAR(mine.includedAlongFreePathDiameter, theirs.includedAlongFreePathDiameter, 1.0e-9) << which;
}

}  // namespace

TEST(voronoi_network, a_nearly_symmetric_crystal_has_the_same_network_however_it_is_turned)
{
  const UnitCell upright(double3x3(double3(edge, 0.0, 0.0), double3(0.0, edge, 0.0), double3(0.0, 0.0, edge)));
  const UnitCell rotated(double3x3(turned(double3(edge, 0.0, 0.0)), turned(double3(0.0, edge, 0.0)),
                                   turned(double3(0.0, 0.0, edge))));

  for (const double amplitude : {0.004, 0.008, 0.012, 0.016, 0.02, 0.03})
  {
    for (std::size_t seed = 0; seed < 40; ++seed)
    {
      // The fractional positions are the same in both cells, so the two are one crystal and nothing else.
      const std::vector<double3> fractionalPositions = nearlyCubic(amplitude, seed);
      const std::vector<double> radii(fractionalPositions.size(), 1.0);
      expectTheSameNetwork(VoronoiNetwork::create(upright, fractionalPositions, radii),
                           VoronoiNetwork::create(rotated, fractionalPositions, radii),
                           std::format("amplitude {}, seed {}", amplitude, seed));
    }
  }
}

TEST(voronoi_network, a_nearly_symmetric_crystal_has_the_same_network_whatever_order_its_atoms_are_in)
{
  const UnitCell box(edge, edge, edge);

  for (const double amplitude : {0.004, 0.008, 0.012, 0.016, 0.02, 0.03})
  {
    for (std::size_t seed = 0; seed < 40; ++seed)
    {
      const std::vector<double3> fractionalPositions = nearlyCubic(amplitude, seed);
      const std::vector<double> radii(fractionalPositions.size(), 1.0);

      // Back to front, with two of them changed over besides.
      std::vector<double3> reordered(fractionalPositions.rbegin(), fractionalPositions.rend());
      std::swap(reordered[1], reordered[5]);

      expectTheSameNetwork(VoronoiNetwork::create(box, fractionalPositions, radii),
                           VoronoiNetwork::create(box, reordered, radii),
                           std::format("amplitude {}, seed {}", amplitude, seed));
    }
  }
}
