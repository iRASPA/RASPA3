#include <gtest/gtest.h>

import std;

import double3;
import double3x3;
import unit_cell;
import randomnumbers;
import skvoronoi;
import voronoi_network;

// What the network keeps per atom is the corners of that atom's own cell, each as the node it is and the
// vector from the atom to it. In a crystal one node is often several corners of the same cell: the corners
// are a lattice vector apart, so they are one point of the home cell, but they stand in different places
// round the atom. The line-of-sight test of the accessibility needs every one of them, since it asks which
// corner lies beyond a sample point as seen from the atom, and a corner on the far side of the cell is no
// stand-in for the one on the near side. So a node is to be kept once for every corner it is, not once.

namespace
{

// The lattice vector nearest to the given vector, and how far the vector is from it.
double distanceFromLattice(const UnitCell& box, const double3& vector)
{
  double3 fractional = box.inverseCell * vector;
  fractional.x -= std::round(fractional.x);
  fractional.y -= std::round(fractional.y);
  fractional.z -= std::round(fractional.z);
  return (box.cell * fractional).length();
}

}  // namespace

TEST(voronoi_network, a_single_atom_in_a_cube_keeps_all_eight_corners_of_its_cell)
{
  const double a = 6.0;
  UnitCell box{a, a, a};
  VoronoiNetwork network = VoronoiNetwork::create(box, {double3(0.3, 0.6, 0.1)}, {1.0});

  // Every corner of the cube is the same point of the lattice, so there is one node, held eight times.
  ASSERT_EQ(network.nodes.size(), 1uz);
  ASSERT_EQ(network.atomNodeVectors.size(), 1uz);
  const auto& corners = network.atomNodeVectors[0];
  ASSERT_EQ(corners.size(), 8uz);

  std::array<int, 8> timesSeen{};
  for (const auto& [node, vector] : corners)
  {
    EXPECT_EQ(node, 0uz);
    EXPECT_NEAR(std::abs(vector.x), 0.5 * a, 1.0e-9);
    EXPECT_NEAR(std::abs(vector.y), 0.5 * a, 1.0e-9);
    EXPECT_NEAR(std::abs(vector.z), 0.5 * a, 1.0e-9);
    ++timesSeen[(vector.x > 0.0 ? 1uz : 0uz) + (vector.y > 0.0 ? 2uz : 0uz) + (vector.z > 0.0 ? 4uz : 0uz)];
  }
  for (std::size_t corner = 0; corner < 8; ++corner)
  {
    EXPECT_EQ(timesSeen[corner], 1) << "corner " << corner;
  }
}

TEST(voronoi_network, a_node_that_is_two_corners_of_one_cell_is_kept_with_both_vectors)
{
  // Body-centred cubic: each cell is a truncated octahedron of 24 corners, and the lattice has 12 such
  // points to a cube, so each atom meets every node twice, on opposite sides of itself.
  const double a = 5.0;
  UnitCell box{a, a, a};
  VoronoiNetwork network =
      VoronoiNetwork::create(box, {double3(0.1, 0.2, 0.3), double3(0.6, 0.7, 0.8)}, {0.8, 0.8});

  ASSERT_EQ(network.nodes.size(), 12uz);
  for (std::size_t atom = 0; atom < 2; ++atom)
  {
    const auto& corners = network.atomNodeVectors[atom];
    ASSERT_EQ(corners.size(), 24uz) << "atom " << atom;

    std::vector<std::vector<double3>> vectorsOfNode(network.nodes.size());
    for (const auto& [node, vector] : corners) vectorsOfNode[node].push_back(vector);

    for (std::size_t node = 0; node < network.nodes.size(); ++node)
    {
      ASSERT_EQ(vectorsOfNode[node].size(), 2uz) << "atom " << atom << ", node " << node;
      const double3 apart = vectorsOfNode[node][1] - vectorsOfNode[node][0];
      EXPECT_GT(apart.length(), 0.5 * a) << "atom " << atom << ", node " << node;
      EXPECT_LT(distanceFromLattice(box, apart), 1.0e-6) << "atom " << atom << ", node " << node;
    }
  }
}

TEST(voronoi_network, every_corner_of_every_cell_is_kept_and_stands_where_its_node_is)
{
  // A skewed cell with few atoms of unequal radii, so that cells reach round the cell and meet their own
  // nodes more than once, held against the cells the diagram itself gives.
  UnitCell box{8.0, 9.0, 10.0, 75.0 * std::numbers::pi / 180.0, 100.0 * std::numbers::pi / 180.0,
               65.0 * std::numbers::pi / 180.0};
  RandomNumber random{std::optional<std::size_t>(11)};
  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  for (std::size_t i = 0; i < 7; ++i)
  {
    fractionalPositions.push_back(double3(random.uniform(), random.uniform(), random.uniform()));
    radii.push_back(0.6 + 0.8 * random.uniform());
  }

  VoronoiNetwork network = VoronoiNetwork::create(box, fractionalPositions, radii);
  SKVoronoi diagram(box.cell, fractionalPositions, radii);
  std::vector<SKVoronoiCell> cells = diagram.computeAllCells();
  ASSERT_EQ(network.atomNodeVectors.size(), cells.size());

  std::size_t nodesMetMoreThanOnce = 0;
  for (std::size_t atom = 0; atom < cells.size(); ++atom)
  {
    const auto& corners = network.atomNodeVectors[atom];
    ASSERT_EQ(corners.size(), cells[atom].verticesCartesian.size()) << "atom " << atom;

    const double3 atomPosition = box.cell * network.atomPositionsFractional[atom];
    std::vector<std::size_t> timesMet(network.nodes.size(), 0);
    for (std::size_t v = 0; v < corners.size(); ++v)
    {
      const auto& [node, vector] = corners[v];
      ASSERT_LT(node, network.nodes.size());
      EXPECT_LT((vector - cells[atom].verticesCartesian[v]).length(), 1.0e-12)
          << "atom " << atom << ", corner " << v;
      // The node is the corner brought into the home cell, to within the distance nodes are merged over.
      EXPECT_LT(distanceFromLattice(box, atomPosition + vector - network.nodes[node].position), 0.05)
          << "atom " << atom << ", corner " << v;
      ++timesMet[node];
    }
    for (std::size_t times : timesMet)
    {
      if (times > 1) ++nodesMetMoreThanOnce;
    }
  }
  // Without this the test would say nothing about nodes that are several corners of a cell.
  EXPECT_GT(nodesMetMoreThanOnce, 0uz);
}
