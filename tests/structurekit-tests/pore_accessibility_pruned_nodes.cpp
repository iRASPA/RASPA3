#include <gtest/gtest.h>

import std;

import double3;
import int3;
import double3x3;
import unit_cell;
import randomnumbers;
import voronoi_network;
import channel_analysis;
import pore_accessibility;

// A point outside every atom that the classifier decides is decided by a node, and the answer is that
// node's pore. A node with no room for the probe has no pore, so it must not be the one that decides: the
// answer would be void that is in no pore, counted as inaccessible and belonging to no pocket.
//
// On a radical diagram geometry nearly rules it out by itself. The line-of-sight test takes a vertex of
// the nearest atom's cell lying beyond the point, such a vertex is outside that atom, and a vertex of an
// atom's cell that is outside its own atom is outside them all. The exception is a node that several
// vertices were merged into: it has the least room of them, and can have none while a cell's own vertex,
// a little way off, still has some. That takes vertices closer together than the merging distance and
// atoms that reach them to within the same distance, which is what the structure below is made to do: a
// simple cubic lattice, eight cells of it, each atom moved a few thousandths of an Ångström so that the
// one vertex at every cube corner comes apart into several, and the radius taken in small steps through
// the distance from an atom to its corners.

namespace
{

struct Lattice
{
  UnitCell box{8.0, 8.0, 8.0};
  std::vector<double3> fractionalPositions;
};

Lattice nearlyCubicLattice()
{
  RandomNumber random{std::optional<std::size_t>(0)};
  Lattice lattice;
  for (int x = 0; x < 2; ++x)
  {
    for (int y = 0; y < 2; ++y)
    {
      for (int z = 0; z < 2; ++z)
      {
        lattice.fractionalPositions.push_back(double3((x + 0.25 + 0.001 * (random.uniform() - 0.5)) / 2.0,
                                                      (y + 0.25 + 0.001 * (random.uniform() - 0.5)) / 2.0,
                                                      (z + 0.25 + 0.001 * (random.uniform() - 0.5)) / 2.0));
      }
    }
  }
  return lattice;
}

// Void, decided, and in no pore: the answer there is no place for.
bool voidInNoPore(const PointClassification& classification)
{
  return !classification.inside && !classification.resample && classification.poreId < 0;
}

}  // namespace

TEST(pore_accessibility, a_node_with_no_room_for_the_probe_does_not_decide_a_point)
{
  const Lattice lattice = nearlyCubicLattice();
  const double cornerDistance = 0.5 * std::sqrt(3.0) * 4.0;

  std::size_t pointsPut = 0;
  for (int step = -100; step <= 100; ++step)
  {
    const std::vector<double> radii(lattice.fractionalPositions.size(), cornerDistance + 1.0e-4 * step);
    const PoreAccessibility accessibility =
        PoreAccessibility::create(lattice.box, lattice.fractionalPositions, radii, 0.0);

    for (std::size_t atom = 0; atom < accessibility.atomPositions.size(); ++atom)
    {
      const double radius = accessibility.atomRadii[atom];
      for (const auto& [node, vertex] : accessibility.network.atomNodeVectors[atom])
      {
        // A node of no pore whose vertex on this cell is nevertheless outside the atom: the points on the
        // way out to that vertex see it beyond them and nothing nearer.
        if (accessibility.channels.nodePoreId[node] >= 0) continue;
        const double gap = vertex.length() - radius;
        if (gap <= 0.0) continue;

        const double3 outward = vertex * (1.0 / vertex.length());
        for (const double beyond : {1.0e-9, 0.5 * gap})
        {
          const double3 point = accessibility.atomPositions[atom] + outward * (radius + beyond);
          const PointClassification classification = accessibility.classify(point);
          ++pointsPut;
          EXPECT_FALSE(voidInNoPore(classification))
              << "step " << step << ", atom " << atom << ", node " << node << ", " << gap << " Å outside the atom";
          if (classification.poreId >= 0)
          {
            EXPECT_LT(static_cast<std::size_t>(classification.poreId), accessibility.channels.pores.size());
          }
        }
      }
    }
  }

  // Without such points the test would say nothing, and the structure is there to make them.
  EXPECT_GT(pointsPut, 1000uz);
}

TEST(pore_accessibility, void_that_is_decided_is_always_in_a_pore)
{
  // The same over points drawn anywhere, on the lattice above and on atoms of two sizes thrown at random,
  // where a good many nodes of the radical diagram stand inside atoms.
  const Lattice lattice = nearlyCubicLattice();
  const double cornerDistance = 0.5 * std::sqrt(3.0) * 4.0;

  std::vector<PoreAccessibility> structures;
  for (int step : {-40, -25, -10, 0, 10})
  {
    const std::vector<double> radii(lattice.fractionalPositions.size(), cornerDistance + 1.0e-4 * step);
    structures.push_back(PoreAccessibility::create(lattice.box, lattice.fractionalPositions, radii, 0.0));
  }
  for (std::size_t seed = 0; seed < 5; ++seed)
  {
    RandomNumber random{std::optional<std::size_t>(seed)};
    std::vector<double3> positions;
    std::vector<double> radii;
    for (std::size_t i = 0; i < 60; ++i)
    {
      positions.push_back(double3(random.uniform(), random.uniform(), random.uniform()));
      radii.push_back(random.uniform() < 0.5 ? 0.5 : 2.0);
    }
    structures.push_back(PoreAccessibility::create(UnitCell{12.0, 12.0, 12.0}, positions, radii, 0.4));
  }

  RandomNumber random{std::optional<std::size_t>(3)};
  for (std::size_t s = 0; s < structures.size(); ++s)
  {
    const PoreAccessibility& accessibility = structures[s];
    std::size_t inNoPore = 0;
    for (std::size_t i = 0; i < 20000; ++i)
    {
      const double3 point = accessibility.unitCell.cell * double3(random.uniform(), random.uniform(), random.uniform());
      if (voidInNoPore(accessibility.classify(point))) ++inNoPore;
    }
    // On the surfaces as well, which is where a sweep for area puts every one of its points.
    for (std::size_t atom = 0; atom < accessibility.atomPositions.size(); ++atom)
    {
      for (std::size_t i = 0; i < 2000; ++i)
      {
        const double3 point =
            accessibility.atomPositions[atom] + accessibility.atomRadii[atom] * random.randomVectorOnUnitSphere();
        if (accessibility.overlapsAtom(point, atom)) continue;
        if (voidInNoPore(accessibility.classify(point))) ++inNoPore;
      }
    }
    EXPECT_EQ(inNoPore, 0uz) << "structure " << s;
  }
}

// A point whose fractional coordinate is just below zero, from -2^-54, which is where a point on a face of the
// cell lands after rounding: the cell it was wrapped out of and the place it was wrapped to add up to the point.
// They used to be a whole cell apart there, so that the images of the neighbours were named one cell off and a
// point on a sphere was buried by that sphere itself.
TEST(pore_accessibility, a_point_just_below_a_cell_face_is_wrapped_consistently)
{
  const UnitCell box{10.0, 10.0, 10.0};
  const std::vector<double3> fractionalPositions{double3(0.0, 0.0, 0.0), double3(0.5, 0.5, 0.5)};
  const std::vector<double> radii{1.5, 1.5};
  const PoreAccessibility accessibility = PoreAccessibility::create(box, fractionalPositions, radii, 0.0);

  for (const double below : {-1.0e-16, -3.0e-16, -5.0e-16})
  {
    // inside the sphere of the atom at the origin, on the far side of the face x = 0 by a rounding error
    const double3 point(below, 1.2, 0.3);
    ASSERT_LT(box.inverseCell.ax * point.x, 0.0);

    // drawn on the atom at the origin, the point is not buried by that atom, and no other atom reaches it
    EXPECT_FALSE(accessibility.overlapsAtom(point, 0)) << below;
    EXPECT_TRUE(accessibility.overlapsAtom(point, 1)) << below;

    const std::vector<NeighbourImage> neighbours = accessibility.neighbourAtomImages(point, 6.0);
    ASSERT_FALSE(neighbours.empty());
    for (const NeighbourImage& neighbour : neighbours)
    {
      const double3 centre =
          accessibility.atomPositions[neighbour.index] +
          box.cell * double3(static_cast<double>(neighbour.image.x), static_cast<double>(neighbour.image.y),
                             static_cast<double>(neighbour.image.z));
      EXPECT_LT((centre - point - neighbour.delta).length(), 1.0e-9)
          << below << ": atom " << neighbour.index << " image " << neighbour.image.x << " " << neighbour.image.y
          << " " << neighbour.image.z;
    }
  }
}

// A cell without volume is refused: the walks over shells of bins never ended on it.
TEST(pore_accessibility, a_cell_without_volume_is_refused)
{
  VoronoiNetwork network;
  network.unitCell = UnitCell(double3x3(double3(5.0, 0.0, 0.0), double3(0.0, 5.0, 0.0), double3(5.0, 5.0, 0.0)));
  network.atomPositionsFractional = {double3(0.1, 0.2, 0.3)};
  network.atomRadii = {1.0};
  EXPECT_THROW(PoreAccessibility::createFromNetwork(network, PoreAccessibility::Metric::Power, 0.0),
               std::invalid_argument);
}
