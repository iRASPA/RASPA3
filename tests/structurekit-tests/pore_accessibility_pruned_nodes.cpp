#include <gtest/gtest.h>

import std;

import double3;
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
