#include <gtest/gtest.h>

import std;

import double3;
import unit_cell;
import randomnumbers;
import pore_accessibility;
import apollonius_accessibility;
import voronoi_network;

// The radical network's nodes walk from their vertices up to the peaks of the clearance, and the walk has to
// get there. The Apollonius network has the peaks exactly, each a vertex of its own, so the two can be held to
// each other: wherever a radical vertex stands inside the free ball of a true peak it is in that peak's own
// pocket, and the walk from it, or from one of the others in the same ball, must reach a clearance at least
// the peak's. On unequal radii the vertices sit well off the peaks, often near a wall with little room, and a
// walk bounded by that room or blind to the atoms the next sphere touches stops short.

namespace
{

struct Structure
{
  UnitCell box{12.0, 12.0, 12.0};
  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
};

Structure randomStructure(std::size_t seed, std::size_t count, bool bimodal)
{
  RandomNumber random{std::optional<std::size_t>(seed)};
  Structure structure;
  for (std::size_t i = 0; i < count; ++i)
  {
    structure.fractionalPositions.push_back(double3(random.uniform(), random.uniform(), random.uniform()));
    double u = random.uniform();
    structure.radii.push_back(bimodal ? (u < 0.5 ? 0.5 : 2.5) : 0.5 + 2.0 * u);
  }
  return structure;
}

// Whether every small step from the point lowers the clearance.
bool isPeak(const PoreAccessibility& accessibility, const double3& point, RandomNumber& random)
{
  const double here = accessibility.clearance(point);
  for (int k = 0; k < 60; ++k)
  {
    double3 direction(random.uniform() - 0.5, random.uniform() - 0.5, random.uniform() - 0.5);
    direction = direction * (1.0e-3 / direction.length());
    if (accessibility.clearance(point + direction) > here + 1.0e-12) return false;
  }
  return true;
}

}  // namespace


TEST(voronoi_network, the_ascent_from_a_vertex_in_a_pockets_free_ball_reaches_the_pockets_peak)
{
  RandomNumber directions{std::optional<std::size_t>(7)};
  std::size_t peaks = 0;
  std::size_t checked = 0;
  for (std::size_t seed = 1; seed <= 24; ++seed)
  {
    const bool bimodal = (seed % 2 == 0);
    Structure structure = randomStructure(seed, bimodal ? 40 : 24, bimodal);
    PoreAccessibility radical;
    PoreAccessibility apollonius;
    try
    {
      radical = PoreAccessibility::create(structure.box, structure.fractionalPositions, structure.radii, 0.0);
      apollonius = ApolloniusAccessibility::create(structure.box, structure.fractionalPositions, structure.radii, 0.0)
                       .accessibility;
    }
    catch (const std::exception&)
    {
      continue;  // a degenerate diagram; another seed serves as well
    }

    for (const VoronoiNode& peak : apollonius.network.nodes)
    {
      if (peak.radius <= 1.0e-3 || !isPeak(radical, peak.position, directions)) continue;
      ++peaks;

      double reached = -1.0;
      bool anyInside = false;
      for (const VoronoiNode& node : radical.network.nodes)
      {
        if (node.radius <= 0.0) continue;
        double3 delta = structure.box.applyPeriodicBoundaryConditions(node.position - peak.position);
        if (delta.length() >= peak.radius) continue;
        anyInside = true;
        reached = std::max(reached, node.maximalRadius);
      }
      if (!anyInside) continue;
      ++checked;
      EXPECT_GE(reached, peak.radius - 1.0e-6)
          << "seed " << seed << ": a peak of clearance " << peak.radius << " at (" << peak.position.x << ", "
          << peak.position.y << ", " << peak.position.z << ") has a vertex in its ball that climbed only to "
          << reached;
    }
  }
  EXPECT_GT(checked, 500uz) << "of " << peaks << " peaks";
}


// The largest sphere is the one the walk is for, and the Apollonius network has it exactly.
TEST(voronoi_network, the_largest_included_sphere_matches_the_apollonius_network)
{
  for (std::size_t seed = 1; seed <= 12; ++seed)
  {
    Structure structure = randomStructure(seed + 50, 24, seed % 2 == 0);
    double radical = 0.0;
    double apollonius = 0.0;
    try
    {
      radical = PoreAccessibility::create(structure.box, structure.fractionalPositions, structure.radii, 0.0)
                    .network.largestIncludedSphereDiameter();
      apollonius = ApolloniusAccessibility::create(structure.box, structure.fractionalPositions, structure.radii, 0.0)
                       .accessibility.network.largestIncludedSphereDiameter();
    }
    catch (const std::exception&)
    {
      continue;
    }
    EXPECT_NEAR(radical, apollonius, 1.0e-9) << "seed " << seed;
  }
}
