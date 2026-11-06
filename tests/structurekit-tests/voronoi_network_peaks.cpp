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


// The tangent spheres the walk steps between, in the two arrangements its closed form is weakest in.

namespace
{

// By how much a sphere fails to touch the worst of the four it is said to be tangent to.
double worstTangency(const std::array<double3, 4>& centres, const std::array<double, 4>& radii,
                     const ApolloniusSphere& sphere)
{
  double worst = 0.0;
  for (std::size_t i = 0; i < 4; ++i)
  {
    worst = std::max(worst, std::abs((sphere.centre - centres[i]).length() - (radii[i] + sphere.radius)));
  }
  return worst;
}

}  // namespace

TEST(voronoi_network, four_centres_in_a_plane_with_unequal_radii_have_a_pair_of_tangent_spheres)
{
  // The centre cannot be written as a function of the radius here, there being nothing to tell one side of the
  // plane from the other, and the solve that writes it so used to return nothing. But the spheres are there:
  // the radii are taken off one of them, a sphere of radius one standing 1.5 above the plane, and its mirror
  // image below is the other.
  const std::array<double3, 4> centres{double3(0.0, 0.0, 0.0), double3(4.0, 0.3, 0.0), double3(0.2, 4.1, 0.0),
                                       double3(4.3, 4.4, 0.0)};
  const double3 above(2.1, 2.0, 1.5);
  std::array<double, 4> radii{};
  for (std::size_t i = 0; i < 4; ++i) radii[i] = (above - centres[i]).length() - 1.0;

  const std::vector<ApolloniusSphere> spheres = apolloniusTangentSpheres(centres, radii);
  ASSERT_EQ(spheres.size(), 2uz);

  std::size_t aboveThePlane = 0;
  for (const ApolloniusSphere& sphere : spheres)
  {
    EXPECT_NEAR(sphere.radius, 1.0, 1.0e-9);
    EXPECT_NEAR(sphere.centre.x, above.x, 1.0e-9);
    EXPECT_NEAR(sphere.centre.y, above.y, 1.0e-9);
    EXPECT_NEAR(std::abs(sphere.centre.z), above.z, 1.0e-9);
    EXPECT_LT(worstTangency(centres, radii, sphere), 1.0e-9);
    if (sphere.centre.z > 0.0) ++aboveThePlane;
  }
  EXPECT_EQ(aboveThePlane, 1uz);

  // Radii all alike in a plane are another matter, with a whole family of tangent spheres or none: still nothing.
  EXPECT_TRUE(apolloniusTangentSpheres(centres, {1.0, 1.0, 1.0, 1.0}).empty());
}

TEST(voronoi_network, the_tangent_sphere_is_as_exact_when_its_twin_is_nearly_a_plane)
{
  // Four spheres standing on a table have the table for one of their two tangent spheres, a sphere of no
  // curvature, and the one nestled above them for the other. Lift one of the four off the table by a little
  // and the table becomes a sphere of enormous radius, the quadratic's leading coefficient being what little
  // the lift made it. The sphere above is then the small root, which the formula as it is usually written gets
  // by taking two nearly equal numbers apart: at a lift of 1e-10 it missed the four by a millionth of an
  // Ångström, and at 1e-13 by a thousandth.
  for (const double lift :
       {1.0e-2, 1.0e-4, 1.0e-6, 1.0e-8, 1.0e-10, 1.0e-11, 1.0e-12, 1.0e-13, 0.0, -1.0e-10, -1.0e-8})
  {
    const std::array<double, 4> radii{1.0, 1.3, 0.8, 1.1};
    const std::array<double3, 4> centres{double3(0.0, 0.0, radii[0]), double3(3.1, 0.2, radii[1]),
                                         double3(0.4, 2.9, radii[2]), double3(3.3, 3.4, radii[3] + lift)};

    // The one above the four, about 21.7 in radius; the other, where there is one, is thousands of times that.
    const std::vector<ApolloniusSphere> spheres = apolloniusTangentSpheres(centres, radii);
    std::size_t nestled = 0;
    for (const ApolloniusSphere& sphere : spheres)
    {
      if (sphere.radius > 100.0) continue;
      ++nestled;
      EXPECT_NEAR(sphere.radius, 21.7437, 2.3) << "lift " << lift;
      EXPECT_LT(worstTangency(centres, radii, sphere), 1.0e-12) << "lift " << lift;
    }
    EXPECT_EQ(nestled, 1uz) << "lift " << lift;
  }
}


// Nodes that symmetry makes alike walk to peaks that symmetry makes alike. The walk climbs with the atoms
// nearest the point it stands at, and at a symmetric point those come in shells, several at exactly one
// clearance. Cut off at a fixed number, the last shell is taken in part, and which part is the rounding of
// distances that are equal: the four atoms the peak touches were then among the candidates or not by accident.
//
// An icosahedron of twelve atoms in a cell otherwise empty: the open space round it has its widest point at the
// corner of the cell, and six vertices on the faces of the cell that stand six atoms away at one clearance and
// six more at the next. All six have to reach the corner, however the cage is listed.
TEST(voronoi_network, nodes_that_symmetry_makes_alike_reach_the_same_peak)
{
  const double edge = 14.0;
  const double shellRadius = 3.0;
  const double atomRadius = 1.8;
  const UnitCell box(edge, edge, edge);

  std::vector<double3> cage;
  const double golden = 0.5 * (1.0 + std::sqrt(5.0));
  const double scale = shellRadius / (std::sqrt(1.0 + golden * golden) * edge);
  for (double first : {-1.0, 1.0})
  {
    for (double second : {-golden, golden})
    {
      cage.push_back(double3(0.5, 0.5 + first * scale, 0.5 + second * scale));
      cage.push_back(double3(0.5 + first * scale, 0.5 + second * scale, 0.5));
      cage.push_back(double3(0.5 + second * scale, 0.5, 0.5 + first * scale));
    }
  }
  const std::vector<double> radii(cage.size(), atomRadius);

  // The room inside the cage, and the room at the corner of the cell.
  const double inside = shellRadius - atomRadius;
  double atTheCorner = std::numeric_limits<double>::max();
  for (const double3& fractional : cage)
  {
    atTheCorner = std::min(atTheCorner, (box.cell * fractional).length() - atomRadius);
  }

  // As listed, back to front, and from the middle round.
  std::vector<std::vector<double3>> listings{cage, std::vector<double3>(cage.rbegin(), cage.rend()), cage};
  std::rotate(listings[2].begin(), listings[2].begin() + 5, listings[2].end());

  for (std::size_t listing = 0; listing < listings.size(); ++listing)
  {
    const VoronoiNetwork network = VoronoiNetwork::create(box, listings[listing], radii);

    std::size_t insideTheCage = 0;
    std::size_t reachedTheCorner = 0;
    for (const VoronoiNode& node : network.nodes)
    {
      if (std::abs(node.maximalRadius - inside) < 1.0e-9)
        ++insideTheCage;
      else if (std::abs(node.maximalRadius - atTheCorner) < 1.0e-9)
        ++reachedTheCorner;
      else
        ADD_FAILURE() << "listing " << listing << ": a node stopped at " << node.maximalRadius << " Å, from "
                      << node.radius << " Å";
    }
    EXPECT_EQ(insideTheCage, 1uz) << "listing " << listing;
    EXPECT_EQ(reachedTheCorner + insideTheCage, network.nodes.size()) << "listing " << listing;
    EXPECT_NEAR(network.largestIncludedSphereDiameter(), 2.0 * atTheCorner, 1.0e-9) << "listing " << listing;
  }
}
