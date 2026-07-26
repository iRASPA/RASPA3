#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <numbers>
#include <optional>
#include <random>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

import double3;
import double3x3;
import int3;
import skapolloniusdiagram;
import skvoronoi;

namespace
{
// Clearance min_j(|x - x_j| - r_j) over every site and its 27 nearest images. The Apollonius
// diagram is defined entirely by this field, so it is the reference every test measures against.
double clearanceAt(const double3x3& cell, const std::vector<double3>& fractionalPositions,
                   const std::vector<double>& radii, const double3& fractional)
{
  double clearance = std::numeric_limits<double>::max();
  for (std::size_t j = 0; j < fractionalPositions.size(); ++j)
  {
    double3 delta = fractionalPositions[j] - fractional;
    delta -= double3(std::round(delta.x), std::round(delta.y), std::round(delta.z));
    clearance = std::min(clearance, (cell * delta).length() - radii[j]);
  }
  return clearance;
}

// Clearance without the minimum-image shortcut. Rounding the fractional difference finds the nearest
// image only in a cell that is not far from orthogonal, so in a skewed cell the images are scanned
// outright. Slow, and meant for checking a skewed cell against something that cannot itself be wrong.
double clearanceByScan(const double3x3& cell, const std::vector<double3>& fractionalPositions,
                       const std::vector<double>& radii, const double3& fractional, int reach)
{
  double clearance = std::numeric_limits<double>::max();
  for (std::size_t j = 0; j < fractionalPositions.size(); ++j)
    for (int dx = -reach; dx <= reach; ++dx)
      for (int dy = -reach; dy <= reach; ++dy)
        for (int dz = -reach; dz <= reach; ++dz)
        {
          double3 delta = fractionalPositions[j] + double3(dx, dy, dz) - fractional;
          clearance = std::min(clearance, (cell * delta).length() - radii[j]);
        }
  return clearance;
}

struct RandomSites
{
  double3x3 cell;
  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
};

// Three sites on one line, the middle one then pushed across that line by `displacement`, so that zero gives
// exactly collinear centres and a small value gives nearly collinear ones. Their trisector is a circle about
// the line, or nearly one, and two further sites are laid against it to cut known arcs out, which leaves four
// vertices on it and two edges. The rest of the cell is filled in so that the diagram about the circle is a
// diagram and not a few sites in a void; no filler may come within reach of the circle, or it would cut it too.
//
// The circle is known in closed form when the centres are collinear. With the sites at -a, 0 and +a along the
// line and radii r, s and r, the outer two being alike puts the trisector in the plane through the middle
// site across the line, and equating the clearance to an outer site and to the middle one there fixes the
// radius as (a^2 - (r - s)^2) / 2(r - s).
struct CollinearSites
{
  double3x3 cell;
  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  double3 origin;  // the centre of the circle
  double circleRadius;
  double clearance;          // the clearance all along the circle
  std::vector<double> ends;  // azimuths of the four vertices, sorted
};

CollinearSites makeNearlyCollinearSites(double displacement)
{
  CollinearSites sites;
  double edge = 30.0;
  sites.cell = double3x3(double3(edge, 0.0, 0.0), double3(0.0, edge, 0.0), double3(0.0, 0.0, edge));
  sites.origin = double3(15.0, 15.0, 15.0);

  double a = 3.0;
  double outer = 2.0;
  double middle = 1.0;
  sites.circleRadius = (a * a - (outer - middle) * (outer - middle)) / (2.0 * (outer - middle));
  sites.clearance = sites.circleRadius - middle;

  std::vector<double3> positions{sites.origin + double3(0.0, 0.0, -a), sites.origin + double3(0.0, displacement, 0.0),
                                 sites.origin + double3(0.0, 0.0, a)};
  sites.radii = {outer, middle, outer};

  // A site off the axis at azimuth `azimuth`, `reach` from the axis and `height` above the circle's plane,
  // stands nearer than the circle's own clearance over an arc about its azimuth, so it cuts that arc out and
  // leaves the two ends of it as vertices. Solving |x(angle) - site| = clearance + radius gives them.
  auto addCuttingSite = [&](double azimuth, double reach, double height, double radius)
  {
    positions.push_back(sites.origin + double3(reach * std::cos(azimuth), reach * std::sin(azimuth), height));
    sites.radii.push_back(radius);
    double tangent = sites.clearance + radius;
    double meeting = (sites.circleRadius * sites.circleRadius + reach * reach + height * height - tangent * tangent) /
                     (2.0 * sites.circleRadius);
    double half = std::acos(meeting / reach);
    sites.ends.push_back(azimuth - half);
    sites.ends.push_back(azimuth + half);
  };
  addCuttingSite(0.0, 6.0, 2.0, 1.0);
  addCuttingSite(std::numbers::pi, 6.5, -1.5, 1.1);
  std::sort(sites.ends.begin(), sites.ends.end());

  auto distanceToCircle = [&](const double3& point)
  {
    double smallest = std::numeric_limits<double>::max();
    for (int dx = -1; dx <= 1; ++dx)
      for (int dy = -1; dy <= 1; ++dy)
        for (int dz = -1; dz <= 1; ++dz)
        {
          double3 offset = point + sites.cell * double3(dx, dy, dz) - sites.origin;
          double across = std::sqrt(offset.x * offset.x + offset.y * offset.y) - sites.circleRadius;
          smallest = std::min(smallest, std::sqrt(across * across + offset.z * offset.z));
        }
    return smallest;
  };

  std::mt19937 engine(20260726u);
  std::uniform_real_distribution<double> coordinate(0.0, edge);
  std::uniform_real_distribution<double> size(0.8, 1.6);
  while (sites.radii.size() < 24)
  {
    double3 candidate(coordinate(engine), coordinate(engine), coordinate(engine));
    double radius = size(engine);
    if (distanceToCircle(candidate) < sites.clearance + radius + 0.3) continue;
    bool clear = true;
    for (std::size_t j = 0; j < sites.radii.size(); ++j)
    {
      double3 separation = candidate - positions[j];
      for (std::size_t k = 0; k < 3; ++k) separation[k] -= edge * std::round(separation[k] / edge);
      if (separation.length() < radius + sites.radii[j] + 0.3) clear = false;
    }
    if (!clear) continue;
    positions.push_back(candidate);
    sites.radii.push_back(radius);
  }

  for (const double3& position : positions) sites.fractionalPositions.push_back(sites.cell.inverse() * position);
  return sites;
}

// The points of a trisector at one clearance, solved independently of the module: two of the tangency
// conditions are subtracted to give planes, and the remaining sphere is intersected with their line.
// Used to confirm from outside that a reported ring really is a curve of equal clearance.
std::vector<double3> trisectorAtClearance(const std::array<double3, 3>& centres, const std::array<double, 3>& radii,
                                          double clearance)
{
  double reference = radii[0] + clearance;
  if (reference < 0.0) return {};

  std::array<double3, 2> normals;
  std::array<double, 2> offsets;
  for (std::size_t i = 0; i < 2; ++i)
  {
    normals[i] = 2.0 * (centres[i + 1] - centres[0]);
    offsets[i] = (double3::dot(centres[i + 1], centres[i + 1]) - radii[i + 1] * radii[i + 1]) -
                 (double3::dot(centres[0], centres[0]) - radii[0] * radii[0]) -
                 2.0 * (radii[i + 1] - radii[0]) * clearance;
  }

  double3 direction = double3::cross(normals[0], normals[1]);
  if (direction.length() < 1.0e-12) return {};
  double3x3 matrix(double3(normals[0].x, normals[1].x, direction.x), double3(normals[0].y, normals[1].y, direction.y),
                   double3(normals[0].z, normals[1].z, direction.z));
  double3 base = matrix.inverse() * double3(offsets[0], offsets[1], 0.0);

  double3 delta = base - centres[0];
  double a = double3::dot(direction, direction);
  double b = 2.0 * double3::dot(delta, direction);
  double c = double3::dot(delta, delta) - reference * reference;
  double discriminant = b * b - 4.0 * a * c;
  if (discriminant < 0.0) return {};

  double root = std::sqrt(discriminant);
  return {base + ((-b + root) / (2.0 * a)) * direction, base + ((-b - root) / (2.0 * a)) * direction};
}

// The bisector of two sites is one sheet of a hyperboloid of revolution about the axis through their
// centres, so it is parametrised by the shared clearance and an angle about that axis: at each clearance
// the two tangent spheres meet in a circle. This returns a point of that circle, or nothing where the
// clearance is below the apex and the circle does not exist.
std::optional<double3> bisectorPoint(const double3& firstCentre, double firstRadius, const double3& secondCentre,
                                     double secondRadius, double clearance, double angle)
{
  double firstDistance = firstRadius + clearance;
  double secondDistance = secondRadius + clearance;
  if (firstDistance < 0.0 || secondDistance < 0.0) return std::nullopt;

  double3 axis = secondCentre - firstCentre;
  double separation = axis.length();
  if (separation < 1.0e-12) return std::nullopt;
  axis = (1.0 / separation) * axis;

  double along = (separation * separation + firstDistance * firstDistance - secondDistance * secondDistance) /
                 (2.0 * separation);
  double radialSquared = firstDistance * firstDistance - along * along;
  if (radialSquared < 0.0) return std::nullopt;
  double radial = std::sqrt(radialSquared);

  double3 reference = std::abs(axis.x) < 0.9 ? double3(1.0, 0.0, 0.0) : double3(0.0, 1.0, 0.0);
  double3 first = double3::cross(axis, reference);
  first = (1.0 / first.length()) * first;
  double3 second = double3::cross(axis, first);

  return firstCentre + along * axis + (radial * std::cos(angle)) * first + (radial * std::sin(angle)) * second;
}

// How many separate regions of a pair's bisector sheet belong to the diagram, counted on the sheet itself
// rather than from the edges the diagram found. The sheet is swept in clearance and angle; a sample
// belongs to the diagram when no third site is nearer than the clearance the pair shares there. The
// regions are the connected components of those samples, with the angle wrapping and the whole apex row
// treated as the single point it is.
std::size_t bisectorRegions(const double3x3& cell, const std::vector<double3>& fractionalPositions,
                            const std::vector<double>& radii, std::size_t first, std::size_t second,
                            const int3& secondImage, double clearanceSpan, std::size_t clearanceSteps,
                            std::size_t angleSteps, bool& reachedSpanEnd)
{
  double3 firstCentre = cell * fractionalPositions[first];
  double3 secondCentre =
      cell * (fractionalPositions[second] + double3(secondImage.x, secondImage.y, secondImage.z));

  // The apex, where the circle shrinks to a point, is the lowest clearance the sheet reaches.
  double low = -std::min(radii[first], radii[second]);
  double high = low + clearanceSpan;
  if (!bisectorPoint(firstCentre, radii[first], secondCentre, radii[second], high, 0.0).has_value()) return 0;
  for (std::size_t step = 0; step < 60; ++step)
  {
    double middle = 0.5 * (low + high);
    if (bisectorPoint(firstCentre, radii[first], secondCentre, radii[second], middle, 0.0).has_value())
      high = middle;
    else
      low = middle;
  }
  double apex = high;

  std::vector<char> inside(clearanceSteps * angleSteps, 0);
  for (std::size_t row = 0; row < clearanceSteps; ++row)
  {
    double clearance = apex + clearanceSpan * static_cast<double>(row) / static_cast<double>(clearanceSteps - 1);
    for (std::size_t column = 0; column < angleSteps; ++column)
    {
      double angle = 2.0 * std::numbers::pi * static_cast<double>(column) / static_cast<double>(angleSteps);
      std::optional<double3> point =
          bisectorPoint(firstCentre, radii[first], secondCentre, radii[second], clearance, angle);
      if (!point.has_value()) continue;
      double3 fractional = cell.inverse() * point.value();
      if (clearanceAt(cell, fractionalPositions, radii, fractional) > clearance - 1.0e-9)
        inside[row * angleSteps + column] = 1;
    }
  }

  // If the sweep still finds the diagram on its last row it has not covered the whole sheet, and any
  // count taken from it would be meaningless.
  reachedSpanEnd = false;
  for (std::size_t column = 0; column < angleSteps; ++column)
    if (inside[(clearanceSteps - 1) * angleSteps + column]) reachedSpanEnd = true;

  std::vector<std::size_t> parent(inside.size());
  for (std::size_t node = 0; node < parent.size(); ++node) parent[node] = node;
  std::function<std::size_t(std::size_t)> find = [&](std::size_t node)
  {
    while (parent[node] != node) node = parent[node] = parent[parent[node]];
    return node;
  };
  auto join = [&](std::size_t left, std::size_t right)
  {
    if (!inside[left] || !inside[right]) return;
    std::size_t a = find(left);
    std::size_t b = find(right);
    if (a != b) parent[a] = b;
  };

  for (std::size_t row = 0; row < clearanceSteps; ++row)
    for (std::size_t column = 0; column < angleSteps; ++column)
    {
      std::size_t node = row * angleSteps + column;
      join(node, row * angleSteps + (column + 1) % angleSteps);
      if (row + 1 < clearanceSteps) join(node, (row + 1) * angleSteps + column);
      if (row == 0) join(node, 0);  // the apex row is one point
    }

  // A sample or two clinging to the rim of a region is the grid failing to resolve its boundary, not a
  // region of its own: such specks come and go as the grid is refined, while a real region grows with it.
  // Anything this small is therefore not counted.
  constexpr std::size_t smallestRegion = 4;
  std::map<std::size_t, std::size_t> sizes;
  for (std::size_t node = 0; node < inside.size(); ++node)
    if (inside[node]) ++sizes[find(node)];
  std::size_t regions = 0;
  for (const auto& [root, size] : sizes)
    if (size >= smallestRegion) ++regions;
  return regions;
}

RandomSites makeRandomSites(double edge, std::size_t count, unsigned seed, double minimumRadius, double maximumRadius)
{
  RandomSites sites;
  sites.cell = double3x3(double3(edge, 0.0, 0.0), double3(0.0, edge, 0.0), double3(0.0, 0.0, edge));
  std::mt19937 generator(seed);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  std::uniform_real_distribution<double> radius(minimumRadius, maximumRadius);
  for (std::size_t i = 0; i < count; ++i)
  {
    sites.fractionalPositions.push_back(double3(uniform(generator), uniform(generator), uniform(generator)));
    sites.radii.push_back(radius(generator));
  }
  return sites;
}
}  // namespace

// The tangent-sphere primitive, against an analytic reference.
TEST(apollonius, tangent_sphere_closed_form)
{
  // Regular tetrahedron of equal spheres: the tangent sphere sits at the circumcentre, here the
  // origin, with radius the circumradius less the common atom radius.
  std::array<double3, 4> centres{double3(1.0, 1.0, 1.0), double3(1.0, -1.0, -1.0), double3(-1.0, 1.0, -1.0),
                                 double3(-1.0, -1.0, 1.0)};
  std::array<double, 4> equal{0.5, 0.5, 0.5, 0.5};

  std::vector<SKApolloniusTangentSphere> spheres = skApolloniusTangentSpheres(centres, equal);
  ASSERT_EQ(spheres.size(), 1u);
  EXPECT_NEAR(spheres[0].radius, std::sqrt(3.0) - 0.5, 1.0e-12);
  EXPECT_NEAR(spheres[0].centre.length(), 0.0, 1.0e-12);

  // Unequal radii have no simple reference, so assert the defining property: the solution must
  // touch all four spheres from outside.
  std::array<double, 4> unequal{0.9, 0.4, 0.6, 0.25};
  std::vector<SKApolloniusTangentSphere> weighted = skApolloniusTangentSpheres(centres, unequal);
  EXPECT_FALSE(weighted.empty());
  for (const SKApolloniusTangentSphere& sphere : weighted)
    for (std::size_t i = 0; i < 4; ++i)
      EXPECT_NEAR((sphere.centre - centres[i]).length(), unequal[i] + sphere.radius, 1.0e-9);

  // Coplanar centres leave the tangent sphere undetermined along the plane normal.
  std::array<double3, 4> coplanar{double3(0.0, 0.0, 0.0), double3(1.0, 0.0, 0.0), double3(0.0, 1.0, 0.0),
                                  double3(1.0, 1.0, 0.0)};
  EXPECT_TRUE(skApolloniusTangentSpheres(coplanar, equal).empty());
}

// Points sampled along an edge must lie on the trisector curve of its three sites, meaning the
// sphere there touches all three. This is what makes edge geometry exact rather than interpolated.
TEST(apollonius, trisector_points_are_tangent_to_three_sites)
{
  std::array<double3, 3> centres{double3(0.0, 0.0, 0.0), double3(4.0, 0.0, 0.0), double3(1.5, 3.5, 0.0)};
  std::array<double, 3> radii{1.2, 0.7, 0.9};
  double3 from(1.8, 1.0, -2.5);
  double3 to(1.8, 1.0, 2.5);

  std::size_t sampled = 0;
  for (std::size_t s = 0; s <= 10; ++s)
  {
    double parameter = static_cast<double>(s) / 10.0;
    std::optional<SKApolloniusTangentSphere> point =
        skApolloniusTrisectorPoint(centres, radii, from, to, parameter);
    if (!point.has_value()) continue;
    ++sampled;
    for (std::size_t i = 0; i < 3; ++i)
      EXPECT_NEAR((point->centre - centres[i]).length(), radii[i] + point->radius, 1.0e-9);
  }
  EXPECT_GT(sampled, 8u);
}

// Every vertex must be exactly what it claims: tangent to its four sites, and empty of all others.
// Emptiness is the certification the construction promises, so it is checked here independently.
TEST(apollonius, every_vertex_is_tangent_and_empty)
{
  RandomSites sites = makeRandomSites(14.0, 40, 5, 0.4, 2.2);
  SKApolloniusDiagram diagram =
      SKApolloniusDiagram::create(sites.cell, sites.fractionalPositions, sites.radii);
  ASSERT_FALSE(diagram.vertices.empty());

  for (const SKApolloniusVertex& vertex : diagram.vertices)
  {
    // Four sites, or more where a sphere happens to touch more, and the tangency must hold for all of them.
    EXPECT_GE(vertex.siteIndices.size(), 4u);
    for (std::size_t s = 0; s < vertex.siteIndices.size(); ++s)
    {
      double3 image = sites.cell * (double3::fract(sites.fractionalPositions[vertex.siteIndices[s]]) +
                                    double3(vertex.siteImages[s].x, vertex.siteImages[s].y, vertex.siteImages[s].z));
      EXPECT_NEAR((vertex.position - image).length(), sites.radii[vertex.siteIndices[s]] + vertex.radius, 1.0e-7);
    }

    // The clearance at the vertex equals its radius exactly when nothing intrudes.
    double3 fractional = sites.cell.inverse() * vertex.position;
    EXPECT_NEAR(clearanceAt(sites.cell, sites.fractionalPositions, sites.radii, fractional), vertex.radius, 1.0e-7);
  }
}

// With equal radii the clearance and power orderings agree, so the Apollonius diagram is the
// ordinary Voronoi diagram. Its deepest vertex must then match the radical construction exactly,
// which ties the new code to the existing one at the point where they must agree.
TEST(apollonius, equal_radii_reduce_to_ordinary_voronoi)
{
  RandomSites sites = makeRandomSites(12.0, 30, 9, 1.0, 1.0);
  SKApolloniusDiagram diagram =
      SKApolloniusDiagram::create(sites.cell, sites.fractionalPositions, sites.radii);
  ASSERT_FALSE(diagram.vertices.empty());

  SKVoronoi radical(sites.cell, sites.fractionalPositions, sites.radii);
  std::vector<SKVoronoiCell> cells = radical.computeAllCells();

  double radicalDeepest = 0.0;
  for (const SKVoronoiCell& cell : cells)
    for (const double3& vertex : cell.verticesCartesian)
      radicalDeepest = std::max(radicalDeepest, vertex.length() - sites.radii[cell.siteIndex]);

  EXPECT_NEAR(diagram.largestEmptySphereRadius(), radicalDeepest, 1.0e-6);
}

// The complex closes as far as the free region permits.
//
// Three causes of earlier failures were found by following the literature, and are worth recording
// because each was a wrong model of the structure rather than a coding slip.
//
// First, a vertex is not identified by its four sites. A quadruple admits two tangent spheres, the
// two roots of the quadratic, and both can be empty, giving two distinct vertices: Kamarianakis
// ("Predicates of the 3D Apollonius Diagram", arXiv:2007.06658) writes them v_ijkl and v_ikjl and
// separates them by the orientation of the tetrahedron of tangency points, and Wang et al. ("Robust
// Computation of 3D Apollonius Diagrams", CGF 39(7) 2020) call the same thing a Type II quadruple in
// their Theorem 4. Keying on the sites alone collided the two and dropped one.
//
// Second, a triple does not carry exactly two vertices, so pairing them off arbitrarily is wrong.
// Kamarianakis analyses trisectors bounding three vertices (Section 2.7.6, Case C) and orders them
// with an Order predicate. Here the vertices are ordered along the trisector, which Wang et al.
// Theorem 5 shows is a planar curve, either an open branch or a closed loop, and consecutive ones are
// paired; which arcs are real is then settled by the direction each edge leaves its vertex, from the
// tangent of their Theorem 5 and the marching test of their Equations 26 to 27.
//
// Third, and the reason this test asserts a weaker identity than Euler's formula: only tangency
// solutions of non-negative radius are kept, since a probe cannot occupy a sphere of negative radius.
// That restricts the diagram to the free region and clips it against the boundary of the union of the
// spheres. On this input exactly one quadruple, whose two roots are -15.2 and -28.8, has its vertex
// inside the spheres; its four triples are therefore each left with one end, and the four vertices at
// the other ends are legitimately three-valent. Euler's formula counts that vertex and those arcs and
// so does not apply. `truncatedTriples` records the clipping.
TEST(apollonius, complex_closes_on_well_behaved_sites)
{
  RandomSites sites = makeRandomSites(14.0, 40, 5, 0.4, 2.2);
  SKApolloniusDiagram diagram =
      SKApolloniusDiagram::create(sites.cell, sites.fractionalPositions, sites.radii);

  EXPECT_EQ(diagram.verification.unpairedTriples, 0u);
  EXPECT_EQ(diagram.verification.overpairedTriples, 0u);
  EXPECT_TRUE(diagram.verification.isComplete());

  // Every vertex is four-valent except where a clipped arc leaves it one short, and each clipped arc
  // accounts for exactly one such shortfall.
  EXPECT_EQ(diagram.verification.verticesOfFullValence + diagram.verification.truncatedTriples,
            diagram.verification.vertexCount);

  // An open face can only be one bordering a clipped arc, and an arc borders the three bisector
  // patches of the three pairs among its sites.
  EXPECT_LE(diagram.verification.unclosedFaces, 3 * diagram.verification.truncatedTriples);
}

// Over all of space nothing is clipped, so the identities that the free-region diagram can only
// satisfy up to its truncation must hold outright: no triple is left with one end, every vertex is
// four-valent, and Euler's formula closes on the 3-torus.
TEST(apollonius, entire_space_diagram_closes_exactly)
{
  RandomSites sites = makeRandomSites(14.0, 40, 5, 0.4, 2.2);
  SKApolloniusDiagram diagram = SKApolloniusDiagram::create(sites.cell, sites.fractionalPositions, sites.radii, 1,
                                                           SKApolloniusRegion::EntireSpace);

  EXPECT_EQ(diagram.verification.truncatedTriples, 0u);
  EXPECT_EQ(diagram.verification.unpairedTriples, 0u);
  EXPECT_EQ(diagram.verification.overpairedTriples, 0u);
  EXPECT_EQ(diagram.verification.unclosedFaces, 0u);
  EXPECT_EQ(diagram.verification.verticesOfFullValence, diagram.verification.vertexCount);
  EXPECT_TRUE(diagram.verification.isComplete());

  // V - E + F - C = 0 for a cell complex on the 3-torus. Faces are stored once per owning site, so
  // each geometric patch is counted twice. These sites need no ring, so every edge here runs between two
  // vertices; the count holds with a ring too, as the test on those sites shows.
  EXPECT_EQ(diagram.verification.vertexlessLoops, 0u);
  std::size_t v = diagram.vertices.size();
  std::size_t e = diagram.edges.size();
  std::size_t f = diagram.faces.size() / 2;
  std::size_t c = 0;
  for (const SKApolloniusCell& cell : diagram.cells)
    if (!cell.isEmpty) ++c;
  // Four-valence everywhere forces exactly two edge-ends per vertex per pair, so this is exact.
  EXPECT_EQ(e, 2 * v);

  // A bisector patch carrying no vertex at all cannot appear in a structure assembled from vertices.
  // Wang et al. call these face humps and closed vertexless edges (Sections 4.2 and 5.5) and report
  // that their own sweep misses them as well. Their number is measured here independently of the
  // diagram, by sampling: the two sites of least clearance at a point are Apollonius neighbours there
  // and so must share a patch. Euler's formula then has to close once they are counted back in.
  // A patch is named by its two sites and the lattice shift between them, so that a patch between a
  // site and its own image counts too, and so that the naming does not depend on which end is read
  // first.
  using PatchKey = std::tuple<std::size_t, std::size_t, int, int, int>;
  auto canonicalPatch = [](std::size_t i, std::size_t j, int dx, int dy, int dz)
  {
    if (i < j) return PatchKey{i, j, dx, dy, dz};
    if (j < i) return PatchKey{j, i, -dx, -dy, -dz};
    if (std::make_tuple(dx, dy, dz) < std::make_tuple(-dx, -dy, -dz)) return PatchKey{i, i, -dx, -dy, -dz};
    return PatchKey{i, i, dx, dy, dz};
  };

  std::set<PatchKey> sampledPatches;
  std::mt19937 engine(7);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  for (std::size_t sample = 0; sample < 20000; ++sample)
  {
    double3 point = sites.cell * double3(uniform(engine), uniform(engine), uniform(engine));
    double best = std::numeric_limits<double>::max();
    double second = std::numeric_limits<double>::max();
    std::size_t bestSite = 0;
    std::size_t secondSite = 0;
    int3 bestImage(0, 0, 0);
    int3 secondImage(0, 0, 0);
    for (std::size_t i = 0; i < sites.radii.size(); ++i)
      for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
          for (int dz = -1; dz <= 1; ++dz)
          {
            double3 centre = sites.cell * (sites.fractionalPositions[i] + double3(dx, dy, dz));
            double clearance = (point - centre).length() - sites.radii[i];
            if (clearance < best)
            {
              second = best;
              secondSite = bestSite;
              secondImage = bestImage;
              best = clearance;
              bestSite = i;
              bestImage = int3(dx, dy, dz);
            }
            else if (clearance < second)
            {
              second = clearance;
              secondSite = i;
              secondImage = int3(dx, dy, dz);
            }
          }
    sampledPatches.insert(canonicalPatch(bestSite, secondSite, secondImage.x - bestImage.x,
                                         secondImage.y - bestImage.y, secondImage.z - bestImage.z));
  }

  std::set<PatchKey> diagramPatches;
  for (const SKApolloniusFace& face : diagram.faces)
    diagramPatches.insert(
        canonicalPatch(face.site1, face.site2, face.site2Image.x, face.site2Image.y, face.site2Image.z));

  std::size_t vertexlessPatches = 0;
  for (const PatchKey& patch : sampledPatches)
    if (!diagramPatches.contains(patch)) ++vertexlessPatches;

  std::int64_t characteristic = static_cast<std::int64_t>(v) - static_cast<std::int64_t>(e) +
                                static_cast<std::int64_t>(f) - static_cast<std::int64_t>(c);
  EXPECT_EQ(vertexlessPatches, 0u);
  EXPECT_EQ(characteristic, 0);
}

// A trisector on which no fourth site ever intrudes is a closed curve carrying no vertex at all, so
// nothing that starts from vertices can reach it. These sites hold one: three overlapping spheres whose
// shared clearance runs from -1.7083 to -1.0515, a ring buried inside the spheres and therefore part of
// the diagram over the entire space but not of the free-region one.
TEST(apollonius, closed_trisector_without_vertices_is_found)
{
  double edge = 12.0;
  double3x3 cell(double3(edge, 0.0, 0.0), double3(0.0, edge, 0.0), double3(0.0, 0.0, edge));
  std::vector<double3> positions{double3(9.971800, 4.331360, 8.432872),  double3(10.321425, 7.695810, 6.580359),
                                 double3(9.147761, 8.595773, 5.606025),  double3(6.869507, 8.955868, 0.762767),
                                 double3(7.762513, 8.832686, 4.785106),  double3(6.086190, 2.744648, 7.802330)};
  std::vector<double> radii{5.842063, 2.142922, 3.045569, 5.403876, 3.532705, 2.818056};

  std::vector<double3> fractionalPositions;
  for (const double3& position : positions) fractionalPositions.push_back((1.0 / edge) * position);

  SKApolloniusDiagram diagram =
      SKApolloniusDiagram::create(cell, fractionalPositions, radii, 1, SKApolloniusRegion::EntireSpace);

  EXPECT_EQ(diagram.verification.vertexlessLoops, 1u);
  EXPECT_EQ(diagram.verification.ringsOfUncertainFace, 0u);

  std::vector<std::size_t> loops;
  for (std::size_t e = 0; e < diagram.edges.size(); ++e)
    if (diagram.edges[e].isLoop) loops.push_back(e);
  ASSERT_EQ(loops.size(), 1u);

  const SKApolloniusEdge& ring = diagram.edges[loops.front()];
  EXPECT_EQ(ring.from, std::numeric_limits<std::size_t>::max());
  EXPECT_EQ(ring.to, std::numeric_limits<std::size_t>::max());
  EXPECT_GT(ring.length, 0.0);
  EXPECT_NEAR(ring.bottleneckRadius, -1.708260, 1.0e-6);

  // The ring is buried, so its clearance is negative throughout and no probe travels it. That is why it
  // belongs to the diagram over the entire space alone.
  EXPECT_LT(ring.bottleneckRadius, 0.0);

  std::array<double3, 3> centres;
  std::array<double, 3> ringRadii;
  for (std::size_t s = 0; s < 3; ++s)
  {
    centres[s] = positions[ring.siteIndices[s]] + cell * double3(ring.siteImages[s].x, ring.siteImages[s].y,
                                                                ring.siteImages[s].z);
    ringRadii[s] = radii[ring.siteIndices[s]];
  }

  // Confirm from outside that the curve is what the diagram claims: at its narrowest the two branches
  // have met, above that it has two points, and at every one of them the three sites share the clearance
  // and no other site undercuts it. That last part is what makes the ring vertexless.
  EXPECT_TRUE(trisectorAtClearance(centres, ringRadii, ring.bottleneckRadius - 1.0e-4).empty());

  // The ring is a hole in the region it borders, not a region of its own, and once it is recorded that way
  // the counts close: V - E + F - C = 0 as for any other cell complex here. Recording it as a face of its
  // own instead leaves the count one over, which is the check that caught that error.
  std::size_t v = diagram.vertices.size();
  std::size_t e = diagram.edges.size();
  std::size_t f = diagram.faces.size() / 2;
  std::size_t c = 0;
  for (const SKApolloniusCell& diagramCell : diagram.cells)
    if (!diagramCell.isEmpty) ++c;
  EXPECT_EQ(static_cast<std::int64_t>(v) - static_cast<std::int64_t>(e) + static_cast<std::int64_t>(f) -
                static_cast<std::int64_t>(c),
            0);
  EXPECT_EQ(e, 2 * v + 1);  // every vertex four-valent, plus the ring, which has no vertex

  std::size_t sampled = 0;
  for (std::size_t sample = 1; sample <= 24; ++sample)
  {
    double clearance = ring.bottleneckRadius + 0.6 * static_cast<double>(sample) / 24.0;
    for (const double3& point : trisectorAtClearance(centres, ringRadii, clearance))
    {
      for (std::size_t s = 0; s < 3; ++s) EXPECT_NEAR((point - centres[s]).length() - ringRadii[s], clearance, 1.0e-8);
      double3 fractional = cell.inverse() * point;
      EXPECT_GT(clearanceAt(cell, fractionalPositions, radii, fractional), clearance - 1.0e-7);
      ++sampled;
    }
  }
  EXPECT_GT(sampled, 20u);
}

// Faces are assembled from edges, so the question is whether a region of a bisector sheet can exist that
// no edge borders and that would therefore be missed. It cannot: a sheet is one branch of a hyperboloid
// of revolution and so is topologically a plane, other sites dominate it far from the pair, and the
// boundary of a bounded region on it consists of points equally near a third site, which is an edge. What
// can happen is that the only bounding edge is a closed trisector carrying no vertex, which is the ring
// handled above. This checks that reasoning against the sheets themselves: every sheet is swept in
// clearance and angle, its regions counted on the surface, and the count compared with the faces the
// diagram built from edges. The sites are the ones carrying the ring, the worst case available.
TEST(apollonius, every_bisector_region_is_a_face)
{
  double edge = 12.0;
  double3x3 cell(double3(edge, 0.0, 0.0), double3(0.0, edge, 0.0), double3(0.0, 0.0, edge));
  std::vector<double3> positions{double3(9.971800, 4.331360, 8.432872),  double3(10.321425, 7.695810, 6.580359),
                                 double3(9.147761, 8.595773, 5.606025),  double3(6.869507, 8.955868, 0.762767),
                                 double3(7.762513, 8.832686, 4.785106),  double3(6.086190, 2.744648, 7.802330)};
  std::vector<double> radii{5.842063, 2.142922, 3.045569, 5.403876, 3.532705, 2.818056};

  std::vector<double3> fractionalPositions;
  for (const double3& position : positions) fractionalPositions.push_back((1.0 / edge) * position);

  SKApolloniusDiagram diagram =
      SKApolloniusDiagram::create(cell, fractionalPositions, radii, 1, SKApolloniusRegion::EntireSpace);

  // Faces are stored once per owning site, so a geometric region appears twice.
  using PairKey = std::tuple<std::size_t, std::size_t, int, int, int>;
  std::map<PairKey, std::size_t> faceCount;
  for (const SKApolloniusFace& face : diagram.faces)
    ++faceCount[PairKey{face.site1, face.site2, face.site2Image.x, face.site2Image.y, face.site2Image.z}];

  // Counting regions on a sampled surface is delicate: too coarse a grid can split one region across a
  // narrow neck, and can also bridge two across a narrow gap. So a disagreement is not taken at face value
  // but re-examined on a finer grid, and only a disagreement that survives refinement is a failure.
  std::size_t comparedPairs = 0;
  std::size_t regionsFound = 0;
  for (std::size_t first = 0; first < radii.size(); ++first)
    for (std::size_t second = 0; second < radii.size(); ++second)
      for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
          for (int dz = -1; dz <= 1; ++dz)
          {
            if (first == second && dx == 0 && dy == 0 && dz == 0) continue;

            std::size_t faces = faceCount[PairKey{first, second, dx, dy, dz}];
            std::size_t regions = 0;
            for (std::size_t resolution : {240u, 480u, 960u})
            {
              bool reachedSpanEnd = false;
              regions = bisectorRegions(cell, fractionalPositions, radii, first, second, int3(dx, dy, dz), 10.0,
                                        resolution, resolution, reachedSpanEnd);
              // A sweep that still sees the diagram at its far end has not covered the sheet.
              ASSERT_FALSE(reachedSpanEnd) << "pair " << first << "/" << second;
              if (regions == faces) break;
            }

            EXPECT_EQ(regions, faces) << "pair " << first << "/" << second << " image " << dx << "," << dy << "," << dz;
            regionsFound += regions;
            if (regions > 0 || faces > 0) ++comparedPairs;
          }

  EXPECT_GT(comparedPairs, 20u);
  EXPECT_GT(regionsFound, 20u);
}

// The diagram belongs to the crystal, not to the box drawn around it. Replacing the lattice basis by
// another basis of the same lattice, here a shear that leaves no cell angle at ninety degrees, describes
// the identical set of spheres in a triclinic cell, so every quantity that does not depend on which cell
// the points were wrapped into has to come out the same. Nothing else in these tests is triclinic, so this
// is what covers skewed cells: a formula that silently assumed orthogonality, such as reading the lattice
// vectors out of the cell matrix transposed, changes the answer here and nowhere else.
TEST(apollonius, diagram_is_independent_of_the_lattice_basis)
{
  RandomSites sites = makeRandomSites(14.0, 40, 5, 0.4, 2.2);

  // Columns a, a + b, a + b + c: an integer basis change of determinant one, so the same lattice seen in a
  // triclinic cell. Stronger shears are possible but describe a cell no one would use, and cost far more
  // to build without covering anything this does not.
  double3x3 shear(double3(1.0, 0.0, 0.0), double3(1.0, 1.0, 0.0), double3(1.0, 1.0, 1.0));
  double3x3 shearedCell = sites.cell * shear;
  double3x3 inverseShear = shear.inverse();
  std::vector<double3> shearedPositions;
  for (const double3& position : sites.fractionalPositions) shearedPositions.push_back(inverseShear * position);

  // Confirm the new cell really is triclinic before drawing conclusions from it.
  double3 a = shearedCell[0];
  double3 b = shearedCell[1];
  double3 c = shearedCell[2];
  EXPECT_GT(std::abs(double3::dot(a, b)), 1.0e-6);
  EXPECT_GT(std::abs(double3::dot(b, c)), 1.0e-6);
  EXPECT_GT(std::abs(double3::dot(a, c)), 1.0e-6);
  EXPECT_NEAR(std::abs(shearedCell.determinant()), std::abs(sites.cell.determinant()), 1.0e-9);

  for (SKApolloniusRegion region : {SKApolloniusRegion::FreeSpace, SKApolloniusRegion::EntireSpace})
  {
    SKApolloniusDiagram upright = SKApolloniusDiagram::create(sites.cell, sites.fractionalPositions, sites.radii, 1,
                                                             region);
    SKApolloniusDiagram skewed = SKApolloniusDiagram::create(shearedCell, shearedPositions, sites.radii, 1, region);

    ASSERT_EQ(skewed.vertices.size(), upright.vertices.size());
    EXPECT_EQ(skewed.edges.size(), upright.edges.size());
    EXPECT_EQ(skewed.faces.size(), upright.faces.size());
    EXPECT_EQ(skewed.verification.verticesOfFullValence, upright.verification.verticesOfFullValence);
    EXPECT_EQ(skewed.verification.unpairedTriples, upright.verification.unpairedTriples);
    EXPECT_EQ(skewed.verification.overpairedTriples, upright.verification.overpairedTriples);
    EXPECT_EQ(skewed.verification.vertexlessLoops, upright.verification.vertexlessLoops);
    EXPECT_NEAR(skewed.largestEmptySphereRadius(), upright.largestEmptySphereRadius(), 1.0e-9);

    // Wrapping puts the vertices in different places in the two cells, so the comparison is of the
    // clearances themselves, which the choice of cell cannot touch.
    std::vector<double> uprightRadii;
    std::vector<double> skewedRadii;
    for (const SKApolloniusVertex& vertex : upright.vertices) uprightRadii.push_back(vertex.radius);
    for (const SKApolloniusVertex& vertex : skewed.vertices) skewedRadii.push_back(vertex.radius);
    std::sort(uprightRadii.begin(), uprightRadii.end());
    std::sort(skewedRadii.begin(), skewedRadii.end());
    for (std::size_t i = 0; i < uprightRadii.size(); ++i) EXPECT_NEAR(skewedRadii[i], uprightRadii[i], 1.0e-9);

    // The neighbour search bins sites by fractional coordinate and sizes its reach from the perpendicular
    // widths of the cell, which are the one place a skewed cell is easy to get wrong. Reading those widths
    // off the wrong axis makes the search too short and lets a sphere that is not actually empty through,
    // so every vertex is checked against a direct scan over images in the skewed cell.
    for (const SKApolloniusVertex& vertex : skewed.vertices)
    {
      double3 fractional = shearedCell.inverse() * vertex.position;
      EXPECT_GT(clearanceByScan(shearedCell, shearedPositions, sites.radii, fractional, 3), vertex.radius - 1.0e-6);
    }

    std::vector<double> uprightBottlenecks;
    std::vector<double> skewedBottlenecks;
    for (const SKApolloniusEdge& edge : upright.edges) uprightBottlenecks.push_back(edge.bottleneckRadius);
    for (const SKApolloniusEdge& edge : skewed.edges) skewedBottlenecks.push_back(edge.bottleneckRadius);
    std::sort(uprightBottlenecks.begin(), uprightBottlenecks.end());
    std::sort(skewedBottlenecks.begin(), skewedBottlenecks.end());
    ASSERT_EQ(skewedBottlenecks.size(), uprightBottlenecks.size());
    for (std::size_t i = 0; i < uprightBottlenecks.size(); ++i)
      EXPECT_NEAR(skewedBottlenecks[i], uprightBottlenecks[i], 1.0e-6);
  }
}

// A vertex belongs to four sites, but nothing stops a sphere from touching five, and then all five meet
// at that point in a vertex of degree higher than four.
//
// It has to come out as one vertex. Each of the five quadruples satisfies the tangency, so the vertex is
// constructed five times over, and left that way it would be five vertices at one point, each pairing
// along triples the others also claim, which wrecks the complex around it. So the test is that the
// diagram holds one vertex there carrying all five sites, that it has the six edges the geometry gives
// it, and that the complex still closes exactly as it does on a generic input.
//
// The degeneracy is implanted rather than searched for. A configuration known to close is taken, and one
// more site is laid exactly against the tangent sphere of one of its vertices, in a direction where it
// overlaps nothing. That sphere is still empty and still touches all four sites it did, and now touches
// five, so the vertex becomes degenerate while everything else stays as well conditioned as it was and
// the closure of the whole diagram remains a fair thing to demand.
TEST(apollonius, degenerate_vertex_is_one_vertex_of_higher_degree)
{
  RandomSites sites = makeRandomSites(14.0, 40, 5, 0.4, 2.2);
  SKApolloniusDiagram base = SKApolloniusDiagram::create(sites.cell, sites.fractionalPositions, sites.radii, 1,
                                                        SKApolloniusRegion::EntireSpace);
  ASSERT_TRUE(base.verification.isComplete());
  ASSERT_EQ(base.verification.degenerateVertices, 0u);

  double newRadius = 1.0;
  std::size_t chosen = base.vertices.size();
  double3 newCentre(0.0, 0.0, 0.0);
  for (std::size_t v = 0; v < base.vertices.size() && chosen == base.vertices.size(); ++v)
  {
    // A middling vertex, so that the new site lands neither on top of the four already touching the
    // sphere nor far across the cell.
    if (base.vertices[v].radius < 0.5 || base.vertices[v].radius > 2.5) continue;

    std::mt19937 generator(11 + v);
    std::normal_distribution<double> gaussian(0.0, 1.0);
    for (std::size_t attempt = 0; attempt < 200; ++attempt)
    {
      double3 direction(gaussian(generator), gaussian(generator), gaussian(generator));
      double length = direction.length();
      if (length < 1.0e-6) continue;
      double3 candidate = base.vertices[v].position + ((base.vertices[v].radius + newRadius) / length) * direction;
      if (clearanceAt(sites.cell, sites.fractionalPositions, sites.radii, sites.cell.inverse() * candidate) <
          newRadius + 0.15)
        continue;  // the new site would overlap, or come close enough to overlapping, one already there
      chosen = v;
      newCentre = candidate;
      break;
    }
  }
  ASSERT_LT(chosen, base.vertices.size());
  double3 degeneratePosition = base.vertices[chosen].position;
  double degenerateRadius = base.vertices[chosen].radius;

  std::vector<double3> fractionalPositions = sites.fractionalPositions;
  std::vector<double> radii = sites.radii;
  fractionalPositions.push_back(sites.cell.inverse() * newCentre);
  radii.push_back(newRadius);
  std::vector<double3> positions;
  for (const double3& fractional : fractionalPositions) positions.push_back(sites.cell * fractional);

  SKApolloniusDiagram diagram =
      SKApolloniusDiagram::create(sites.cell, fractionalPositions, radii, 1, SKApolloniusRegion::EntireSpace);

  // The sweep meets the degeneracy and one vertex comes of it rather than a cluster of copies. It takes no
  // special effort over it: five exactly cotangent sites cannot be separated however far a box around them
  // is refined, and the sweep does not try to, finishing a box as soon as few sites can be nearest in it and
  // solving those. So it never refines to the floor, which is what the count of degenerate boxes records.
  EXPECT_EQ(diagram.verification.degenerateSweepBoxes, 0u);
  EXPECT_EQ(diagram.verification.degenerateVertices, 1u);
  EXPECT_EQ(diagram.verification.coincidentVertices, 0u);

  std::size_t atVertex = 0;
  for (const SKApolloniusVertex& vertex : diagram.vertices)
  {
    if ((vertex.position - degeneratePosition).length() > 1.0e-6) continue;
    ++atVertex;

    // The sphere has not moved, and all five sites are on it.
    EXPECT_NEAR(vertex.radius, degenerateRadius, 1.0e-9);
    EXPECT_EQ(vertex.siteIndices.size(), 5u);
    for (std::size_t s = 0; s < vertex.siteIndices.size(); ++s)
    {
      const int3& image = vertex.siteImages[s];
      double3 centre = positions[vertex.siteIndices[s]] + sites.cell * double3(image.x, image.y, image.z);
      EXPECT_NEAR((vertex.position - centre).length() - radii[vertex.siteIndices[s]], vertex.radius, 1.0e-7);
    }

    // The tangency points of n sites in general position on the tangent sphere have a simplicial convex
    // hull, which by Euler has 2n - 4 triangles, and each is an edge leaving the vertex. So the degree is
    // 6 here, where it is 4 for the four sites of the general case.
    EXPECT_EQ(vertex.branches.size(), 2 * vertex.siteIndices.size() - 4);

    // Each branch heads off along the trisector of its three sites, which requires those three to stay
    // equally near while the others fall away.
    for (const SKApolloniusBranch& branch : vertex.branches)
    {
      double3 point = vertex.position + 1.0e-4 * branch.direction;
      std::vector<double> onBranch;
      std::vector<double> offBranch;
      for (std::size_t s = 0; s < vertex.siteIndices.size(); ++s)
      {
        const int3& image = vertex.siteImages[s];
        double3 centre = positions[vertex.siteIndices[s]] + sites.cell * double3(image.x, image.y, image.z);
        double clearance = (point - centre).length() - radii[vertex.siteIndices[s]];
        bool inTriple = std::find(branch.sites.begin(), branch.sites.end(), s) != branch.sites.end();
        (inTriple ? onBranch : offBranch).push_back(clearance);
      }
      ASSERT_EQ(onBranch.size(), 3u);
      for (double clearance : onBranch) EXPECT_NEAR(clearance, onBranch.front(), 1.0e-8);
      for (double clearance : offBranch) EXPECT_GT(clearance, onBranch.front() + 1.0e-9);
    }
  }
  EXPECT_EQ(atVertex, 1u);

  // Held as one vertex the complex closes just as it does without any degeneracy: every branch carries an
  // edge, every triple pairs, and Euler's formula gives zero on the 3-torus.
  EXPECT_EQ(diagram.verification.ambiguousBranches, 0u);
  EXPECT_EQ(diagram.verification.unpairedTriples, 0u);
  EXPECT_EQ(diagram.verification.overpairedTriples, 0u);
  EXPECT_EQ(diagram.verification.truncatedTriples, 0u);
  EXPECT_EQ(diagram.verification.unclosedFaces, 0u);
  EXPECT_EQ(diagram.verification.verticesOfFullValence, diagram.verification.vertexCount);
  EXPECT_TRUE(diagram.verification.isComplete());

  std::size_t branchCount = 0;
  for (const SKApolloniusVertex& vertex : diagram.vertices) branchCount += vertex.branches.size();
  std::size_t arcs = 0;
  for (const SKApolloniusEdge& edge : diagram.edges)
    if (!edge.isLoop) ++arcs;
  EXPECT_EQ(branchCount, 2 * arcs);  // each arc consumes a branch at either end

  std::size_t v = diagram.vertices.size();
  std::size_t e = diagram.edges.size();
  std::size_t f = diagram.faces.size() / 2;  // faces are stored once per owning site
  std::size_t c = 0;
  for (const SKApolloniusCell& cell : diagram.cells)
    if (!cell.isEmpty) ++c;
  EXPECT_EQ(v + f, e + c);
}

// The full diagram must contain the free-region one: restricting to non-negative clearance may only
// remove vertices, never move or invent them, so the deepest vertex is the same in both.
TEST(apollonius, entire_space_diagram_extends_the_free_region_one)
{
  RandomSites sites = makeRandomSites(14.0, 40, 5, 0.4, 2.2);
  SKApolloniusDiagram free =
      SKApolloniusDiagram::create(sites.cell, sites.fractionalPositions, sites.radii, 1, SKApolloniusRegion::FreeSpace);
  SKApolloniusDiagram entire = SKApolloniusDiagram::create(sites.cell, sites.fractionalPositions, sites.radii, 1,
                                                          SKApolloniusRegion::EntireSpace);

  EXPECT_GT(entire.vertices.size(), free.vertices.size());
  EXPECT_NEAR(entire.largestEmptySphereRadius(), free.largestEmptySphereRadius(), 1.0e-9);

  // The extra vertices are exactly those of negative clearance, and the free diagram has none.
  std::size_t negative = 0;
  for (const SKApolloniusVertex& vertex : entire.vertices)
    if (vertex.radius < 0.0) ++negative;
  EXPECT_EQ(entire.vertices.size() - negative, free.vertices.size());
  for (const SKApolloniusVertex& vertex : free.vertices) EXPECT_GE(vertex.radius, 0.0);
}

// The diagram must report its own state rather than presenting a clipped complex as unconditionally
// finished, or a broken one as sound.
TEST(apollonius, verification_reports_incomplete_rather_than_hiding_it)
{
  RandomSites sites = makeRandomSites(14.0, 40, 5, 0.4, 2.2);
  SKApolloniusDiagram diagram =
      SKApolloniusDiagram::create(sites.cell, sites.fractionalPositions, sites.radii);

  EXPECT_GT(diagram.verification.vertexCount, 0u);
  EXPECT_EQ(diagram.verification.isComplete(),
            diagram.verification.verticesOfFullValence + diagram.verification.truncatedTriples >=
                    diagram.verification.vertexCount &&
                diagram.verification.unpairedTriples == 0 && diagram.verification.overpairedTriples == 0 &&
                diagram.verification.unclosedFaces <= 3 * diagram.verification.truncatedTriples);
}

// The deepest vertex of the diagram is the largest sphere that fits anywhere. Validated against a
// brute-force scan of the clearance field with local refinement from many starts.
TEST(apollonius, largest_empty_sphere_matches_brute_force)
{
  RandomSites sites = makeRandomSites(14.0, 40, 11, 0.4, 2.2);

  std::size_t resolution = 60;
  std::vector<std::pair<double, double3>> scan;
  for (std::size_t ix = 0; ix < resolution; ++ix)
    for (std::size_t iy = 0; iy < resolution; ++iy)
      for (std::size_t iz = 0; iz < resolution; ++iz)
      {
        double3 fractional(static_cast<double>(ix) / resolution, static_cast<double>(iy) / resolution,
                           static_cast<double>(iz) / resolution);
        scan.emplace_back(clearanceAt(sites.cell, sites.fractionalPositions, sites.radii, fractional), fractional);
      }

  std::size_t starts = 64;
  std::partial_sort(scan.begin(), scan.begin() + static_cast<std::ptrdiff_t>(starts), scan.end(),
                    [](const auto& lhs, const auto& rhs) { return lhs.first > rhs.first; });

  double reference = -std::numeric_limits<double>::max();
  for (std::size_t start = 0; start < starts; ++start)
  {
    double3 current = scan[start].second;
    double currentClearance = scan[start].first;
    for (double step = 1.0 / static_cast<double>(resolution); step > 1.0e-9; step *= 0.5)
    {
      bool improved = true;
      while (improved)
      {
        improved = false;
        for (int ix = -1; ix <= 1; ++ix)
          for (int iy = -1; iy <= 1; ++iy)
            for (int iz = -1; iz <= 1; ++iz)
            {
              double3 trial = current + step * double3(ix, iy, iz);
              double clearance = clearanceAt(sites.cell, sites.fractionalPositions, sites.radii, trial);
              if (clearance > currentClearance)
              {
                currentClearance = clearance;
                current = trial;
                improved = true;
              }
            }
      }
    }
    reference = std::max(reference, currentClearance);
  }

  SKApolloniusDiagram diagram =
      SKApolloniusDiagram::create(sites.cell, sites.fractionalPositions, sites.radii);

  // The diagram must be at least as deep as the scan. It is routinely deeper: a coordinate search
  // stalls on this field, which is a minimum of cones and so not smooth at its peak, whereas the
  // diagram solves for the peak exactly. So the binding check is not agreement with the scan but
  // that the reported radius is genuinely realised at the reported centre, which proves it is a
  // sphere that really fits rather than a number that merely came out large.
  EXPECT_GE(diagram.largestEmptySphereRadius(), reference - 1.0e-6);

  double3 fractional = sites.cell.inverse() * diagram.largestEmptySpherePosition();
  EXPECT_NEAR(clearanceAt(sites.cell, sites.fractionalPositions, sites.radii, fractional),
              diagram.largestEmptySphereRadius(), 1.0e-7);
}

// Equal radii are a case of their own, not merely an easy one.
//
// The bisector of two spheres is one sheet of a hyperboloid of revolution, and it degenerates to a plane
// exactly when the two radii match, so with equal radii every trisector is a straight line. Any ordering
// of the vertices along a trisector that appeals to the shape of the curve, such as the angle they subtend
// at a point inside it, has nothing to work with there. What makes that dangerous is that it breaks the
// complex while leaving every individual vertex perfectly correct, so a test that only checks vertices
// passes. And equal radii are no exotic input: it is the case that reduces to the ordinary Voronoi diagram.
//
// The second configuration is a small cluster in a large void, where the widest empty sphere is wider than
// the cell is long, so the arcs are long compared with the cell and the sites meeting at a vertex come from
// several different periodic images of the cluster.
TEST(apollonius, equal_radii_close_the_complex_too)
{
  std::vector<RandomSites> configurations{makeRandomSites(14.0, 40, 5, 1.0, 1.0)};
  {
    RandomSites clustered;
    double edge = 20.0;
    clustered.cell = double3x3(double3(edge, 0.0, 0.0), double3(0.0, edge, 0.0), double3(0.0, 0.0, edge));
    double3 centre(10.0, 10.0, 10.0);
    std::vector<double3> directions{double3(1.0, 0.0, 0.0), double3(0.3, 0.95, 0.1), double3(-0.8, 0.4, 0.45),
                                    double3(-0.5, -0.7, 0.5), double3(0.2, -0.4, -0.9)};
    for (const double3& direction : directions)
    {
      clustered.fractionalPositions.push_back(clustered.cell.inverse() *
                                              (centre + (4.0 / direction.length()) * direction));
      clustered.radii.push_back(1.0);
    }
    configurations.push_back(clustered);
  }

  for (const RandomSites& sites : configurations)
  {
    SKApolloniusDiagram diagram = SKApolloniusDiagram::create(sites.cell, sites.fractionalPositions, sites.radii, 1,
                                                             SKApolloniusRegion::EntireSpace);
    ASSERT_FALSE(diagram.edges.empty());

    EXPECT_EQ(diagram.verification.unpairedTriples, 0u);
    EXPECT_EQ(diagram.verification.overpairedTriples, 0u);
    EXPECT_EQ(diagram.verification.truncatedTriples, 0u);
    EXPECT_EQ(diagram.verification.unclosedFaces, 0u);
    EXPECT_EQ(diagram.verification.verticesOfFullValence, diagram.verification.vertexCount);
    EXPECT_TRUE(diagram.verification.isComplete());

    std::size_t branchCount = 0;
    for (const SKApolloniusVertex& vertex : diagram.vertices) branchCount += vertex.branches.size();
    std::size_t arcs = 0;
    for (const SKApolloniusEdge& edge : diagram.edges)
      if (!edge.isLoop) ++arcs;
    EXPECT_EQ(branchCount, 2 * arcs);

    std::size_t v = diagram.vertices.size();
    std::size_t e = diagram.edges.size();
    std::size_t f = diagram.faces.size() / 2;
    std::size_t c = 0;
    for (const SKApolloniusCell& cell : diagram.cells)
      if (!cell.isEmpty) ++c;
    EXPECT_EQ(v + f, e + c);

    // The arcs really are straight here, which is what makes this configuration the one it claims to be:
    // the midpoint of each arc sits on the chord between its endpoints rather than bulging off it.
    for (const SKApolloniusEdge& edge : diagram.edges)
    {
      if (edge.isLoop) continue;
      std::array<double3, 3> centres;
      std::array<double, 3> edgeRadii;
      for (std::size_t s = 0; s < 3; ++s)
      {
        centres[s] = sites.cell * (sites.fractionalPositions[edge.siteIndices[s]] +
                                   double3(edge.siteImages[s].x, edge.siteImages[s].y, edge.siteImages[s].z));
        edgeRadii[s] = sites.radii[edge.siteIndices[s]];
      }
      double3 from = diagram.vertices[edge.from].position;
      double3 to =
          diagram.vertices[edge.to].position + sites.cell * double3(edge.toImage.x, edge.toImage.y, edge.toImage.z);
      std::optional<SKApolloniusTangentSphere> middle = skApolloniusTrisectorPoint(centres, edgeRadii, from, to, 0.5);
      ASSERT_TRUE(middle.has_value());
      EXPECT_NEAR((middle->centre - 0.5 * (from + to)).length(), 0.0, 1.0e-9);
      EXPECT_NEAR(edge.length, (to - from).length(), 1.0e-9);
    }
  }
}

// Close packing is the configuration that is degenerate everywhere at once: every octahedral hole of a
// face-centred lattice of equal spheres has six of them cotangent and every tetrahedral hole four, and the
// holes sit on the faces and edges of the cell, so each is found many times over and in several wrappings.
// It is where both the sweep and the gathering of the copies are worked hardest, and the point of the test
// is that neither has to work hard: the sweep never refines to its floor, and the copies come together into
// as many vertices as the lattice has holes and no more.
TEST(apollonius, close_packing_is_degenerate_everywhere_and_costs_nothing_extra)
{
  double edge = 6.0;
  double3x3 latticeCell(double3(edge, 0.0, 0.0), double3(0.0, edge, 0.0), double3(0.0, 0.0, edge));
  std::vector<double3> fractionalPositions{double3(0.0, 0.0, 0.0), double3(0.0, 0.5, 0.5),
                                           double3(0.5, 0.0, 0.5), double3(0.5, 0.5, 0.0)};
  std::vector<double> radii(fractionalPositions.size(), 1.0);

  SKApolloniusDiagram diagram =
      SKApolloniusDiagram::create(latticeCell, fractionalPositions, radii, 1, SKApolloniusRegion::EntireSpace);

  // Four octahedral holes and eight tetrahedral ones to the cell, which is what the lattice has.
  EXPECT_EQ(diagram.vertices.size(), 12u);
  std::size_t sixfold = 0;
  for (const SKApolloniusVertex& vertex : diagram.vertices)
    if (vertex.siteIndices.size() == 6) ++sixfold;
  EXPECT_EQ(sixfold, 4u);
  EXPECT_EQ(diagram.verification.degenerateVertices, 4u);

  // Nothing was left behind by the gathering, and nothing was left to the refinement floor.
  EXPECT_EQ(diagram.verification.coincidentVertices, 0u);
  EXPECT_EQ(diagram.verification.degenerateSweepBoxes, 0u);
  EXPECT_TRUE(diagram.verification.isComplete());

  std::size_t v = diagram.vertices.size();
  std::size_t e = diagram.edges.size();
  std::size_t f = diagram.faces.size() / 2;
  std::size_t c = 0;
  for (const SKApolloniusCell& cell : diagram.cells)
    if (!cell.isEmpty) ++c;
  EXPECT_EQ(v + f, e + c);
}

// An edge is the narrowest passage between its endpoints, so its bottleneck cannot exceed either
// endpoint, and a sphere of that size must fit everywhere along the arc.
TEST(apollonius, edge_bottleneck_is_the_narrowest_point_of_the_arc)
{
  RandomSites sites = makeRandomSites(14.0, 40, 5, 0.4, 2.2);
  SKApolloniusDiagram diagram =
      SKApolloniusDiagram::create(sites.cell, sites.fractionalPositions, sites.radii);
  ASSERT_FALSE(diagram.edges.empty());

  for (const SKApolloniusEdge& edge : diagram.edges)
  {
    EXPECT_LE(edge.bottleneckRadius, diagram.vertices[edge.from].radius + 1.0e-9);
    EXPECT_LE(edge.bottleneckRadius, diagram.vertices[edge.to].radius + 1.0e-9);
    EXPECT_GT(edge.length, 0.0);

    // The straight chord between the endpoints is a lower bound on the arc length.
    double3 from = diagram.vertices[edge.from].position;
    double3 to = diagram.vertices[edge.to].position +
                 sites.cell * double3(edge.toImage.x, edge.toImage.y, edge.toImage.z);
    EXPECT_GE(edge.length, (to - from).length() - 1.0e-9);

    // The bottleneck is a place on the arc and not only a number, which is what lets the window across
    // the passage be measured there. It is a point of the trisector at that clearance, so all three
    // sites of the arc stand at the bottleneck radius from it, and the direction kept with it is the
    // unit tangent of the arc.
    for (std::size_t s = 0; s < 3; ++s)
    {
      double3 centre = sites.cell * (sites.fractionalPositions[edge.siteIndices[s]] +
                                     double3(edge.siteImages[s].x, edge.siteImages[s].y, edge.siteImages[s].z));
      EXPECT_NEAR((edge.bottleneckPosition - centre).length() - sites.radii[edge.siteIndices[s]],
                  edge.bottleneckRadius, 1.0e-6);
    }
    EXPECT_NEAR(edge.bottleneckDirection.length(), 1.0, 1.0e-9);
  }
}

// Three collinear centres of unequal radii, which is the arrangement in which the trisector is neither a
// branch nor a line but a circle about the axis through the centres, of one clearance the whole way round.
// It defeats both halves of the machinery at once: the tangency solve cannot write the centre of the tangent
// sphere as a function of its radius, because two of its rows have parallel coefficients, and the clearance
// cannot order points of a curve it takes one value on. The configuration is laid out so that the circle is
// known in closed form and two sites cut known arcs from it, leaving four vertices and two edges to find.
TEST(apollonius, collinear_centres_meet_in_a_circle)
{
  CollinearSites sites = makeNearlyCollinearSites(0.0);
  const double3x3& cell = sites.cell;
  const std::vector<double3>& fractionalPositions = sites.fractionalPositions;
  const std::vector<double>& radii = sites.radii;
  const std::vector<double>& ends = sites.ends;
  double circleRadius = sites.circleRadius;
  double clearance = sites.clearance;
  ASSERT_NEAR(circleRadius, 4.0, 1.0e-12);
  ASSERT_NEAR(clearance, 3.0, 1.0e-12);
  ASSERT_EQ(ends.size(), 4u);

  auto onCircle = [&](double angle)
  { return sites.origin + circleRadius * double3(std::cos(angle), std::sin(angle), 0.0); };

  // The configuration is what it is claimed to be before anything is asked of the diagram: the ends of the
  // cut arcs stand at the circle's clearance, the arcs between them stay there, and the cut arcs fall below.
  for (double angle : ends)
    ASSERT_NEAR(clearanceByScan(cell, fractionalPositions, radii, cell.inverse() * onCircle(angle), 1), clearance,
                1.0e-9);
  ASSERT_NEAR(clearanceByScan(cell, fractionalPositions, radii, cell.inverse() * onCircle(0.5 * std::numbers::pi), 1),
              clearance, 1.0e-9);
  ASSERT_NEAR(clearanceByScan(cell, fractionalPositions, radii, cell.inverse() * onCircle(1.5 * std::numbers::pi), 1),
              clearance, 1.0e-9);
  ASSERT_LT(clearanceByScan(cell, fractionalPositions, radii, cell.inverse() * onCircle(0.0), 1), clearance - 0.1);
  ASSERT_LT(clearanceByScan(cell, fractionalPositions, radii, cell.inverse() * onCircle(std::numbers::pi), 1),
            clearance - 0.1);

  for (SKApolloniusRegion region : {SKApolloniusRegion::FreeSpace, SKApolloniusRegion::EntireSpace})
  {
    SKApolloniusDiagram diagram = SKApolloniusDiagram::create(cell, fractionalPositions, radii, 1, region);

    // Every end of a cut arc is a vertex, and each is tangent to all three collinear sites.
    for (double angle : ends)
    {
      std::size_t found = 0;
      for (const SKApolloniusVertex& vertex : diagram.vertices)
      {
        if ((vertex.position - onCircle(angle)).length() > 1.0e-7) continue;
        ++found;
        EXPECT_NEAR(vertex.radius, clearance, 1.0e-7);
        for (std::size_t site = 0; site < 3; ++site)
          EXPECT_NE(std::find(vertex.siteIndices.begin(), vertex.siteIndices.end(), site), vertex.siteIndices.end());
      }
      EXPECT_EQ(found, 1u);
    }

    // The circle contributes exactly the two arcs left over from the two cuts. Each runs at the circle's
    // clearance the whole way, so that is its bottleneck, and its length is the arc it subtends.
    std::vector<double> lengths;
    for (const SKApolloniusEdge& e : diagram.edges)
    {
      if (e.siteIndices[0] > 2 || e.siteIndices[1] > 2 || e.siteIndices[2] > 2) continue;
      bool home = true;
      for (const int3& image : e.siteImages)
        if (image.x != 0 || image.y != 0 || image.z != 0) home = false;
      if (!home) continue;

      EXPECT_FALSE(e.isLoop);
      EXPECT_NEAR(e.bottleneckRadius, clearance, 1.0e-7);
      lengths.push_back(e.length);

      // Both ends are ends of cut arcs, and they are the two ends of one and the same arc.
      double3 from = diagram.vertices[e.from].position;
      double3 to = diagram.vertices[e.to].position + cell * double3(e.toImage.x, e.toImage.y, e.toImage.z);
      std::size_t fromEnd = ends.size();
      std::size_t toEnd = ends.size();
      for (std::size_t k = 0; k < ends.size(); ++k)
      {
        if ((from - onCircle(ends[k])).length() < 1.0e-7) fromEnd = k;
        if ((to - onCircle(ends[k])).length() < 1.0e-7) toEnd = k;
      }
      ASSERT_LT(fromEnd, ends.size());
      ASSERT_LT(toEnd, ends.size());

      // Sorted by azimuth the ends run as the first cut's two, then the second cut's two, so the arcs the
      // cuts leave over are the ones from an odd end to the next end round the circle. Joining any other
      // pair would mean an arc that a cutting site has been shown to occupy.
      std::size_t lower = (fromEnd % 2u == 1u) ? fromEnd : toEnd;
      std::size_t upper = (fromEnd % 2u == 1u) ? toEnd : fromEnd;
      ASSERT_EQ(lower % 2u, 1u);
      ASSERT_EQ(upper, (lower + 1u) % ends.size());

      // Sampling the arc as a polyline falls a little short of the true length, and never over it.
      double sweep = ends[upper] - ends[lower] + ((upper < lower) ? 2.0 * std::numbers::pi : 0.0);
      double subtended = circleRadius * sweep;
      EXPECT_LE(e.length, subtended + 1.0e-9);
      EXPECT_GT(e.length, 0.999 * subtended);
    }
    EXPECT_EQ(lengths.size(), 2u);

    // And the diagram as a whole closes, so nothing was traded away to get the circle right.
    EXPECT_EQ(diagram.verification.unpairedTriples, 0u);
    EXPECT_EQ(diagram.verification.unclosedFaces, 0u);
    EXPECT_EQ(diagram.verification.verticesOfFullValence, diagram.vertices.size());
    EXPECT_TRUE(diagram.verification.isComplete());
  }
}

// Centres that are nearly but not exactly collinear, which is the harder half of the same problem. Their
// trisector is a closed curve very close to a circle, so the clearance does vary along it, but by so little
// that it carries no ordering: at the smaller displacements below, the four vertices on the curve stand at
// clearances agreeing to the last bits of a double, and any ordering of them by clearance is a coin toss.
// Nothing may turn on how nearly collinear the centres are, so the whole band is walked, from exactly
// collinear up to plainly not, and the complex is required to close at every step.
TEST(apollonius, nearly_collinear_centres_close_the_complex_too)
{
  for (double displacement : {0.0, 1.0e-12, 1.0e-9, 1.0e-6, 1.0e-4, 1.0e-2, 1.0e-1})
  {
    CollinearSites sites = makeNearlyCollinearSites(displacement);
    SKApolloniusDiagram diagram = SKApolloniusDiagram::create(sites.cell, sites.fractionalPositions, sites.radii, 1,
                                                              SKApolloniusRegion::EntireSpace);

    // The three near-collinear sites still meet in one closed curve with the two cut arcs taken out of it,
    // however little their centres depart from a line.
    std::vector<double> clearances;
    std::size_t arcs = 0;
    for (const SKApolloniusEdge& e : diagram.edges)
    {
      if (e.siteIndices[0] > 2 || e.siteIndices[1] > 2 || e.siteIndices[2] > 2) continue;
      bool home = true;
      for (const int3& image : e.siteImages)
        if (image.x != 0 || image.y != 0 || image.z != 0) home = false;
      if (!home) continue;
      ++arcs;
      clearances.push_back(diagram.vertices[e.from].radius);
      clearances.push_back(diagram.vertices[e.to].radius);
    }
    EXPECT_EQ(arcs, 2u) << "at displacement " << displacement;

    // How little the clearance has to work with here. It spreads across the four vertices by about four times
    // the displacement, so with the centres exactly on a line it is the same at all four to the last bit, and
    // a displacement of 1.0e-12 leaves it the same to a part in 1.0e12, which is as closely as the vertices
    // themselves are placed. Ordering by clearance is guesswork at this end of the band; the ordering has to
    // come from the shape of the curve instead.
    ASSERT_FALSE(clearances.empty());
    double spread = *std::max_element(clearances.begin(), clearances.end()) -
                    *std::min_element(clearances.begin(), clearances.end());
    if (displacement == 0.0) EXPECT_LT(spread, 1.0e-13) << "at displacement " << displacement;
    if (displacement <= 1.0e-12) EXPECT_LT(spread, 1.0e-11) << "at displacement " << displacement;

    EXPECT_EQ(diagram.verification.unpairedTriples, 0u) << "at displacement " << displacement;
    EXPECT_EQ(diagram.verification.unclosedFaces, 0u) << "at displacement " << displacement;
    EXPECT_EQ(diagram.verification.verticesOfFullValence, diagram.vertices.size())
        << "at displacement " << displacement;
    EXPECT_TRUE(diagram.verification.isComplete()) << "at displacement " << displacement;

    std::size_t cells = 0;
    for (const SKApolloniusCell& cell : diagram.cells)
      if (!cell.isEmpty) ++cells;
    EXPECT_EQ(static_cast<long long>(diagram.vertices.size()) - static_cast<long long>(diagram.edges.size()) +
                  static_cast<long long>(diagram.faces.size() / 2) - static_cast<long long>(cells),
              0)
        << "at displacement " << displacement;
  }
}
