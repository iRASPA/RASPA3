#include <gtest/gtest.h>

import std;

import double3;
import simulationbox;
import pore_accessibility;
import exact_union_volume;
import exact_surface_patches;
import exact_solvent_excluded;

// The solvent excluded surface of an arrangement of balls, the volume it encloses, and the pore-size
// distribution that is the derivative of that volume.
//
// The cases below are of three kinds. A few are arrangements whose excluded volume is a ball or a sum of
// balls, so the answer is known before the code is run. One is the pair of atoms, which is a solid of
// revolution and so has a closed form of its own: the excluded region is bounded by the sphere of each atom
// where the probe touches it and by the torus the probe sweeps in the crease between them, and revolving that
// profile gives an integral that is elementary. It is the only closed form here that is derived without any of
// the machinery under test, and it covers both shapes the torus takes: the ordinary one, and the spindle it
// becomes when the crease is tighter than the probe and the surface folds back and meets itself at a cusp.
//
// The rest hold the volume and the distribution against each other. They are computed by routes with almost
// nothing in common -- the volume by subtracting the shell between the accessible and the excluded surface
// from the union of the inflated atoms, the distribution as an integral of the normal speed over the reentrant
// surface alone -- so a difference of the one reproducing the other is a real check on both.

namespace
{

double ballVolume(double radius) { return 4.0 / 3.0 * std::numbers::pi * radius * radius * radius; }


PoreAccessibility inflatedGeometry(const SimulationBox& box, const std::vector<double3>& positions,
                                      const std::vector<double>& bareRadii, double probeRadius)
{
  std::vector<double3> fractionalPositions;
  for (const double3& position : positions)
  {
    fractionalPositions.push_back(double3::fract(box.inverseCell * position));
  }
  // The bare radii and the probe go in separately, so that the structure carries both: the boundary is taken on
  // the inflated spheres, and the excluded surface behind it on the bare ones.
  return PoreAccessibility::create(box, fractionalPositions, bareRadii, probeRadius);
}


// The excluded volume of two equal atoms on an axis, from the profile of the region revolved about it.
//
// At heights near the far side of either atom the profile is that atom's own sphere of radius `bare`. Nearer
// the middle it is the torus the probe sweeps, whose profile is s - sqrt(r^2 - z^2) with s the radius of the
// circle the probe's centre rolls on. The two meet where the probe touches the atom, and the profile is
// continuous there. Where s is smaller than the probe the torus reaches the axis before the middle and stops:
// that is the cusp, and below it the two atoms' regions are separate.
double excludedVolumeOfTwoAtoms(double bare, double separation, double probe)
{
  const double inflated = bare + probe;
  const double half = 0.5 * separation;
  const double rolling = std::sqrt(inflated * inflated - half * half);  // s

  const double join = probe * half / inflated;                                              // where the two meet
  const double cusp = (rolling < probe) ? std::sqrt(probe * probe - rolling * rolling) : 0.0;

  // The torus part, twice: the integral of pi (s - sqrt(r^2 - z^2))^2 from the cusp up to the join.
  auto arcsine = [&](double z)
  { return 0.5 * (z * std::sqrt(probe * probe - z * z) + probe * probe * std::asin(z / probe)); };
  auto torus = [&](double z)
  {
    return (rolling * rolling + probe * probe) * z - z * z * z / 3.0 - 2.0 * rolling * arcsine(z);
  };
  double crease = 2.0 * std::numbers::pi * (torus(join) - torus(cusp));

  // The spherical part, twice: the integral of pi (bare^2 - u^2) over u from the join to the far pole.
  double lower = join - half;
  double caps = 2.0 * std::numbers::pi *
                (2.0 / 3.0 * bare * bare * bare - bare * bare * lower + lower * lower * lower / 3.0);

  return crease + caps;
}


// The area of the excluded surface of two atoms on an axis, by the same revolution as the volume above.
//
// The surface has two kinds of piece here. Away from the crease each atom is seen bare, and what hides the rest
// of it is the probe: the probe touches the atom along the direction to its own centre, so the patch is
// everything outside the cap of half angle arccos(h / R) about the direction to the neighbour, whose area is
// 2 pi a^2 (1 + h / R). In the crease the surface is the torus the probe sweeps, at a distance s - r cos psi
// from the axis, and its area is that distance integrated over the turn and round the axis. Where the rolling
// circle is smaller than the probe the turn is cut at the cusp, the surface having folded through the axis.
double excludedAreaOfTwoAtoms(double bare, double separation, double probe)
{
  const double pi = std::numbers::pi;
  const double inflated = bare + probe;
  const double half = 0.5 * separation;
  const double rolling = std::sqrt(inflated * inflated - half * half);  // s

  const double caps = 2.0 * (2.0 * pi * bare * bare * (1.0 + half / inflated));

  // sin psi = h / R at the tangency, and the cusp is where the distance to the axis reaches zero.
  const double outer = std::asin(std::clamp(half / inflated, -1.0, 1.0));
  const double inner = (rolling < probe) ? std::acos(std::clamp(rolling / probe, -1.0, 1.0)) : 0.0;
  const double band =
      (outer > inner) ? 2.0 * (2.0 * pi * probe *
                               (rolling * (outer - inner) - probe * (std::sin(outer) - std::sin(inner))))
                      : 0.0;

  return caps + band;
}


// The derivative of the excluded volume in the probe radius, by a central difference of the routine's own
// volume, with one Richardson step to take out the truncation error of the step.
//
// The volume has a corner in the probe radius wherever the surface changes shape: a patch appearing, a cusp
// opening on a crease, a neighbouring probe beginning to reach into a concave patch. A difference taken across
// a corner is not an estimate of the slope on either side of it, so the shape is read at every point of the
// stencil and the comparison is refused where it is not the same throughout. That refusal is not a weakening
// of the test: an equilateral triangle passes through such a corner at a particular probe radius, and asking
// for the slope there is asking for something that does not exist.
struct DifferencedSlope
{
  double slope{0.0};
  bool comparable{false};
};

DifferencedSlope differencedDistribution(const SimulationBox& box, const std::vector<double3>& positions,
                                         const std::vector<double>& bareRadii, double probeRadius, double step)
{
  auto shapeAt = [&](double probe)
  {
    SolventExcludedGeometry geometry =
        solventExcludedGeometry(inflatedGeometry(box, positions, bareRadii, probe), probe);
    const SolventExcludedGeometry::Diagnostics& counts = geometry.diagnostics;
    std::array<std::size_t, 4> shape = {counts.numberOfArcs, counts.cuspedArcs, counts.numberOfVertices,
                                        counts.clippedVertices};
    return std::make_pair(geometry.excludedVolume, shape);
  };

  std::array<std::size_t, 4> shape = shapeAt(probeRadius).second;
  std::array<double, 4> volumes{};
  const std::array<double, 4> offsets = {-2.0 * step, -step, step, 2.0 * step};
  DifferencedSlope answer;
  for (std::size_t i = 0; i < offsets.size(); ++i)
  {
    auto [volume, other] = shapeAt(probeRadius + offsets[i]);
    if (other != shape) return answer;
    volumes[i] = volume;
  }

  double near = (volumes[2] - volumes[1]) / (2.0 * step);
  double far = (volumes[3] - volumes[0]) / (4.0 * step);
  answer.slope = (4.0 * near - far) / 3.0;
  answer.comparable = true;
  return answer;
}

}  // namespace


// The region of a sphere left by a set of caps, in the cases where the answer is a fraction of the sphere.
TEST(exact_solvent_excluded, a_region_of_a_sphere_cut_out_by_caps)
{
  const double pi = std::numbers::pi;

  // No cap at all leaves the whole sphere, whose moment vanishes by symmetry.
  SphericalRegion whole = regionOutsideCaps({});
  EXPECT_NEAR(whole.solidAngle, 4.0 * pi, 1.0e-12);
  EXPECT_NEAR(whole.moment.length(), 0.0, 1.0e-12);

  // One cap of half angle gamma leaves 2 pi (1 + cos gamma), and the moment of what is left points away
  // from the cap: the integral of the direction over a cap is pi sin^2(gamma) along its axis.
  for (double gamma : {0.3, 1.0, std::numbers::pi / 2.0, 2.4})
  {
    SphericalRegion left = regionOutsideCaps({{double3(0.0, 0.0, 1.0), std::cos(gamma)}});
    EXPECT_NEAR(left.solidAngle, 2.0 * pi * (1.0 + std::cos(gamma)), 1.0e-9) << "half angle " << gamma;
    EXPECT_NEAR(left.moment.z, -pi * std::sin(gamma) * std::sin(gamma), 1.0e-9) << "half angle " << gamma;
    EXPECT_NEAR(std::hypot(left.moment.x, left.moment.y), 0.0, 1.0e-9);
  }

  // Three great circles through the centre, mutually perpendicular, leave an octant. The caps here are the
  // complements of the half spaces kept, which is how a concave patch is cut out.
  SphericalRegion octant = regionOutsideCaps({{double3(-1.0, 0.0, 0.0), 0.0},
                                              {double3(0.0, -1.0, 0.0), 0.0},
                                              {double3(0.0, 0.0, -1.0), 0.0}});
  EXPECT_NEAR(octant.solidAngle, 0.5 * pi, 1.0e-9);
  // Over an octant the integral of each component of the direction is pi/4.
  EXPECT_NEAR(octant.moment.x, 0.25 * pi, 1.0e-9);
  EXPECT_NEAR(octant.moment.y, 0.25 * pi, 1.0e-9);
  EXPECT_NEAR(octant.moment.z, 0.25 * pi, 1.0e-9);

  // A cap that covers the sphere leaves nothing, and one that reaches nothing leaves all of it.
  EXPECT_NEAR(regionOutsideCaps({{double3(0.0, 0.0, 1.0), -1.0}}).solidAngle, 0.0, 1.0e-12);
  EXPECT_NEAR(regionOutsideCaps({{double3(0.0, 0.0, 1.0), 1.0}}).solidAngle, 4.0 * pi, 1.0e-12);
}


// One atom with room around it. Nothing is reentrant, so the excluded region is the atom itself whatever the
// probe, and the distribution is zero: a convex solid has no pore size at all.
TEST(exact_solvent_excluded, a_lone_atom_excludes_itself)
{
  double a = 40.0;
  SimulationBox box(a, a, a);
  const double bare = 1.7;

  for (double probe : {0.0, 0.5, 1.4, 3.0})
  {
    SolventExcludedGeometry geometry =
        solventExcludedGeometry(inflatedGeometry(box, {double3(20.0, 20.0, 20.0)}, {bare}, probe), probe);

    EXPECT_NEAR(geometry.excludedVolume, ballVolume(bare), 1.0e-9 * ballVolume(bare)) << "probe " << probe;
    EXPECT_NEAR(geometry.distribution, 0.0, 1.0e-9) << "probe " << probe;
    EXPECT_EQ(geometry.diagnostics.numberOfArcs, 0uz);
    EXPECT_EQ(geometry.diagnostics.numberOfVertices, 0uz);
  }
}


// Two atoms whose inflated spheres do not reach one another. The probe can pass between them, so there is no
// crease to fill and the two atoms are excluded separately.
TEST(exact_solvent_excluded, two_atoms_out_of_reach_of_one_another)
{
  double a = 60.0;
  SimulationBox box(a, a, a);
  const double bare = 1.7;
  const double probe = 1.0;

  SolventExcludedGeometry geometry = solventExcludedGeometry(
      inflatedGeometry(box, {double3(20.0, 30.0, 30.0), double3(30.0, 30.0, 30.0)}, {bare, bare}, probe), probe);

  EXPECT_NEAR(geometry.excludedVolume, 2.0 * ballVolume(bare), 1.0e-9 * ballVolume(bare));
  EXPECT_NEAR(geometry.distribution, 0.0, 1.0e-9);
}


// Two atoms in the crease of one another, against the closed form for the solid of revolution. The last two
// separations are tight enough that the rolling circle is smaller than the probe, which puts the cusp on the
// toroidal patch and is the case a naive integral over the surface gets wrong.
TEST(exact_solvent_excluded, two_atoms_against_the_solid_of_revolution)
{
  double a = 60.0;
  SimulationBox box(a, a, a);

  // bare radius, separation, probe
  const std::array<std::array<double, 3>, 7> cases = {{{1.7, 4.0, 1.0},
                                                       {1.7, 4.4, 1.4},
                                                       {1.7, 5.2, 1.4},
                                                       {1.2, 3.0, 0.8},
                                                       {2.0, 4.6, 2.0},
                                                       {1.7, 4.2, 1.8},
                                                       {1.5, 3.4, 1.6}}};

  for (const std::array<double, 3>& one : cases)
  {
    const double bare = one[0];
    const double separation = one[1];
    const double probe = one[2];

    std::vector<double3> positions = {double3(30.0 - 0.5 * separation, 30.0, 30.0),
                                      double3(30.0 + 0.5 * separation, 30.0, 30.0)};
    SolventExcludedGeometry geometry =
        solventExcludedGeometry(inflatedGeometry(box, positions, {bare, bare}, probe), probe);

    double expected = excludedVolumeOfTwoAtoms(bare, separation, probe);
    EXPECT_NEAR(geometry.excludedVolume, expected, 1.0e-8 * expected)
        << "bare " << bare << " separation " << separation << " probe " << probe;

    // The pair has one circle of intersection and no vertex, so this is the toroidal term on its own.
    EXPECT_EQ(geometry.diagnostics.numberOfVertices, 0uz);
    EXPECT_GT(geometry.diagnostics.numberOfArcs, 0uz);
  }
}


// The same pair, now with the distribution held against a difference of the volume. It is the toroidal term
// alone that has to reproduce the slope, and both shapes of it are covered.
TEST(exact_solvent_excluded, the_distribution_of_a_pair_is_the_slope_of_its_volume)
{
  double a = 60.0;
  SimulationBox box(a, a, a);
  const double bare = 1.7;

  std::size_t compared = 0;
  for (double separation : {4.0, 4.6, 5.4, 6.2})
  {
    for (double probe : {0.9, 1.4, 2.0})
    {
      std::vector<double3> positions = {double3(30.0 - 0.5 * separation, 30.0, 30.0),
                                        double3(30.0 + 0.5 * separation, 30.0, 30.0)};
      SolventExcludedGeometry geometry =
          solventExcludedGeometry(inflatedGeometry(box, positions, {bare, bare}, probe), probe);
      if (geometry.diagnostics.numberOfArcs == 0) continue;  // out of reach at this probe, and both sides are zero

      DifferencedSlope differenced = differencedDistribution(box, positions, {bare, bare}, probe, 0.01);
      if (!differenced.comparable) continue;

      ++compared;
      EXPECT_NEAR(geometry.distribution, differenced.slope, 2.0e-4 * std::max(1.0, std::abs(differenced.slope)))
          << "separation " << separation << " probe " << probe;
    }
  }
  EXPECT_GE(compared, 8uz);
}


// Three atoms in a triangle, which is the smallest arrangement with a vertex in it: the probe touches all
// three at once, and the concave patch it leaves enters both the volume and the distribution. There are two
// such patches, one on either side of the plane, and the tighter triangles bring them into one another so that
// each has to be clipped against the other; the last cases, where the probe is larger than the atoms, have the
// creases folding back on themselves as well, so that the clipped concave patch and the cusped torus are both
// under test at once.
TEST(exact_solvent_excluded, the_distribution_of_a_triangle_is_the_slope_of_its_volume)
{
  double a = 60.0;
  SimulationBox box(a, a, a);

  // bare radius, side of the triangle, probe. A vertex needs the circumradius within reach of the inflated
  // spheres and a cusp needs the triangle wider than twice sqrt(bare^2 + 2 bare probe), and the two together
  // ask for a probe larger than the atoms.
  const std::array<std::array<double, 3>, 8> cases = {{{1.7, 4.4, 0.9},
                                                       {1.7, 4.4, 1.4},
                                                       {1.7, 4.8, 1.4},
                                                       {1.7, 4.4, 1.8},
                                                       {1.7, 4.8, 1.8},
                                                       {1.7, 5.6, 1.8},
                                                       {1.0, 4.6, 2.0},
                                                       {1.0, 4.9, 2.0}}};

  std::size_t compared = 0;
  std::size_t withClipping = 0;
  std::size_t withCusps = 0;
  for (const std::array<double, 3>& one : cases)
  {
    const double bare = one[0];
    const double side = one[1];
    const double probe = one[2];

    const double reach = side / std::sqrt(3.0);
    std::vector<double3> positions;
    for (std::size_t k = 0; k < 3; ++k)
    {
      double angle = 2.0 * std::numbers::pi * static_cast<double>(k) / 3.0;
      positions.push_back(double3(30.0 + reach * std::cos(angle), 30.0 + reach * std::sin(angle), 30.0));
    }
    std::vector<double> radii = {bare, bare, bare};

    SolventExcludedGeometry geometry = solventExcludedGeometry(inflatedGeometry(box, positions, radii, probe), probe);
    ASSERT_GT(geometry.diagnostics.numberOfVertices, 0uz) << "bare " << bare << " side " << side << " probe " << probe;

    DifferencedSlope differenced = differencedDistribution(box, positions, radii, probe, 0.01);
    if (!differenced.comparable) continue;

    ++compared;
    if (geometry.diagnostics.clippedVertices > 0) ++withClipping;
    if (geometry.diagnostics.cuspedArcs > 0) ++withCusps;
    EXPECT_NEAR(geometry.distribution, differenced.slope, 2.0e-4 * std::max(1.0, std::abs(differenced.slope)))
        << "bare " << bare << " side " << side << " probe " << probe << ", vertices "
        << geometry.diagnostics.numberOfVertices << " of them " << geometry.diagnostics.clippedVertices
        << " clipped, cusped arcs " << geometry.diagnostics.cuspedArcs;
  }

  // Both of the things the reentrant surface has to cope with have to have been met along the way, or the
  // comparison has only been made where it is easy.
  EXPECT_GE(compared, 6uz);
  EXPECT_GE(withClipping, 1uz);
  EXPECT_GE(withCusps, 1uz);
}


// At zero probe radius the excluded surface is the surface of the bare atoms, so the excluded volume is the
// volume of their union and the shell between the two surfaces is nothing at all. It is the one probe radius
// at which the answer is already known by another route in the same code, and it has to agree exactly.
TEST(exact_solvent_excluded, a_vanishing_probe_leaves_the_union_of_the_atoms)
{
  double a = 12.0;
  SimulationBox box(a, a, a);

  std::vector<double3> positions;
  std::vector<double> radii;
  const std::size_t side = 3;
  for (std::size_t i = 0; i < side; ++i)
  {
    for (std::size_t j = 0; j < side; ++j)
    {
      for (std::size_t k = 0; k < side; ++k)
      {
        positions.push_back(double3(static_cast<double>(i) * 4.0 + 0.7, static_cast<double>(j) * 4.0 + 1.3,
                                    static_cast<double>(k) * 4.0 + 2.1));
        radii.push_back(1.5 + 0.1 * static_cast<double>((i + j + k) % 3));
      }
    }
  }

  PoreAccessibility geometry = inflatedGeometry(box, positions, radii, 0.0);
  SolventExcludedGeometry excluded = solventExcludedGeometry(geometry, 0.0);

  double expected = unionOfBallsVolume(geometry);
  EXPECT_NEAR(excluded.excludedVolume, expected, 1.0e-9 * expected);
  EXPECT_NEAR(excluded.shellVolume, 0.0, 1.0e-9 * expected);
  EXPECT_NEAR(excluded.distribution, 0.0, 1.0e-9);
}


// The three kinds of patch on the excluded surface, one arrangement at a time, each against a closed form
// derived without any of the machinery under test.
//
// A lone atom is convex and nothing else, whatever the probe: there is no crease to roll into and nowhere to
// wedge. It is the only case where the area is known without any integration at all.
TEST(exact_solvent_excluded, a_lone_atom_is_all_convex)
{
  double a = 40.0;
  SimulationBox box(a, a, a);
  const double bare = 1.7;
  const double sphere = 4.0 * std::numbers::pi * bare * bare;

  for (double probe : {0.0, 0.5, 1.4, 3.0})
  {
    SolventExcludedGeometry geometry =
        solventExcludedGeometry(inflatedGeometry(box, {double3(20.0, 20.0, 20.0)}, {bare}, probe), probe);

    EXPECT_NEAR(geometry.area.convex, sphere, 1.0e-9 * sphere) << "probe " << probe;
    EXPECT_NEAR(geometry.area.saddle, 0.0, 1.0e-12) << "probe " << probe;
    EXPECT_NEAR(geometry.area.concave, 0.0, 1.0e-12) << "probe " << probe;
    EXPECT_NEAR(geometry.area.convexFraction(), 1.0, 1.0e-12) << "probe " << probe;
  }
}


// Two atoms in the crease of one another: two convex caps and the band of the torus between them, against the
// closed form for the same surface of revolution the volume was checked against. The last cases are tight
// enough that the band folds back on itself at a cusp, which is where the turn has to be cut and where an
// integral over the whole torus would return more area than the surface has.
TEST(exact_solvent_excluded, two_atoms_against_the_surface_of_revolution)
{
  double a = 60.0;
  SimulationBox box(a, a, a);

  // bare radius, separation, probe. The rolling circle is smaller than the probe, and so the band cusped, once
  // the atoms are further apart than twice sqrt(a^2 + 2 a r): the last two are, and the rest are not.
  const std::array<std::array<double, 3>, 9> cases = {{{1.7, 4.0, 1.0},
                                                       {1.7, 4.4, 1.4},
                                                       {1.7, 5.2, 1.4},
                                                       {1.2, 3.0, 0.8},
                                                       {2.0, 4.6, 2.0},
                                                       {1.7, 4.2, 1.8},
                                                       {1.5, 3.4, 1.6},
                                                       {1.7, 6.4, 1.8},
                                                       {1.2, 5.0, 1.6}}};

  std::size_t withCusps = 0;
  for (const std::array<double, 3>& one : cases)
  {
    const double bare = one[0];
    const double separation = one[1];
    const double probe = one[2];

    std::vector<double3> positions = {double3(30.0 - 0.5 * separation, 30.0, 30.0),
                                      double3(30.0 + 0.5 * separation, 30.0, 30.0)};
    SolventExcludedGeometry geometry =
        solventExcludedGeometry(inflatedGeometry(box, positions, {bare, bare}, probe), probe);

    const double expected = excludedAreaOfTwoAtoms(bare, separation, probe);
    const double caps = 2.0 * (2.0 * std::numbers::pi * bare * bare * (1.0 + 0.5 * separation / (bare + probe)));

    if (geometry.diagnostics.cuspedArcs > 0) ++withCusps;
    EXPECT_NEAR(geometry.area.total(), expected, 1.0e-9 * expected)
        << "bare " << bare << " separation " << separation << " probe " << probe;
    EXPECT_NEAR(geometry.area.convex, caps, 1.0e-9 * caps)
        << "bare " << bare << " separation " << separation << " probe " << probe;
    EXPECT_NEAR(geometry.area.saddle, expected - caps, 1.0e-9 * expected)
        << "bare " << bare << " separation " << separation << " probe " << probe;
    EXPECT_NEAR(geometry.area.concave, 0.0, 1.0e-12);
  }
  EXPECT_GE(withCusps, 2uz);
}


// Three atoms in a triangle, for the concave patch. The probe wedges into the triangle from either side and
// stops, and the surface there is the piece of its own sphere facing the void: the spherical triangle whose
// corners are the three directions it touches the atoms in, its edges the great circles it leaves the tori
// along. Girard's theorem gives its area from the angles alone, and for an equilateral triangle those follow
// from the side. Sides are chosen so that the two patches are out of reach of one another, since a clipped
// patch is no longer a triangle and Girard says nothing about it.
TEST(exact_solvent_excluded, a_triangle_leaves_two_spherical_triangles)
{
  double a = 60.0;
  SimulationBox box(a, a, a);

  // bare radius, side of the triangle, probe. The two patches are out of reach of one another as long as the
  // probe's rest is further above the plane than the probe's own radius, which asks for 3 (a + r)^2 > side^2 + 3 r^2.
  const std::array<std::array<double, 3>, 4> cases = {
      {{1.7, 4.0, 1.0}, {1.7, 4.2, 1.0}, {1.7, 4.0, 1.2}, {1.2, 3.2, 0.9}}};

  for (const std::array<double, 3>& one : cases)
  {
    const double bare = one[0];
    const double side = one[1];
    const double probe = one[2];
    const double inflated = bare + probe;

    const double reach = side / std::sqrt(3.0);
    std::vector<double3> positions;
    for (std::size_t k = 0; k < 3; ++k)
    {
      double angle = 2.0 * std::numbers::pi * static_cast<double>(k) / 3.0;
      positions.push_back(double3(30.0 + reach * std::cos(angle), 30.0 + reach * std::sin(angle), 30.0));
    }

    SolventExcludedGeometry geometry =
        solventExcludedGeometry(inflatedGeometry(box, positions, {bare, bare, bare}, probe), probe);
    ASSERT_EQ(geometry.diagnostics.numberOfVertices, 2uz) << "bare " << bare << " side " << side << " probe " << probe;
    ASSERT_EQ(geometry.diagnostics.clippedVertices, 0uz) << "bare " << bare << " side " << side << " probe " << probe;

    // The angular side of the spherical triangle, from the two directions subtending the atoms' separation, and
    // the angle at a corner of an equilateral one.
    const double cosineSide = 1.0 - 0.5 * side * side / (inflated * inflated);
    const double corner = std::acos(std::clamp(cosineSide / (1.0 + cosineSide), -1.0, 1.0));
    const double expected = 2.0 * probe * probe * (3.0 * corner - std::numbers::pi);

    // To the quadrature of the sweep, which is what the region of the sphere is measured by rather than by
    // Girard's theorem: the two agree here because the region happens to be a triangle, and only here.
    EXPECT_NEAR(geometry.area.concave, expected, 1.0e-7 * expected)
        << "bare " << bare << " side " << side << " probe " << probe;
  }
}


// At zero probe radius the excluded surface is the surface of the bare atoms themselves, so it is all convex
// and its area is the area of the union, which the same sweep returns by the route the surface area uses. As
// the probe grows from there the creases and the wedges take over, and what they take is continuous: they eat
// into the convex area and give back their own, so the total leaves the union's area smoothly.
TEST(exact_solvent_excluded, a_vanishing_probe_leaves_the_surface_of_the_atoms)
{
  double a = 12.0;
  SimulationBox box(a, a, a);

  std::vector<double3> positions;
  std::vector<double> radii;
  for (std::size_t i = 0; i < 3; ++i)
  {
    for (std::size_t j = 0; j < 3; ++j)
    {
      for (std::size_t k = 0; k < 3; ++k)
      {
        positions.push_back(double3(static_cast<double>(i) * 4.0 + 0.7, static_cast<double>(j) * 4.0 + 1.3,
                                    static_cast<double>(k) * 4.0 + 2.1));
        radii.push_back(1.5 + 0.1 * static_cast<double>((i + j + k) % 3));
      }
    }
  }

  const double union0 = exactSurfaceArea(inflatedGeometry(box, positions, radii, 0.0), 1).area;

  SolventExcludedGeometry bare = solventExcludedGeometry(inflatedGeometry(box, positions, radii, 0.0), 0.0);
  EXPECT_NEAR(bare.area.convex, union0, 1.0e-9 * union0);
  EXPECT_NEAR(bare.area.reentrant(), 0.0, 1.0e-12);

  // The reentrant area appears at first order in the probe, the creases being lines and the wedges points, so
  // it goes with the probe rather than with its square: halving the probe halves the fraction.
  double previous = 0.0;
  for (double probe : {0.08, 0.04, 0.02})
  {
    SolventExcludedGeometry geometry = solventExcludedGeometry(inflatedGeometry(box, positions, radii, probe), probe);
    EXPECT_LT(std::abs(geometry.area.total() - union0), 0.1 * union0) << "probe " << probe;
    if (previous > 0.0) EXPECT_LT(geometry.area.reentrantFraction(), previous) << "probe " << probe;
    previous = geometry.area.reentrantFraction();
  }
}


// The three kinds account for the wall and no more, and so do the sides, whatever the arrangement. Both are
// bookkeeping rather than geometry, and both are the sort of thing that a patch counted on two sides or left off
// one would break silently.
TEST(exact_solvent_excluded, the_kinds_and_the_sides_each_account_for_the_whole_wall)
{
  double a = 11.0;
  SimulationBox box(a, a, a);

  // A lattice with one atom left out, which leaves a pocket the probe cannot reach from outside, so that the
  // sides being added up is a test of something rather than of one side being everything.
  std::vector<double3> positions;
  std::vector<double> radii;
  for (std::size_t i = 0; i < 4; ++i)
  {
    for (std::size_t j = 0; j < 4; ++j)
    {
      for (std::size_t k = 0; k < 4; ++k)
      {
        if (i == 1 && j == 1 && k == 1) continue;
        positions.push_back(double3(static_cast<double>(i), static_cast<double>(j), static_cast<double>(k)) * 2.75);
        radii.push_back(1.4);
      }
    }
  }

  for (double probe : {0.3, 0.7, 1.1})
  {
    SolventExcludedGeometry geometry = solventExcludedGeometry(inflatedGeometry(box, positions, radii, probe), probe);
    const double total = geometry.area.total();
    ASSERT_GT(total, 0.0) << "probe " << probe;

    EXPECT_NEAR(geometry.area.convexFraction() + geometry.area.saddleFraction() + geometry.area.concaveFraction(), 1.0,
                1.0e-12)
        << "probe " << probe;

    const double bySide =
        geometry.accessibleArea.total() + geometry.inaccessibleArea.total() + geometry.undecidedArea.total();
    EXPECT_NEAR(bySide, total, 1.0e-9 * total) << "probe " << probe;

    for (double ExcludedSurfaceAreas::*kind :
         {&ExcludedSurfaceAreas::convex, &ExcludedSurfaceAreas::saddle, &ExcludedSurfaceAreas::concave})
    {
      double sum = geometry.accessibleArea.*kind + geometry.inaccessibleArea.*kind + geometry.undecidedArea.*kind;
      EXPECT_NEAR(sum, geometry.area.*kind, 1.0e-9 * total) << "probe " << probe;
    }

    // Nothing negative, and the concave patches are on the probe's own sphere so they cannot exceed one whole
    // sphere apiece.
    EXPECT_GE(geometry.area.convex, 0.0);
    EXPECT_GE(geometry.area.saddle, 0.0);
    EXPECT_GE(geometry.area.concave, 0.0);
    EXPECT_LE(geometry.area.concave, 4.0 * std::numbers::pi * probe * probe *
                                         static_cast<double>(geometry.diagnostics.numberOfVertices) + 1.0e-9);
  }
}


// A shell of overlapping atoms with room inside it, which is a cage nothing can get into and out of: twelve
// atoms at the vertices of an icosahedron, large enough that the shell closes over the faces as well as the
// edges. Whatever the probe, the room inside is a pore of its own, walled by one closed surface and sealed off
// from the void outside.
std::pair<std::vector<double3>, std::vector<double>> icosahedralCage(const double3& centre, double shellRadius,
                                                                    double atomRadius)
{
  const double golden = 0.5 * (1.0 + std::sqrt(5.0));
  const double scale = shellRadius / std::sqrt(1.0 + golden * golden);

  std::vector<double3> positions;
  std::vector<double> radii;
  for (double first : {-1.0, 1.0})
  {
    for (double second : {-golden, golden})
    {
      positions.push_back(centre + double3(0.0, first, second) * scale);
      positions.push_back(centre + double3(first, second, 0.0) * scale);
      positions.push_back(centre + double3(second, 0.0, first) * scale);
    }
  }
  radii.assign(positions.size(), atomRadius);
  return {positions, radii};
}


// What a sealed cage holds is its own share of the pore volume, and it has to behave as a pore volume does:
// fall as the probe grows, and go to nothing when the probe no longer fits inside.
//
// The cage's share is not measured the way the total is. The total is the cell less the excluded volume; the
// cage's is what its own closed surface encloses, by the divergence theorem over the surface's own arcs, plus
// the shell standing over that surface, which is a sum over the patches, creases and corners belonging to it.
// So this holds a decomposition against a property of the thing decomposed, and it is what a piece of the
// shell put behind the wrong surface breaks: the wall of a cage has patches facing the cage and patches facing
// the void outside, and a corner attributed across it moves volume from the one pore to the other. That shows
// up here and nowhere else, the total being right either way.
TEST(exact_solvent_excluded, a_sealed_cage_holds_a_pore_volume_that_only_falls)
{
  const double a = 10.0;
  SimulationBox box(a, a, a);

  // Room for a probe of radius 0.9 inside, and sealed to any probe at all: the nearest a face of the shell
  // comes to being open is 2.5 sin(37.377 degrees) = 1.518, which is inside the atoms.
  const double shellRadius = 2.5;
  const double atomRadius = 1.6;
  auto [positions, radii] = icosahedralCage(double3(0.5 * a, 0.5 * a, 0.5 * a), shellRadius, atomRadius);

  double previous = std::numeric_limits<double>::max();
  for (double probe = 0.1; probe < 0.9; probe += 0.05)
  {
    PoreAccessibility accessibility = inflatedGeometry(box, positions, radii, probe);
    SolventExcludedGeometry geometry = solventExcludedGeometry(accessibility, probe);

    // The cage is the only pore shut off from the void running past it outside, and the room in it is worth
    // something at every probe up to the one that fills it.
    EXPECT_GT(geometry.inaccessiblePoreVolume, 0.0) << "probe " << probe;
    EXPECT_NEAR(geometry.undecidedPoreVolume, 0.0, 1.0e-9) << "probe " << probe;
    EXPECT_LT(geometry.inaccessiblePoreVolume, geometry.poreVolume) << "probe " << probe;

    // The room the cage opens to a larger probe is room it opens to a smaller one, so this cannot rise. A
    // corner of the shell attributed through the wall moves several cubic Angstrom at a step and is caught
    // here by a margin of orders of magnitude rather than by a tolerance.
    EXPECT_LE(geometry.inaccessiblePoreVolume, previous + 1.0e-8) << "probe " << probe;
    previous = geometry.inaccessiblePoreVolume;

    // Every piece of the shell belongs to one of the surfaces, so the pieces add up to the whole of it.
    double byComponent = 0.0;
    for (double piece : geometry.componentShellVolume) byComponent += piece;
    EXPECT_NEAR(byComponent, geometry.shellVolume, 1.0e-8 * std::max(1.0, geometry.shellVolume)) << "probe " << probe;
  }

  // Past the room inside there is nowhere in the cage for the probe's centre, and the cage holds nothing.
  PoreAccessibility filled = inflatedGeometry(box, positions, radii, shellRadius - atomRadius + 0.05);
  SolventExcludedGeometry sealed = solventExcludedGeometry(filled, shellRadius - atomRadius + 0.05);
  EXPECT_NEAR(sealed.inaccessiblePoreVolume, 0.0, 1.0e-8);
}


// A probe too large to fit anywhere leaves no pore volume at all, and none of it is negative on the way
// there: the volume a probe opens up cannot grow as the probe does.
TEST(exact_solvent_excluded, the_pore_volume_falls_to_nothing_and_never_rises)
{
  double a = 9.0;
  SimulationBox box(a, a, a);

  std::vector<double3> positions;
  std::vector<double> radii;
  for (std::size_t i = 0; i < 3; ++i)
  {
    for (std::size_t j = 0; j < 3; ++j)
    {
      for (std::size_t k = 0; k < 3; ++k)
      {
        positions.push_back(
            double3(static_cast<double>(i) * 3.0, static_cast<double>(j) * 3.0, static_cast<double>(k) * 3.0));
        radii.push_back(1.2);
      }
    }
  }

  double previous = std::numeric_limits<double>::max();
  for (double probe = 0.0; probe < 1.6; probe += 0.2)
  {
    SolventExcludedGeometry excluded =
        solventExcludedGeometry(inflatedGeometry(box, positions, radii, probe), probe);
    EXPECT_GE(excluded.poreVolume, -1.0e-8) << "probe " << probe;
    EXPECT_LE(excluded.poreVolume, previous + 1.0e-8) << "probe " << probe;
    EXPECT_GE(excluded.distribution, -1.0e-8) << "probe " << probe;
    previous = excluded.poreVolume;
  }

  // Well past the point where the atoms of this lattice touch one another there is nowhere left for the
  // probe's centre, so the whole cell is excluded and nothing is open.
  SolventExcludedGeometry sealed = solventExcludedGeometry(inflatedGeometry(box, positions, radii, 2.5), 2.5);
  EXPECT_NEAR(sealed.poreVolume, 0.0, 1.0e-8);
  EXPECT_NEAR(sealed.excludedVolume, box.volume, 1.0e-8 * box.volume);
}
