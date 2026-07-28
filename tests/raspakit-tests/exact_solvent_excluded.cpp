#include <gtest/gtest.h>

import std;

import double3;
import simulationbox;
import voronoi_accessibility;
import exact_union_volume;
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


VoronoiAccessibility inflatedGeometry(const SimulationBox& box, const std::vector<double3>& positions,
                                      const std::vector<double>& bareRadii, double probeRadius)
{
  std::vector<double3> fractionalPositions;
  std::vector<double> inflated;
  for (std::size_t i = 0; i < positions.size(); ++i)
  {
    fractionalPositions.push_back(double3::fract(box.inverseCell * positions[i]));
    inflated.push_back(bareRadii[i] + probeRadius);
  }
  // Inflated once, here, so that the sphere the boundary is taken on is the same one the excluded volume is
  // measured against; the routine takes the probe radius separately and reads the bare radii back from it.
  return VoronoiAccessibility::create(box, fractionalPositions, inflated, 0.0);
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
    std::array<std::size_t, 4> shape = {geometry.numberOfArcs, geometry.cuspedArcs, geometry.numberOfVertices,
                                        geometry.clippedVertices};
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
    EXPECT_EQ(geometry.numberOfArcs, 0uz);
    EXPECT_EQ(geometry.numberOfVertices, 0uz);
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
    EXPECT_EQ(geometry.numberOfVertices, 0uz);
    EXPECT_GT(geometry.numberOfArcs, 0uz);
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
      if (geometry.numberOfArcs == 0) continue;  // out of reach at this probe, and both sides are zero

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
    ASSERT_GT(geometry.numberOfVertices, 0uz) << "bare " << bare << " side " << side << " probe " << probe;

    DifferencedSlope differenced = differencedDistribution(box, positions, radii, probe, 0.01);
    if (!differenced.comparable) continue;

    ++compared;
    if (geometry.clippedVertices > 0) ++withClipping;
    if (geometry.cuspedArcs > 0) ++withCusps;
    EXPECT_NEAR(geometry.distribution, differenced.slope, 2.0e-4 * std::max(1.0, std::abs(differenced.slope)))
        << "bare " << bare << " side " << side << " probe " << probe << ", vertices " << geometry.numberOfVertices
        << " of them " << geometry.clippedVertices << " clipped, cusped arcs " << geometry.cuspedArcs;
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

  VoronoiAccessibility geometry = inflatedGeometry(box, positions, radii, 0.0);
  SolventExcludedGeometry excluded = solventExcludedGeometry(geometry, 0.0);

  double expected = unionOfBallsVolume(geometry);
  EXPECT_NEAR(excluded.excludedVolume, expected, 1.0e-9 * expected);
  EXPECT_NEAR(excluded.shellVolume, 0.0, 1.0e-9 * expected);
  EXPECT_NEAR(excluded.distribution, 0.0, 1.0e-9);
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
