#include <gtest/gtest.h>

import std;

import double3;
import unit_cell;
import randomnumbers;
import brute_force_structure;
import brute_force_surface_area;
import brute_force_solvent_excluded;

// The solvent-excluded surface as the brute force measures it: the concave patches, where a probe wedged
// against three atoms has part of its own surface on the boundary.

TEST(brute_force_solvent_excluded, a_concave_patch_reaching_past_a_right_angle_from_its_middle_is_all_sampled)
{
  // A probe wedged against three atoms, two of them further than a right angle from the mean direction to the
  // three. The patch is the spherical triangle between the three contact directions, and its side between those
  // two swings out beyond them, further from the middle than any contact is. The probe's mirror position on the
  // other side of the atoms is wedged just the same, and the two probes each cover part of the other's patch;
  // that is the only burial there is here, the surface sample being left empty, and it is worked out below
  // from the geometry alone. A draw confined to the cap that reaches only the contacts misses the swing.
  const double probe = 1.0;
  const double atomRadius = 1.5;
  const double reach = atomRadius + probe;
  const UnitCell box(30.0, 30.0, 30.0);
  const double3 centre(15.0, 15.0, 15.0);

  const std::array<double3, 3> contacts{double3(-0.6, 0.0, 0.8), double3(0.3, 0.933, -0.2), double3(0.3, -0.933, -0.2)};
  std::vector<double3> positions;
  std::vector<double> bare;
  std::vector<double> inflatedRadii;
  for (const double3 &contact : contacts)
  {
    positions.push_back(centre + reach * ((1.0 / contact.length()) * contact));
    bare.push_back(atomRadius);
    inflatedRadii.push_back(reach);
  }
  BruteForceStructure structure = BruteForceStructure::create("corner", box, positions, bare);
  BruteForceStructure inflated = BruteForceStructure::create("corner inflated", box, positions, inflatedRadii);
  BruteForceSurfaceArea surface;  // no rolling probe positions: only the two wedged ones

  BruteForceSolventExcluded excluded =
      BruteForceSolventExcluded::compute(structure, inflated, surface, probe, 64, 256, 400000);
  ASSERT_EQ(excluded.numberOfCorners, 2uz);

  // The patch of the probe at `centre`: inside the three great circles, each oriented to hold the third
  // contact, and not within a probe radius of the mirror position. The mirror probe's patch is its reflection.
  std::array<double3, 3> unit;
  for (std::size_t i = 0; i < 3; ++i) unit[i] = (1.0 / contacts[i].length()) * contacts[i];
  std::array<double3, 3> edges;
  for (std::size_t i = 0; i < 3; ++i)
  {
    double3 normal = double3::cross(unit[i], unit[(i + 1) % 3]);
    normal = (1.0 / normal.length()) * normal;
    if (double3::dot(normal, unit[(i + 2) % 3]) < 0.0) normal = -1.0 * normal;
    edges[i] = normal;
  }
  double3 planeNormal = double3::cross(positions[1] - positions[0], positions[2] - positions[0]);
  planeNormal = (1.0 / planeNormal.length()) * planeNormal;
  double height = double3::dot(centre - positions[0], planeNormal);
  double3 mirror = centre - 2.0 * height * planeNormal;

  RandomNumber random{std::optional<std::size_t>(3)};
  const std::size_t draws = 4000000;
  std::size_t inside = 0;
  for (std::size_t draw = 0; draw < draws; ++draw)
  {
    double z = 2.0 * random.uniform() - 1.0;
    double ring = std::sqrt(std::max(0.0, 1.0 - z * z));
    double angle = 2.0 * std::numbers::pi * random.uniform();
    double3 direction(ring * std::cos(angle), ring * std::sin(angle), z);
    bool inPatch = std::ranges::all_of(edges, [&](const double3 &normal) { return double3::dot(normal, direction) >= 0.0; });
    if (!inPatch) continue;
    if ((centre + probe * direction - mirror).length() < probe - 1.0e-7) continue;
    ++inside;
  }
  const double expected = 2.0 * 4.0 * std::numbers::pi * probe * probe * static_cast<double>(inside) /
                          static_cast<double>(draws);

  EXPECT_NEAR(excluded.concaveArea, expected, 0.01 * expected);
}
