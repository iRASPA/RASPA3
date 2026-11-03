#include <gtest/gtest.h>

import std;

import int3;
import double3;
import unit_cell;
import brute_force_structure;
import brute_force_voxels;
import brute_force_pore_spikes;

// The pore-size spikes read off the bare void: the maximal spheres the walk uphill finds, and which families
// of them are reported.

namespace
{

// A tube of rings of six atoms about an axis through the middle of the cell. `radiusAt` gives the ring's
// radius a distance `t` along the axis, for t from `from` to `to` in steps of `spacing`; the rings' atoms all
// sit at the same six azimuths, so the axis is the middle of the tube everywhere along it.
struct Tube
{
  UnitCell box;
  double3 middle;
  double3 axis;
  std::vector<double3> positions;
  std::vector<double> radii;

  Tube(double side, double3 along) : box(side, side, side), middle(0.5 * side, 0.5 * side, 0.5 * side)
  {
    this->axis = (1.0 / along.length()) * along;
  }

  void rings(double from, double to, double spacing, double atomRadius, auto &&radiusAt)
  {
    double3 u = double3::cross(this->axis, double3(0.0, 0.0, 1.0));
    if (u.length() < 0.5) u = double3::cross(this->axis, double3(1.0, 0.0, 0.0));
    u = (1.0 / u.length()) * u;
    double3 v = double3::cross(this->axis, u);
    for (double t = from; t <= to + 1.0e-9; t += spacing)
    {
      double radius = radiusAt(t);
      for (int k = 0; k < 6; ++k)
      {
        double azimuth = k * std::numbers::pi / 3.0;
        this->positions.push_back(this->middle + t * this->axis +
                                  (radius * std::cos(azimuth)) * u + (radius * std::sin(azimuth)) * v);
        this->radii.push_back(atomRadius);
      }
    }
  }

  BruteForceStructure structure() const
  {
    return BruteForceStructure::create("tube", this->box, this->positions, this->radii);
  }
};

// The most room anywhere on the axis between `from` and `to`, by a scan refined three times over.
double mostRoomOnTheAxis(const BruteForceStructure &structure, const Tube &tube, double from, double to)
{
  double bestT = from;
  double best = -1.0;
  for (double step : {1.0e-3, 1.0e-6, 1.0e-9})
  {
    double low = (best < 0.0) ? from : bestT - 1000.0 * step;
    double high = (best < 0.0) ? to : bestT + 1000.0 * step;
    for (double t = low; t <= high; t += step)
    {
      double room = structure.clearance(tube.middle + t * tube.axis);
      if (room > best)
      {
        best = room;
        bestT = t;
      }
    }
  }
  return best;
}

// A sealed cage of `count` spheres on a shell about `centre`, placed on the Fibonacci spiral.
std::vector<double3> fibonacciShell(const double3 &centre, double shellRadius, std::size_t count)
{
  std::vector<double3> positions;
  const double golden = std::numbers::pi * (1.0 + std::sqrt(5.0));
  for (std::size_t i = 0; i < count; ++i)
  {
    double fraction = (static_cast<double>(i) + 0.5) / static_cast<double>(count);
    double cosinePolar = 1.0 - 2.0 * fraction;
    double sinePolar = std::sqrt(std::max(0.0, 1.0 - cosinePolar * cosinePolar));
    double azimuth = golden * (static_cast<double>(i) + 0.5);
    positions.push_back(centre + double3(shellRadius * sinePolar * std::cos(azimuth),
                                         shellRadius * sinePolar * std::sin(azimuth), shellRadius * cosinePolar));
  }
  return positions;
}

}  // namespace


TEST(brute_force_pore_spikes, the_walk_climbs_a_crease_that_runs_between_the_lattice_directions)
{
  // A tube that widens slowly along an axis pointing nowhere in particular, closed off by a narrower ring at
  // its wide end. The rings are close enough that the room on the axis rises all the way to the closing ring
  // and nowhere in between, so the only maximal sphere in the tube sits where the widening meets the closing.
  // Along the axis the room rises by five parts in a hundred, and any step off the axis costs its full length
  // in room, so the only way up is along the axis to within a few degrees. No lattice direction lies that
  // close to it, and a walk confined to the lattice directions stops on the axis wherever it first reaches it.
  Tube tube(20.0, double3(1.0, 2.0, 3.0));
  tube.rings(-6.0, 6.0, 0.5, 1.5, [](double t) { return 2.5 + 0.05 * (t + 6.0); });
  tube.rings(6.5, 7.0, 0.5, 1.5, [](double) { return 2.0; });
  BruteForceStructure structure = tube.structure();
  BruteForceVoxels voxels = BruteForceVoxels::build(structure, 0.5);

  const double room = mostRoomOnTheAxis(structure, tube, -6.0, 7.0);
  ASSERT_GT(room, 1.0);

  BruteForcePoreSpikes spikes = BruteForcePoreSpikes::compute(structure, voxels, 100000, 8);

  bool found = false;
  for (const BruteForceSpikeFamily &family : spikes.families)
  {
    if (std::abs(family.diameter - 2.0 * room) <= 1.0e-6) found = true;
  }
  EXPECT_TRUE(found) << "no family has the diameter " << 2.0 * room << " of the sphere that fits at the top of the tube";
}


TEST(brute_force_pore_spikes, families_beyond_the_dozen_widest_are_counted_and_the_widest_above_the_floor_reported)
{
  // A tube of rings two steps apart, each a little wider than the last: the room on the axis bulges between
  // every pair of rings, and each bulge is a maximal sphere of its own, so the tube alone gives a dozen
  // families of distinct diameters. Beside it, eight sealed cages, all alike and narrower than any bulge, which
  // are one family. There are more families than are weighed in full without a first look, so the first look
  // is taken; every family is still counted, and the ones reported are the widest with any void to their name.
  Tube tube(30.0, double3(1.0, 0.0, 0.0));
  tube.rings(-11.0, 11.0, 2.0, 1.5, [](double t) { return 3.1 + 0.05 * (t + 11.0); });
  for (double3 corner : {double3(4.0, 4.0, 4.0), double3(26.0, 4.0, 4.0), double3(4.0, 26.0, 4.0),
                         double3(26.0, 26.0, 4.0), double3(4.0, 4.0, 26.0), double3(26.0, 4.0, 26.0),
                         double3(4.0, 26.0, 26.0), double3(26.0, 26.0, 26.0)})
  {
    for (const double3 &position : fibonacciShell(corner, 3.0, 20))
    {
      tube.positions.push_back(position);
      tube.radii.push_back(1.5);
    }
  }
  BruteForceStructure structure = tube.structure();
  BruteForceVoxels voxels = BruteForceVoxels::build(structure, 0.75);

  BruteForcePoreSpikes spikes = BruteForcePoreSpikes::compute(structure, voxels, 400000, 3);
  ASSERT_GT(spikes.numberOfFamilies, 12uz) << "the bulges alone were to give a dozen families";
  ASSERT_EQ(spikes.families.size(), 3uz);
  for (std::size_t f = 0; f + 1 < spikes.families.size(); ++f)
  {
    EXPECT_GT(spikes.families[f].diameter, spikes.families[f + 1].diameter);
  }
  for (const BruteForceSpikeFamily &family : spikes.families) EXPECT_GT(family.weight, 1.0e-5);

  // Asked for every family, the cages come too: eight balls of radius one and a half, one family of diameter
  // three, the narrowest of all.
  BruteForcePoreSpikes all = BruteForcePoreSpikes::compute(structure, voxels, 400000, 32);
  ASSERT_EQ(all.families.size(), all.numberOfFamilies);
  const BruteForceSpikeFamily &cages = all.families.back();
  EXPECT_NEAR(cages.diameter, 3.0, 1.0e-6);
  EXPECT_EQ(cages.centres.size(), 8uz);
  EXPECT_GT(cages.weight, 1.0e-3);
}
