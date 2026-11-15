#include <gtest/gtest.h>

import std;

import int3;
import double3;
import double3x3;
import simulationbox;
import randomnumbers;
import dctst_barrier;
import move_statistics;
import mc_moves_move_types;
import mc_moves_statistics;

// Regression tests for the geometry of the dcTST dividing surface and cage lattice, and for the bookkeeping of
// Monte Carlo moves that are undone after acceptance.

// An oblique normal: the reaction coordinate is the distance from the plane wrapped on the repeat distance of the
// box along the normal. A component-wise minimum image followed by a dot product jumped by O(L) where the vector to
// the barrier point crossed a face of the box.
TEST(DCTST_GEOMETRY_FIXES, oblique_projection_is_continuous_across_the_box_faces)
{
  const double L = 20.0;
  const SimulationBox box(L, L, L);
  const DcTSTBarrier barrier(box, double3(0.5, 0.5, 0.5), double3(1.0, 1.0, 0.0), L / std::sqrt(2.0));
  EXPECT_EQ(barrier.millerIndices, int3(1, 1, 0));
  EXPECT_NEAR(barrier.period, L / std::sqrt(2.0), 1e-12);

  const double q1 = barrier.projection(box, double3(19.99, 2.0, 10.0));
  const double q2 = barrier.projection(box, double3(20.01, 2.0, 10.0));
  EXPECT_NEAR(q2 - q1, 0.02 / std::sqrt(2.0), 1e-9);
}

// The projection of a uniform point of the box onto an oblique plane is uniform on the periodic plane (it was a
// triangular 'tent' density, twice the mean in the middle and zero at the ends, without a Hastings correction).
TEST(DCTST_GEOMETRY_FIXES, random_points_on_an_oblique_plane_are_uniform)
{
  const double L = 20.0;
  const SimulationBox box(L, L, L);
  const DcTSTBarrier barrier(box, double3(0.5, 0.5, 0.5), double3(1.0, 1.0, 0.0), L / std::sqrt(2.0));
  RandomNumber random(7);
  const double3 inPlane = double3::normalize(double3(1.0, -1.0, 0.0));
  const double period = L / std::sqrt(2.0);  // the in-plane repeat along (1, -1, 0)
  constexpr std::size_t numberOfBins = 10;
  constexpr std::size_t numberOfSamples = 200000;
  std::vector<double> histogram(numberOfBins, 0.0);
  for (std::size_t i = 0; i < numberOfSamples; ++i)
  {
    const double3 p = barrier.randomPointOnPlane(random, box);
    ASSERT_NEAR(barrier.projection(box, p), 0.0, 1e-9);
    double t = double3::dot(p - barrier.point, inPlane);
    t -= period * std::floor(t / period);
    histogram[std::min(numberOfBins - 1, static_cast<std::size_t>(t / period * numberOfBins))] += 1.0;
  }
  for (const double count : histogram)
  {
    EXPECT_NEAR(count / (static_cast<double>(numberOfSamples) / numberOfBins), 1.0, 0.05);
  }
}

// A normal that is not the normal of a lattice plane, or a box that does not repeat a whole number of lattice
// distances along it, has no periodic dividing surface.
TEST(DCTST_GEOMETRY_FIXES, barrier_refuses_a_non_periodic_geometry)
{
  const SimulationBox box(30.0, 30.0, 30.0);
  EXPECT_THROW(DcTSTBarrier(box, double3(0.5, 0.5, 0.5), double3(1.0, 0.0, 0.0), 7.0), std::invalid_argument);
  EXPECT_THROW(DcTSTBarrier(box, double3(0.5, 0.5, 0.5), double3(1.0, 0.37, 0.0), 10.0), std::invalid_argument);
  EXPECT_NO_THROW(DcTSTBarrier(box, double3(0.5, 0.5, 0.5), double3(1.0, 0.0, 0.0), 10.0));
}

// Hexagonal cell (gamma = 120 degrees) with the barrier normal to the (100) planes: the neighbor across the window
// is one cell along a, not the diagonal cell (the y-component of the shifted center rounded up to +1).
TEST(DCTST_GEOMETRY_FIXES, hexagonal_cage_neighbor_steps_along_one_axis)
{
  const double a = 10.0;
  const SimulationBox unitCell(a, a, 12.0, 0.5 * std::numbers::pi, 0.5 * std::numbers::pi, 2.0 * std::numbers::pi / 3.0);
  const SimulationBox supercell(4.0 * a, 4.0 * a, 12.0, 0.5 * std::numbers::pi, 0.5 * std::numbers::pi,
                                2.0 * std::numbers::pi / 3.0);
  const double3 normal = double3::cross(unitCell.cell[1], unitCell.cell[2]);
  const double spacing = a * std::sin(std::numbers::pi / 3.0);
  const DcTSTBarrier barrier(supercell, double3(0.5, 0.5, 0.5), normal, spacing);
  EXPECT_NEAR(barrier.period, 4.0 * spacing, 1e-9);

  const DcTSTCageLattice lattice(unitCell, int3(4, 4, 1), barrier.point + 0.5 * barrier.lambda * barrier.normal);
  EXPECT_FALSE(lattice.inconsistencyWith(supercell, barrier).has_value());
  EXPECT_NEAR(lattice.planeSpacing(0), spacing, 1e-9);
  for (std::size_t cage = 0; cage < lattice.numberOfCages(); ++cage)
  {
    const int3 from = lattice.indices(cage);
    const int3 to = lattice.indices(lattice.neighbor(cage, barrier.lambda * barrier.normal));
    EXPECT_EQ(to.y, from.y);
    EXPECT_EQ(to.z, from.z);
    EXPECT_EQ((to.x - from.x + 4) % 4, 1);
  }
}

// A unit cell with several cages along the normal (the conventional LTA cell holds 2 x 2 x 2 alpha cages): the
// lattice of one cage per cell merged two cages into one occupancy counter and is now reported as inconsistent;
// 'CagesPerUnitCell' subdivides the cell.
TEST(DCTST_GEOMETRY_FIXES, cages_per_unit_cell_subdivide_the_lattice)
{
  const double a = 24.0;
  const SimulationBox unitCell(a, a, a);
  const SimulationBox supercell(2.0 * a, 2.0 * a, 2.0 * a);
  const DcTSTBarrier barrier(supercell, double3(0.0, 0.25, 0.25), double3(1.0, 0.0, 0.0), 0.5 * a);
  const double3 origin = barrier.point + 0.5 * barrier.lambda * barrier.normal;

  const DcTSTCageLattice single(unitCell, int3(2, 2, 2), origin);
  EXPECT_TRUE(single.inconsistencyWith(supercell, barrier).has_value());

  const DcTSTCageLattice subdivided(unitCell, int3(2, 2, 2), origin, int3(2, 2, 2));
  EXPECT_FALSE(subdivided.inconsistencyWith(supercell, barrier).has_value());
  EXPECT_EQ(subdivided.numberOfCages(), 64uz);
  EXPECT_NEAR(subdivided.planeSpacing(0), 0.5 * a, 1e-12);
}

// A move accepted and then undone by a constraint (the cage occupancies of the dcTST drivers) counts as a rejected
// trial, not as an accepted one that inflates the step size.
TEST(DCTST_GEOMETRY_FIXES, revoked_acceptance_counts_as_a_rejected_trial)
{
  MCMoveStatistics statistics;
  const Move::Types move = Move::Types::Translation;
  statistics.addTrial(move, 1);
  const auto before = statistics.stats[std::to_underlying(move)];
  const SubMoveStatistics subBefore = statistics.subMoveStatistics(move);
  statistics.addConstructed(move, 1);
  statistics.addAccepted(move, 1);
  statistics.revokeAccepted(move, before, subBefore);

  const MoveStatistics<double3>& translation = std::get<MoveStatistics<double3>>(statistics.stats[std::to_underlying(move)]);
  EXPECT_EQ(translation.counts.y, 1.0);
  EXPECT_EQ(translation.accepted.y, 0.0);
  EXPECT_EQ(translation.totalAccepted.y, 0.0);
}

// The two cages on either side of the dividing plane a bead is at share the lateral indices of the bead. A step
// along the Cartesian normal changes the lateral fractional coordinates of an oblique cell: for a bead near a
// lateral cage wall it landed in the cages beside the column.
TEST(DCTST_GEOMETRY_FIXES, cages_across_the_plane_stay_in_the_column_of_the_bead)
{
  const double a = 10.0;
  const double gamma = 2.0 * std::numbers::pi / 3.0;
  const SimulationBox unitCell(a, a, 12.0, 0.5 * std::numbers::pi, 0.5 * std::numbers::pi, gamma);
  const SimulationBox supercell(4.0 * a, 4.0 * a, 12.0, 0.5 * std::numbers::pi, 0.5 * std::numbers::pi, gamma);
  const double3 normal = double3::cross(unitCell.cell[1], unitCell.cell[2]);
  const double spacing = a * std::sin(std::numbers::pi / 3.0);
  const DcTSTBarrier barrier(supercell, double3(0.5, 0.5, 0.5), normal, spacing);
  const DcTSTCageLattice lattice(unitCell, int3(4, 4, 1), barrier.point + 0.5 * barrier.lambda * barrier.normal);
  ASSERT_FALSE(lattice.inconsistencyWith(supercell, barrier).has_value());

  const std::size_t cage = lattice.index(int3(1, 2, 0));
  const std::size_t ahead = lattice.neighborAlongAxis(cage, 0, 1);
  const std::size_t behind = lattice.neighborAlongAxis(cage, 0, -1);

  // on the plane between 'cage' and the cage ahead, close to the lateral wall of the column
  const double3 onPlane = lattice.center(cage) + lattice.cell * double3(0.5, 0.45, 0.0);
  EXPECT_NEAR(barrier.reducedCoordinate(barrier.projection(supercell, onPlane)), 0.0, 1e-9);
  const auto pair = lattice.cagesAcrossNearestPlane(onPlane, barrier.normal);
  EXPECT_TRUE((pair == std::pair{cage, ahead}));
  // the quarter-cage steps along the normal leave the column there
  const double3 step = 0.25 * barrier.lambda * barrier.normal;
  EXPECT_NE(lattice.indices(lattice.cageOf(onPlane + step)).y, 2);

  // inside the cage: the nearest plane is the one on the side the bead is on
  EXPECT_TRUE((lattice.cagesAcrossNearestPlane(lattice.center(cage) + lattice.cell * double3(0.4, 0.45, 0.0),
                                               barrier.normal) == std::pair{cage, ahead}));
  EXPECT_TRUE((lattice.cagesAcrossNearestPlane(lattice.center(cage) + lattice.cell * double3(-0.4, -0.45, 0.0),
                                               barrier.normal) == std::pair{behind, cage}));
  // against the normal the roles of the two cages exchange
  EXPECT_TRUE((lattice.cagesAcrossNearestPlane(onPlane, -1.0 * barrier.normal) == std::pair{ahead, cage}));
}

// A plane slightly inclined to the cage planes passes a cosine test ((2 1 0) in a cell with b >> a here, |cos| =
// 0.99955; (12 1 0) with b = 2a is the same case) but cuts across the cages beside the column: only a plane with a
// single non-zero Miller index is a cage boundary.
TEST(DCTST_GEOMETRY_FIXES, an_inclined_plane_is_not_a_cage_boundary)
{
  const SimulationBox box(6.0, 100.0, 6.0);
  const double3 reciprocal(2.0 / 6.0, 1.0 / 100.0, 0.0);
  const DcTSTBarrier barrier(box, double3(0.5, 0.5, 0.5), reciprocal, 1.0 / reciprocal.length());
  ASSERT_EQ(barrier.millerIndices, int3(2, 1, 0));

  const DcTSTCageLattice lattice(box, int3(1, 1, 1), barrier.point + 0.5 * barrier.lambda * barrier.normal,
                                 int3(2, 1, 1));
  ASSERT_GT(lattice.axisAlongDirection(barrier.normal).cosine, 0.999);
  ASSERT_LT(std::abs(lattice.planeSpacing(0) - barrier.lambda), 1e-3 * barrier.lambda);
  const std::optional<std::string> inconsistency = lattice.inconsistencyWith(box, barrier);
  ASSERT_TRUE(inconsistency.has_value());
  EXPECT_NE(inconsistency->find("(2 1 0)"), std::string::npos) << *inconsistency;

  // the plane (1 0 0) of the same box is one
  const DcTSTBarrier axial(box, double3(0.5, 0.5, 0.5), double3(1.0, 0.0, 0.0), 3.0);
  const DcTSTCageLattice axialLattice(box, int3(1, 1, 1), axial.point + 0.5 * axial.lambda * axial.normal,
                                      int3(2, 1, 1));
  EXPECT_FALSE(axialLattice.inconsistencyWith(box, axial).has_value());
}
