#include <gtest/gtest.h>

import std;

import uint3;
import double3;
import unit_cell;
import crystal;
import units;
import pair_interactions;
import energy_opencl_void_fraction;
import energy_opencl_surface_area;
import opencl;

// The single-site device routes for the void fraction and the surface area. The device works over a range
// padded out to whole work-groups, so the grid is followed by points that are not on it; those must weigh
// nothing in the average. And a grid point can lie on an atom's centre, where the pair's separation is zero;
// an atom with no strength of its own has to count for nothing there rather than poison the whole field.

namespace
{

constexpr double edge = 10.0;
constexpr double sigma = 3.0;
constexpr double epsilon = 100.0;
constexpr double cutOff = 4.0;
constexpr double temperature = 298.0;

// Three types: the probe, an atom with strength, and one without.
PairInteractions probeAndTwoAtoms()
{
  PairInteractions interactions;
  interactions.numberOfTypes = 3;
  interactions.names = {"He", "O", "H"};
  interactions.charges = {0.0, 0.0, 0.0};
  interactions.parameters = {PairParameters{sigma, epsilon, 0.0}, PairParameters{sigma, epsilon, 0.0},
                             PairParameters{2.5, 0.0, 0.0},       PairParameters{sigma, epsilon, 0.0},
                             PairParameters{sigma, epsilon, 0.0}, PairParameters{2.5, 0.0, 0.0},
                             PairParameters{2.5, 0.0, 0.0},       PairParameters{2.5, 0.0, 0.0},
                             PairParameters{2.5, 0.0, 0.0}};
  interactions.cutOffVDW = cutOff;
  interactions.cutOffCoulomb = cutOff;
  return interactions;
}

// One atom with strength in the middle of the cell. The origin, which is the first grid point and the
// position every padding point is given, is a cutoff and more away from it.
Crystal atomInTheMiddle(const std::string &name)
{
  Crystal framework;
  framework.name = name;
  framework.unitCell = UnitCell(edge, edge, edge);
  framework.mass = 16.0;
  framework.atoms.push_back(CrystalAtom{double3(0.5 * edge, 0.5 * edge, 0.5 * edge), 0.0, 1});
  framework.fractionalPositions.push_back(double3(0.5, 0.5, 0.5));
  return framework;
}

// The same, with a strengthless atom on the origin, exactly on the first grid point.
Crystal withStrengthlessAtomOnTheOrigin(const std::string &name)
{
  Crystal framework = atomInTheMiddle(name);
  framework.atoms.push_back(CrystalAtom{double3(0.0, 0.0, 0.0), 0.0, 2});
  framework.fractionalPositions.push_back(double3(0.0, 0.0, 0.0));
  return framework;
}

// The average of exp(-U/kT) over the grid, summed in double on the host from the one atom with strength.
double expectedVoidFraction(const uint3 &gridSize)
{
  const double beta = 1.0 / (Units::KB * temperature);
  double sum = 0.0;
  for (std::size_t k = 0; k < gridSize.z; ++k)
  {
    for (std::size_t j = 0; j < gridSize.y; ++j)
    {
      for (std::size_t i = 0; i < gridSize.x; ++i)
      {
        double3 ds(static_cast<double>(i) / static_cast<double>(gridSize.x) - 0.5,
                   static_cast<double>(j) / static_cast<double>(gridSize.y) - 0.5,
                   static_cast<double>(k) / static_cast<double>(gridSize.z) - 0.5);
        ds.x -= std::rint(ds.x);
        ds.y -= std::rint(ds.y);
        ds.z -= std::rint(ds.z);
        double rr = edge * edge * double3::dot(ds, ds);
        double energy = 0.0;
        if (rr < cutOff * cutOff)
        {
          double ratio3 = std::pow(sigma * sigma / rr, 3.0);
          energy = std::min(4.0 * epsilon * ratio3 * (ratio3 - 1.0), 1.0e7);
        }
        sum += std::exp(-beta * energy);
      }
    }
  }
  return sum / (static_cast<double>(gridSize.x) * static_cast<double>(gridSize.y) * static_cast<double>(gridSize.z));
}

bool noDevice()
{
  OpenCL::initialize();
  return !OpenCL::clContext.has_value();
}

}  // namespace


// An odd number of grid points is not a multiple of any even work-group size, so the range is padded.
TEST(energy_legacy_void_fraction, padding_points_weigh_nothing_in_the_average)
{
  if (noDevice()) GTEST_SKIP() << "no OpenCL device";
  const uint3 gridSize{9, 9, 9};
  PairInteractions interactions = probeAndTwoAtoms();
  Crystal framework = atomInTheMiddle("legacy-void-fraction-padding");

  EnergyOpenCLVoidFraction vf;
  vf.run(interactions, framework, "He", gridSize, temperature);
  std::filesystem::remove(framework.name + ".energy.vf.gpu.txt");

  const double expected = expectedVoidFraction(gridSize);
  EXPECT_NEAR(vf.voidFraction, expected, 1.0e-4 * expected);
}


TEST(energy_legacy_void_fraction, a_strengthless_atom_on_a_grid_point_counts_for_nothing)
{
  if (noDevice()) GTEST_SKIP() << "no OpenCL device";
  const uint3 gridSize{8, 8, 8};
  PairInteractions interactions = probeAndTwoAtoms();
  Crystal framework = withStrengthlessAtomOnTheOrigin("legacy-void-fraction-origin");

  EnergyOpenCLVoidFraction vf;
  vf.run(interactions, framework, "He", gridSize, temperature);
  std::filesystem::remove(framework.name + ".energy.vf.gpu.txt");

  const double expected = expectedVoidFraction(gridSize);
  ASSERT_TRUE(std::isfinite(vf.voidFraction));
  EXPECT_NEAR(vf.voidFraction, expected, 1.0e-4 * expected);
}


TEST(energy_legacy_void_fraction, a_strengthless_atom_on_a_grid_point_leaves_the_surface_alone)
{
  if (noDevice()) GTEST_SKIP() << "no OpenCL device";
  const uint3 gridSize{32, 32, 32};
  // On the repulsive side of the well, so that the surface is one sphere about the atom with strength and
  // every point beyond the cutoff, the origin among them, lies below it.
  const double isoValue = 10.0;
  PairInteractions interactions = probeAndTwoAtoms();
  Crystal alone = atomInTheMiddle("legacy-surface-area-alone");
  Crystal withStrengthless = withStrengthlessAtomOnTheOrigin("legacy-surface-area-origin");

  EnergyOpenCLSurfaceArea without;
  without.run(interactions, alone, isoValue, "He", gridSize);
  std::filesystem::remove(alone.name + ".energy.sa.gpu.txt");

  EnergyOpenCLSurfaceArea with;
  with.run(interactions, withStrengthless, isoValue, "He", gridSize);
  std::filesystem::remove(withStrengthless.name + ".energy.sa.gpu.txt");

  ASSERT_TRUE(std::isfinite(without.area));
  ASSERT_GT(without.area, 0.0);
  EXPECT_TRUE(std::isfinite(with.area));
  EXPECT_DOUBLE_EQ(with.area, without.area);
}
