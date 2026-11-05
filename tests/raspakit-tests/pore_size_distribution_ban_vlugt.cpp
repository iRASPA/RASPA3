#include <gtest/gtest.h>

import std;

import int3;
import uint3;
import double3;
import double2;
import atom;
import forcefield;
import framework;
import simulationbox;
import pore_size_distribution_ban_vlugt;
import opencl;

// The pore size of Ban and Vlugt at a grid point is the radius of the largest sphere that holds the point
// and no framework atom, the atoms being hard spheres of half their Lennard-Jones size. The sphere is grown
// from every grid point out to its nearest atom and laid over every grid point within it, wrapping round the
// cell as the structure does. The same thing, done in double on the host over the same grid points, is the
// reference.

namespace
{

constexpr double edge = 10.0;
const uint3 gridSize{10, 10, 10};

bool noDevice()
{
  OpenCL::initialize();
  return !OpenCL::clContext.has_value();
}

// One oxygen near a corner of a cubic cell, off every mirror plane of the cell: an atom on a corner would
// let every sphere that has to wrap round the cell be matched by its mirror image that does not, and the
// wrapping would go untested.
const double3 atomFractional(0.15, 0.1, 0.05);

Framework oneAtom(const ForceField &forceField)
{
  std::size_t typeO = forceField.findPseudoAtom("O").value();
  std::vector<Atom> atoms{Atom(atomFractional, 0.0, 1.0, 0, static_cast<std::uint16_t>(typeO), 0, false, false)};
  return Framework(forceField, "one-oxygen", SimulationBox(edge, edge, edge), 1, atoms, int3(1, 1, 1));
}

double3 wrapped(double3 d)
{
  d.x -= std::rint(d.x);
  d.y -= std::rint(d.y);
  d.z -= std::rint(d.z);
  return d;
}

std::vector<double> reference(double atomRadius)
{
  const std::size_t n = static_cast<std::size_t>(gridSize.x) * gridSize.y * gridSize.z;
  std::vector<double3> points(n);
  std::vector<double> clearance(n);
  for (std::size_t k = 0; k < gridSize.z; ++k)
  {
    for (std::size_t j = 0; j < gridSize.y; ++j)
    {
      for (std::size_t i = 0; i < gridSize.x; ++i)
      {
        std::size_t v = (k * gridSize.y + j) * gridSize.x + i;
        points[v] = double3(static_cast<double>(i) / static_cast<double>(gridSize.x),
                            static_cast<double>(j) / static_cast<double>(gridSize.y),
                            static_cast<double>(k) / static_cast<double>(gridSize.z));
        clearance[v] = edge * wrapped(points[v] - atomFractional).length() - atomRadius;
      }
    }
  }

  std::vector<double> largest(n, 0.0);
  for (std::size_t centre = 0; centre < n; ++centre)
  {
    if (!(clearance[centre] > 0.0)) continue;
    for (std::size_t p = 0; p < n; ++p)
    {
      double distance = edge * wrapped(points[p] - points[centre]).length();
      if (distance <= clearance[centre]) largest[p] = std::max(largest[p], clearance[centre]);
    }
  }
  return largest;
}

}  // namespace


TEST(pore_size_distribution_ban_vlugt, the_largest_sphere_is_grown_from_the_atom_surface_and_wraps_round_the_cell)
{
  if (noDevice()) GTEST_SKIP() << "no OpenCL device";
  ForceField forceField = ForceField::makeZeoliteForceField(12.0, true, false, true);
  Framework framework = oneAtom(forceField);
  const double sigma = framework.atomUnitCellLennardJonesPotentialParameters(forceField)[0].y;
  ASSERT_GT(sigma, 0.0);

  BanVlugtPoreSizeDistribution psd(gridSize);
  psd.run(forceField, framework);
  std::filesystem::remove(framework.name + ".PSD_ban_vlugt.gpu.txt");

  std::vector<double> expected = reference(0.5 * sigma);
  ASSERT_EQ(psd.poreSizes.size(), expected.size());

  // The grid point across the cell from the atom holds a sphere wider than half the cell.
  const std::size_t farCorner = (5 * gridSize.y + 6) * gridSize.x + 6;
  EXPECT_GT(expected[farCorner], 0.5 * edge);
  for (std::size_t v = 0; v < expected.size(); ++v)
  {
    EXPECT_NEAR(static_cast<double>(psd.poreSizes[v]), expected[v], 1.0e-4 * std::max(1.0, expected[v]))
        << "voxel " << v;
  }
}
