#include <gtest/gtest.h>

import std;

import uint3;
import double3;
import forcefield;
import simulationbox;
import interpolation_energy_grid;

namespace
{
// A smooth periodic test function of the fractional coordinate.
double periodicFunction(const double3& s)
{
  return std::cos(2.0 * std::numbers::pi * s.x) + 0.5 * std::sin(2.0 * std::numbers::pi * s.y) +
         0.25 * std::cos(4.0 * std::numbers::pi * s.z) + 0.1 * std::sin(2.0 * std::numbers::pi * (s.x + s.y + s.z));
}

// Tabulates the function on the polynomial (energy-only) grid, whose last plane duplicates the first one
// (numberOfGridPoints = numberOfCells + 1).
InterpolationEnergyGrid makePolynomialGrid(const SimulationBox& box, std::size_t cells)
{
  InterpolationEnergyGrid grid(box, double3(0.0, 0.0, 0.0), uint3(cells + 1, cells + 1, cells + 1),
                               ForceFieldSettings::InterpolationScheme::Polynomial);
  const std::size_t points = cells + 1;
  for (std::size_t k = 0; k < points; ++k)
  {
    for (std::size_t j = 0; j < points; ++j)
    {
      for (std::size_t i = 0; i < points; ++i)
      {
        const double3 s(static_cast<double>(i) / static_cast<double>(cells),
                        static_cast<double>(j) / static_cast<double>(cells),
                        static_cast<double>(k) / static_cast<double>(cells));
        grid.data[i + points * (j + points * k)] = periodicFunction(s);
      }
    }
  }
  return grid;
}
}  // namespace

// The stencil of the polynomial scheme reaches across the cell boundaries. The grid stores numberOfCells + 1
// planes per axis with the last plane a copy of the first (s = 0 and s = 1), so a periodic wrap of a stencil
// index must be modulo numberOfCells: wrapping modulo numberOfGridPoints shifts the stencil values by one
// plane on either side of the unit cell. Positions in the first and last cell must be interpolated as
// accurately as positions in the interior.
TEST(interpolation_energy_grid, polynomial_stencil_wraps_periodically)
{
  const SimulationBox box(10.0, 10.0, 10.0);
  const std::size_t cells = 24;
  const InterpolationEnergyGrid grid = makePolynomialGrid(box, cells);

  double maximumInteriorError = 0.0;
  for (double t : {0.37, 0.5, 0.63})
  {
    const double3 s(t, t + 0.05, t - 0.05);
    maximumInteriorError =
        std::max(maximumInteriorError, std::abs(grid.interpolate(box.cell * s) - periodicFunction(s)));
  }
  EXPECT_LT(maximumInteriorError, 1e-4);

  // the boundary cells on both sides, along every axis
  for (double edge : {0.003, 0.02, 0.98, 0.997})
  {
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      double3 s(0.43, 0.51, 0.58);
      s[axis] = edge;
      const double error = std::abs(grid.interpolate(box.cell * s) - periodicFunction(s));
      EXPECT_LT(error, 10.0 * maximumInteriorError + 1e-6) << "axis " << axis << " s = " << edge;
    }
  }

  // the two copies of the boundary plane interpolate to the same value
  const double3 sLow(0.0, 0.3, 0.7);
  const double3 sHigh(1.0 - 1e-12, 0.3, 0.7);
  EXPECT_NEAR(grid.interpolate(box.cell * sLow), grid.interpolate(box.cell * sHigh), 1e-8);
}
