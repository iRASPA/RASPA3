#include <gtest/gtest.h>

import std;

import double3;
import double3x3;
import simulationbox;
import property_simulationbox;
import property_block_average;

TEST(SIMULATIONBOX_BLOCK_AVERAGE, cell_error_is_element_wise_for_a_triclinic_cell)
{
  // The block-average error squares the deviation of each block from the mean and takes the square root of the
  // summed squares, so the product of two boxes must be element-wise; the matrix product mixes the cell elements
  // of a triclinic cell (and the square root of a negative mixed element is NaN).
  const std::size_t numberOfBlocks = 5;
  PropertySimulationBox boxAverage(numberOfBlocks);
  std::array<BlockAverage<double>, 9> elementAverages;
  elementAverages.fill(BlockAverage<double>(numberOfBlocks));

  for (std::size_t block = 0; block < numberOfBlocks; ++block)
  {
    const double shift = 0.3 * static_cast<double>(block * block) - 0.7 * static_cast<double>(block);
    const SimulationBox box(20.0 + shift, 22.0 - 0.5 * shift, 25.0 + 0.2 * shift,
                            (80.0 + shift) * std::numbers::pi / 180.0, (100.0 - shift) * std::numbers::pi / 180.0,
                            (110.0 + 0.5 * shift) * std::numbers::pi / 180.0);
    boxAverage.addSample(block, box, 1.0);
    const std::array<double, 9> elements{box.cell.ax, box.cell.ay, box.cell.az, box.cell.bx, box.cell.by,
                                         box.cell.bz, box.cell.cx, box.cell.cy, box.cell.cz};
    for (std::size_t i = 0; i < 9; ++i)
    {
      elementAverages[i].addSample(block, elements[i], 1.0);
    }
  }

  const auto [average, error] = boxAverage.average();
  const std::array<double, 9> cellErrors{error.cell.ax, error.cell.ay, error.cell.az, error.cell.bx, error.cell.by,
                                         error.cell.bz, error.cell.cx, error.cell.cy, error.cell.cz};
  for (std::size_t i = 0; i < 9; ++i)
  {
    const double expected = elementAverages[i].average().second;
    EXPECT_NEAR(cellErrors[i], expected, 1.0e-10 * std::max(1.0, expected)) << "cell element " << i;
  }
}

// A cell with right angles that is rotated away from the axes is not rectangular: the rectangular minimum image
// reads only the diagonal of the cell, and with a rotated cell it moved a lattice vector to somewhere else than zero.
TEST(SIMULATIONBOX_FIXES, rotated_orthogonal_cell_is_triclinic_and_wraps_lattice_vectors_to_zero)
{
  const double angle = 0.3;
  const double c = std::cos(angle);
  const double s = std::sin(angle);
  const double3 a(20.0 * c, 20.0 * s, 0.0);
  const double3 b(-22.0 * s, 22.0 * c, 0.0);
  const double3 cVector(0.0, 0.0, 24.0);
  const SimulationBox rotated(double3x3(a, b, cVector));
  EXPECT_EQ(rotated.type, SimulationBox::Type::Triclinic);
  for (const double3& lattice : {a, b, cVector, a + b, a - b + cVector})
  {
    const double3 wrapped = rotated.applyPeriodicBoundaryConditions(lattice);
    EXPECT_NEAR(wrapped.length(), 0.0, 1.0e-9);
  }

  const SimulationBox alongTheAxes(
      double3x3(double3(20.0, 0.0, 0.0), double3(0.0, 22.0, 0.0), double3(0.0, 0.0, 24.0)));
  EXPECT_EQ(alongTheAxes.type, SimulationBox::Type::Rectangular);
}
