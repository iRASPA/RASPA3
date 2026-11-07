#include <gtest/gtest.h>

import std;

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
