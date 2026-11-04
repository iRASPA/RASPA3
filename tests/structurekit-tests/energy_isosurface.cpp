#include <gtest/gtest.h>

import std;

import uint3;
import double3;
import energy_isosurface;

// The processor extractor's vertex normals. The field is periodic and the extractor repeats the first plane
// at the far end, so the normals on the two planes that meet there have to be central differences read round
// the join: a one-sided difference stopping at the join tilts every normal on a cell face.

namespace
{

constexpr std::size_t n = 32;
const uint3 gridSize{n, n, n};

// cos(2πx) + cos(2πy)/2, which has a maximum on the plane x = 0 and so no gradient along x there. Its level
// set at one crosses that plane, so there are vertices on it to ask about.
std::vector<float> ridge()
{
  std::vector<float> field(n * n * n);
  for (std::size_t k = 0; k < n; ++k)
  {
    for (std::size_t j = 0; j < n; ++j)
    {
      for (std::size_t i = 0; i < n; ++i)
      {
        double x = static_cast<double>(i) / static_cast<double>(n);
        double y = static_cast<double>(j) / static_cast<double>(n);
        field[(k * n + j) * n + i] = static_cast<float>(std::cos(2.0 * std::numbers::pi * x) +
                                                        0.5 * std::cos(2.0 * std::numbers::pi * y));
      }
    }
  }
  return field;
}

}  // namespace


TEST(energy_isosurface, normals_on_the_cell_face_read_round_the_join)
{
  std::vector<float> field = ridge();
  std::vector<double3> gradients;
  std::vector<double3> corners = EnergyIsosurface::trianglesOfIsosurface(field, gridSize, 1.0, &gradients);
  ASSERT_EQ(gradients.size(), corners.size());

  std::size_t onTheFace = 0;
  double worstTilt = 0.0;
  for (std::size_t v = 0; v < corners.size(); ++v)
  {
    const double x = corners[v].x;
    if (!(x < 1.0e-9 || x > 1.0 - 1.0e-9)) continue;
    ++onTheFace;
    const double length = gradients[v].length();
    ASSERT_GT(length, 0.0);
    worstTilt = std::max(worstTilt, std::abs(gradients[v].x) / length);
  }

  ASSERT_GT(onTheFace, 0u);
  // The field is even in x about the face, so the exact gradient there has no x at all.
  EXPECT_LT(worstTilt, 1.0e-3);
}
