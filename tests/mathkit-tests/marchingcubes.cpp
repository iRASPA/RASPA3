#include <gtest/gtest.h>

import std;

import marching_cubes;

TEST(marchingcubes, Test_a_second_run_does_not_append_to_the_first)
{
  // the distance from the centre of a 12^3 grid, so isosurfaces are spheres
  constexpr std::size_t n = 12;
  MarchingCubes cubes(static_cast<int>(n), static_cast<int>(n), static_cast<int>(n));
  cubes.init_all();
  for (std::size_t k = 0; k < n; ++k)
    for (std::size_t j = 0; j < n; ++j)
      for (std::size_t i = 0; i < n; ++i)
      {
        double x = static_cast<double>(i) - 5.5;
        double y = static_cast<double>(j) - 5.5;
        double z = static_cast<double>(k) - 5.5;
        cubes.set_data(std::sqrt(x * x + y * y + z * z), i, j, k);
      }

  cubes.run(3.0);
  std::size_t vertices = cubes.nverts();
  std::size_t triangles = cubes.ntrigs();
  ASSERT_GT(triangles, 0uz);

  cubes.run(3.0);
  EXPECT_EQ(cubes.nverts(), vertices);
  EXPECT_EQ(cubes.ntrigs(), triangles);

  // a smaller sphere on the same instance has fewer triangles, and every triangle refers to a current vertex
  cubes.run(2.0);
  EXPECT_LT(cubes.ntrigs(), triangles);
  for (std::size_t t = 0; t < cubes.ntrigs(); ++t)
  {
    EXPECT_LT(static_cast<std::size_t>(cubes.triangles()[t].v1), cubes.nverts());
    EXPECT_LT(static_cast<std::size_t>(cubes.triangles()[t].v2), cubes.nverts());
    EXPECT_LT(static_cast<std::size_t>(cubes.triangles()[t].v3), cubes.nverts());
  }
}
