#include <gtest/gtest.h>

import std;

import bool3;
import double3;
import double4;
import float4;
import double3x3;
import double4x4;
import int3x3;
import matrix;
import simd_quatd;
import ringmatrix;
import cubic;
import randomnumbers;
import marching_cubes;

TEST(MATHKIT_FIXES, normalise_divides_by_the_four_dimensional_length)
{
  double4 d(1.0, 2.0, 2.0, 4.0);
  d.normalise();
  EXPECT_NEAR(d.x * d.x + d.y * d.y + d.z * d.z + d.w * d.w, 1.0, 1e-12);
  EXPECT_NEAR(d.x, 0.2, 1e-12);

  float4 f(1.0f, 2.0f, 2.0f, 4.0f);
  f.normalise();
  EXPECT_NEAR(f.x * f.x + f.y * f.y + f.z * f.z + f.w * f.w, 1.0f, 1e-6f);
}

TEST(MATHKIT_FIXES, float4_subtraction_subtracts_every_component)
{
  float4 a(5.0f, 6.0f, 7.0f, 8.0f);
  a -= float4(1.0f, 2.0f, 3.0f, 4.0f);
  EXPECT_EQ(a.w, 4.0f);
  EXPECT_EQ(a.x, 4.0f);
}

TEST(MATHKIT_FIXES, double4x4_transpose_swaps_every_element)
{
  double4x4 m(11.0, 21.0, 31.0, 41.0, 12.0, 22.0, 32.0, 42.0, 13.0, 23.0, 33.0, 43.0, 14.0, 24.0, 34.0, 44.0);
  const double4x4 t = m.transpose();
  EXPECT_EQ(t.m32, 23.0);
  EXPECT_EQ(t.m23, 32.0);
  EXPECT_EQ(t.m42, 24.0);
  EXPECT_EQ(t.m24, 42.0);
}

TEST(MATHKIT_FIXES, symmetric_eigen_system_pairs_each_eigenvalue_with_its_eigenvector)
{
  double3x3 tensor(double3(4.0, 1.0, 0.5), double3(1.0, 3.0, 0.2), double3(0.5, 0.2, 1.0));
  double3 eigenvalues{};
  double3x3 eigenvectors{};
  double3x3 copy = tensor;
  copy.EigenSystemSymmetric(eigenvalues, eigenvectors);

  EXPECT_NEAR(eigenvectors.determinant(), 1.0, 1e-12);
  for (std::size_t i = 0; i < 3; ++i)
  {
    const double3 v = eigenvectors[i];
    const double3 Av = tensor * v;
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
      EXPECT_NEAR(Av[axis], eigenvalues[i] * v[axis], 1e-10) << "eigenpair " << i;
    }
  }
}

TEST(MATHKIT_FIXES, euler_angles_return_the_attitude_angle)
{
  // a rotation by 0.6 about z: the attitude (second) angle is 0.6
  const double angle = 0.6;
  simd_quatd q(0.0, 0.0, std::sin(0.5 * angle), std::cos(0.5 * angle));
  EXPECT_NEAR(q.EulerAngles().y, angle, 1e-12);
}

TEST(MATHKIT_FIXES, cubic_returns_a_double_root_twice)
{
  // (x - 1)^2 (x + 2) = x^3 - 3x + 2: discriminant exactly zero
  double coefficients[4] = {2.0, -3.0, 0.0, 1.0};
  double roots[3] = {0.0, 0.0, 0.0};
  int numberOfRoots = 0;
  cubic(coefficients, roots, &numberOfRoots);
  ASSERT_EQ(numberOfRoots, 3);
  EXPECT_NEAR(roots[0], -2.0, 1e-12);
  EXPECT_NEAR(roots[1], 1.0, 1e-12);
  EXPECT_NEAR(roots[2], 1.0, 1e-12);
}

TEST(MATHKIT_FIXES, cubic_newton_step_leaves_a_double_root_in_place)
{
  // (x - 5)^2 (x - 2) = x^3 - 12x^2 + 45x - 50: the trigonometric form lands a few ulps from the double root,
  // where the derivative is rounding noise; the Newton step used to jump that root from 5 to 7
  double coefficients[4] = {-50.0, 45.0, -12.0, 1.0};
  double roots[3] = {0.0, 0.0, 0.0};
  int numberOfRoots = 0;
  cubic(coefficients, roots, &numberOfRoots);
  ASSERT_EQ(numberOfRoots, 3);
  EXPECT_NEAR(roots[0], 2.0, 1e-12);
  EXPECT_NEAR(roots[1], 5.0, 1e-12);
  EXPECT_NEAR(roots[2], 5.0, 1e-12);
  EXPECT_LE(roots[0], roots[1]);
  EXPECT_LE(roots[1], roots[2]);
}

TEST(MATHKIT_FIXES, cubic_finds_a_double_root_whose_discriminant_rounds_positive)
{
  // (x - 3)^2 (x + 1) = x^3 - 5x^2 + 3x + 9: A[2] / 3 rounds, the discriminant comes out +1e-14 and the
  // solver used to return the single root -1
  double coefficients[4] = {9.0, 3.0, -5.0, 1.0};
  double roots[3] = {0.0, 0.0, 0.0};
  int numberOfRoots = 0;
  cubic(coefficients, roots, &numberOfRoots);
  ASSERT_EQ(numberOfRoots, 3);
  EXPECT_NEAR(roots[0], -1.0, 1e-12);
  EXPECT_NEAR(roots[1], 3.0, 1e-12);
  EXPECT_NEAR(roots[2], 3.0, 1e-12);
}

TEST(MATHKIT_FIXES, cubic_returns_a_triple_root_three_times)
{
  // (x - 1)^3 and (3x - 1)^3, the latter with non-dyadic root and rounded intermediates
  double dyadic[4] = {-1.0, 3.0, -3.0, 1.0};
  double roots[3] = {0.0, 0.0, 0.0};
  int numberOfRoots = 0;
  cubic(dyadic, roots, &numberOfRoots);
  ASSERT_EQ(numberOfRoots, 3);
  for (std::size_t i = 0; i < 3; ++i) EXPECT_DOUBLE_EQ(roots[i], 1.0) << "root " << i;

  double nonDyadic[4] = {-1.0, 9.0, -27.0, 27.0};
  cubic(nonDyadic, roots, &numberOfRoots);
  ASSERT_EQ(numberOfRoots, 3);
  for (std::size_t i = 0; i < 3; ++i) EXPECT_NEAR(roots[i], 1.0 / 3.0, 1e-15) << "root " << i;
}

TEST(MATHKIT_FIXES, cubic_keeps_three_distinct_and_single_real_roots)
{
  // (x + 2)(x - 1)(x - 4) = x^3 - 3x^2 - 6x + 8
  double distinct[4] = {8.0, -6.0, -3.0, 1.0};
  double roots[3] = {0.0, 0.0, 0.0};
  int numberOfRoots = 0;
  cubic(distinct, roots, &numberOfRoots);
  ASSERT_EQ(numberOfRoots, 3);
  EXPECT_NEAR(roots[0], -2.0, 1e-12);
  EXPECT_NEAR(roots[1], 1.0, 1e-12);
  EXPECT_NEAR(roots[2], 4.0, 1e-12);

  // (x - 2)(x^2 + 1) = x^3 - 2x^2 + x - 2: one real root
  double single[4] = {-2.0, 1.0, -2.0, 1.0};
  cubic(single, roots, &numberOfRoots);
  ASSERT_EQ(numberOfRoots, 1);
  EXPECT_NEAR(roots[0], 2.0, 1e-12);
}

TEST(MATHKIT_FIXES, random_vector_on_cone_redraws_a_direction_parallel_to_the_axis)
{
  // the axis is the unit vector the seeded generator draws first, so the sampler's first draw is exactly
  // (anti)parallel to it and its perpendicular component vanishes; it used to return +-axis
  for (const double sign : {1.0, -1.0})
  {
    RandomNumber draw(7);
    const double3 axis = sign * draw.randomVectorOnUnitSphere();
    RandomNumber random(7);
    const double angle = 0.4;
    const double3 w = random.randomVectorOnCone(axis, angle);
    EXPECT_NEAR(double3::angle(axis, w), angle, 1e-10) << "sign " << sign;
    EXPECT_NEAR(w.length(), 1.0, 1e-12) << "sign " << sign;
  }
}

TEST(MATHKIT_FIXES, marching_cubes_degenerate_interior_test_keeps_the_vertices_finite)
{
  // a single cube with opposite corners 0 and 6 positive (case 4) whose interior-test quadratic has a == b == 0,
  // so the stationary point -b / (2a) is NaN; the components stay separated: two triangles
  MarchingCubes cubes(2, 2, 2);
  cubes.init_all();
  const double corner[8] = {1.0, -1.0, -1.0, -3.0, -1.0, -3.0, 1.0, -1.0};
  for (std::size_t p = 0; p < 8; ++p)
  {
    cubes.set_data(corner[p], (p ^ (p >> 1)) & 1, (p >> 1) & 1, p >> 2);
  }
  cubes.run(0.0);
  EXPECT_EQ(cubes.ntrigs(), 2uz);
  ASSERT_GT(cubes.nverts(), 0uz);
  for (std::size_t v = 0; v < cubes.nverts(); ++v)
  {
    const Vertex& vertex = cubes.vertices()[v];
    EXPECT_TRUE(std::isfinite(vertex.x) && std::isfinite(vertex.y) && std::isfinite(vertex.z)) << "vertex " << v;
    EXPECT_TRUE(std::isfinite(vertex.nx) && std::isfinite(vertex.ny) && std::isfinite(vertex.nz)) << "vertex " << v;
  }
}

TEST(MATHKIT_FIXES, marching_cubes_periodic_single_plane_has_a_finite_gradient)
{
  // a periodic grid of one plane along x: its gradient along x used to take a remainder modulo zero; the field
  // changes sign along y, so the intersection pass evaluates gradients at every point of the plane
  MarchingCubes cubes(1, 4, 4);
  cubes.init_all();
  cubes.set_periodic(true);
  for (std::size_t k = 0; k < 4; ++k)
    for (std::size_t j = 0; j < 4; ++j) cubes.set_data(static_cast<double>(j) - 1.5, 0, j, k);
  cubes.run(0.0);
  ASSERT_GT(cubes.nverts(), 0uz);
  for (std::size_t v = 0; v < cubes.nverts(); ++v)
  {
    const Vertex& vertex = cubes.vertices()[v];
    EXPECT_TRUE(std::isfinite(vertex.x) && std::isfinite(vertex.y) && std::isfinite(vertex.z)) << "vertex " << v;
    EXPECT_TRUE(std::isfinite(vertex.nx) && std::isfinite(vertex.ny) && std::isfinite(vertex.nz)) << "vertex " << v;
    EXPECT_NEAR(vertex.ny, 1.0, 1e-12) << "vertex " << v;
  }
}

TEST(MATHKIT_FIXES, int3x3_from_a_value_is_diagonal)
{
  const int3x3 m(5);
  for (std::size_t i = 0; i < 9; ++i) EXPECT_EQ(m.m[i], (i % 4 == 0) ? 5 : 0) << "element " << i;
}

TEST(MATHKIT_FIXES, hermite_and_smith_normal_forms_transform_the_matrix)
{
  // a rank-deficient (rank one) and a full-rank matrix with more rows than the rank or equal to it
  for (const RingMatrix& input : {RingMatrix(3, 2, {{2, 4}, {1, 2}, {3, 6}}), RingMatrix(2, 2, {{2, 1}, {1, 3}}),
                                  RingMatrix(3, 3, {{2, 0, 1}, {4, 1, 3}, {6, 1, 4}})})
  {
    RingMatrix A = input;
    const auto [T, H, rankProfile] = A.HermiteNormalForm();
    const RingMatrix TA = T * input;
    for (std::size_t i = 0; i < input.rows(); ++i)
    {
      for (std::size_t j = 0; j < input.columns(); ++j)
      {
        EXPECT_EQ(TA(i, j), H(i, j)) << "Hermite: T A != H at " << i << " " << j;
      }
    }

    RingMatrix B = input;
    const auto [U, V, S] = B.SmithNormalForm();
    const RingMatrix UAV = U * input * V;
    for (std::size_t i = 0; i < input.rows(); ++i)
    {
      for (std::size_t j = 0; j < input.columns(); ++j)
      {
        EXPECT_EQ(UAV(i, j), S(i, j)) << "Smith: U A V != S at " << i << " " << j;
        if (i != j) EXPECT_EQ(S(i, j), 0) << "Smith form not diagonal at " << i << " " << j;
      }
    }
  }
}

TEST(MATHKIT_FIXES, random_ring_matrix_spans_both_signs)
{
  std::srand(7);
  const RingMatrix m = RingMatrix::createRandomRingMatrix(20, 20, 3);
  bool negative = false;
  bool positive = false;
  for (std::size_t i = 0; i < m.rows(); ++i)
  {
    for (std::size_t j = 0; j < m.columns(); ++j)
    {
      EXPECT_GE(m(i, j), -3);
      EXPECT_LE(m(i, j), 3);
      negative = negative || m(i, j) < 0;
      positive = positive || m(i, j) > 0;
    }
  }
  EXPECT_TRUE(negative);
  EXPECT_TRUE(positive);
}

TEST(MATHKIT_FIXES, smith_normal_form_of_random_matrices)
{
  std::srand(11);
  for (std::size_t trial = 0; trial < 50; ++trial)
  {
    const RingMatrix input = RingMatrix::createRandomRingMatrix(3, 3, 4);
    RingMatrix B = input;
    const auto [U, V, S] = B.SmithNormalForm();
    const RingMatrix UAV = U * input * V;
    for (std::size_t i = 0; i < 3; ++i)
    {
      for (std::size_t j = 0; j < 3; ++j)
      {
        EXPECT_EQ(UAV(i, j), S(i, j)) << "trial " << trial << ": U A V != S at " << i << " " << j;
        if (i != j) EXPECT_EQ(S(i, j), 0) << "trial " << trial << ": Smith form not diagonal at " << i << " " << j;
      }
    }
  }
}

TEST(MATHKIT_FIXES, bool3_indices_are_its_components)
{
  const bool3 flags(true, false, true);
  EXPECT_TRUE(flags[0]);
  EXPECT_FALSE(flags[1]);
  EXPECT_TRUE(flags[2]);
  const bool3 cleared;
  EXPECT_FALSE(cleared.x || cleared.y || cleared.z);
  EXPECT_EQ(sizeof(bool3), 3uz);
}

TEST(MATHKIT_FIXES, matrix_inverse_inverts_and_refuses_singular_and_non_square_matrices)
{
  Matrix a(3, 3, 0.0);
  const double values[3][3] = {{4.0, 1.0, 2.0}, {1.0, 3.0, 0.5}, {2.0, 0.5, 5.0}};
  for (std::size_t i = 0; i < 3; ++i)
    for (std::size_t j = 0; j < 3; ++j) a(i, j) = values[i][j];
  Matrix inverse = a;
  ASSERT_TRUE(inverse.inverse());
  const Matrix product = a * inverse;
  for (std::size_t i = 0; i < 3; ++i)
    for (std::size_t j = 0; j < 3; ++j) EXPECT_NEAR(product(i, j), i == j ? 1.0 : 0.0, 1e-12);

  Matrix singular(2, 2, 1.0);
  EXPECT_FALSE(singular.inverse());

  Matrix nonSquare(2, 3, 1.0);
  EXPECT_FALSE(nonSquare.inverse());
}

TEST(MATHKIT_FIXES, cubic_linear_and_repeated_quadratic_roots)
{
  // 2 + 4 x = 0
  double linear[4] = {2.0, 4.0, 0.0, 0.0};
  double roots[3] = {0.0, 0.0, 0.0};
  int numberOfRoots = 0;
  cubic(linear, roots, &numberOfRoots);
  ASSERT_EQ(numberOfRoots, 1);
  EXPECT_DOUBLE_EQ(roots[0], -0.5);

  // x^2 - 2 x + 1 = 0: the repeated root 1
  double quadratic[4] = {1.0, -2.0, 1.0, 0.0};
  cubic(quadratic, roots, &numberOfRoots);
  ASSERT_EQ(numberOfRoots, 2);
  EXPECT_DOUBLE_EQ(roots[0], 1.0);
  EXPECT_DOUBLE_EQ(roots[1], 1.0);
}
