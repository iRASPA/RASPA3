#include <gtest/gtest.h>

import std;

import symmetric_eigensolver;

TEST(symmetric_eigensolver, ascending_orthonormal_eigenpairs_have_small_residuals)
{
  constexpr std::size_t n = 3;
  const std::array<double, n * n> matrix = {4.0, 1.0, -2.0, 1.0, 2.0, 0.5, -2.0, 0.5, 3.0};
  const SymmetricEigenSystem result = diagonalizeSymmetric(matrix, n);

  ASSERT_EQ(result.eigenvalues.size(), n);
  EXPECT_LE(result.eigenvalues[0], result.eigenvalues[1]);
  EXPECT_LE(result.eigenvalues[1], result.eigenvalues[2]);
  for (std::size_t mode = 0; mode < n; ++mode)
  {
    double norm = 0.0;
    for (std::size_t row = 0; row < n; ++row)
    {
      norm += result.eigenvector(row, mode) * result.eigenvector(row, mode);
      double residual = -result.eigenvalues[mode] * result.eigenvector(row, mode);
      for (std::size_t column = 0; column < n; ++column)
      {
        residual += matrix[row * n + column] * result.eigenvector(column, mode);
      }
      EXPECT_NEAR(residual, 0.0, 1e-12);
    }
    EXPECT_NEAR(norm, 1.0, 1e-12);
  }
}
