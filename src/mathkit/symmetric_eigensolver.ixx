module;

export module symmetric_eigensolver;

import std;

export struct SymmetricEigenSystem
{
  std::size_t size{};
  std::vector<double> eigenvalues;
  // Column-major: eigenvectors[row + size * mode].
  std::vector<double> eigenvectors;

  double eigenvector(std::size_t row, std::size_t mode) const
  {
    return eigenvectors[row + size * mode];
  }
};

/**
 * Diagonalize a real symmetric matrix supplied in row-major order.
 * Eigenvalues are ascending and eigenvectors are returned as columns.
 */
export SymmetricEigenSystem diagonalizeSymmetric(std::span<const double> matrix, std::size_t size);
