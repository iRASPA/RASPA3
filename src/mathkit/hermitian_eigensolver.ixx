module;

export module hermitian_eigensolver;

import std;

export struct HermitianEigenSystem
{
  std::size_t size{};
  std::vector<double> eigenvalues;
  // Column-major complex eigenvectors: eigenvectors[row + size * mode].
  std::vector<std::complex<double>> eigenvectors;

  std::complex<double> eigenvector(std::size_t row, std::size_t mode) const { return eigenvectors[row + size * mode]; }
};

/**
 * Diagonalize a complex Hermitian matrix supplied in row-major order.
 * Eigenvalues are ascending and (complex) eigenvectors are returned as columns.
 */
export HermitianEigenSystem diagonalizeHermitian(std::span<const std::complex<double>> matrix, std::size_t size);
