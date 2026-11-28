module;

#ifdef BLAS_ILP64
using blas_int = long long;
#else
using blas_int = int;
#endif

extern "C"
{
  void dsyevd_(char *jobz, char *uplo, blas_int *n, double *a, blas_int *lda, double *w, double *work,
               blas_int *lwork, blas_int *iwork, blas_int *liwork, blas_int *info);
}

module symmetric_eigensolver;

import std;

import finite_value;

SymmetricEigenSystem diagonalizeSymmetric(std::span<const double> matrix, std::size_t size)
{
  if (matrix.size() != size * size)
  {
    throw std::invalid_argument("diagonalizeSymmetric: matrix extent does not match size");
  }
  if (size == 0)
  {
    return {};
  }
  if (size > static_cast<std::size_t>(std::numeric_limits<blas_int>::max()))
  {
    throw std::overflow_error("diagonalizeSymmetric: matrix is too large for the configured LAPACK integer type");
  }

  SymmetricEigenSystem result{.size = size, .eigenvalues = std::vector<double>(size),
                              .eigenvectors = std::vector<double>(size * size)};
  // A symmetric matrix has the same values under row-major/column-major reinterpretation.
  // Copy explicitly while symmetrizing to avoid propagating small assembly asymmetries.
  // LAPACK reads the upper triangle (row <= column) of the column-major array: only that half is filled, column
  // by column at stride one; the symmetrized value takes both triangles of the input in.
  for (std::size_t column = 0; column < size; ++column)
  {
    for (std::size_t row = 0; row <= column; ++row)
    {
      const double value = 0.5 * (matrix[row * size + column] + matrix[column * size + row]);
      if (!isFiniteValue(value))
      {
        throw std::runtime_error("diagonalizeSymmetric: matrix contains a non-finite value");
      }
      result.eigenvectors[row + size * column] = value;
    }
  }

  char jobz = 'V';
  char uplo = 'U';
  // The workspace of the eigenvectors, 2 n^2 + 6 n + 1, has to be counted in the LAPACK integer type, which LAPACK
  // does itself in the size query (from n = 32768 on with 32-bit integers it overflowed)
  const double nn = static_cast<double>(size);
  if (2.0 * nn * nn + 6.0 * nn + 1.0 > static_cast<double>(std::numeric_limits<blas_int>::max()))
  {
    throw std::overflow_error("diagonalizeSymmetric: the workspace is too large for the configured LAPACK integer type");
  }
  blas_int n = static_cast<blas_int>(size);
  blas_int lda = n;
  blas_int info = 0;
  blas_int lwork = -1;
  blas_int liwork = -1;
  double workQuery = 0.0;
  blas_int iworkQuery = 0;
  dsyevd_(&jobz, &uplo, &n, result.eigenvectors.data(), &lda, result.eigenvalues.data(), &workQuery, &lwork,
          &iworkQuery, &liwork, &info);
  if (info != 0)
  {
    throw std::runtime_error(std::format("LAPACK dsyevd workspace query failed with info={}", info));
  }

  lwork = static_cast<blas_int>(workQuery);
  liwork = iworkQuery;
  std::vector<double> work(static_cast<std::size_t>(lwork));
  std::vector<blas_int> iwork(static_cast<std::size_t>(liwork));
  dsyevd_(&jobz, &uplo, &n, result.eigenvectors.data(), &lda, result.eigenvalues.data(), work.data(), &lwork,
          iwork.data(), &liwork, &info);
  if (info < 0)
  {
    throw std::invalid_argument(std::format("LAPACK dsyevd received an invalid argument at position {}", -info));
  }
  if (info > 0)
  {
    throw std::runtime_error(std::format("LAPACK dsyevd failed to converge (info={})", info));
  }
  return result;
}
