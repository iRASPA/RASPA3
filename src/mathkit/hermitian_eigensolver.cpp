module;

#ifdef BLAS_ILP64
using blas_int = long long;
#else
using blas_int = int;
#endif

extern "C"
{
  void zheevd_(char* jobz, char* uplo, blas_int* n, double* a, blas_int* lda, double* w, double* work, blas_int* lwork,
               double* rwork, blas_int* lrwork, blas_int* iwork, blas_int* liwork, blas_int* info);
}

module hermitian_eigensolver;

import std;

HermitianEigenSystem diagonalizeHermitian(std::span<const std::complex<double>> matrix, std::size_t size,
                                          bool computeEigenvectors)
{
  if (matrix.size() != size * size)
  {
    throw std::invalid_argument("diagonalizeHermitian: matrix extent does not match size");
  }
  if (size == 0)
  {
    return {};
  }
  if (size > static_cast<std::size_t>(std::numeric_limits<blas_int>::max()))
  {
    throw std::overflow_error("diagonalizeHermitian: matrix is too large for the configured LAPACK integer type");
  }

  HermitianEigenSystem result{
      .size = size,
      .eigenvalues = std::vector<double>(size),
      .eigenvectors = std::vector<std::complex<double>>(computeEigenvectors ? size * size : std::size_t{0})};

  // In eigenvalues-only mode LAPACK still overwrites the input matrix, so use scratch storage instead of
  // the (empty) eigenvector array.
  std::vector<std::complex<double>> scratch;
  std::complex<double>* storage = result.eigenvectors.data();
  if (!computeEigenvectors)
  {
    scratch.resize(size * size);
    storage = scratch.data();
  }

  // Hermitize into column-major storage: a[row + size*column] = 0.5 (M[row][column] + conj(M[column][row])).
  for (std::size_t row = 0; row < size; ++row)
  {
    for (std::size_t column = 0; column < size; ++column)
    {
      const std::complex<double> value = 0.5 * (matrix[row * size + column] + std::conj(matrix[column * size + row]));
      if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
      {
        throw std::runtime_error("diagonalizeHermitian: matrix contains a non-finite value");
      }
      storage[row + size * column] = value;
    }
  }

  char jobz = computeEigenvectors ? 'V' : 'N';
  char uplo = 'U';
  blas_int n = static_cast<blas_int>(size);
  blas_int lda = n;
  blas_int info = 0;
  blas_int lwork = -1;
  blas_int lrwork = -1;
  blas_int liwork = -1;
  double workQuery[2] = {0.0, 0.0};
  double rworkQuery = 0.0;
  blas_int iworkQuery = 0;
  double* a = reinterpret_cast<double*>(storage);
  zheevd_(&jobz, &uplo, &n, a, &lda, result.eigenvalues.data(), workQuery, &lwork, &rworkQuery, &lrwork, &iworkQuery,
          &liwork, &info);
  if (info != 0)
  {
    throw std::runtime_error(std::format("LAPACK zheevd workspace query failed with info={}", info));
  }

  lwork = static_cast<blas_int>(workQuery[0]);
  lrwork = static_cast<blas_int>(rworkQuery);
  liwork = iworkQuery;
  std::vector<std::complex<double>> work(static_cast<std::size_t>(lwork));
  std::vector<double> rwork(static_cast<std::size_t>(lrwork));
  std::vector<blas_int> iwork(static_cast<std::size_t>(liwork));
  zheevd_(&jobz, &uplo, &n, a, &lda, result.eigenvalues.data(), reinterpret_cast<double*>(work.data()), &lwork,
          rwork.data(), &lrwork, iwork.data(), &liwork, &info);
  if (info < 0)
  {
    throw std::invalid_argument(std::format("LAPACK zheevd received an invalid argument at position {}", -info));
  }
  if (info > 0)
  {
    throw std::runtime_error(std::format("LAPACK zheevd failed to converge (info={})", info));
  }
  return result;
}
