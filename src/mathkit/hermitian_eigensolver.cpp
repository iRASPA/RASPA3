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

namespace
{
// Whether a number is finite, read off its bits: the library is compiled with finite mathematics only, under which
// 'std::isfinite' is true of everything. (Taken by reference and read as bytes, so that the test is not folded
// away with it.)
bool isFiniteNumber(const double& value)
{
  std::uint64_t bits{};
  std::memcpy(&bits, &value, sizeof(bits));
  return (bits & 0x7ff0000000000000ull) != 0x7ff0000000000000ull;
}
}  // namespace

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
  // LAPACK reads the upper triangle (row <= column): only that half is filled, column by column at stride one.
  for (std::size_t column = 0; column < size; ++column)
  {
    for (std::size_t row = 0; row <= column; ++row)
    {
      const std::complex<double> value = 0.5 * (matrix[row * size + column] + std::conj(matrix[column * size + row]));
      if (!isFiniteNumber(value.real()) || !isFiniteNumber(value.imag()))
      {
        throw std::runtime_error("diagonalizeHermitian: matrix contains a non-finite value");
      }
      storage[row + size * column] = value;
    }
  }

  // The workspace of the eigenvectors, 2 n^2 + 5 n + 1 reals, has to be counted in the LAPACK integer type, which
  // LAPACK does itself in the size query (from n = 32768 on with 32-bit integers it overflowed)
  const double nn = static_cast<double>(size);
  if (computeEigenvectors &&
      2.0 * nn * nn + 5.0 * nn + 1.0 > static_cast<double>(std::numeric_limits<blas_int>::max()))
  {
    throw std::overflow_error("diagonalizeHermitian: the workspace is too large for the configured LAPACK integer type");
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
