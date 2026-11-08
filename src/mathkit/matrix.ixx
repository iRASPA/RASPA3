module;

// the LAPACK integer: 64 bits for an ILP64 build of BLAS/LAPACK, 32 bits otherwise (as in double3x3.cpp)
#ifdef BLAS_ILP64
typedef long long blas_int;
#else
typedef int blas_int;
#endif

export module matrix;

import std;

import ringmatrix;

extern "C"
{
  // LU decomoposition of a general matrix
  void dgetrf_(blas_int* M, blas_int* N, double* A, blas_int* lda, blas_int* IPIV, blas_int* INFO);

  // generate inverse of a matrix given its LU decomposition
  void dgetri_(blas_int* N, double* A, blas_int* lda, blas_int* IPIV, double* WORK, blas_int* lwork, blas_int* INFO);
}

export class Matrix
{
 public:
  Matrix(std::size_t rows, std::size_t columns, double initialValue);

  std::size_t rows() const { return _rows; }
  std::size_t columns() const { return _columns; }

  double& operator()(std::size_t row, std::size_t col) { return _grid[row * _columns + col]; }
  double operator()(std::size_t row, std::size_t col) const { return _grid[row * _columns + col]; }

  friend Matrix operator*(const Matrix& v1, const Matrix& v2);

  // Almost all LAPACK routines return an info::BlasInt argument which is one of :
  //
  // 0 if everything went smoothly,
  //     a negative number - n specifying that the nth argument was invalid,
  //     a positive number whose meaning depends on the routine being called.
  //     For some routines, additional error information is encoded in arrays, often in allocated work, lwork or rwork
  //     arrays.
  /// Inverts the matrix in place. Returns false when the matrix is singular (a zero pivot in the LU
  /// factorization): the contents are then not an inverse (the LU factors, or a zero for a 1 x 1 matrix)
  /// and must not be used. A matrix that is not square has no inverse either.
  [[nodiscard]] inline bool inverse()
  {
    if (_rows != _columns) return false;
    if (_rows == 0) return true;

    blas_int N = static_cast<blas_int>(_rows);
    blas_int LWORK = N * N;
    std::vector<blas_int> IPIV(_rows);
    std::vector<double> WORK(_rows * _rows);
    blas_int INFO = 0;

    dgetrf_(&N, &N, _grid.data(), &N, IPIV.data(), &INFO);
    bool invertible = (INFO == 0);
    if (invertible)
    {
      dgetri_(&N, _grid.data(), &N, IPIV.data(), WORK.data(), &LWORK, &INFO);
      invertible = (INFO == 0);
    }
    return invertible;
  }

 private:
  std::size_t _rows;
  std::size_t _columns;
  std::vector<double> _grid;
};

export inline Matrix operator*(const RingMatrix& a, const Matrix& b)
{
  if (a.columns() != b.rows())
  {
    throw std::invalid_argument(
        std::format("Matrix product: {} columns of the left factor do not match {} rows of the right factor",
                    a.columns(), b.rows()));
  }
  Matrix results = Matrix(a.rows(), b.columns(), 0);
  for (std::size_t i = 0; i < a.rows(); i++)
  {
    for (std::size_t j = 0; j < b.columns(); j++)
    {
      double v = 0;
      for (std::size_t k = 0; k < a.columns(); k++)
      {
        v += a(i, k) * b(k, j);
      }
      results(i, j) = v;
    }
  }
  return results;
}

export inline Matrix operator*(const Matrix& a, const RingMatrix& b)
{
  if (a.columns() != b.rows())
  {
    throw std::invalid_argument(
        std::format("Matrix product: {} columns of the left factor do not match {} rows of the right factor",
                    a.columns(), b.rows()));
  }
  Matrix results = Matrix(a.rows(), b.columns(), 0);
  for (std::size_t i = 0; i < a.rows(); i++)
  {
    for (std::size_t j = 0; j < b.columns(); j++)
    {
      double v = 0;
      for (std::size_t k = 0; k < a.columns(); k++)
      {
        v += a(i, k) * b(k, j);
      }
      results(i, j) = v;
    }
  }
  return results;
}
