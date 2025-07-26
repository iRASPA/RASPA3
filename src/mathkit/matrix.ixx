module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <iostream>
#include <optional>
#include <vector>
#endif

export module matrix;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import ringmatrix;

extern "C"
{
  // LU decomoposition of a general matrix
  void dgetrf_(size_t* M, size_t* N, double* A, size_t* lda, int* IPIV, int* INFO);

  // generate inverse of a matrix given its LU decomposition
  void dgetri_(size_t* N, double* A, size_t* lda, int* IPIV, double* WORK, size_t* lwork, int* INFO);
}

export class Matrix
{
 public:
  Matrix(size_t rows, size_t columns, double initialValue);

  size_t rows() const { return _rows; }
  size_t columns() const { return _columns; }

  double& operator()(size_t row, size_t col) { return _grid[row * _columns + col]; }
  double operator()(size_t row, size_t col) const { return _grid[row * _columns + col]; }

  friend Matrix operator*(const Matrix& v1, const Matrix& v2);

  // Almost all LAPACK routines return an info::BlasInt argument which is one of :
  //
  // 0 if everything went smoothly,
  //     a negative number - n specifying that the nth argument was invalid,
  //     a positive number whose meaning depends on the routine being called.
  //     For some routines, additional error information is encoded in arrays, often in allocated work, lwork or rwork
  //     arrays.
  inline void inverse()
  {
    size_t N = _rows;
    if (N > 0)
    {
      int* IPIV = new int[N];
      size_t LWORK = N * N;
      double* WORK = new double[LWORK];
      int INFO;

      dgetrf_(&N, &N, _grid.data(), &N, IPIV, &INFO);
      dgetri_(&N, _grid.data(), &N, IPIV, WORK, &LWORK, &INFO);

      delete[] IPIV;
      delete[] WORK;
    }
  }

 private:
  size_t _rows;
  size_t _columns;
  std::vector<double> _grid;
};

export inline Matrix operator*(const RingMatrix& a, const Matrix& b)
{
  Matrix results = Matrix(a.rows(), b.columns(), 0);
  for (size_t i = 0; i < a.rows(); i++)
  {
    for (size_t j = 0; j < b.columns(); j++)
    {
      double v = 0;
      for (size_t k = 0; k < a.columns(); k++)
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
  Matrix results = Matrix(a.rows(), b.columns(), 0);
  for (size_t i = 0; i < a.rows(); i++)
  {
    for (size_t j = 0; j < b.columns(); j++)
    {
      double v = 0;
      for (size_t k = 0; k < a.columns(); k++)
      {
        v += a(i, k) * b(k, j);
      }
      results(i, j) = v;
    }
  }
  return results;
}
