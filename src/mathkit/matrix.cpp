module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cctype>
#include <cstddef>
#include <vector>
#endif

module matrix;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

Matrix::Matrix(std::size_t rows, std::size_t columns, double initialValue)
    : _rows(rows), _columns(columns), _grid(rows * columns)
{
  std::fill(_grid.begin(), _grid.end(), initialValue);
}

Matrix operator*(const Matrix& lhs, const Matrix& rhs)
{
  Matrix results = Matrix(lhs._rows, rhs._columns, 0);
  for (std::size_t i = 0; i < lhs._rows; i++)
  {
    for (std::size_t j = 0; j < rhs._columns; j++)
    {
      double v = 0;
      for (std::size_t k = 0; k < lhs.columns(); k++)
      {
        v += lhs(i, k) * rhs(k, j);
      }
      results(i, j) = v;
    }
  }
  return results;
}
