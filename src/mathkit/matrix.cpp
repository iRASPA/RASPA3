module;

module matrix;

import std;

Matrix::Matrix(std::size_t rows, std::size_t columns, double initialValue)
    : _rows(rows), _columns(columns), _grid(rows * columns, initialValue)
{
}

Matrix operator*(const Matrix& lhs, const Matrix& rhs)
{
  // the inner sum reads rhs(k, j) for k below the columns of lhs: past its storage when the shapes do not match
  if (lhs._columns != rhs._rows)
  {
    throw std::invalid_argument(
        std::format("Matrix product: {} columns of the left factor do not match {} rows of the right factor",
                    lhs._columns, rhs._rows));
  }
  Matrix results = Matrix(lhs._rows, rhs._columns, 0);
  // i, k, j: the inner loop runs along a row of rhs and the same row of the result at stride one (the sum over k
  // of every entry is still taken in ascending k, so the values are the ones of the i, j, k order)
  for (std::size_t i = 0; i < lhs._rows; i++)
  {
    for (std::size_t k = 0; k < lhs.columns(); k++)
    {
      const double aik = lhs(i, k);
      for (std::size_t j = 0; j < rhs._columns; j++)
      {
        results(i, j) += aik * rhs(k, j);
      }
    }
  }
  return results;
}
