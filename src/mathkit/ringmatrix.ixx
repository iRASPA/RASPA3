module;

#ifdef USE_LEGACY_HEADERS
#include <any>
#include <cstddef>
#include <optional>
#include <string>
#include <variant>
#include <vector>
#endif

export module ringmatrix;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import int3x3;

export class RingMatrix
{
 public:
  RingMatrix(std::size_t rows, std::size_t columns, int initialValue);
  RingMatrix(std::size_t rows, std::size_t columns, std::vector<std::vector<int>> data);
  RingMatrix(int3x3 t1, int3x3 t2, int3x3 t3);

  std::size_t rows() const { return _rows; }
  std::size_t columns() const { return _columns; }

  static RingMatrix createRandomRingMatrix(std::size_t rows, std::size_t columns, int limit);

  RingMatrix transposed();

  int& operator()(std::size_t row, std::size_t col) { return _grid.at(row * _columns + col); }
  int operator()(std::size_t row, std::size_t col) const { return _grid.at(row * _columns + col); }
  static RingMatrix identity(std::size_t size);
  static RingMatrix zeros(std::size_t rows, std::size_t columns);
  static RingMatrix ones(std::size_t rows, std::size_t columns);

  inline bool operator==(const RingMatrix& b) const
  {
    if (this->_rows != b._rows) return false;
    if (this->_columns != b._columns) return false;

    for (std::size_t i = 0; i < this->_rows; i++)
    {
      for (std::size_t j = 0; j < b._columns; j++)
      {
        if ((*this)(i, j) != b(i, j)) return false;
      }
    }
    return true;
  }

  std::tuple<RingMatrix, RingMatrix, std::vector<std::size_t>> HermiteNormalForm();
  std::tuple<RingMatrix, RingMatrix, RingMatrix> SmithNormalForm();

 private:
  std::size_t _rows;
  std::size_t _columns;
  std::vector<int> _grid;

  const static int HNF_C_Iwaniec;

  void swapColumns(std::size_t a, std::size_t b);

  RingMatrix submatrix(std::size_t startRow, std::size_t startColumn, std::size_t numberOfRows,
                       std::size_t numberOfColumns);
  void assignSubmatrix(std::size_t startRow, std::size_t startColumn, RingMatrix replacement);
  std::tuple<RingMatrix, RingMatrix, RingMatrix> ColumnReduction(RingMatrix A1, std::size_t col_1, std::size_t col_2,
                                                                 std::size_t row_start);
  std::pair<std::vector<int>, std::vector<int>> Conditioning(RingMatrix A, std::size_t col_1, std::size_t col_2,
                                                             std::size_t row_start);
  std::variant<std::vector<int>, int> Algorithm_6_14(int a, int b, int N, std::vector<int> Nfact);
  std::vector<int> RemovedDuplicates(std::vector<int> array);

  void Smith_Theorem5(RingMatrix& A, RingMatrix& U, RingMatrix& V, std::size_t col);
  void Smith_Theorem8(RingMatrix& A, RingMatrix& U, RingMatrix& V, std::size_t row, std::size_t r);
  std::vector<int> Algorithm_6_15(int a, std::vector<int> bi, int N);
};

export inline RingMatrix operator*(const RingMatrix& a, const RingMatrix& b)
{
  RingMatrix results = RingMatrix(a.rows(), b.columns(), 0);
  for (std::size_t i = 0; i < a.rows(); i++)
  {
    for (std::size_t j = 0; j < b.columns(); j++)
    {
      int v = 0;
      for (std::size_t k = 0; k < a.columns(); k++)
      {
        v += a(i, k) * b(k, j);
      }
      results(i, j) = v;
    }
  }
  return results;
}
