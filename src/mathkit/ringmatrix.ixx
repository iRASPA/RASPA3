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
import <vector>;
import <optional>;
import <any>;
import <variant>;
import <string>;
#endif

import int3x3;

export class RingMatrix
{
 public:
  RingMatrix(size_t rows, size_t columns, int initialValue);
  RingMatrix(size_t rows, size_t columns, std::vector<std::vector<int>> data);
  RingMatrix(int3x3 t1, int3x3 t2, int3x3 t3);

  size_t rows() const { return _rows; }
  size_t columns() const { return _columns; }

  static RingMatrix createRandomRingMatrix(size_t rows, size_t columns, int limit);

  RingMatrix transposed();

  int& operator()(size_t row, size_t col) { return _grid.at(row * _columns + col); }
  int operator()(size_t row, size_t col) const { return _grid.at(row * _columns + col); }
  static RingMatrix identity(size_t size);
  static RingMatrix zeros(size_t rows, size_t columns);
  static RingMatrix ones(size_t rows, size_t columns);

  inline bool operator==(const RingMatrix& b) const
  {
    if (this->_rows != b._rows) return false;
    if (this->_columns != b._columns) return false;

    for (size_t i = 0; i < this->_rows; i++)
    {
      for (size_t j = 0; j < b._columns; j++)
      {
        if ((*this)(i, j) != b(i, j)) return false;
      }
    }
    return true;
  }

  std::tuple<RingMatrix, RingMatrix, std::vector<size_t>> HermiteNormalForm();
  std::tuple<RingMatrix, RingMatrix, RingMatrix> SmithNormalForm();

 private:
  size_t _rows;
  size_t _columns;
  std::vector<int> _grid;

  const static int HNF_C_Iwaniec;

  void swapColumns(size_t a, size_t b);

  RingMatrix submatrix(size_t startRow, size_t startColumn, size_t numberOfRows, size_t numberOfColumns);
  void assignSubmatrix(size_t startRow, size_t startColumn, RingMatrix replacement);
  std::tuple<RingMatrix, RingMatrix, RingMatrix> ColumnReduction(RingMatrix A1, size_t col_1, size_t col_2,
                                                                 size_t row_start);
  std::pair<std::vector<int>, std::vector<int>> Conditioning(RingMatrix A, size_t col_1, size_t col_2,
                                                             size_t row_start);
  std::variant<std::vector<int>, int> Algorithm_6_14(int a, int b, int N, std::vector<int> Nfact);
  std::vector<int> RemovedDuplicates(std::vector<int> array);

  void Smith_Theorem5(RingMatrix& A, RingMatrix& U, RingMatrix& V, size_t col);
  void Smith_Theorem8(RingMatrix& A, RingMatrix& U, RingMatrix& V, size_t row, size_t r);
  std::vector<int> Algorithm_6_15(int a, std::vector<int> bi, int N);
};

export inline RingMatrix operator*(const RingMatrix& a, const RingMatrix& b)
{
  RingMatrix results = RingMatrix(a.rows(), b.columns(), 0);
  for (size_t i = 0; i < a.rows(); i++)
  {
    for (size_t j = 0; j < b.columns(); j++)
    {
      int v = 0;
      for (size_t k = 0; k < a.columns(); k++)
      {
        v += a(i, k) * b(k, j);
      }
      results(i, j) = v;
    }
  }
  return results;
}
