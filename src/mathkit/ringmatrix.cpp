module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <variant>
#include <vector>
#endif

module ringmatrix;

#ifndef USE_LEGACY_HEADERS
import <algorithm>;
import <cmath>;
import <limits>;
import <vector>;
import <variant>;
#endif

import ring;
import int3x3;

template <typename T>
int sign(T val)
{
  return (T(0) < val) - (val < T(0));
}

RingMatrix::RingMatrix(size_t rows, size_t columns, int initialValue)
    : _rows(rows), _columns(columns), _grid(rows * columns)
{
  std::fill(_grid.begin(), _grid.end(), initialValue);
}

RingMatrix RingMatrix::createRandomRingMatrix(size_t rows, size_t columns, int limit)
{
  RingMatrix m = RingMatrix(rows, columns, 0);

  for (size_t i = 0; i < rows; i++)
  {
    for (size_t j = 0; j < columns; j++)
    {
      m(i, j) = (rand() % limit) - limit;
    }
  }
  return m;
}

RingMatrix::RingMatrix(size_t rows, size_t columns, std::vector<std::vector<int>> data)
    : _rows(rows), _columns(columns), _grid(rows * columns)
{
  std::fill(_grid.begin(), _grid.end(), 0);
  for (size_t i = 0; i < data.size(); i++)
  {
    for (size_t j = 0; j < data[i].size(); j++)
    {
      (*this)(i, j) = data[i][j];
    }
  }
}

RingMatrix::RingMatrix(int3x3 t1, int3x3 t2, int3x3 t3)
{
  _columns = 3;
  _rows = 9;
  _grid = std::vector<int>(_rows * _columns);

  for (size_t row = 0; row < 3; row++)
  {
    for (size_t column = 0; column < 3; column++)
    {
      (*this)(row + 0 * 3, column) = t1.mm[column][row];
    }
  }
  for (size_t row = 0; row < 3; row++)
  {
    for (size_t column = 0; column < 3; column++)
    {
      (*this)(row + 1 * 3, column) = t2.mm[column][row];
    }
  }
  for (size_t row = 0; row < 3; row++)
  {
    for (size_t column = 0; column < 3; column++)
    {
      (*this)(row + 2 * 3, column) = t3.mm[column][row];
    }
  }
}

RingMatrix RingMatrix::transposed()
{
  RingMatrix result = RingMatrix(this->_columns, this->_rows, 0);
  for (size_t i = 0; i < this->_rows; i++)
  {
    for (size_t j = 0; j < this->_columns; j++)
    {
      result(j, i) = (*this)(i, j);
    }
  }
  return result;
}

void RingMatrix::swapColumns(size_t a, size_t b)
{
  for (size_t i = 0; i < this->_rows; i++)
  {
    int temp = (*this)(i, a);
    (*this)(i, a) = (*this)(i, b);
    (*this)(i, b) = temp;
  }
}

// MARK: - HermiteNormalForm

RingMatrix RingMatrix::submatrix(size_t startRow, size_t startColumn, size_t numberOfRows, size_t numberOfColumns)
{
  RingMatrix matrix = RingMatrix(numberOfRows, numberOfColumns, 0);
  for (size_t i = 0; i < numberOfRows; i++)
  {
    for (size_t j = 0; j < numberOfColumns; j++)
    {
      matrix(i, j) = (*this)(i + startRow, j + startColumn);
    }
  }

  return matrix;
}

void RingMatrix::assignSubmatrix(size_t startRow, size_t startColumn, RingMatrix replacement)
{
  for (size_t i = 0; i < replacement._rows; i++)
  {
    for (size_t j = 0; j < replacement._columns; j++)
    {
      (*this)(startRow + i, startColumn + j) = replacement(i, j);
    }
  }
}

RingMatrix RingMatrix::zeros(size_t rows, size_t columns) { return RingMatrix(rows, columns, 0); }

RingMatrix RingMatrix::ones(size_t rows, size_t columns) { return RingMatrix(rows, columns, 1); }

/* Creates a (square) identity RingMatrix. */
RingMatrix RingMatrix::identity(size_t size)
{
  RingMatrix m = RingMatrix::zeros(size, size);
  for (size_t i = 0; i < size; i++)
  {
    m(i, i) = 1;
  }
  return m;
}

// Hermite Normal Form
// ====================================================================================

std::tuple<RingMatrix, RingMatrix, std::vector<size_t>> RingMatrix::HermiteNormalForm()
{
  // Create larger matrix
  RingMatrix Apad = RingMatrix(this->_rows + 2, this->_columns + 2, 0);
  Apad(0, 0) = 1;
  Apad(this->_rows + 1, this->_columns + 1) = 1;
  Apad.assignSubmatrix(1, 1, *this);

  std::vector<size_t> rp = {0};
  size_t r = 0;

  // Create transformation matrix
  RingMatrix QQ = RingMatrix::identity(this->_rows + 2);
  RingMatrix CC = RingMatrix::identity(this->_rows + 2);

  // Begin computing the HNF
  for (size_t j = 1; j < Apad._columns; j++)
  {
    // Search for next rank increase
    bool found = false;
    for (size_t k = (r + 1); k < Apad._rows; k++)
    {
      if (Apad(r, rp[r]) * Apad(k, j) != Apad(r, j) * Apad(k, rp[r]))
      {
        found = true;
      }
    }

    // Found something?
    if (found)
    {
      // Increase rank
      rp.push_back(j);
      r += 1;

      // Do column reduction
      // (Q: RingMatrix, C: RingMatrix, Apad: RingMatrix)
      std::tuple<RingMatrix, RingMatrix, RingMatrix> columnReduction = ColumnReduction(Apad, rp[r - 1], rp[r], r - 1);
      Apad = std::get<2>(columnReduction);

      // Update CC
      for (size_t i = (r + 1); i < std::get<1>(columnReduction)._columns; i++)
      {
        CC(r, i) = std::get<1>(columnReduction)(r, i);
      }
      // Update QQ to QQ * C^{-1}
      for (size_t l = (r + 1); l < std::get<1>(columnReduction)._columns; l++)
      {
        if (std::get<1>(columnReduction)(r, l) != 0)
        {
          for (size_t i = 0; i < QQ._rows; i++)
          {
            QQ(i, l) -= std::get<1>(columnReduction)(r, l) * QQ(i, r);
          }
        }
      }
      // Update QQ to C * QQ
      for (size_t i = (r + 1); i < std::get<1>(columnReduction)._columns; i++)
      {
        if (std::get<1>(columnReduction)(r, i) != 0)
        {
          for (size_t k = 0; k < QQ._rows; k++)
          {
            QQ(r, k) += std::get<1>(columnReduction)(r, i) * QQ(i, k);
          }
        }
      }
      // Update QQ to Q * QQ
      for (size_t i = 0; i < QQ._rows; i++)
      {
        if (i != r - 1 && i != r)
        {
          for (size_t l = 0; l < QQ._columns; l++)
          {
            QQ(i, l) = QQ(i, l) + std::get<0>(columnReduction)(i, r - 1) * QQ(r - 1, l) +
                       std::get<0>(columnReduction)(i, r) * QQ(r, l);
          }
        }
      }
      int a = std::get<0>(columnReduction)(r - 1, r - 1);
      int b = std::get<0>(columnReduction)(r - 1, r);
      int c = std::get<0>(columnReduction)(r, r - 1);
      int d = std::get<0>(columnReduction)(r, r);
      for (size_t l = 0; l < QQ._columns; l++)
      {
        int temp1 = a * QQ(r - 1, l) + b * QQ(r, l);
        int temp2 = c * QQ(r - 1, l) + d * QQ(r, l);
        QQ(r - 1, l) = temp1;
        QQ(r, l) = temp2;
      }
    }
  }

  // Compute the transformation matrix
  RingMatrix T = QQ * CC;

  // Extract the necessary matrices
  RingMatrix TT = RingMatrix(_rows, _rows, 0);

  if (r > 1)
  {
    TT.assignSubmatrix(0, 0, T.submatrix(1, 1, r - 2, _rows));
    TT.assignSubmatrix(r - 2, 0, T.submatrix(r - 1, 1, _rows - r + 2, _rows));
  }
  RingMatrix AA = RingMatrix(Apad.submatrix(1, 1, _rows, _columns));

  // Extract rank profile  FIX
  rp = std::vector<size_t>(rp.begin() + 1,
                           rp.begin() + static_cast<std::vector<size_t>::difference_type>(r));  // rp = Array(rp[1..<r])

  for (size_t i = 0; i < rp.size(); i++)
  {
    rp[i] -= 1;
  }
  return std::make_tuple(TT, AA, rp);
}

std::tuple<RingMatrix, RingMatrix, RingMatrix> RingMatrix::ColumnReduction(RingMatrix A1, size_t col_1, size_t col_2,
                                                                           size_t row_start)
{
  RingMatrix A = A1;
  size_t n = A._rows;
  size_t m = A._columns;

  // Apply conditioning subroutine
  std::pair<std::vector<int>, std::vector<int>> conditioning = Conditioning(A, col_1, col_2, row_start);

  // Initialize C
  RingMatrix C = RingMatrix::identity(n);

  for (size_t j = 0; j < conditioning.second.size(); j++)
  {
    C(row_start + 1, row_start + 2 + j) = conditioning.second[j];
  }
  // Transform A
  for (size_t j = col_1; j < m; j++)
  {
    int v = A(row_start + 1, j);
    for (size_t i = 0; i < conditioning.second.size(); i++)
    {
      v += conditioning.second[i] * A(row_start + 2 + i, j);
    }
    A(row_start + 1, j) = v;
  }

  // Compute Q transform
  // (t1: Int, m1: Int, m2: Int)
  std::tuple<int, int, int> extendedGCD =
      Ring::extendedGreatestCommonDivisor(A(row_start, col_1), A(row_start + 1, col_1));

  int s = sign(A(row_start, col_1) * A(row_start + 1, col_2) - A(row_start, col_2) * A(row_start + 1, col_1));
  RingMatrix Q = RingMatrix::identity(n);
  int q1 = -s * A(row_start + 1, col_1) / std::get<0>(extendedGCD);
  int q2 = s * A(row_start, col_1) / std::get<0>(extendedGCD);
  Q(row_start, row_start) = std::get<1>(extendedGCD);
  Q(row_start, row_start + 1) = std::get<2>(extendedGCD);
  Q(row_start + 1, row_start) = q1;
  Q(row_start + 1, row_start + 1) = q2;

  // Transform A
  int v = std::get<1>(extendedGCD) * A(row_start, col_2) + std::get<2>(extendedGCD) * A(row_start + 1, col_2);
  int t2 = q1 * A(row_start, col_2) + q2 * A(row_start + 1, col_2);
  A(row_start, col_1) = std::get<0>(extendedGCD);
  A(row_start, col_2) = v;
  A(row_start + 1, col_1) = 0;
  A(row_start + 1, col_2) = t2;
  for (size_t j = (col_1 + 1); j < m; j++)
  {
    if (j != col_2)
    {
      int v1 = std::get<1>(extendedGCD) * A(row_start, j) + std::get<2>(extendedGCD) * A(row_start + 1, j);
      int v2 = q1 * A(row_start, j) + q2 * A(row_start + 1, j);
      A(row_start, j) = v1;
      A(row_start + 1, j) = v2;
    }
  }

  // Clean up above
  for (size_t i = 0; i < row_start; i++)
  {
    int s1 = -Ring::floorDivision(A(i, col_1), std::get<0>(extendedGCD));
    for (size_t j = col_1; j < m; j++)
    {
      A(i, j) = A(i, j) + s1 * A(row_start, j);
    }
    int s2 = -Ring::floorDivision(A(i, col_2), t2);
    for (size_t j = col_1; j < m; j++)
    {
      A(i, j) = A(i, j) + s2 * A(row_start + 1, j);
    }
    for (size_t j = 0; j < n; j++)
    {
      int temp = Q(i, j) + s1 * Q(row_start, j) + s2 * Q(row_start + 1, j);
      Q(i, j) = temp;
    }
  }

  // Clean up below
  for (size_t i = (row_start + 2); i < n; i++)
  {
    // assert A[i, col_1] % t1 == 0
    int s1 = -Ring::floorDivision(A(i, col_1), std::get<0>(extendedGCD));
    for (size_t j = col_1; j < m; j++)
    {
      A(i, j) += s1 * A(row_start, j);
    }
    int s2 = -Ring::floorDivision(A(i, col_2), t2);
    for (size_t j = col_1; j < m; j++)
    {
      A(i, j) += s2 * A(row_start + 1, j);
    }
    for (size_t j = 0; j < n; j++)
    {
      Q(i, j) += s1 * Q(row_start, j) + s2 * Q(row_start + 1, j);
    }
  }

  return std::make_tuple(Q, C, A);
}

std::pair<std::vector<int>, std::vector<int>> RingMatrix::Conditioning(RingMatrix A, size_t col_1, size_t col_2,
                                                                       size_t row_start)
{
  size_t k = A._rows - row_start - 2;
  int d11 = A(row_start, col_1);
  int d12 = A(row_start, col_2);
  int d21 = A(row_start + 1, col_1);
  int d22 = A(row_start + 1, col_2);
  std::vector<int> ci = std::vector<int>(k, 0);  //[Int](repeating: 0, count : k)
  if (d11 * d22 == d12 * d21)
  {
    for (size_t s = (row_start + 2); s < A._rows; s++)
    {
      if (d11 * A(s, col_2) != d12 * A(s, col_1))
      {
        ci[s - row_start - 2] = 1;
        d21 += A(s, col_1);
        d22 += A(s, col_2);
        break;
      }
    }
  }

  // We now have  det( d11 & d12 \\ d21 & d22 ) \neq 0.
  if (d11 > 1)
  {
    // Perform a modified Algorithm 6.15:
    std::vector<int> F = {d11};
    int ahat = d21;
    size_t i = 0;
    bool has_gi = false;
    bool neg = false;

    int biprime = 0;
    int ahatprime = 0;
    int bi = 0;
    int bip = 0;
    while (i < k)
    {
      if (!has_gi)
      {
        bi = A(row_start + 2 + i, col_1);
        bip = A(row_start + 2 + i, col_2);
        int gi = Ring::greatestCommonDivisor(ahat, bi);
        if (gi == 0)
        {
          i += 1;
          continue;
        }
        ahatprime = Ring::modulo(ahat / gi, d11);
        biprime = Ring::modulo(bi / gi, d11);
        neg = false;
        if (sign(d11 * d22 - d12 * d21) != sign(d11 * bip - d12 * bi))
        {
          biprime = -biprime;
          neg = true;
        }
        has_gi = true;
      }

      std::variant<std::vector<int>, int> res = Algorithm_6_14(ahatprime, biprime, d11, F);

      if (res.index() == 0)
      {
        F = std::get<0>(res);
        std::sort(F.begin(), F.end());
        F = RemovedDuplicates(F);
      }
      else
      {
        int convertedRes = std::get<1>(res);
        if (neg)
        {
          ci[i] -= convertedRes;
          d21 -= convertedRes * bi;
          d22 -= convertedRes * bip;
          ahat = Ring::modulo((ahat - convertedRes * bi), d11);
        }
        else
        {
          ci[i] += convertedRes;
          d21 += convertedRes * bi;
          d22 += convertedRes * bip;
          ahat = Ring::modulo((ahat + convertedRes * bi), d11);
        }
        has_gi = false;
        i += 1;
      }
    }
  }

  return std::make_pair(std::vector<int>{d21, d22}, ci);
}

const int RingMatrix::HNF_C_Iwaniec = 3;

std::variant<std::vector<int>, int> RingMatrix::Algorithm_6_14(int a, int b, int N, std::vector<int> Nfact)
{
  size_t k = 0;
  int HNF_C_IwaniecLocal = RingMatrix::HNF_C_Iwaniec;

  if (N <= 1)
  {
    // qDebug() << "STOP ERROR";
    // assert(false);
    return {};
  }

  while (true)
  {
    if (N == 2)
    {
      k = 1;
    }
    else
    {
      double temp = log(double(N)) / log(2.0);
      k = static_cast<size_t>(static_cast<double>(HNF_C_IwaniecLocal) * temp * (std::pow(std::log(temp), 2)));
    }

    // Prepare B
    std::vector<bool> B = std::vector<bool>(k + 1, true);  //   [Bool](repeating: true, count: k+1)

    // Compute residues
    size_t t = Nfact.size();
    std::vector<int> ai = std::vector<int>(t, 0);  // [Int](repeating: 0, count: t)
    std::vector<int> bi = std::vector<int>(t, 0);  // [Int](repeating: 0, count: t)

    for (size_t i = 0; i < t; i++)
    {
      ai[i] = Ring::modulo(a, Nfact[i]);
      bi[i] = Ring::modulo(b, Nfact[i]);
    }

    // Compute extended GCDs
    std::vector<int> xi = std::vector<int>(t, 0);  //[Int](repeating: 0, count: t)
    for (size_t i = 0; i < t; i++)
    {
      // (gi: Int, xi: Int, yi: Int)
      std::tuple<int, int, int> extendedGCD = Ring::extendedGreatestCommonDivisor(bi[i], Nfact[i]);
      xi[i] = std::get<1>(extendedGCD);
      if ((1 < std::get<0>(extendedGCD)) && (std::get<0>(extendedGCD) < Nfact[i]))
      {
        std::vector<int> res{};
        for (size_t l = 0; l < i; l++) res.push_back(Nfact[l]);
        res.push_back(std::get<0>(extendedGCD));
        res.push_back(Nfact[i] / std::get<0>(extendedGCD));
        for (size_t l = i + 1; l < t; l++) res.push_back(Nfact[l]);
        return res;
      }
    }

    // Do sieving
    for (size_t i = 0; i < t; i++)
    {
      if (bi[i] != 0)
      {
        int si = Ring::modulo(-ai[i] * xi[i], Nfact[i]);
        size_t idx = static_cast<size_t>(si);
        while (idx <= k)
        {
          B[idx] = false;
          idx += static_cast<size_t>(Nfact[i]);
        }
      }
    }
    // Find result

    for (size_t c = 0; c < (k + 1); c++)
    {
      if (B[c] == true)
      {
        for (size_t i = 0; i < t; i++)
        {
          int gi = Ring::greatestCommonDivisor(ai[i] + static_cast<int>(c) * bi[i], Nfact[i]);

          if (gi > 1)
          {
            std::vector<int> res{};
            for (size_t l = 0; l < i; l++) res.push_back(Nfact[l]);
            res.push_back(gi);
            res.push_back(Nfact[i] / gi);
            for (size_t l = (i + 1); l < t; l++) res.push_back(Nfact[l]);

            return res;
          }
        }
        return int(c);
      }
    }
    HNF_C_IwaniecLocal *= 2;
  }
  // qDebug() << "Should not get here";
  // assert(false);
}

std::vector<int> RingMatrix::RemovedDuplicates(std::vector<int> array)
{
  std::vector<int> res{};

  if (array.size() > 0)
  {
    res.push_back(array[0]);
    for (size_t i = 1; i < array.size(); i++)
    {
      if (array[i] != res[res.size() - 1])
      {
        res.push_back(array[i]);
      }
    }
  }
  return res;
}

// Smith Normal Form
// ====================================================================================

std::tuple<RingMatrix, RingMatrix, RingMatrix> RingMatrix::SmithNormalForm()
{
  size_t n = this->_rows;
  size_t m = this->_columns;

  // (U: RingMatrix, A: RingMatrix, rp: [Int])
  std::tuple<RingMatrix, RingMatrix, std::vector<size_t>> hnf = HermiteNormalForm();

  RingMatrix U = std::get<0>(hnf);
  RingMatrix A = std::get<1>(hnf);
  std::vector rp = std::get<2>(hnf);
  size_t r = rp.size();

  RingMatrix V = RingMatrix::identity(m);

  // Transform A via V so that the left r x r block of A is invertible
  for (size_t i = 0; i < r; i++)
  {
    if (rp[i] > i)
    {
      A.swapColumns(i, rp[i]);
      V.swapColumns(i, rp[i]);
    }
  }

  // Phase one
  for (size_t i = 0; i < r; i++)
  {
    Smith_Theorem5(A, U, V, i);
  }

  size_t beg = 0;
  while (beg < r && A(beg, beg) == 1)
  {
    beg += 1;
  }

  // Phase two
  if (beg < r && r < m)
  {
    for (size_t i = beg; i < r; i++)
    {
      Smith_Theorem8(A, U, V, i, r);
    }

    // Run transposed Phase One
    RingMatrix AA = A.submatrix(beg, beg, r - beg, r - beg).transposed();
    RingMatrix UU = RingMatrix::identity(r - beg);
    RingMatrix VV = RingMatrix::identity(r - beg);
    // Check if it is actually not a diagonal matrix
    for (size_t i = 0; i < (r - beg); i++)
    {
      Smith_Theorem5(AA, UU, VV, i);
    }

    // Insert AA
    AA = AA.transposed();
    A.assignSubmatrix(beg, beg, AA.transposed());

    // Insert transformations
    RingMatrix temp = UU;
    UU = VV.transposed();
    VV = temp.transposed();

    U.assignSubmatrix(beg, 0, UU * U.submatrix(beg, 0, r - beg, n));
    V.assignSubmatrix(0, beg, V.submatrix(0, beg, m, r - beg) * VV);
  }

  // V.denominator = self.denominator
  return std::make_tuple(U, V, A);
}

void RingMatrix::Smith_Theorem5(RingMatrix& A, RingMatrix& U, RingMatrix& V, size_t col)
{
  size_t n = A._rows;
  size_t m = A._columns;

  // Lemma 6:
  for (size_t i = col; i-- > 0;)  // note: starts at col-1
  {
    // Compute ci[0] such that GCD(A[i, col] + ci[0] * A[i + 1, col], A[i, i])
    // equals GCD(A[i, col], A[i + 1, col], A[i, i])
    std::vector<int> ci = Algorithm_6_15(A(i, col), {A(i + 1, col)}, A(i, i));

    // Add ci[0] times the (i+1)-th row to the i-th row
    for (size_t j = 0; j < m; j++)
    {
      A(i, j) = A(i, j) + ci[0] * A(i + 1, j);
    }
    for (size_t j = 0; j < n; j++)
    {
      U(i, j) = U(i, j) + ci[0] * U(i + 1, j);
    }

    // Reduce i-th row modulo A[i, i]
    for (size_t j = (i + 1); j < m; j++)
    {
      std::pair<int, int> divmod = Ring::divisionModulo(A(i, j), A(i, i));
      if (divmod.first != 0)
      {
        // Subtract d times the i-th column from the j-th column
        A(i, j) = divmod.second;
        for (size_t k = 0; k < m; k++)
        {
          V(k, j) = V(k, j) - divmod.first * V(k, i);
        }
      }
    }
  }

  // Lemma 7
  for (size_t j = 0; j < col; j++)
  {
    // Apply lemma 7 to submatrix starting at (j, j)
    std::tuple<int, int, int> extendedGCD = Ring::extendedGreatestCommonDivisor(A(j, j), A(j, col));
    int ss = -A(j, col) / std::get<0>(extendedGCD);
    int tt = A(j, j) / std::get<0>(extendedGCD);
    // Transform columns j and col by a 2x2 matrix
    A(j, j) = std::get<0>(extendedGCD);
    A(j, col) = 0;
    for (size_t i = (j + 1); i < n; i++)
    {
      int temp = A(i, j);
      A(i, j) = std::get<1>(extendedGCD) * A(i, j) + std::get<2>(extendedGCD) * A(i, col);
      A(i, col) = ss * temp + tt * A(i, col);
    }
    for (size_t i = 0; i < m; i++)
    {
      int temp = V(i, j);
      V(i, j) = std::get<1>(extendedGCD) * V(i, j) + std::get<2>(extendedGCD) * V(i, col);
      V(i, col) = ss * temp + tt * V(i, col);
    }

    // Clear column j in rows below
    for (size_t i = (j + 1); i < n; i++)
    {
      int mul = A(i, j) / A(j, j);
      if (mul != 0)
      {
        for (size_t jj = 0; jj < m; jj++)
        {
          A(i, jj) = A(i, jj) - mul * A(j, jj);
        }
        for (size_t jj = 0; jj < n; jj++)
        {
          U(i, jj) = U(i, jj) - mul * U(j, jj);
        }
      }
    }

    // Reduce j-th row modulo A[j, j]
    for (size_t jj = (j + 1); jj < m; jj++)
    {
      std::pair<int, int> divmod = Ring::divisionModulo(A(j, jj), A(j, j));
      if (divmod.first != 0)
      {
        // Subtract d times the i-th column from the j-th column
        A(j, jj) = divmod.second;
        for (size_t k = 0; k < m; k++)
        {
          V(k, jj) = V(k, jj) - divmod.first * V(k, j);
        }
      }
    }
  }

  // Make A[col, col] positive
  if (A(col, col) < 0)
  {
    for (size_t jj = col; jj < m; jj++)
    {
      A(col, jj) = -A(col, jj);
    }
    for (size_t jj = 0; jj < n; jj++)
    {
      U(col, jj) = -U(col, jj);
    }
  }

  // Reduce col-th row modulo A[col, col]
  for (size_t j = (col + 1); j < m; j++)
  {
    std::pair<int, int> divmod = Ring::divisionModulo(A(col, j), A(col, col));
    if (divmod.first != 0)
    {
      // Subtract d times the col-th column from the j-th column
      A(col, j) = divmod.second;
      for (size_t k = 0; k < m; k++)
      {
        V(k, j) = V(k, j) - divmod.first * V(k, col);
      }
    }
  }
}

void RingMatrix::Smith_Theorem8(RingMatrix& A, RingMatrix& U, RingMatrix& V, size_t row, size_t r)
{
  size_t n = A._rows;
  size_t m = A._columns;

  for (size_t j = r; j < m; j++)
  {
    if (A(row, j) != 0)
    {
      std::tuple<int, int, int> extendedGCD = Ring::extendedGreatestCommonDivisor(A(row, row), A(row, j));
      int ss = -A(row, j) / std::get<0>(extendedGCD);
      int tt = A(row, row) / std::get<0>(extendedGCD);
      // Transform columns row and j by a 2x2 matrix
      A(row, row) = std::get<0>(extendedGCD);
      A(row, j) = 0;

      for (size_t i = (row + 1); i < n; i++)
      {
        int temp = A(i, row);
        A(i, row) = std::get<1>(extendedGCD) * A(i, row) + std::get<2>(extendedGCD) * A(i, j);
        A(i, j) = ss * temp + tt * A(i, j);
      }
      for (size_t i = 0; i < m; i++)
      {
        int temp = V(i, row);
        V(i, row) = std::get<1>(extendedGCD) * V(i, row) + std::get<2>(extendedGCD) * V(i, j);
        V(i, j) = ss * temp + tt * V(i, j);
      }

      // Reduce column row
      for (size_t i = (row + 1); i < n; i++)
      {
        int d = Ring::floorDivision(A(i, row), A(row, row));
        if (d != 0)
        {
          for (size_t jj = 0; jj < m; jj++)
          {
            A(i, jj) = A(i, jj) - d * A(row, jj);
          }
          for (size_t jj = 0; jj < n; jj++)
          {
            U(i, jj) = U(i, jj) - d * U(row, jj);
          }
        }
      }
      // Reduce column row
      for (size_t i = (row + 1); i < n; i++)
      {
        int d = Ring::floorDivision(A(i, row), A(row, row));
        if (d != 0)
        {
          for (size_t jj = 0; jj < m; jj++)
          {
            A(i, jj) = A(i, jj) - d * A(row, jj);
          }
          for (size_t jj = 0; jj < n; jj++)
          {
            U(i, jj) = U(i, jj) - d * U(row, jj);
          }
        }
      }
    }
  }
}

std::vector<int> RingMatrix::Algorithm_6_15(int a, std::vector<int> bi, int N)
{
  if (N == 1)
  {
    return std::vector<int>(bi.size(), 0);
  }
  std::vector<int> F = {N};
  int ahat = a;
  size_t i = 0;
  size_t n = bi.size();
  bool has_gi = false;
  std::vector<int> ci(n, 0);

  int ahatprime = 0;
  int biprime = 0;

  while (i < n)
  {
    if (!has_gi)
    {
      int gi = Ring::greatestCommonDivisor(ahat, bi[i]);
      if (gi == 0)
      {
        // both ahat and bi[i] are zero: take 0 as ci[i] and continue with next entry
        ci[i] = 0;
        i += 1;
        continue;
      }
      ahatprime = (ahat / gi) % N;
      biprime = (bi[i] / gi) % N;
      has_gi = true;
    }
    std::variant<std::vector<int>, int> res = Algorithm_6_14(ahatprime, biprime, N, F);

    if (res.index() == 0)
    {
      F = std::get<0>(res);
      std::sort(F.begin(), F.end());
      F = RemovedDuplicates(F);
    }
    else
    {
      int convertedRes = std::get<1>(res);
      ci[i] = convertedRes;
      ahat = (ahat + convertedRes * bi[i]) % N;
      i += 1;
      has_gi = false;
    }
  }
  return ci;
}
