module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <set>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <vector>
#endif

export module skintegersymmetryoperationset;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import int3;
import int3x3;
import double3;
import double3x3;
import skdefinitions;
import skseitzintegermatrix;

export struct SKIntegerSymmetryOperationSet
{
  std::unordered_set<SKSeitzIntegerMatrix, SKSeitzIntegerMatrix::hashFunction> operations;
  Centring centring;

  SKIntegerSymmetryOperationSet();
  // SKIntegerSymmetryOperationSet(std::unordered_set<SKSeitzIntegerMatrix> &operations);
  SKIntegerSymmetryOperationSet(std::vector<SKSeitzIntegerMatrix> operations);

  inline std::size_t size() { return operations.size(); }
  SKIntegerSymmetryOperationSet fullSeitzMatrices();

  std::vector<std::tuple<double3, std::size_t, double>> symmetrize(
      double3x3 lattice, std::vector<std::tuple<double3, std::size_t, double>> atoms, double symmetryPrecision);
  std::vector<std::tuple<double3, std::size_t, double>> asymmetricAtoms(
      std::size_t HallNumber, std::vector<std::tuple<double3, std::size_t, double>>& atoms, double3x3 lattice,
      bool allowPartialOccupancies, double symmetryPrecision);
};
