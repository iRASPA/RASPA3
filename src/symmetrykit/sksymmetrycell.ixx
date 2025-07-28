module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <vector>
#endif

export module sksymmetrycell;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import int3;
import int3x3;
import double3;
import double3x3;
import skdefinitions;
import skrotationmatrix;
import sktransformationmatrix;
import skpointsymmetryset;

export class SKSymmetryCell
{
 public:
  SKSymmetryCell();
  SKSymmetryCell(double a, double b, double c, double alpha, double beta, double gamma);
  static SKSymmetryCell createFromMetricTensor(double3x3 metricTensor);
  static SKSymmetryCell createFromUnitCell(double3x3 unitCell);
  SKSymmetryCell idealized(std::size_t pointGroupNumber, std::string qualifier);
  double3x3 unitCell() const;
  double3x3 metricTensor();
  double volume();
  static double3x3 findSmallestPrimitiveCell(std::vector<std::tuple<double3, std::size_t, double>> reducedAtoms,
                                             std::vector<std::tuple<double3, std::size_t, double>> atoms,
                                             double3x3 unitCell, bool allowPartialOccupancies,
                                             double symmetryPrecision);
  static bool testTranslationalSymmetry(double3 translationVector,
                                        std::vector<std::tuple<double3, std::size_t, double>> atoms, double3x3 unitCell,
                                        bool allowPartialOccupancies, double precision);
  static bool testSymmetry(double3 translationVector, SKRotationMatrix rotationMatrix,
                           std::vector<std::tuple<double3, std::size_t, double>> atoms, double3x3 unitCell,
                           bool allowPartialOccupancies, double precision);
  static std::vector<double3> primitiveTranslationVectors(
      double3x3 unitCell, std::vector<std::tuple<double3, std::size_t, double>> reducedAtoms,
      std::vector<std::tuple<double3, std::size_t, double>> atoms, SKRotationMatrix rotationMatrix,
      bool allowPartialOccupancies, double symmetryPrecision);
  static std::optional<double3x3> computeDelaunayReducedCell(double3x3 unitCell, double symmetryPrecision);
  static std::optional<double3x3> computeDelaunayReducedCell2D(double3x3 unitCell, double symmetryPrecision);
  static SKPointSymmetrySet findLatticeSymmetry(double3x3 reducedLattice, double symmetryPrecision);
  static bool checkMetricSimilarity(double3x3 transformedMetricTensor, double3x3 metricTensor,
                                    double symmetryPrecision);
  static std::vector<std::tuple<double3, std::size_t, double>> trim(
      std::vector<std::tuple<double3, std::size_t, double>> atoms, double3x3 from, double3x3 to,
      bool allowPartialOccupancies, double symmetryPrecision);
  std::optional<std::pair<SKSymmetryCell, SKTransformationMatrix>> computeReducedNiggliCellAndChangeOfBasisMatrix();
  static bool isOverlap(double3 a, double3 b, double3x3 lattice, double symmetryPrecision);

  double a() const { return _a; }
  double b() const { return _b; }
  double c() const { return _c; }
  double alpha() const { return _alpha; }
  double beta() const { return _beta; }
  double gamma() const { return _gamma; }

 private:
  double _a;
  double _b;
  double _c;
  double _alpha;
  double _beta;
  double _gamma;

  static bool isSmallerThen(double x, double y);
  static bool isLargerThen(double x, double y);
  static bool isSmallerEqualThen(double x, double y);
  static bool isLargerEqualThen(double x, double y);
  static bool isEqualTo(double x, double y);
  static bool isLargerThenZeroXiEtaZeta(double xi, double eta, double zeta);
};
