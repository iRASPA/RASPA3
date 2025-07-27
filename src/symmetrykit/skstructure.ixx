module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <vector>
#endif

export module skstructure;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import int3;
import double3;
import skcell;
import skatom;

export class SKStructure
{
 public:
  enum class DataType
  {
    Uint8,
    Int8,
    Uint16,
    Int16,
    Uint32,
    Int32,
    Uint64,
    Int64,
    Float,
    Double
  };

  SKStructure() : cell(std::make_shared<SKCell>()) {};

  enum class Kind : std::int64_t
  {
    none = -1,
    object = 0,
    structure = 1,
    crystal = 2,
    molecularCrystal = 3,
    molecule = 4,
    protein = 5,
    proteinCrystal = 6,
    proteinCrystalSolvent = 7,
    crystalSolvent = 8,
    molecularCrystalSolvent = 9,
    crystalEllipsoidPrimitive = 10,
    crystalCylinderPrimitive = 11,
    crystalPolygonalPrismPrimitive = 12,
    ellipsoidPrimitive = 13,
    cylinderPrimitive = 14,
    polygonalPrismPrimitive = 15,
    gridVolume = 16,
    RASPADensityVolume = 17,
    VTKDensityVolume = 18,
    VASPDensityVolume = 19,
    GaussianCubeVolume = 20
  };

  Kind kind = Kind::crystal;
  std::vector<std::shared_ptr<SKAsymmetricAtom>> atoms;
  std::set<std::string> unknownAtoms;

  std::optional<std::string> displayName;
  std::shared_ptr<SKCell> cell;
  std::optional<int> spaceGroupHallNumber;
  bool drawUnitCell = false;
  bool periodic = false;

  std::optional<std::string> creationDate;
  std::optional<std::string> creationMethod;
  std::optional<std::string> chemicalFormulaSum;
  std::optional<std::string> chemicalFormulaStructural;
  std::optional<int> cellFormulaUnitsZ;

  std::optional<int> numberOfChannels;
  std::optional<int> numberOfPockets;
  std::optional<int> dimensionality;
  std::optional<double> Di;
  std::optional<double> Df;
  std::optional<double> Dif;

  int3 dimensions;
  double3 origin;
  double3 spacing;
  DataType dataType;
  std::pair<double, double> range;
  double average;
  double variance;
  std::vector<char> byteData;
};
