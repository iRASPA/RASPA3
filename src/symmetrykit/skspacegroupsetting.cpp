module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <iostream>
#include <vector>
#endif

module skspacegroupsetting;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import int3;
import int3x3;
import double3;
import double3x3;

import skdefinitions;
import skasymmetricunit;
import skrotationmatrix;
import sktransformationmatrix;
import skintegersymmetryoperationset;
import skseitzintegermatrix;
import skpointgroup;

SKSpaceGroupSetting::SKSpaceGroupSetting(std::size_t number, std::size_t spaceGroupNumber, std::size_t order, char ext,
                                         std::string qualifier, std::string HM, std::string Hall,
                                         bool inversionAtOrigin, int3 inversionCenter, Symmorphicity symmorphicity,
                                         bool standard, Centring centring, std::vector<int3> latticeTranslations,
                                         std::size_t pointGroupNumber, std::string schoenflies, std::string generators,
                                         std::string encoding, SKAsymmetricUnit asymmetricUnit,
                                         SKTransformationMatrix transformationMatrix)
{
  _HallNumber = number;
  _spaceGroupNumber = spaceGroupNumber;
  _order = order;
  _ext = ext;
  _qualifier = qualifier;
  _HMString = HM;
  _HallString = Hall;
  _encodedGenerators = generators;
  _encodedSeitz = encoding;
  _inversionAtOrigin = inversionAtOrigin;
  _inversionCenter = inversionCenter;
  _standard = standard;
  _symmorphicity = symmorphicity;
  _centring = centring;
  _latticeTranslations = latticeTranslations;
  _schoenflies = schoenflies;
  _pointGroupNumber = pointGroupNumber;
  _asymmetricUnit = asymmetricUnit;
  _transformationMatrix = transformationMatrix;
}

std::string SKSpaceGroupSetting::symmorphicityString() const
{
  switch (_symmorphicity)
  {
    case Symmorphicity::asymmorphic:
      return "asymmorphic";
    case Symmorphicity::symmorphic:
      return "symmorphic";
    case Symmorphicity::hemisymmorphic:
      return "hemisymmorphic";
  }
  return std::string();
}

std::string SKSpaceGroupSetting::centringString() const
{
  switch (_centring)
  {
    case Centring::none:
      return "none";
    case Centring::primitive:
      return "primitive";
    case Centring::body:
      return "body";
    case Centring::a_face:
      return "a";
    case Centring::b_face:
      return "b";
    case Centring::c_face:
      return "c";
    case Centring::face:
      return "face";
    case Centring::base:
      return "base";
    case Centring::r:
      return "r";
    case Centring::h:
      return "h";
    case Centring::d:
      return "d";
  }
  return std::string();
}

std::ostream& operator<<(std::ostream& os, const SKSpaceGroupSetting& setting)
{
  os << "SpaceGroupSetting: " << setting._HallNumber << '/' << setting._spaceGroupNumber << '/'
     << setting._encodedSeitz;
  return os;
}

SKIntegerSymmetryOperationSet SKSpaceGroupSetting::fullSeitzMatrices() const
{
  // assert(_encodedSeitz.size() % 3 == 0);
  // assert(_encodedSeitz.size() > 0);

  bool centrosymmetric = SKPointGroup::pointGroupData[static_cast<std::size_t>(_pointGroupNumber)].centrosymmetric();
  std::size_t m = _encodedSeitz.size() / 3;

  std::size_t size = centrosymmetric ? 2 * m : m;
  std::vector<int3> translationVectors = _latticeTranslations;
  std::vector<SKSeitzIntegerMatrix> matrices = std::vector<SKSeitzIntegerMatrix>();
  matrices.resize(size * translationVectors.size());

  for (std::size_t i = 0; i < m; i++)
  {
    char x = _encodedSeitz[3 * i];
    char y = _encodedSeitz[3 * i + 1];
    char z = _encodedSeitz[3 * i + 2];

    matrices[i] = SKSeitzIntegerMatrix(x, y, z);
  }

  if (centrosymmetric)
  {
    for (std::size_t i = 0; i < m; i++)
    {
      char x = _encodedSeitz[3 * i];
      char y = _encodedSeitz[3 * i + 1];
      char z = _encodedSeitz[3 * i + 2];

      SKSeitzIntegerMatrix seitz = SKSeitzIntegerMatrix(x, y, z);

      int3 translation = seitz.translation + seitz.rotation * _inversionCenter;
      matrices[m + i] = SKSeitzIntegerMatrix(-seitz.rotation, translation);
    }
  }

  // use the translation vectors on all Seitz matrices
  for (std::size_t k = 1; k < translationVectors.size(); k++)
  {
    for (std::size_t i = 0; i < size; i++)
    {
      matrices[size * k + i] = matrices[i];
      matrices[size * k + i].translation = matrices[size * k + i].translation + translationVectors[k];
    }
  }

  return SKIntegerSymmetryOperationSet(matrices);
}

std::vector<SKSeitzIntegerMatrix> SKSpaceGroupSetting::SeitzMatricesWithoutTranslation() const
{
  // assert(_encodedSeitz.size() % 3 == 0);
  // assert(_encodedSeitz.size() > 0);

  bool centrosymmetric = SKPointGroup::pointGroupData[static_cast<std::size_t>(_pointGroupNumber)].centrosymmetric();
  std::size_t m = _encodedSeitz.size() / 3;

  std::size_t size = centrosymmetric ? 2 * m : m;
  std::vector<int3> translationVectors = _latticeTranslations;
  std::vector<SKSeitzIntegerMatrix> matrices = std::vector<SKSeitzIntegerMatrix>();
  matrices.resize(size);

  for (std::size_t i = 0; i < m; i++)
  {
    char x = _encodedSeitz[3 * i];
    char y = _encodedSeitz[3 * i + 1];
    char z = _encodedSeitz[3 * i + 2];

    matrices[i] = SKSeitzIntegerMatrix(x, y, z);
  }

  if (centrosymmetric)
  {
    for (std::size_t i = 0; i < m; i++)
    {
      char x = _encodedSeitz[3 * i];
      char y = _encodedSeitz[3 * i + 1];
      char z = _encodedSeitz[3 * i + 2];

      SKSeitzIntegerMatrix seitz = SKSeitzIntegerMatrix(x, y, z);

      int3 translation = seitz.translation + seitz.rotation * _inversionCenter;
      matrices[m + i] = SKSeitzIntegerMatrix(-seitz.rotation, translation);
    }
  }

  return matrices;
}
