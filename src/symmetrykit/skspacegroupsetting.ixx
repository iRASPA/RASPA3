module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <cstdint>
#include <ostream>
#include <string>
#include <vector>
#endif

export module skspacegroupsetting;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import int3;
import int3x3;
import double3;
import double3x3;
import skdefinitions;
import skasymmetricunit;
import sktransformationmatrix;
import skintegersymmetryoperationset;
import skseitzintegermatrix;

export class SKSpaceGroupSetting
{
 public:
  // SKSpaceGroupSetting() {}
  SKSpaceGroupSetting(std::size_t number, std::size_t spaceGroupNumber, std::size_t order, char ext,
                      std::string qualifier, std::string HM, std::string Hall, bool inversionAtOrigin,
                      int3 inversionCenter, Symmorphicity symmorphicity, bool standard, Centring centring,
                      std::vector<int3> latticeTranslations, std::size_t pointGroupNumber, std::string schoenflies,
                      std::string generators, std::string encoding, SKAsymmetricUnit asymmetricUnit,
                      SKTransformationMatrix transformationMatrix);

  SKIntegerSymmetryOperationSet fullSeitzMatrices() const;
  std::vector<SKSeitzIntegerMatrix> SeitzMatricesWithoutTranslation() const;

  std::size_t number() const { return _spaceGroupNumber; }
  std::size_t HallNumber() const { return _HallNumber; }
  std::string HallString() const { return _HallString; }
  std::string HMString() const { return _HMString; }
  std::size_t pointGroupNumber() const { return _pointGroupNumber; }
  std::string qualifier() const { return _qualifier; }
  Symmorphicity symmorphicity() const { return _symmorphicity; }
  std::string symmorphicityString() const;
  std::string centringString() const;

  const std::string encodedGenerators() const { return _encodedGenerators; }

  bool inversionAtOrigin() const { return _inversionAtOrigin; }
  int3 inversionCenter() const { return _inversionCenter; }

  const std::vector<int3> latticeTranslations() const { return _latticeTranslations; }
  Centring centring() const { return _centring; }

  // check
  SKAsymmetricUnit asymmetricUnit() const { return _asymmetricUnit; }

  friend std::ostream& operator<<(std::ostream& os, const SKSpaceGroupSetting& setting);

 private:
  std::size_t _HallNumber = 1;
  std::size_t _spaceGroupNumber = 1;  // space group number (1-230)
  std::size_t _order;
  char _ext;                       // '1', '2', 'H', 'R' or '\0'
  std::string _qualifier;          // e.g. "-cba" or "b1"
  std::string _HMString;           // H-M symbol; nul-terminated string
  std::string _HallString;         // Hall symbol; nul-terminated string
  std::string _encodedGenerators;  // encoded seitz matrix-generators
  std::string _encodedSeitz;       // encoded seitz matrix
  bool _inversionAtOrigin;
  int3 _inversionCenter;
  bool _standard = false;
  Symmorphicity _symmorphicity = Symmorphicity::asymmorphic;
  Centring _centring = Centring::primitive;
  std::vector<int3> _latticeTranslations;
  std::string _schoenflies;
  std::size_t _pointGroupNumber;
  SKAsymmetricUnit _asymmetricUnit;              // = {{10,0}, {20,1}, {30,2}};
  SKTransformationMatrix _transformationMatrix;  // the inverse of the transformation to "standard" setting (so:
                                                 // standard to unconventional setting)
};
