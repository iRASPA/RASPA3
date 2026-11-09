module;

module skspacegroupsetting;

import std;

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
  return SKIntegerSymmetryOperationSet(SKSeitzIntegerMatrix::SeitzMatrices(_encodedSeitz));
}

// The rotation parts of the setting's operations, each once and with a zero translation: the point-group
// operations in the setting's basis (48 for F d -3 m, whose 192 operations are these times the 4 centring
// translations). It used to return the full operations, translations included, under this name.
std::vector<SKSeitzIntegerMatrix> SKSpaceGroupSetting::SeitzMatricesWithoutTranslation() const
{
  std::vector<SKSeitzIntegerMatrix> rotations{};
  for (const SKSeitzIntegerMatrix& operation : SKSeitzIntegerMatrix::SeitzMatrices(_encodedSeitz))
  {
    const SKSeitzIntegerMatrix rotation(operation.rotation, int3(0, 0, 0));
    if (std::find(rotations.begin(), rotations.end(), rotation) == rotations.end()) rotations.push_back(rotation);
  }
  return rotations;
}
