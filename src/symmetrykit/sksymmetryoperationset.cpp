module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <vector>
#endif

module sksymmetryoperationset;

#ifndef USE_LEGACY_HEADERS
import <algorithm>;
import <vector>;
import <iterator>;
#endif

import double3;
import int3x3;

import skdefinitions;
import skseitzmatrix;
import skrotationmatrix;
import sktransformationmatrix;

SKSymmetryOperationSet::SKSymmetryOperationSet() {}

SKSymmetryOperationSet::SKSymmetryOperationSet(std::vector<SKSeitzMatrix> operations) { this->operations = operations; }

const std::vector<SKRotationMatrix> SKSymmetryOperationSet::rotations() const
{
  std::vector<SKRotationMatrix> rotationMatrices{};

  std::transform(this->operations.begin(), this->operations.end(), std::back_inserter(rotationMatrices),
                 [](SKSeitzMatrix m) { return m.rotation; });

  return rotationMatrices;
}

const std::vector<SKRotationMatrix> SKSymmetryOperationSet::properRotations() const
{
  std::vector<SKRotationMatrix> rotationMatrices{};

  std::transform(this->operations.begin(), this->operations.end(), std::back_inserter(rotationMatrices),
                 [](SKSeitzMatrix m) { return m.rotation.proper(); });

  return rotationMatrices;
}
// the transformationMatrix does not have a translational part
const SKSymmetryOperationSet SKSymmetryOperationSet::changedBasis(SKTransformationMatrix transformationMatrix) const
{
  std::vector<SKSeitzMatrix> newSet{};

  for (const SKSeitzMatrix& seitzMatrix : this->operations)
  {
    SKTransformationMatrix inverseTransformation = transformationMatrix.adjugate();
    SKTransformationMatrix rotation =
        (inverseTransformation * seitzMatrix.rotation * transformationMatrix) / transformationMatrix.determinant();
    double3 translation = inverseTransformation * seitzMatrix.translation;
    SKSeitzMatrix transformedSeitzMatrix = SKSeitzMatrix(SKRotationMatrix(rotation.transformation),
                                                         translation / double(transformationMatrix.determinant()));
    newSet.push_back(transformedSeitzMatrix);
  }
  return SKSymmetryOperationSet(newSet);
}

const SKSymmetryOperationSet SKSymmetryOperationSet::addingCenteringOperations(Centring centering) const
{
  std::vector<double3> shifts{};

  switch (centering)
  {
    case Centring::none:
    case Centring::primitive:
      shifts = {double3(0, 0, 0)};
      break;
    case Centring::face:
      shifts = {double3(0, 0, 0), double3(0, 0.5, 0.5), double3(0.5, 0, 0.5), double3(0.5, 0.5, 0)};
      break;
    case Centring::r:
      shifts = {double3(0, 0, 0), double3(8.0 / 12.0, 4.0 / 12.0, 4.0 / 12.0),
                double3(4.0 / 12.0, 8.0 / 12.0, 8.0 / 12.0)};
      break;
    case Centring::h:
      shifts = {double3(0, 0, 0), double3(8.0 / 12.0, 4.0 / 12.0, 0), double3(0, 8.0 / 12.0, 4.0 / 12.0)};
      break;
    case Centring::d:
      shifts = {double3(0, 0, 0), double3(4.0 / 12.0, 4.0 / 12.0, 4.0 / 12.0),
                double3(8.0 / 12.0, 8.0 / 12.0, 8.0 / 12.0)};
      break;
    case Centring::body:
      shifts = {double3(0, 0, 0), double3(0.5, 0.5, 0.5)};
      break;
    case Centring::a_face:
      shifts = {double3(0, 0, 0), double3(0, 0.5, 0.5)};
      break;
    case Centring::b_face:
      shifts = {double3(0, 0, 0), double3(0.5, 0, 0.5)};
      break;
    case Centring::c_face:
      shifts = {double3(0, 0, 0), double3(0.5, 0.5, 0)};
      break;
    default:
      shifts = {double3(0, 0, 0)};
      break;
  }

  std::vector<SKSeitzMatrix> symmetry{};

  for (const SKSeitzMatrix& seitzMatrix : this->operations)
  {
    for (const double3& shift : shifts)
    {
      symmetry.push_back(SKSeitzMatrix(seitzMatrix.rotation, seitzMatrix.translation + shift));
    }
  }
  return SKSymmetryOperationSet(symmetry);
}
