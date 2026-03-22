module;

export module sksymmetryoperationset;

import std;

import skdefinitions;
import skseitzmatrix;
import skrotationmatrix;
import sktransformationmatrix;
import sksymmetrycell;

export struct SKSymmetryOperationSet
{
  std::vector<SKSeitzMatrix> operations;
  Centring centring;

  SKSymmetryOperationSet();
  SKSymmetryOperationSet(std::vector<SKSeitzMatrix> operations);

  inline std::size_t size() { return this->operations.size(); }

  SKSymmetryOperationSet fullSeitzMatrices();
  const std::vector<SKRotationMatrix> rotations() const;
  const std::vector<SKRotationMatrix> properRotations() const;

  const SKSymmetryOperationSet changedBasis(SKTransformationMatrix transformationMatrix) const;
  const SKSymmetryOperationSet addingCenteringOperations(Centring centering) const;
};
