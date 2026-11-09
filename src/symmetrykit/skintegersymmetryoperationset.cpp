module;

module skintegersymmetryoperationset;

import std;

import int3;
import int3x3;
import double3;
import double3x3;

import sksymmetrycell;
import skrotationmatrix;
import skseitzintegermatrix;
import skasymmetricunit;

SKIntegerSymmetryOperationSet::SKIntegerSymmetryOperationSet() {}

// SKIntegerSymmetryOperationSet::SKIntegerSymmetryOperationSet(std::unordered_set<SKSeitzIntegerMatrix> &operations)
//{
//     this->operations = operations;
//     this->centring = Centring::primitive;
// }

SKIntegerSymmetryOperationSet::SKIntegerSymmetryOperationSet(std::vector<SKSeitzIntegerMatrix> operations)
{
  this->operations = std::unordered_set<SKSeitzIntegerMatrix, SKSeitzIntegerMatrix::hashFunction>(operations.begin(),
                                                                                                  operations.end());
}

// A set holds all of its operations already (the setting's version of this decodes the whole group from the
// database); this returned an empty set
SKIntegerSymmetryOperationSet SKIntegerSymmetryOperationSet::fullSeitzMatrices() { return *this; }

std::vector<std::tuple<double3, std::size_t, double>> SKIntegerSymmetryOperationSet::symmetrize(
    double3x3 lattice, std::vector<std::tuple<double3, std::size_t, double>> atoms, double symmetryPrecision = 1e-2)
{
  std::vector<std::tuple<double3, std::size_t, double>> symmetrizedAtoms{};
  symmetrizedAtoms.reserve(atoms.size());

  for (std::size_t i = 0; i < atoms.size(); i++)
  {
    SKRotationMatrix averageRotation = SKRotationMatrix::zero;
    double3 averageTranslation = double3(0.0, 0.0, 0.0);
    int count = 0;

    for (const SKSeitzIntegerMatrix& operation : this->operations)
    {
      double3 translation = double3(double(operation.translation.x) / 24.0, double(operation.translation.y) / 24.0,
                                    double(operation.translation.z) / 24.0);
      double3 position = operation.rotation * std::get<0>(atoms[i]) + translation;

      if (SKSymmetryCell::isOverlap(position, std::get<0>(atoms[i]), lattice, symmetryPrecision))
      {
        averageRotation = averageRotation + operation.rotation;
        averageTranslation += translation - double3::rint(position - std::get<0>(atoms[i]));
        count = count + 1;
      }
    }

    double3x3 averagedRotation = double3x3(averageRotation.int3x3_m) / double(count);
    double3 averagedTranslation =
        double3(double(averageTranslation.x), double(averageTranslation.y), double(averageTranslation.z)) /
        double(count);

    std::tuple<double3, std::size_t, double> symmetrizedAtom = std::make_tuple(
        averagedRotation * std::get<0>(atoms[i]) + averagedTranslation, std::get<1>(atoms[i]), std::get<2>(atoms[i]));
    symmetrizedAtoms.push_back(symmetrizedAtom);
  }

  return symmetrizedAtoms;
}

std::vector<std::tuple<double3, std::size_t, double>> SKIntegerSymmetryOperationSet::asymmetricAtoms(
    const SKAsymmetricUnit& asymmetricUnit, std::vector<std::tuple<double3, std::size_t, double>>& atoms,
    double3x3 lattice, bool allowPartialOccupancies, double symmetryPrecision = 1e-2)
{
  std::vector<std::tuple<double3, std::size_t, double, std::make_signed_t<std::size_t>>> atomData{};
  std::transform(atoms.begin(), atoms.end(), std::back_inserter(atomData),
                 [](const std::tuple<double3, std::size_t, double>& atom)
                 { return std::make_tuple(std::get<0>(atom), std::get<1>(atom), std::get<2>(atom), -1); });

  if (atoms.empty()) return {};

  std::vector<std::tuple<double3, std::size_t, double>> asymmetricAtoms = {};
// std::get<3>(atomData[0]) = 0;

// loop over all atoms
loop:
  for (std::size_t i = 0; i < atoms.size(); i++)
  {
    // skip if already tagged
    if (std::get<3>(atomData[i]) == -1)
    {
      // loop over all current asymmetric atoms, and see if one of the symmetry-copies matches with an asymmetric
      // atom
      for (std::size_t j = 0; j < asymmetricAtoms.size(); j++)
      {
        // with partial occupancies the types of a shared site may differ (the comparison used to apply regardless)
        if (allowPartialOccupancies || std::get<1>(atomData[i]) == std::get<1>(asymmetricAtoms[j]))
        {
          for (const SKSeitzIntegerMatrix& operation : this->operations)
          {
            double3 position = operation.rotation * std::get<0>(atomData[i]) +
                               double3(double(operation.translation.x) / 24.0, double(operation.translation.y) / 24.0,
                                       double(operation.translation.z) / 24.0);
            if (SKSymmetryCell::isOverlap(position, std::get<0>(asymmetricAtoms[j]), lattice, symmetryPrecision))
            {
              // overlap and of the same type: the atom is therefore a copy of the asymmetric atom 'j'
              std::get<3>(atomData[i]) = static_cast<std::make_signed_t<std::size_t>>(j);
              break;
            }
          }
        }
      }

      // not typed yet
      if (std::get<3>(atomData[i]) == -1)
      {
        asymmetricAtoms.push_back(
            std::make_tuple(std::get<0>(atomData[i]), std::get<1>(atomData[i]), std::get<2>(atomData[i])));
        std::get<3>(atomData[i]) = static_cast<std::make_signed_t<std::size_t>>(asymmetricAtoms.size() - 1);
        goto loop;
      }
    }
  }

  // move each representative to the image of its orbit that lies in the asymmetric unit of the setting. The check
  // had been disabled, which left the loop computing images and discarding them; the unit is passed in, since the
  // space-group modules depend on this one, and the containment test was repaired (open bounds at 1, images below
  // -1/2).
  for (std::size_t i = 0; i < asymmetricAtoms.size(); i++)
  {
    for (const SKSeitzIntegerMatrix& operation : operations)
    {
      double3 position = operation.rotation * std::get<0>(asymmetricAtoms[i]) +
                         double3(double(operation.translation.x) / 24.0, double(operation.translation.y) / 24.0,
                                 double(operation.translation.z) / 24.0);
      if (asymmetricUnit.contains(position))
      {
        std::get<0>(asymmetricAtoms[i]) = double3::fract(position);
        break;
      }
    }
  }

  return asymmetricAtoms;
}
