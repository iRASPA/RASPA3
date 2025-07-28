module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>
#endif

module skspacegroup;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import int3;
import int3x3;
import double3;
import double3x3;
import ringmatrix;
import matrix;

import skdefinitions;
import skrotationmatrix;
import sktransformationmatrix;
import skseitzmatrix;
import skpointsymmetryset;
import sksymmetryoperationset;
import sksymmetrycell;
import skintegerchangeofbasis;
import skpointgroup;
import skspacegroupdatabase;
import skspacegroupsetting;
import skrotationalchangeofbasis;
import skseitzintegermatrix;
import skintegersymmetryoperationset;

std::string simplified(std::string a) { return a; };

std::string toLower(std::string a) { return a; }

SKSpaceGroup::SKSpaceGroup(std::size_t HallNumber)
{
  _spaceGroupSetting = SKSpaceGroupDataBase::spaceGroupData[HallNumber];
}

bool SKSpaceGroup::matchSpacegroup(std::string spaceSearchGroupString, std::string storedSpaceGroupString)
{
  if (storedSpaceGroupString == spaceSearchGroupString)
  {
    return true;
  }

  if ("'" + storedSpaceGroupString + "'" == spaceSearchGroupString)
  {
    return true;
  }

  if ("\"" + storedSpaceGroupString + "\"" == spaceSearchGroupString)
  {
    return true;
  }
  return false;
}

std::optional<std::size_t> SKSpaceGroup::HallNumber(std::string string)
{
  std::string spaceSearchGroupString = toLower(simplified(string));

  for (std::size_t i = 0; i <= 530; i++)
  {
    std::string storedSpaceGroupString = toLower(simplified(SKSpaceGroupDataBase::spaceGroupData[i].HallString()));

    if (SKSpaceGroup::matchSpacegroup(spaceSearchGroupString, storedSpaceGroupString))
    {
      return i;
    }

    storedSpaceGroupString =
        toLower(simplified(SKSpaceGroupDataBase::spaceGroupData[i].HallString()));  // .remove(' ');
    if (SKSpaceGroup::matchSpacegroup(spaceSearchGroupString, storedSpaceGroupString))
    {
      return i;
    }
  }

  return std::nullopt;
}

std::optional<std::size_t> SKSpaceGroup::HallNumberFromHMString(std::string string)
{
  std::string spaceSearchGroupString = toLower(simplified(string));

  for (std::size_t i = 0; i <= 530; i++)
  {
    std::string storedSpaceGroupString = toLower(simplified(SKSpaceGroupDataBase::spaceGroupData[i].HMString()));

    if (SKSpaceGroup::matchSpacegroup(spaceSearchGroupString, storedSpaceGroupString))
    {
      return i;
    }

    if (SKSpaceGroupDataBase::spaceGroupData[i].qualifier() == "H")
    {
      storedSpaceGroupString = storedSpaceGroupString.replace(0, 1, "h");
      if (SKSpaceGroup::matchSpacegroup(spaceSearchGroupString, storedSpaceGroupString))
      {
        return i;
      }
    }

    storedSpaceGroupString = toLower(simplified(SKSpaceGroupDataBase::spaceGroupData[i].HMString()));  //.remove(' ');
    if (SKSpaceGroup::matchSpacegroup(spaceSearchGroupString, storedSpaceGroupString))
    {
      return i;
    }
  }

  return std::nullopt;
}

std::optional<std::size_t> SKSpaceGroup::HallNumberFromSpaceGroupNumber([[maybe_unused]] std::size_t number)
{
  if (number > 0 && number <= 230)
  {
    std::vector<std::size_t> hall_numbers = SKSpaceGroupDataBase::spaceGroupHallData[number];
    if (hall_numbers.size() > 0)
    {
      return hall_numbers.front();
    }
  }
  return std::nullopt;
}

std::vector<double3> SKSpaceGroup::listOfSymmetricPositions(double3 pos)
{
  std::unordered_set<SKSeitzIntegerMatrix, SKSeitzIntegerMatrix::hashFunction> seitzMatrices =
      _spaceGroupSetting.fullSeitzMatrices().operations;
  std::size_t m = seitzMatrices.size();

  std::vector<double3> positions = std::vector<double3>{};
  positions.reserve(m);

  for (const SKSeitzIntegerMatrix& elem : seitzMatrices)
  {
    positions.push_back(elem * pos);
  }
  return positions;
}

std::vector<std::string> SKSpaceGroup::latticeTranslationStrings(std::size_t HallNumber)
{
  std::vector<std::string> latticeStrings{"", "", "", ""};

  SKSpaceGroupSetting setting = SKSpaceGroupDataBase::spaceGroupData[HallNumber];
  std::vector<int3> latticeVectors = setting.latticeTranslations();

  std::size_t index = 0;
  for (int3 latticeVector : latticeVectors)
  {
    int3 gcd = int3::greatestCommonDivisor(latticeVector, 12);
    std::string latticeStringX =
        latticeVector.x == 0 ? "0" : std::to_string(latticeVector.x / gcd.x) + "/" + std::to_string(12 / gcd.x);
    std::string latticeStringY =
        latticeVector.y == 0 ? "0" : std::to_string(latticeVector.y / gcd.y) + "/" + std::to_string(12 / gcd.y);
    std::string latticeStringZ =
        latticeVector.z == 0 ? "0" : std::to_string(latticeVector.z / gcd.z) + "/" + std::to_string(12 / gcd.z);
    latticeStrings[index] = "(" + latticeStringX + "," + latticeStringY + "," + latticeStringZ + ")";
    index++;
  }

  return latticeStrings;
}

std::string SKSpaceGroup::inversionCenterString(std::size_t HallNumber)
{
  SKSpaceGroupSetting setting = SKSpaceGroupDataBase::spaceGroupData[HallNumber];
  int3 inversionCenter = setting.inversionCenter();
  int3 gcd = int3::greatestCommonDivisor(inversionCenter, 12);
  std::string latticeStringX =
      inversionCenter.x == 0 ? "0" : std::to_string(inversionCenter.x / gcd.x) + "/" + std::to_string(12 / gcd.x);
  std::string latticeStringY =
      inversionCenter.y == 0 ? "0" : std::to_string(inversionCenter.y / gcd.y) + "/" + std::to_string(12 / gcd.y);
  std::string latticeStringZ =
      inversionCenter.z == 0 ? "0" : std::to_string(inversionCenter.z / gcd.z) + "/" + std::to_string(12 / gcd.z);
  return "(" + latticeStringX + "," + latticeStringY + "," + latticeStringZ + ")";
}

SKSymmetryOperationSet SKSpaceGroup::findSpaceGroupSymmetry(
    double3x3 unitCell, std::vector<std::tuple<double3, std::size_t, double>> reducedAtoms,
    std::vector<std::tuple<double3, std::size_t, double>> atoms, SKPointSymmetrySet latticeSymmetries,
    bool allowPartialOccupancies, double symmetryPrecision = 1e-2)
{
  std::vector<SKSeitzMatrix> spaceGroupSymmetries{};

  for (const SKRotationMatrix& rotationMatrix : latticeSymmetries.rotations())
  {
    std::vector<double3> translations = SKSymmetryCell::primitiveTranslationVectors(
        unitCell, reducedAtoms, atoms, rotationMatrix, allowPartialOccupancies, symmetryPrecision);

    for (const double3& translation : translations)
    {
      SKSeitzMatrix matrix = SKSeitzMatrix(rotationMatrix, translation);

      // avoid duplicate Seitz-matrices
      if (!(std::find(spaceGroupSymmetries.begin(), spaceGroupSymmetries.end(), matrix) != spaceGroupSymmetries.end()))
      {
        spaceGroupSymmetries.push_back(matrix);
      }
    }
  }
  return SKSymmetryOperationSet(spaceGroupSymmetries);
}

std::optional<SKSpaceGroup::FoundPrimitiveCellInfo> SKSpaceGroup::SKFindPrimitive(
    double3x3 unitCell, std::vector<std::tuple<double3, std::size_t, double>> atoms, bool allowPartialOccupancies,
    double symmetryPrecision = 1e-2)
{
  std::optional<FoundSpaceGroupInfo> foundSpaceGroup =
      findSpaceGroup(unitCell, atoms, allowPartialOccupancies, symmetryPrecision);
  if (foundSpaceGroup)
  {
    const auto& [HallNumber, origin, cell, changeOfBasis, transformationMatrix, rotationMatrix, atoms1,
                 asymemtricAtoms] = *foundSpaceGroup;

    SKSpaceGroup spaceGroup = SKSpaceGroup(HallNumber);
    Centring centring = spaceGroup.spaceGroupSetting().centring();

    double3x3 transformation = double3x3::identity();
    switch (centring)
    {
      case Centring::primitive:
        transformation = double3x3::identity();
        break;
      case Centring::body:
        transformation = SKTransformationMatrix::bodyCenteredToPrimitive;
        break;
      case Centring::face:
        transformation = SKTransformationMatrix::faceCenteredToPrimitive;
        break;
      case Centring::a_face:
        transformation = SKTransformationMatrix::ACenteredToPrimitive;
        break;
      case Centring::b_face:
        transformation = SKTransformationMatrix::BCenteredToPrimitive;
        break;
      case Centring::c_face:
        transformation = SKTransformationMatrix::CCenteredToPrimitive;
        break;
      case Centring::r:
        transformation = SKTransformationMatrix::rhombohedralToPrimitive;
        break;
      case Centring::h:
        transformation = SKTransformationMatrix::hexagonalToPrimitive;
        break;
      default:
        transformation = double3x3::identity();
        break;
    }

    double3x3 primitiveUnitCell = cell.unitCell() * transformation;

    SKSymmetryCell primitiveCell = SKSymmetryCell::createFromUnitCell(primitiveUnitCell);

    std::vector<std::tuple<double3, std::size_t, double>> positionInPrimitiveCell =
        SKSymmetryCell::trim(atoms1, cell.unitCell(), primitiveUnitCell, allowPartialOccupancies, symmetryPrecision);

    return SKSpaceGroup::FoundPrimitiveCellInfo{primitiveCell, positionInPrimitiveCell};
  }

  return std::nullopt;
}

std::optional<SKSpaceGroup::FoundSpaceGroupInfo> SKSpaceGroup::findSpaceGroup(
    double3x3 unitCell, std::vector<std::tuple<double3, std::size_t, double>> atoms, bool allowPartialOccupancies,
    double symmetryPrecision = 1e-2)
{
  std::vector<std::tuple<double3, std::size_t, double>> reducedAtoms{};

  if (allowPartialOccupancies)
  {
    reducedAtoms = atoms;
  }
  else
  {
    std::map<std::size_t, std::size_t> histogram{};
    for (const std::tuple<double3, std::size_t, double>& atom : atoms)
    {
      histogram[std::get<1>(atom)] = histogram[std::get<1>(atom)] + 1;
    }
    std::map<std::size_t, std::size_t>::iterator index = std::min_element(
        histogram.begin(), histogram.end(), [](const auto& l, const auto& r) { return l.second < r.second; });
    std::size_t leastOccuringAtomType = index->first;

    std::copy_if(atoms.begin(), atoms.end(), std::back_inserter(reducedAtoms),
                 [leastOccuringAtomType](std::tuple<double3, std::size_t, double> a)
                 { return std::get<1>(a) == leastOccuringAtomType; });
  }

  double3x3 smallestUnitCell = SKSymmetryCell::findSmallestPrimitiveCell(reducedAtoms, atoms, unitCell,
                                                                         allowPartialOccupancies, symmetryPrecision);

  std::optional<double3x3> primitiveDelaunayUnitCell =
      SKSymmetryCell::computeDelaunayReducedCell(smallestUnitCell, symmetryPrecision);

  if (primitiveDelaunayUnitCell)
  {
    SKPointSymmetrySet latticeSymmetries =
        SKSymmetryCell::findLatticeSymmetry(*primitiveDelaunayUnitCell, symmetryPrecision);

    std::vector<std::tuple<double3, std::size_t, double>> positionInPrimitiveCell =
        SKSymmetryCell::trim(atoms, unitCell, *primitiveDelaunayUnitCell, allowPartialOccupancies, symmetryPrecision);

    SKSymmetryOperationSet spaceGroupSymmetries =
        SKSpaceGroup::findSpaceGroupSymmetry(unitCell, positionInPrimitiveCell, positionInPrimitiveCell,
                                             latticeSymmetries, allowPartialOccupancies, symmetryPrecision);

    SKPointSymmetrySet pointSymmetry = SKPointSymmetrySet(spaceGroupSymmetries.rotations());

    std::optional<SKPointGroup> pointGroup = SKPointGroup(pointSymmetry);

    if (pointGroup)
    {
      // Use the axes directions of the Laue group-specific symmetry as a new basis
      std::optional<SKTransformationMatrix> Mprime = pointGroup->constructAxes(spaceGroupSymmetries.properRotations());

      if (Mprime)
      {
        // adjustment of (M',0) to (M,0) for certain combination of Laue and centring types
        switch (pointGroup->laue())
        {
          case Laue::laue_1:
          {
            SKSymmetryCell::createFromUnitCell((*primitiveDelaunayUnitCell) * (*Mprime));
            std::optional<std::pair<SKSymmetryCell, SKTransformationMatrix>> symmetryCell =
                SKSymmetryCell::createFromUnitCell((*primitiveDelaunayUnitCell) * (*Mprime))
                    .computeReducedNiggliCellAndChangeOfBasisMatrix();
            if (!symmetryCell)
            {
              return std::nullopt;
            }
            Mprime = std::get<1>(*symmetryCell);
            break;
          }
          case Laue::laue_2m:
          {
            // Change the basis for this monoclinic centrosymmetric point group using Delaunay reduction in
            // 2D (algorithm of Atsushi Togo used) The unique axis is chosen as b, choose shortest a, c
            // lattice vectors
            // (|a| < |c|)
            std::optional<double3x3> computedDelaunayReducedCell2D = SKSymmetryCell::computeDelaunayReducedCell2D(
                (*primitiveDelaunayUnitCell) * (*Mprime), symmetryPrecision);
            if (!computedDelaunayReducedCell2D)
            {
              return std::nullopt;
            }
            double3x3 temp = primitiveDelaunayUnitCell->inverse() * *computedDelaunayReducedCell2D;
            Mprime = SKTransformationMatrix(temp.toInt3x3());
            break;
          }
          default:
            break;
        }

        Centring centering = pointGroup->computeCentering(*Mprime);

        SKTransformationMatrix correctedBasis = pointGroup->computeBasisCorrection(*Mprime, centering);

        double3x3 primitiveLattice = (*primitiveDelaunayUnitCell) * correctedBasis;

        // transform the symmetries (rotation and translation) from the primtive cell to the conventional cell
        // the centering is used to add the additional translations
        SKSymmetryOperationSet symmetryInConventionalCell =
            spaceGroupSymmetries.changedBasis(correctedBasis).addingCenteringOperations(centering);

        for (std::size_t i = 1; i <= 230; i++)
        {
          std::size_t HallNumber = SKSpaceGroupDataBase::spaceGroupHallData[i].front();
          if (SKSpaceGroupDataBase::spaceGroupData[HallNumber].pointGroupNumber() == pointGroup->number())
          {
            std::optional<std::pair<double3, SKRotationalChangeOfBasis>> foundSpaceGroup =
                SKSpaceGroup::matchSpaceGroup(HallNumber, primitiveLattice, centering,
                                              symmetryInConventionalCell.operations, symmetryPrecision);
            if (foundSpaceGroup)
            {
              double3 origin = std::get<0>(*foundSpaceGroup);
              SKRotationalChangeOfBasis changeOfBasis = std::get<1>(*foundSpaceGroup);
              double3x3 conventionalBravaisLattice = primitiveLattice * changeOfBasis.inverseRotationMatrix;

              double3x3 transformationMatrix = conventionalBravaisLattice.inverse() * unitCell;

              const SKSpaceGroup spaceGroup = SKSpaceGroup(HallNumber);

              SKIntegerSymmetryOperationSet dataBaseSpaceGroupSymmetries =
                  spaceGroup.spaceGroupSetting().fullSeitzMatrices();

              double3x3 transform = conventionalBravaisLattice.inverse() * *primitiveDelaunayUnitCell;

              std::vector<std::tuple<double3, std::size_t, double>> atomsInConventionalCell{};
              std::transform(positionInPrimitiveCell.begin(), positionInPrimitiveCell.end(),
                             std::back_inserter(atomsInConventionalCell),
                             [transform, origin](const std::tuple<double3, std::size_t, double>& tuple)
                             {
                               return std::make_tuple(double3::fract(transform * std::get<0>(tuple) + origin),
                                                      std::get<1>(tuple), std::get<2>(tuple));
                             });

              std::vector<std::tuple<double3, std::size_t, double>> symmetrizedAtomsInConventionalCell =
                  dataBaseSpaceGroupSymmetries.symmetrize(conventionalBravaisLattice, atomsInConventionalCell,
                                                          symmetryPrecision);

              std::vector<std::tuple<double3, std::size_t, double>> asymmetricAtoms =
                  dataBaseSpaceGroupSymmetries.asymmetricAtoms(HallNumber, symmetrizedAtomsInConventionalCell,
                                                               conventionalBravaisLattice, allowPartialOccupancies,
                                                               symmetryPrecision);

              SKSymmetryCell temp = SKSymmetryCell::createFromUnitCell(conventionalBravaisLattice);
              SKSymmetryCell cell = temp.idealized(spaceGroup.spaceGroupSetting().pointGroupNumber(),
                                                   spaceGroup.spaceGroupSetting().qualifier());
              double3x3 rotationMatrix = conventionalBravaisLattice * cell.unitCell().inverse();

              // must be a rigid rotation
              // assert((rotationMatrix.determinant() - 1.0) < 1e-5);

              return SKSpaceGroup::FoundSpaceGroupInfo{HallNumber,
                                                       origin,
                                                       cell,
                                                       SKRotationalChangeOfBasis(SKRotationMatrix()),
                                                       transformationMatrix,
                                                       rotationMatrix,
                                                       symmetrizedAtomsInConventionalCell,
                                                       asymmetricAtoms};
            }
          }
        }

        // special cases
        // Gross-Kunstleve: special case Pa-3 (205) hallSymbol 501
        std::size_t HallNumber = SKSpaceGroupDataBase::spaceGroupHallData[205].front();
        if (SKSpaceGroupDataBase::spaceGroupData[HallNumber].pointGroupNumber() == pointGroup->number())
        {
          std::optional<double3> originShift =
              getOriginShift(HallNumber, centering,
                             SKRotationalChangeOfBasis(SKRotationMatrix(int3(0, 0, 1), int3(0, -1, 0), int3(1, 0, 0))),
                             symmetryInConventionalCell.operations, symmetryPrecision);
          if (originShift)
          {
            double3 origin = *originShift;
            SKRotationalChangeOfBasis changeOfBasis =
                SKRotationalChangeOfBasis(SKRotationMatrix(int3(0, 0, 1), int3(0, -1, 0), int3(1, 0, 0)));
            double3x3 conventionalBravaisLattice = primitiveLattice * changeOfBasis.inverseRotationMatrix;

            double3x3 transformationMatrix = conventionalBravaisLattice.inverse() * unitCell;

            SKSpaceGroup spaceGroup = SKSpaceGroup(HallNumber);

            SKIntegerSymmetryOperationSet dataBaseSpaceGroupSymmetries =
                spaceGroup.spaceGroupSetting().fullSeitzMatrices();

            double3x3 transform = conventionalBravaisLattice.inverse() * *primitiveDelaunayUnitCell;

            std::vector<std::tuple<double3, std::size_t, double>> atomsInConventionalCell{};
            std::transform(positionInPrimitiveCell.begin(), positionInPrimitiveCell.end(),
                           std::back_inserter(atomsInConventionalCell),
                           [transform, origin](const std::tuple<double3, std::size_t, double>& tuple)
                           {
                             return std::make_tuple(transform * std::get<0>(tuple) + origin, std::get<1>(tuple),
                                                    std::get<2>(tuple));
                           });

            std::vector<std::tuple<double3, std::size_t, double>> symmetrizedAtomsInConventionalCell =
                dataBaseSpaceGroupSymmetries.symmetrize(conventionalBravaisLattice, atomsInConventionalCell,
                                                        symmetryPrecision);

            std::vector<std::tuple<double3, std::size_t, double>> asymmetricAtoms =
                dataBaseSpaceGroupSymmetries.asymmetricAtoms(HallNumber, symmetrizedAtomsInConventionalCell,
                                                             conventionalBravaisLattice, allowPartialOccupancies,
                                                             symmetryPrecision);

            SKSymmetryCell temp = SKSymmetryCell::createFromUnitCell(conventionalBravaisLattice);
            SKSymmetryCell cell = temp.idealized(spaceGroup.spaceGroupSetting().pointGroupNumber(),
                                                 spaceGroup.spaceGroupSetting().qualifier());
            double3x3 rotationMatrix = conventionalBravaisLattice * cell.unitCell().inverse();

            // must be a rigid rotation
            // assert((rotationMatrix.determinant() - 1.0) < 1e-5);

            return SKSpaceGroup::FoundSpaceGroupInfo{HallNumber,
                                                     origin,
                                                     cell,
                                                     SKRotationalChangeOfBasis(SKRotationMatrix()),
                                                     transformationMatrix,
                                                     rotationMatrix,
                                                     symmetrizedAtomsInConventionalCell,
                                                     asymmetricAtoms};
          }
        }
      }
    }
  }
  return std::nullopt;
}

std::optional<SKSpaceGroup::FoundNiggliCellInfo> SKSpaceGroup::findNiggliCell(
    double3x3 unitCell, std::vector<std::tuple<double3, std::size_t, double>> atoms, bool allowPartialOccupancies,
    double symmetryPrecision = 1e-2)
{
  std::vector<std::tuple<double3, std::size_t, double>> reducedAtoms{};

  std::size_t leastOccuringAtomType;
  if (allowPartialOccupancies)
  {
    reducedAtoms = atoms;
  }
  else
  {
    std::map<std::size_t, std::size_t> histogram{};
    for (const std::tuple<double3, std::size_t, double>& atom : atoms)
    {
      histogram[std::get<1>(atom)] = histogram[std::get<1>(atom)] + 1;
    }
    std::map<std::size_t, std::size_t>::iterator index = std::min_element(
        histogram.begin(), histogram.end(), [](const auto& l, const auto& r) { return l.second < r.second; });
    leastOccuringAtomType = index->first;

    std::copy_if(atoms.begin(), atoms.end(), std::back_inserter(reducedAtoms),
                 [leastOccuringAtomType](std::tuple<double3, std::size_t, double> a)
                 { return std::get<1>(a) == leastOccuringAtomType; });
  }

  double3x3 smallestUnitCell = SKSymmetryCell::findSmallestPrimitiveCell(reducedAtoms, atoms, unitCell,
                                                                         allowPartialOccupancies, symmetryPrecision);

  std::optional<double3x3> primitiveDelaunayUnitCell =
      SKSymmetryCell::computeDelaunayReducedCell(smallestUnitCell, symmetryPrecision);

  if (primitiveDelaunayUnitCell)
  {
    std::vector<std::tuple<double3, std::size_t, double>> positionInPrimitiveDelaunayCell =
        SKSymmetryCell::trim(atoms, unitCell, *primitiveDelaunayUnitCell, allowPartialOccupancies, symmetryPrecision);

    std::optional<std::pair<SKSymmetryCell, SKTransformationMatrix>> NiggliSymmetryCell =
        SKSymmetryCell::createFromUnitCell(*primitiveDelaunayUnitCell).computeReducedNiggliCellAndChangeOfBasisMatrix();

    if (NiggliSymmetryCell)
    {
      double3x3 NiggliUnitCell = std::get<0>(*NiggliSymmetryCell).unitCell();
      SKTransformationMatrix changeOfBasis = std::get<1>(*NiggliSymmetryCell);

      std::vector<std::tuple<double3, std::size_t, double>> positionInNiggliCell{};
      std::transform(positionInPrimitiveDelaunayCell.begin(), positionInPrimitiveDelaunayCell.end(),
                     std::back_inserter(positionInNiggliCell),
                     [changeOfBasis](const std::tuple<double3, std::size_t, double>& tuple)
                     {
                       return std::make_tuple(double3x3(changeOfBasis.transformation).inverse() * std::get<0>(tuple),
                                              std::get<1>(tuple), std::get<2>(tuple));
                     });

      std::vector<std::tuple<double3, std::size_t, double>> reducedPositionsInNiggliCell{};
      if (allowPartialOccupancies)
      {
        std::copy(positionInNiggliCell.begin(), positionInNiggliCell.end(),
                  std::back_inserter(reducedPositionsInNiggliCell));
      }
      else
      {
        std::copy_if(positionInNiggliCell.begin(), positionInNiggliCell.end(),
                     std::back_inserter(reducedPositionsInNiggliCell),
                     [leastOccuringAtomType](const std::tuple<double3, std::size_t, double>& atom)
                     { return std::get<1>(atom) == leastOccuringAtomType; });
      }

      SKPointSymmetrySet latticeSymmetries = SKSymmetryCell::findLatticeSymmetry(NiggliUnitCell, symmetryPrecision);

      SKSymmetryOperationSet spaceGroupSymmetries =
          SKSpaceGroup::findSpaceGroupSymmetry(NiggliUnitCell, reducedPositionsInNiggliCell, positionInNiggliCell,
                                               latticeSymmetries, allowPartialOccupancies, symmetryPrecision);

      for (std::size_t i = 230; i >= 1; i--)
      {
        std::size_t HallNumber = SKSpaceGroupDataBase::spaceGroupHallData[i].front();

        std::optional<std::pair<double3, SKRotationalChangeOfBasis>> foundSpaceGroup = SKSpaceGroup::matchSpaceGroup(
            HallNumber, NiggliUnitCell, Centring::primitive, spaceGroupSymmetries.operations, symmetryPrecision);
        if (foundSpaceGroup)
        {
          double3 origin = std::get<0>(*foundSpaceGroup);

          double3x3 conventionalBravaisLattice = NiggliUnitCell * std::get<1>(*foundSpaceGroup).inverseRotationMatrix;

          SKSpaceGroup spaceGroup = SKSpaceGroup(HallNumber);

          SKIntegerSymmetryOperationSet dataBaseSpaceGroupSymmetries =
              spaceGroup.spaceGroupSetting().fullSeitzMatrices();

          double3x3 transform = conventionalBravaisLattice.inverse() * NiggliUnitCell;

          std::vector<std::tuple<double3, std::size_t, double>> atomsInConventionalCell{};
          std::transform(positionInNiggliCell.begin(), positionInNiggliCell.end(),
                         std::back_inserter(atomsInConventionalCell),
                         [transform, origin](const std::tuple<double3, std::size_t, double>& tuple)
                         {
                           return std::make_tuple(double3::fract(transform * std::get<0>(tuple) + origin),
                                                  std::get<1>(tuple), std::get<2>(tuple));
                         });

          std::vector<std::tuple<double3, std::size_t, double>> symmetrizedAtomsInConventionalCell =
              dataBaseSpaceGroupSymmetries.symmetrize(conventionalBravaisLattice, atomsInConventionalCell,
                                                      symmetryPrecision);

          std::vector<std::tuple<double3, std::size_t, double>> asymmetricAtoms =
              dataBaseSpaceGroupSymmetries.asymmetricAtoms(HallNumber, symmetrizedAtomsInConventionalCell,
                                                           conventionalBravaisLattice, allowPartialOccupancies,
                                                           symmetryPrecision);

          SKSymmetryCell temp = SKSymmetryCell::createFromUnitCell(conventionalBravaisLattice);
          SKSymmetryCell cell = temp.idealized(spaceGroup.spaceGroupSetting().pointGroupNumber(),
                                               spaceGroup.spaceGroupSetting().qualifier());
          return SKSpaceGroup::FoundNiggliCellInfo{HallNumber, cell, asymmetricAtoms};
        }
      }
    }
  }
  return std::nullopt;
}

std::optional<SKPointGroup> SKSpaceGroup::findPointGroup(double3x3 unitCell,
                                                         std::vector<std::tuple<double3, std::size_t, double>> atoms,
                                                         bool allowPartialOccupancies, double symmetryPrecision = 1e-2)
{
  std::vector<std::tuple<double3, std::size_t, double>> reducedAtoms{};

  if (allowPartialOccupancies)
  {
    reducedAtoms = atoms;
  }
  else
  {
    std::map<std::size_t, std::size_t> histogram{};
    for (const std::tuple<double3, std::size_t, double>& atom : atoms)
    {
      histogram[std::get<1>(atom)] = histogram[std::get<1>(atom)] + 1;
    }
    std::map<std::size_t, std::size_t>::iterator index = std::min_element(
        histogram.begin(), histogram.end(), [](const auto& l, const auto& r) { return l.second < r.second; });
    std::size_t leastOccuringAtomType = index->first;

    std::copy_if(atoms.begin(), atoms.end(), std::back_inserter(reducedAtoms),
                 [leastOccuringAtomType](std::tuple<double3, std::size_t, double> a)
                 { return std::get<1>(a) == leastOccuringAtomType; });
  }
  double3x3 smallestUnitCell = SKSymmetryCell::findSmallestPrimitiveCell(reducedAtoms, atoms, unitCell,
                                                                         allowPartialOccupancies, symmetryPrecision);

  std::optional<double3x3> primitiveDelaunayUnitCell =
      SKSymmetryCell::computeDelaunayReducedCell(smallestUnitCell, symmetryPrecision);

  if (primitiveDelaunayUnitCell)
  {
    SKPointSymmetrySet latticeSymmetries =
        SKSymmetryCell::findLatticeSymmetry(*primitiveDelaunayUnitCell, symmetryPrecision);

    std::vector<std::tuple<double3, std::size_t, double>> positionInPrimitiveCell =
        SKSymmetryCell::trim(atoms, unitCell, *primitiveDelaunayUnitCell, allowPartialOccupancies, symmetryPrecision);

    SKSymmetryOperationSet spaceGroupSymmetries =
        SKSpaceGroup::findSpaceGroupSymmetry(unitCell, positionInPrimitiveCell, positionInPrimitiveCell,
                                             latticeSymmetries, allowPartialOccupancies, symmetryPrecision);

    SKPointSymmetrySet pointSymmetry = SKPointSymmetrySet(spaceGroupSymmetries.rotations());

    return SKPointGroup(pointSymmetry);
  }
  return std::nullopt;
}

std::optional<std::pair<double3, SKRotationalChangeOfBasis>> SKSpaceGroup::matchSpaceGroup(
    std::size_t HallNumber, double3x3 lattice, Centring centering, std::vector<SKSeitzMatrix> seitzMatrices,
    double symmetryPrecision = 1e-2)
{
  std::size_t pointGroupNumber = SKSpaceGroupDataBase::spaceGroupData[HallNumber].pointGroupNumber();
  switch (SKPointGroup::pointGroupData[pointGroupNumber].holohedry())
  {
    case Holohedry::none:
      break;
    case Holohedry::triclinic:
    case Holohedry::tetragonal:
    case Holohedry::trigonal:
    case Holohedry::hexagonal:
    case Holohedry::cubic:
    {
      std::optional<double3> originShift =
          getOriginShift(HallNumber, centering, SKRotationalChangeOfBasis(SKRotationMatrix::identity), seitzMatrices,
                         symmetryPrecision);
      if (originShift)
      {
        return std::make_pair(*originShift, SKRotationalChangeOfBasis(SKRotationMatrix::identity));
      }
      return std::nullopt;
    }
    case Holohedry::monoclinic:
    {
      std::vector<std::pair<double3, SKRotationalChangeOfBasis>> solutions{};
      for (const SKRotationalChangeOfBasis& changeOfMonoclinicCentering :
           SKRotationalChangeOfBasis::changeOfMonoclinicCentering)
      {
        std::optional<double3> originShift =
            getOriginShift(HallNumber, centering, changeOfMonoclinicCentering, seitzMatrices, symmetryPrecision);
        if (originShift)
        {
          solutions.push_back(std::make_pair(*originShift, changeOfMonoclinicCentering));
        }
      }
      std::sort(solutions.begin(), solutions.end(),
                [lattice](const std::pair<double3, SKRotationalChangeOfBasis>& a,
                          const std::pair<double3, SKRotationalChangeOfBasis>& b) -> bool
                {
                  double3x3 conventionalBravaisLatticeA = lattice * std::get<1>(a).inverseRotationMatrix;
                  SKSymmetryCell cellA = SKSymmetryCell::createFromUnitCell(conventionalBravaisLatticeA);
                  double3x3 conventionalBravaisLatticeB = lattice * std::get<1>(b).inverseRotationMatrix;
                  SKSymmetryCell cellB = SKSymmetryCell::createFromUnitCell(conventionalBravaisLatticeB);
                  return (cellA.a() + cellA.c()) < (cellB.a() + cellB.c());
                });
      if (!solutions.empty()) return solutions.front();

      return std::nullopt;
    }
    case Holohedry::orthorhombic:
    {
      std::vector<std::pair<double3, SKRotationalChangeOfBasis>> solutions{};
      for (const SKRotationalChangeOfBasis& changeOfOrthorhombicCentering :
           SKRotationalChangeOfBasis::changeOfOrthorhombicCentering)
      {
        std::optional<double3> originShift =
            getOriginShift(HallNumber, centering, changeOfOrthorhombicCentering, seitzMatrices, symmetryPrecision);
        if (originShift)
        {
          solutions.push_back(std::make_pair(*originShift, changeOfOrthorhombicCentering));
        }
      }
      std::sort(solutions.begin(), solutions.end(),
                [lattice](const std::pair<double3, SKRotationalChangeOfBasis>& a,
                          const std::pair<double3, SKRotationalChangeOfBasis>& b) -> bool
                {
                  double3x3 conventionalBravaisLatticeA = lattice * std::get<1>(a).inverseRotationMatrix;
                  SKSymmetryCell cellA = SKSymmetryCell::createFromUnitCell(conventionalBravaisLatticeA);
                  double3x3 conventionalBravaisLatticeB = lattice * std::get<1>(b).inverseRotationMatrix;
                  SKSymmetryCell cellB = SKSymmetryCell::createFromUnitCell(conventionalBravaisLatticeB);
                  return std::pair(cellA.a(), cellA.b()) < std::pair(cellB.a(), cellB.b());
                });
      if (!solutions.empty()) return solutions.front();

      return std::nullopt;
    }
    default:
      break;
  }

  // qDebug() << "SHOULD NOT GET HERE";
  // assert(false);
  return std::nullopt;
}

std::optional<double3> SKSpaceGroup::getOriginShift(std::size_t HallNumber, Centring centering,
                                                    SKRotationalChangeOfBasis changeOfBasis,
                                                    std::vector<SKSeitzMatrix> seitzMatrices,
                                                    double symmetryPrecision = 1e-2)
{
  double3x3 translationsnew = double3x3();

  SKSpaceGroup dataBaseSpaceGroup = SKSpaceGroup(HallNumber);
  std::vector<SKSeitzIntegerMatrix> dataBaseSpaceGroupGenerators =
      SKSeitzIntegerMatrix::SeitzMatrices(dataBaseSpaceGroup.spaceGroupSetting().encodedGenerators());

  // assert(!dataBaseSpaceGroupGenerators.empty());

  // apply change-of-basis to generators
  for (std::size_t i = 0; i < dataBaseSpaceGroupGenerators.size(); i++)
  {
    dataBaseSpaceGroupGenerators[i] = changeOfBasis * dataBaseSpaceGroupGenerators[i];
  }

  // apply change-of-basis to lattice translations
  std::vector<int3> spaceGroupLatticeTranslations = dataBaseSpaceGroup.spaceGroupSetting().latticeTranslations();
  for (std::size_t i = 0; i < spaceGroupLatticeTranslations.size(); i++)
  {
    spaceGroupLatticeTranslations[i] = changeOfBasis * spaceGroupLatticeTranslations[i];
  }

  // apply change-of-basis to centring
  Centring dataBaseSpaceGroupCentering = dataBaseSpaceGroup.spaceGroupSetting().centring();
  switch (dataBaseSpaceGroupCentering)
  {
    case Centring::a_face:
    case Centring::b_face:
    case Centring::c_face:
      if (spaceGroupLatticeTranslations[1].x == 0)
      {
        dataBaseSpaceGroupCentering = Centring::a_face;
      }
      if (spaceGroupLatticeTranslations[1].y == 0)
      {
        dataBaseSpaceGroupCentering = Centring::b_face;
      }
      if (spaceGroupLatticeTranslations[1].z == 0)
      {
        dataBaseSpaceGroupCentering = Centring::c_face;
      }
      break;
    default:
      break;
  }

  // return if the centring is not equal to the spacegroup one
  if (centering != dataBaseSpaceGroupCentering)
  {
    return std::nullopt;
  }

  // apply change-of-basis to the Seitz-matrices
  std::vector<SKSeitzIntegerMatrix> dataBaseSpaceGroupSeitzMatrices =
      dataBaseSpaceGroup.spaceGroupSetting().SeitzMatricesWithoutTranslation();
  for (std::size_t i = 0; i < dataBaseSpaceGroupSeitzMatrices.size(); i++)
  {
    dataBaseSpaceGroupSeitzMatrices[i] = changeOfBasis * dataBaseSpaceGroupSeitzMatrices[i];
  }

  for (std::size_t i = 0; i < dataBaseSpaceGroupGenerators.size(); i++)
  {
    // math the rotional part of the generator with the Seitz-matrices
    SKRotationMatrix toFind = dataBaseSpaceGroupGenerators[i].rotation;
    auto index = std::find_if(seitzMatrices.begin(), seitzMatrices.end(),
                              [&toFind](const SKSeitzMatrix& m) { return m.rotation == toFind; });
    if (index == seitzMatrices.end())
    {
      return std::nullopt;
    }

    // and then take the translation part
    translationsnew[i] = index->translation;
  }

  SKTransformationMatrix transformation = SKTransformationMatrix::identity;
  switch (dataBaseSpaceGroupCentering)
  {
    case Centring::primitive:
      transformation = SKTransformationMatrix::identity;
      break;
    case Centring::body:
      transformation = SKTransformationMatrix::primitiveToBodyCentered;
      break;
    case Centring::face:
      transformation = SKTransformationMatrix::primitiveToFaceCentered;
      break;
    case Centring::a_face:
      transformation = SKTransformationMatrix::primitiveToACentered;
      break;
    case Centring::b_face:
      transformation = SKTransformationMatrix::primitiveToBCentered;
      break;
    case Centring::c_face:
      transformation = SKTransformationMatrix::primitiveToCCentered;
      break;
    case Centring::r:
      transformation = SKTransformationMatrix::primitiveToRhombohedral;
      break;
    case Centring::h:
      transformation = SKTransformationMatrix::primitiveToHexagonal;
      break;
    default:
      break;
  }

  // SKIntegerChangeOfBasis changeToPrimitive = SKIntegerChangeOfBasis(inversionTransformation: transformation)
  SKIntegerChangeOfBasis changeToPrimitive = SKIntegerChangeOfBasis(transformation);

  SKRotationMatrix r1 = dataBaseSpaceGroupGenerators[0].rotation;
  SKRotationMatrix r2 =
      dataBaseSpaceGroupGenerators.size() > 1 ? dataBaseSpaceGroupGenerators[1].rotation : SKRotationMatrix::identity;
  SKRotationMatrix r3 =
      dataBaseSpaceGroupGenerators.size() > 2 ? dataBaseSpaceGroupGenerators[2].rotation : SKRotationMatrix::identity;

  SKTransformationMatrix t1 = changeToPrimitive * SKTransformationMatrix((r1 - SKRotationMatrix::identity).int3x3_m);
  SKTransformationMatrix t2 = changeToPrimitive * SKTransformationMatrix((r2 - SKRotationMatrix::identity).int3x3_m);
  SKTransformationMatrix t3 = changeToPrimitive * SKTransformationMatrix((r3 - SKRotationMatrix::identity).int3x3_m);

  // m is a 9x3 matrix
  RingMatrix m = RingMatrix(t1.transformation, t2.transformation, t3.transformation);

  // The system M * cp = b (mod Z) can be solved by computing the Smith normal form D = PMQ.
  // b is the translation difference, cp the origin shift
  // D is a matrix in diagonal form with diagonal entries d1, . . . , dn.
  // P is square, 9x9, invertible matrix
  // Q is square, 3x3, invertible matrix
  std::tuple<RingMatrix, RingMatrix, RingMatrix> sol = m.SmithNormalForm();
  RingMatrix P = std::get<0>(sol);
  RingMatrix Q = std::get<1>(sol);
  RingMatrix D = std::get<2>(sol);

  Matrix b = Matrix(9, 1, 0.0);
  for (std::size_t i = 0; i < dataBaseSpaceGroupGenerators.size(); i++)
  {
    // let seitzMatrix: SKSeitzIntegerMatrix? = dataBaseSpaceGroupSeitzMatrices.filter{$0.rotation ==
    // dataBaseSpaceGroupGenerators[i].rotation}.first guard seitzMatrix != nil else {return nil}

    double3 transPrimitive = changeToPrimitive * translationsnew[i];

    int3 dataBaseTranslation = changeToPrimitive * dataBaseSpaceGroupGenerators[i].translation;

    double3 translationDifference = double3::fract(transPrimitive - double3(dataBaseTranslation) / 24.0);
    b(3 * i, 0) = translationDifference.x;
    b(3 * i + 1, 0) = translationDifference.y;
    b(3 * i + 2, 0) = translationDifference.z;
  }

  // v (9x1) =  P (9x9) x b (9,1)
  Matrix v = P * b;

  // The system P * b = v, v = [v1,...,vn] has solutions(mod Z) if and only if vi==0 whenever di=0
  if ((D(0, 0) == 0 && std::fabs(v(0, 0) - std::rint(v(0, 0))) > symmetryPrecision) ||
      (D(1, 1) == 0 && std::fabs(v(1, 0) - std::rint(v(1, 0))) > symmetryPrecision) ||
      (D(2, 2) == 0 && std::fabs(v(2, 0) - std::rint(v(2, 0))) > symmetryPrecision))
  {
    return std::nullopt;
  }
  for (std::size_t i = 3; i < 9; i++)
  {
    if (std::fabs(v(i, 0) - std::rint(v(i, 0))) > symmetryPrecision)
    {
      return std::nullopt;
    }
  }

  Matrix Dinv = Matrix(3, 9, 0.0);
  for (std::size_t i = 0; i < 3; i++)
  {
    if (D(i, i) != 0)
    {
      Dinv(i, i) = 1.0 / double(D(i, i));
    }
  }

  // sol.Q (3x3), T (3x9), sol.P (9x9), bm (9x1) -> (3x1)
  Matrix cp = (Q * Dinv * P) * b;

  double3 originShift = double3::fract(double3(cp(0, 0), cp(1, 0), cp(2, 0)));
  SKIntegerChangeOfBasis basis = SKIntegerChangeOfBasis(changeToPrimitive).inverse();
  return double3::fract(changeOfBasis.inverse() * (basis * originShift));
}
