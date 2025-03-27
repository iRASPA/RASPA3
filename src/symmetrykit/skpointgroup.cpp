module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <map>
#include <string>
#include <type_traits>
#include <unordered_set>
#include <vector>
#endif

module skpointgroup;

#ifndef USE_LEGACY_HEADERS
import <type_traits>;
import <algorithm>;
import <iterator>;
import <map>;
import <string>;
import <unordered_set>;
import <vector>;
import <iterator>;
#endif

import int3;
import int3x3;
import double3;
import double3x3;

import skdefinitions;
import sktransformationmatrix;
import skrotationmatrix;
import skrotationaloccurancetable;
import skpointsymmetryset;

SKPointGroup::SKPointGroup(SKRotationalOccuranceTable table, size_t number, std::string symbol, std::string schoenflies,
                           Holohedry holohedry, Laue laue, bool centrosymmetric, bool enantiomorphic)
    : _table(table),
      _number(number),
      _symbol(symbol),
      _schoenflies(schoenflies),
      _holohedry(holohedry),
      _laue(laue),
      _centrosymmetric(centrosymmetric),
      _enantiomorphic(enantiomorphic)
{
}

SKPointGroup::SKPointGroup(SKPointSymmetrySet pointSymmetry) : _table(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
{
  SKRotationalOccuranceTable table = SKRotationalOccuranceTable(0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

  std::vector<SKRotationMatrix> rotationMatrices = pointSymmetry.rotations();

  for (const SKRotationMatrix& rotation : rotationMatrices)
  {
    std::make_signed_t<std::size_t> type =
        static_cast<typename std::underlying_type<SKRotationMatrix::RotationType>::type>(rotation.type());
    table.occurance[type] += 1;
  }

  for (const SKPointGroup& pointGroup : pointGroupData)
  {
    if (pointGroup.table() == table)
    {
      *this = pointGroup;
      return;
    }
  }
}

std::string SKPointGroup::holohedryString() const
{
  switch (_holohedry)
  {
    case Holohedry::none:
      return "None";
    case Holohedry::triclinic:
      return "Triclinic";
    case Holohedry::monoclinic:
      return "Monoclinic";
    case Holohedry::orthorhombic:
      return "Orthorhombic";
    case Holohedry::tetragonal:
      return "Tetragonal";
    case Holohedry::trigonal:
      return "Trigonal";
    case Holohedry::hexagonal:
      return "Hexagonal";
    case Holohedry::cubic:
      return "Cubic";
  }
  return std::string();
}

std::string SKPointGroup::LaueString() const
{
  switch (_laue)
  {
    case Laue::none:
      return "none";
    case Laue::laue_1:
      return "1";
    case Laue::laue_2m:
      return "2m";
    case Laue::laue_mmm:
      return "mmm";
    case Laue::laue_4m:
      return "4m";
    case Laue::laue_4mmm:
      return "4mmm";
    case Laue::laue_3:
      return "3";
    case Laue::laue_3m:
      return "3m";
    case Laue::laue_6m:
      return "6m";
    case Laue::laue_6mmm:
      return "6mmm";
    case Laue::laue_m3:
      return "m3";
    case Laue::laue_m3m:
      return "m3m";
  }
  return std::string();
}

/// Table 3 of spglib lists the conditions to determine the centring types
/// "spacegroup.c" line 2075 of 2312
Centring SKPointGroup::computeCentering(SKTransformationMatrix basis)
{
  // the absolute value of the determinant gives the scale factor by which volume is multiplied under the associated
  // linear transformation, while its sign indicates whether the transformation preserves orientation

  // Number of lattice points per cell (1.2.1 in Hahn 2005 fifth ed.)
  // 1: primitive centred
  // 2: C-face centred, B-face centred, A-face centred, body-centred
  // 3: rhombohedrally centred, hexagonally centred
  // 4: all-face centred

  switch (std::abs(basis.determinant()))
  {
    case 1:
      return Centring::primitive;
    case 2:
    {
      // detect a-center
      for (int i = 0; i < 3; i++)
      {
        // if (1,0,0) is found, then 'a' is detected
        if (std::abs(basis[0][i]) == 1 && basis[1][i] == 0 && basis[2][i] == 0)
        {
          return Centring::a_face;
        }
      }

      // detect b-center
      for (int i = 0; i < 3; i++)
      {
        // if (0,1,0) is found, then 'b' is detected
        if (basis[0][i] == 0 && std::abs(basis[1][i]) == 1 && basis[2][i] == 0)
        {
          return Centring::b_face;
        }
      }

      // detect c-center
      for (int i = 0; i < 3; i++)
      {
        // if (0,0,1) is found, then 'b' is detected
        if (basis[0][i] == 0 && basis[1][i] == 0 && std::abs(basis[2][i]) == 1)
        {
          return Centring::c_face;
        }
      }

      // detect body-center
      int sum = std::abs(basis[0][0]) + std::abs(basis[1][0]) + std::abs(basis[2][0]);
      if (sum == 2)
      {
        return Centring::body;
      }
      return Centring::none;
    }
    case 3:
      return Centring::r;
    case 4:
      return Centring::face;
    default:
      return Centring::none;
  }
}

/// For the convenience in the following steps, the basis vectors are further transformed to have a specific centring
/// type by multiplying a correction matrix M with M′ for the Laue classes of 2/m and mmm and and the rhombohedral
/// system. For the Laue class 2/m, the basis vectors with the I, A, and B centring types are transformed to those with
/// the C centring type. For the Laue class mmm, those with the A, and B centring types are transformed to those with
/// the C centring type. For the rhombohedral system, a rhombohedrally-centred hexagonal cell is obtained by M' in
/// either the obverse or reverse setting.
///     This is transformed to the primitive rhombohedral cell by Mobv if it is the obverse setting or by Mrev if it is
///     the reverse setting. Only one of M′Mobv or M′Mrev has to be an integer matrix, which is chosen as the
///     transformation matrix. By this, it is known whether the rhombohedrally-centred hexagonal cell obtained by M′ is
///     in the obverse or reverse setting.
SKTransformationMatrix SKPointGroup::computeBasisCorrection(SKTransformationMatrix basis, Centring& centering)
{
  int det = std::abs(basis.determinant());
  Laue lau = _laue;

  // the absolute value of the determinant gives the scale factor by which volume is multiplied under the associated
  // linear transformation, while its sign indicates whether the transformation preserves orientation

  // Number of lattice points per cell (1.2.1 in Hahn 2005 fifth ed.)
  // 1: primitive centred (including R-centered description with ‘rhombohedral axes’)
  // 2: C-face centred, B-face centred, A-face centred, body-centred
  // 3: Rhombohedrally centred (description with ‘hexagonal axes’), Hexagonally centred
  // 4: all-face centred

  switch (det)
  {
    case 1:
      return basis;
    case 2:
      switch (centering)
      {
        case Centring::a_face:
          if (lau == Laue::laue_2m)
          {
            // Tranformation monoclinic A-centring to C-centring (preserving b-axis)
            // Axes a and c are swapped, to keep the same handiness b (to keep Beta obtuse) is made negative
            centering = Centring::c_face;
            return basis * SKTransformationMatrix::monoclinicAtoC;
          }
          else
          {
            centering = Centring::c_face;
            return basis * SKTransformationMatrix::AtoC;
          }
        case Centring::b_face:
          centering = Centring::c_face;
          return basis * SKTransformationMatrix::BtoC;
        case Centring::body:
          if (lau == Laue::laue_2m)
          {
            centering = Centring::c_face;
            return basis * SKTransformationMatrix::monoclinicItoC;
          }
          return basis;
        default:
          return basis;
      }
    case 3:
    {
      SKTransformationMatrix m =
          SKTransformationMatrix::primitiveRhombohedralToTripleHexagonalCell_R2 * basis.adjugate();
      if (m.greatestCommonDivisor() == 3)
      {
        // all elements divisable by 3: reverse detected -> change to obverse
        return basis * SKTransformationMatrix(int3(1, 1, 0), int3(-1, 0, 0), int3(0, 0, 1));
      }
      return basis;
    }
    case 4:
      return basis;
    default:
      return basis;
  }
}

const std::optional<SKTransformationMatrix> SKPointGroup::constructAxes(
    std::vector<SKRotationMatrix> properRotations) const
{
  switch (_laue)
  {
    case Laue::laue_1:
      return SKTransformationMatrix::identity;
    case Laue::laue_2m:
    {
      // Look for all proper rotation matrices of rotation type 2
      std::vector<SKRotationMatrix> properRotationMatrices{};
      std::copy_if(properRotations.begin(), properRotations.end(), std::back_inserter(properRotationMatrices),
                   [](SKRotationMatrix m)
                   { return static_cast<typename std::underlying_type<Laue>::type>(m.type()) == 2; });

      if (properRotationMatrices.empty())
      {
        return std::nullopt;
      }

      SKRotationMatrix properRotationmatrix = properRotationMatrices.front();

      SKTransformationMatrix axes = SKTransformationMatrix();

      // set the rotation axis as the first axis
      axes[1] = properRotationmatrix.rotationAxis();

      // possible candidates for the second axis are vectors that are orthogonal to the axes of rotation
      std::vector<int3> orthogonalAxes = properRotationmatrix.orthogonalToAxisDirection(2);

      // sort accoording to smallest length
      std::sort(orthogonalAxes.begin(), orthogonalAxes.end(),
                [](const int3& a, const int3& b) -> bool { return a.length_squared() < b.length_squared(); });

      // assert(!orthogonalAxes.empty());
      if (orthogonalAxes.empty()) return std::nullopt;

      // the second and thirs axis are the shortest orthogonal axes
      axes[0] = orthogonalAxes[0];
      axes[2] = orthogonalAxes[1];

      if (axes.determinant() < 0)
      {
        return SKTransformationMatrix(axes[2], axes[1], axes[0]);
      }
      return axes;
    }
    case Laue::laue_mmm:
    case Laue::laue_m3:
    case Laue::laue_m3m:
    {
      size_t rotationalTypeForBasis = rotationTypeForBasis[this->_laue];

      // look for all proper rotation matrices of the wanted rotation type
      std::vector<SKRotationMatrix> filteredProperRotations{};
      std::copy_if(
          properRotations.begin(), properRotations.end(), std::back_inserter(filteredProperRotations),
          [rotationalTypeForBasis](SKRotationMatrix m)
          { return static_cast<typename std::underlying_type<Laue>::type>(m.type()) == rotationalTypeForBasis; });

      // take their rotation axes (use a set to avoid duplicates)
      std::unordered_set<int3, int3::hashFunction> allAxes;
      std::transform(filteredProperRotations.begin(), filteredProperRotations.end(),
                     std::inserter(allAxes, allAxes.begin()), [](SKRotationMatrix m) { return m.rotationAxis(); });

      std::vector<int3> sortedAxes{allAxes.begin(), allAxes.end()};
      std::sort(sortedAxes.begin(), sortedAxes.end(),
                [](const int3& a, const int3& b) -> bool { return a.length_squared() < b.length_squared(); });

      if (sortedAxes.size() >= 3)
      {
        SKTransformationMatrix axes = SKTransformationMatrix(sortedAxes[0], sortedAxes[1], sortedAxes[2]);

        if (axes.determinant() < 0)
        {
          return SKTransformationMatrix(axes[0], axes[2], axes[1]);
        }

        return axes;
      }
      return std::nullopt;
    }
    case Laue::laue_4m:
    case Laue::laue_4mmm:
    case Laue::laue_3:
    case Laue::laue_3m:
    case Laue::laue_6m:
    case Laue::laue_6mmm:
    {
      size_t rotationalTypeForBasis = rotationTypeForBasis[this->_laue];

      // look for all proper rotation matrices of the wanted rotation type
      std::vector<SKRotationMatrix> properRotationMatrices{};
      std::copy_if(
          properRotations.begin(), properRotations.end(), std::back_inserter(properRotationMatrices),
          [rotationalTypeForBasis](SKRotationMatrix m)
          { return static_cast<typename std::underlying_type<Laue>::type>(m.type()) == rotationalTypeForBasis; });

      if (properRotationMatrices.empty())
      {
        return std::nullopt;
      }

      SKRotationMatrix properRotationmatrix = properRotationMatrices.front();

      SKTransformationMatrix axes{};

      // set the rotation axis as the first axis
      axes[2] = properRotationmatrix.rotationAxis();

      // possible candidates for the second axis are vectors that are orthogonal to the axes of rotation
      std::vector<int3> orthogonalAxes = properRotationmatrix.orthogonalToAxisDirection(rotationalTypeForBasis);

      for (const int3& orthogonalAxis : orthogonalAxes)
      {
        axes[0] = orthogonalAxis;

        int3 axisVector = properRotationmatrix * orthogonalAxis;

        if ((std::find(SKRotationMatrix::allPossibleRotationAxes.begin(),
                       SKRotationMatrix::allPossibleRotationAxes.end(),
                       axisVector) != SKRotationMatrix::allPossibleRotationAxes.end()) ||
            (std::find(SKRotationMatrix::allPossibleRotationAxes.begin(),
                       SKRotationMatrix::allPossibleRotationAxes.end(),
                       -axisVector) != SKRotationMatrix::allPossibleRotationAxes.end()))
        {
          axes[1] = axisVector;

          // to avoid F-center choice det=4
          if (std::fabs(int3x3(axes[0], axes[1], axes[2]).determinant()) < 4)
          {
            if (axes.determinant() < 0)
            {
              return SKTransformationMatrix(axes[1], axes[0], axes[2]);
            }

            return axes;
          }
        }
      }
      return std::nullopt;
    }
    case Laue::none:
      return std::nullopt;
  }
  return std::nullopt;
}

/// The look-up table used for the construction of M'
///    Laue group           n_e           |N|
///    -1                          1                1
///     2/m                     1                2
///     mmm                  3                2, 2, 2
///     4/m                     1                4
///     4/mmm               2                4, 2
///    -3                          1                3
///    -3m                       2                3, 2
///     6/m                     1                3
///     6/mmm               2                3, 2
///     m-3                     2                3, 2
///     m-3m                  2                3, 4
/// Ref: R.W. Grosse-Kunstleve, "Algorithms for deriving crystallographic space-group information", Acta Cryst. A55,
/// 383-395, 1999
std::map<Laue, size_t> SKPointGroup::rotationTypeForBasis = {
    std::pair(Laue::laue_1, 0),  std::pair(Laue::laue_2m, 2),   std::pair(Laue::laue_mmm, 2),
    std::pair(Laue::laue_4m, 4), std::pair(Laue::laue_4mmm, 4), std::pair(Laue::laue_3, 3),
    std::pair(Laue::laue_3m, 3), std::pair(Laue::laue_6m, 3),   std::pair(Laue::laue_6mmm, 3),
    std::pair(Laue::laue_m3, 2), std::pair(Laue::laue_m3m, 4),
};

std::vector<SKPointGroup> SKPointGroup::pointGroupData =
    std::vector<SKPointGroup>{SKPointGroup(SKRotationalOccuranceTable(0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 0, "", "",
                                           Holohedry::none, Laue::none, false, false),
                              SKPointGroup(SKRotationalOccuranceTable(0, 0, 0, 0, 0, 1, 0, 0, 0, 0), 1, "1", "C1",
                                           Holohedry::triclinic, Laue::laue_1, false, true),
                              SKPointGroup(SKRotationalOccuranceTable(0, 0, 0, 0, 1, 1, 0, 0, 0, 0), 2, "-1", "Ci",
                                           Holohedry::triclinic, Laue::laue_1, true, false),
                              SKPointGroup(SKRotationalOccuranceTable(0, 0, 0, 0, 0, 1, 1, 0, 0, 0), 3, "2", "C2",
                                           Holohedry::monoclinic, Laue::laue_2m, false, true),
                              SKPointGroup(SKRotationalOccuranceTable(0, 0, 0, 1, 0, 1, 0, 0, 0, 0), 4, "m", "Cs",
                                           Holohedry::monoclinic, Laue::laue_2m, false, false),
                              SKPointGroup(SKRotationalOccuranceTable(0, 0, 0, 1, 1, 1, 1, 0, 0, 0), 5, "2/m", "C2h",
                                           Holohedry::monoclinic, Laue::laue_2m, true, false),
                              SKPointGroup(SKRotationalOccuranceTable(0, 0, 0, 0, 0, 1, 3, 0, 0, 0), 6, "222", "D2",
                                           Holohedry::orthorhombic, Laue::laue_mmm, false, true),
                              SKPointGroup(SKRotationalOccuranceTable(0, 0, 0, 2, 0, 1, 1, 0, 0, 0), 7, "mm2", "C2v",
                                           Holohedry::orthorhombic, Laue::laue_mmm, false, false),
                              SKPointGroup(SKRotationalOccuranceTable(0, 0, 0, 3, 1, 1, 3, 0, 0, 0), 8, "mmm", "D2h",
                                           Holohedry::orthorhombic, Laue::laue_mmm, true, false),
                              SKPointGroup(SKRotationalOccuranceTable(0, 0, 0, 0, 0, 1, 1, 0, 2, 0), 9, "4", "C4",
                                           Holohedry::tetragonal, Laue::laue_4m, false, true),
                              SKPointGroup(SKRotationalOccuranceTable(0, 2, 0, 0, 0, 1, 1, 0, 0, 0), 10, "-4", "S4",
                                           Holohedry::tetragonal, Laue::laue_4m, false, false),
                              SKPointGroup(SKRotationalOccuranceTable(0, 2, 0, 1, 1, 1, 1, 0, 2, 0), 11, "4/m", "C4h",
                                           Holohedry::tetragonal, Laue::laue_4m, true, false),
                              SKPointGroup(SKRotationalOccuranceTable(0, 0, 0, 0, 0, 1, 5, 0, 2, 0), 12, "422", "D4",
                                           Holohedry::tetragonal, Laue::laue_4mmm, false, true),
                              SKPointGroup(SKRotationalOccuranceTable(0, 0, 0, 4, 0, 1, 1, 0, 2, 0), 13, "4mm", "C4v",
                                           Holohedry::tetragonal, Laue::laue_4mmm, false, false),
                              SKPointGroup(SKRotationalOccuranceTable(0, 2, 0, 2, 0, 1, 3, 0, 0, 0), 14, "-42m", "D2d",
                                           Holohedry::tetragonal, Laue::laue_4mmm, false, false),
                              SKPointGroup(SKRotationalOccuranceTable(0, 2, 0, 5, 1, 1, 5, 0, 2, 0), 15, "4/mmm", "D4h",
                                           Holohedry::tetragonal, Laue::laue_4mmm, true, false),
                              SKPointGroup(SKRotationalOccuranceTable(0, 0, 0, 0, 0, 1, 0, 2, 0, 0), 16, "3", "C3",
                                           Holohedry::trigonal, Laue::laue_3, false, true),
                              SKPointGroup(SKRotationalOccuranceTable(0, 0, 2, 0, 1, 1, 0, 2, 0, 0), 17, "-3", "C3i",
                                           Holohedry::trigonal, Laue::laue_3, true, false),
                              SKPointGroup(SKRotationalOccuranceTable(0, 0, 0, 0, 0, 1, 3, 2, 0, 0), 18, "32", "D3",
                                           Holohedry::trigonal, Laue::laue_3m, false, true),
                              SKPointGroup(SKRotationalOccuranceTable(0, 0, 0, 3, 0, 1, 0, 2, 0, 0), 19, "3m", "C3v",
                                           Holohedry::trigonal, Laue::laue_3m, false, false),
                              SKPointGroup(SKRotationalOccuranceTable(0, 0, 2, 3, 1, 1, 3, 2, 0, 0), 20, "-3m", "D3d",
                                           Holohedry::trigonal, Laue::laue_3m, true, false),
                              SKPointGroup(SKRotationalOccuranceTable(0, 0, 0, 0, 0, 1, 1, 2, 0, 2), 21, "6", "C6",
                                           Holohedry::hexagonal, Laue::laue_6m, false, true),
                              SKPointGroup(SKRotationalOccuranceTable(2, 0, 0, 1, 0, 1, 0, 2, 0, 0), 22, "-6", "C3h",
                                           Holohedry::hexagonal, Laue::laue_6m, false, false),
                              SKPointGroup(SKRotationalOccuranceTable(2, 0, 2, 1, 1, 1, 1, 2, 0, 2), 23, "6/m", "C6h",
                                           Holohedry::hexagonal, Laue::laue_6m, true, false),
                              SKPointGroup(SKRotationalOccuranceTable(0, 0, 0, 0, 0, 1, 7, 2, 0, 2), 24, "622", "D6",
                                           Holohedry::hexagonal, Laue::laue_6mmm, false, true),
                              SKPointGroup(SKRotationalOccuranceTable(0, 0, 0, 6, 0, 1, 1, 2, 0, 2), 25, "6mm", "C6v",
                                           Holohedry::hexagonal, Laue::laue_6mmm, false, false),
                              SKPointGroup(SKRotationalOccuranceTable(2, 0, 0, 4, 0, 1, 3, 2, 0, 0), 26, "-6m", "D3h",
                                           Holohedry::hexagonal, Laue::laue_6mmm, false, false),
                              SKPointGroup(SKRotationalOccuranceTable(2, 0, 2, 7, 1, 1, 7, 2, 0, 2), 27, "6/mmm", "D6h",
                                           Holohedry::hexagonal, Laue::laue_6mmm, true, false),
                              SKPointGroup(SKRotationalOccuranceTable(0, 0, 0, 0, 0, 1, 3, 8, 0, 0), 28, "23", "T",
                                           Holohedry::cubic, Laue::laue_m3, false, true),
                              SKPointGroup(SKRotationalOccuranceTable(0, 0, 8, 3, 1, 1, 3, 8, 0, 0), 29, "m-3", "Th",
                                           Holohedry::cubic, Laue::laue_m3, true, false),
                              SKPointGroup(SKRotationalOccuranceTable(0, 0, 0, 0, 0, 1, 9, 8, 6, 0), 30, "432", "O",
                                           Holohedry::cubic, Laue::laue_m3m, false, true),
                              SKPointGroup(SKRotationalOccuranceTable(0, 6, 0, 6, 0, 1, 3, 8, 0, 0), 31, "-43m", "Td",
                                           Holohedry::cubic, Laue::laue_m3m, false, false),
                              SKPointGroup(SKRotationalOccuranceTable(0, 6, 8, 9, 1, 1, 9, 8, 6, 0), 32, "m-3m", "Oh",
                                           Holohedry::cubic, Laue::laue_m3m, true, false)};
