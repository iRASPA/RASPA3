#include <gtest/gtest.h>

import std;

import int3;
import double3;
import double3x3;
import skdefinitions;
import skcell;
import sksymmetrycell;
import skpointsymmetryset;
import skrotationmatrix;
import skseitzintegermatrix;
import skseitzmatrix;
import sksymmetryoperationset;
import skintegersymmetryoperationset;
import skspacegroupsetting;
import skasymmetricunit;
import skrotationalchangeofbasis;
import skspacegroupdatabase;
import skspacegroup;
import skstructure;
import skparser;
import skposcarparser;
import blocking_spheres;
import energy_shared_energy_backend;
import energy_shared_blocking_mask;

TEST(SYMMETRY_STRUCTURE_FIXES, set_box_divides_the_box_over_the_replicas)
{
  SKCell cell(10.0, 10.0, 10.0, std::numbers::pi / 2.0, std::numbers::pi / 2.0, std::numbers::pi / 2.0);
  cell.setMaximumReplicas(int3(1, 0, 0));  // two replicas along a
  cell.setBox(double3x3(double3(30.0, 0.0, 0.0), double3(0.0, 12.0, 0.0), double3(0.0, 0.0, 14.0)));
  EXPECT_NEAR(cell.unitCell().ax, 15.0, 1e-12);
  EXPECT_NEAR(cell.unitCell().by, 12.0, 1e-12);
  EXPECT_NEAR(cell.inverseUnitCell().ax, 1.0 / 15.0, 1e-12);
}

TEST(SYMMETRY_STRUCTURE_FIXES, angle_setters_keep_the_inverse_unit_cell)
{
  SKCell cell(10.0, 11.0, 12.0, std::numbers::pi / 2.0, std::numbers::pi / 2.0, std::numbers::pi / 2.0);
  cell.setAlphaAngle(80.0 * std::numbers::pi / 180.0);
  const double3x3 product = cell.unitCell() * cell.inverseUnitCell();
  EXPECT_NEAR(product.ax, 1.0, 1e-12);
  EXPECT_NEAR(product.bz, 0.0, 1e-12);
  EXPECT_NEAR(product.cz, 1.0, 1e-12);
}

TEST(SYMMETRY_STRUCTURE_FIXES, equal_seitz_matrices_hash_equally)
{
  const SKSeitzIntegerMatrix a(SKRotationMatrix::identity, int3(1, 2, 3));
  const SKSeitzIntegerMatrix b(SKRotationMatrix::identity, int3(25, 26, 27));
  ASSERT_TRUE(a == b);
  EXPECT_EQ(SKSeitzIntegerMatrix::hashFunction{}(a), SKSeitzIntegerMatrix::hashFunction{}(b));

  // a negative translation (a product with -1 rotation entries) is the same operation as its positive image: the
  // signed remainder left -6 as -6 instead of 18, so these compared unequal and hashed apart
  const SKSeitzIntegerMatrix negative(SKRotationMatrix::identity, int3(-6, -1, 23));
  const SKSeitzIntegerMatrix positive(SKRotationMatrix::identity, int3(18, 23, -1));
  EXPECT_TRUE(negative == positive);
  EXPECT_EQ(SKSeitzIntegerMatrix::hashFunction{}(negative), SKSeitzIntegerMatrix::hashFunction{}(positive));
  EXPECT_EQ(negative.normalizedTranslation(), int3(18, 23, 23));
  EXPECT_FALSE(negative == SKSeitzIntegerMatrix(SKRotationMatrix::identity, int3(6, 1, 23)));
}

TEST(SYMMETRY_STRUCTURE_FIXES, b_centred_settings_are_b_face_centred)
{
  EXPECT_EQ(SKSpaceGroupDataBase::spaceGroupData[13].HallNumber(), 13uz);
  EXPECT_EQ(SKSpaceGroupDataBase::spaceGroupData[13].centring(), Centring::b_face);
  EXPECT_EQ(SKSpaceGroupDataBase::spaceGroupData[15].HallNumber(), 15uz);
  EXPECT_EQ(SKSpaceGroupDataBase::spaceGroupData[15].centring(), Centring::b_face);
}

TEST(SYMMETRY_STRUCTURE_FIXES, space_group_search_without_atoms_finds_nothing)
{
  const double3x3 unitCell(double3(10.0, 0.0, 0.0), double3(0.0, 10.0, 0.0), double3(0.0, 0.0, 10.0));
  EXPECT_FALSE(SKSpaceGroup::findSpaceGroup(unitCell, {}, false, 1e-2).has_value());
}

TEST(SYMMETRY_STRUCTURE_FIXES, poscar_blank_comment_negative_volume_and_scaled_cartesian_coordinates)
{
  // blank comment line, a negative scaling factor (a cell volume of 1000 A^3 from a 2x2x2 cube), Cartesian
  // coordinates scaled like the lattice
  const std::string content =
      "\n"
      "-1000.0\n"
      "2.0 0.0 0.0\n"
      "0.0 2.0 0.0\n"
      "0.0 0.0 2.0\n"
      "Si\n"
      "1\n"
      "Cartesian\n"
      "0.5 1.0 1.5\n";
  SKPOSCARParser parser(content);
  parser.startParsing();
  const double3x3 cell = parser.movies().front().front()->cell->unitCell();
  EXPECT_NEAR(cell.ax, 10.0, 1e-9);
  EXPECT_NEAR(cell.determinant(), 1000.0, 1e-6);

  const std::vector<std::tuple<double3, std::size_t, double>> atoms = parser.firstTestFrame();
  ASSERT_EQ(atoms.size(), 1uz);
  // the Cartesian position (0.5, 1, 1.5) scaled by 5 is (2.5, 5, 7.5): fractional (0.25, 0.5, 0.75)
  EXPECT_NEAR(std::get<0>(atoms[0]).x, 0.25, 1e-12);
  EXPECT_NEAR(std::get<0>(atoms[0]).y, 0.5, 1e-12);
  EXPECT_NEAR(std::get<0>(atoms[0]).z, 0.75, 1e-12);
}

TEST(SYMMETRY_STRUCTURE_FIXES, a_backend_with_blocked_pockets_is_another_backend_for_the_field_cache)
{
  EnergyBackend backend;
  backend.name = "cpu";
  const std::vector<BlockingSphere> first{BlockingSphere{double3(0.1, 0.2, 0.3), 1.5}};
  const std::vector<BlockingSphere> second{BlockingSphere{double3(0.1, 0.2, 0.3), 2.0}};
  const EnergyBackend blocked = withBlockedPockets(backend, first, 1.0, 2.0);
  const EnergyBackend otherBlocked = withBlockedPockets(backend, second, 1.0, 2.0);
  EXPECT_TRUE(backend.blocking.empty());
  EXPECT_FALSE(blocked.blocking.empty());
  EXPECT_NE(blocked.blocking, otherBlocked.blocking);
  EXPECT_EQ(blocked.name, backend.name);
}

TEST(SYMMETRY_STRUCTURE_FIXES, overlap_test_keeps_the_signs_of_the_fractional_difference)
{
  // a sheared cell (symmetric matrix, so the row/column convention does not matter): the Cartesian distance of a
  // fractional difference depends on the relative signs of its components
  const double3x3 lattice(double3(10.0, -10.0, 0.0), double3(-10.0, 20.0, 0.0), double3(0.0, 0.0, 10.0));
  const double precision = 0.15;

  // the difference (0.01, -0.99, 0) wraps to (0.01, 0.01, 0), at 0.1 from each other; the absolute value before
  // the wrap gave (0.01, -0.01, 0), at 0.36, and missed the overlap
  const double3 a(0.5, 0.005, 0.5), b(0.49, 0.995, 0.5);
  EXPECT_TRUE(SKSymmetryCell::isOverlap(a, b, lattice, precision));
  // the difference (0.01, -0.01, 0) is at 0.36; the absolute value made it (0.01, 0.01, 0), a false overlap
  const double3 c(0.5, 0.5, 0.5), d(0.49, 0.51, 0.5);
  EXPECT_FALSE(SKSymmetryCell::isOverlap(c, d, lattice, precision));

  // trim merges the first pair and keeps the second
  using Atom = std::tuple<double3, std::size_t, double>;
  EXPECT_EQ(SKSymmetryCell::trim({Atom{a, 1, 1.0}, Atom{b, 1, 1.0}}, lattice, lattice, false, precision).size(), 1uz);
  EXPECT_EQ(SKSymmetryCell::trim({Atom{c, 1, 1.0}, Atom{d, 1, 1.0}}, lattice, lattice, false, precision).size(), 2uz);
}

TEST(SYMMETRY_STRUCTURE_FIXES, h_centring_adds_the_two_ab_plane_translations)
{
  // the third H-centring shift was (0, 2/3, 1/3), its components shifted by one slot out of the ab-plane
  const SKSymmetryOperationSet identity({SKSeitzMatrix(SKRotationMatrix::identity, double3(0.0, 0.0, 0.0))});
  const SKSymmetryOperationSet centred = identity.addingCenteringOperations(Centring::h);
  ASSERT_EQ(centred.operations.size(), 3uz);
  const double3 expected[3] = {double3(0.0, 0.0, 0.0), double3(2.0 / 3.0, 1.0 / 3.0, 0.0),
                               double3(1.0 / 3.0, 2.0 / 3.0, 0.0)};
  for (std::size_t i = 0; i < 3; ++i)
  {
    EXPECT_NEAR(centred.operations[i].translation.x, expected[i].x, 1e-12) << "shift " << i;
    EXPECT_NEAR(centred.operations[i].translation.y, expected[i].y, 1e-12) << "shift " << i;
    EXPECT_NEAR(centred.operations[i].translation.z, expected[i].z, 1e-12) << "shift " << i;
  }
}

namespace
{
// the operations of a tabulated setting, as a CIF would list them (translations in fractions of the cell)
std::vector<SKSeitzMatrix> operationsOfSetting(std::size_t hallNumber)
{
  std::vector<SKSeitzMatrix> operations;
  for (const SKSeitzIntegerMatrix& m : SKSpaceGroup(hallNumber).spaceGroupSetting().fullSeitzMatrices().operations)
  {
    operations.emplace_back(m.rotation, double3(m.translation.x, m.translation.y, m.translation.z) / 24.0);
  }
  return operations;
}
}  // namespace

TEST(SYMMETRY_STRUCTURE_FIXES, hall_number_search_requires_the_order_and_prefers_the_listed_setting)
{
  const double3x3 cubic(double3(10.0, 0.0, 0.0), double3(0.0, 10.0, 0.0), double3(0.0, 0.0, 10.0));
  const double3x3 orthorhombic(double3(5.0, 0.0, 0.0), double3(0.0, 6.0, 0.0), double3(0.0, 0.0, 7.0));

  // P -1 (Hall 2) without an international number: P 1 has no generators and used to match every primitive group
  const std::vector<SKSeitzMatrix> pBar1 = operationsOfSetting(2);
  ASSERT_EQ(pBar1.size(), 2uz);
  EXPECT_EQ(SKSpaceGroup::HallNumberFromSymmetryOperations(pBar1, cubic, Centring::primitive, std::nullopt), 2uz);

  // F d -3 m: origin choice 1 (Hall 525) and 2 (Hall 526) differ by an origin shift; the listed one is returned,
  // where the standard setting (choice 2) used to win for both
  const std::vector<SKSeitzMatrix> origin1 = operationsOfSetting(525);
  const std::vector<SKSeitzMatrix> origin2 = operationsOfSetting(526);
  ASSERT_EQ(origin1.size(), 192uz);
  EXPECT_EQ(SKSpaceGroup::HallNumberFromSymmetryOperations(origin1, cubic, Centring::face, std::nullopt), 525uz);
  EXPECT_EQ(SKSpaceGroup::HallNumberFromSymmetryOperations(origin2, cubic, Centring::face, std::nullopt), 526uz);
  EXPECT_EQ(SKSpaceGroup::HallNumberFromSymmetryOperations(origin1, cubic, Centring::face, 227), 525uz);

  // the generators of F d -3 (227's subgroup 203) are among the operations of F d -3 m, but the order differs
  EXPECT_FALSE(SKSpaceGroup::HallNumberFromSymmetryOperations(origin2, cubic, Centring::face, 203).has_value());

  // P b n m (Hall 294) is P n m a (Hall 292, standard) in another axis setting: the listed setting is returned
  const std::vector<SKSeitzMatrix> pbnm = operationsOfSetting(294);
  ASSERT_EQ(pbnm.size(), 8uz);
  EXPECT_EQ(SKSpaceGroup::HallNumberFromSymmetryOperations(pbnm, orthorhombic, Centring::primitive, std::nullopt),
            294uz);
  EXPECT_EQ(SKSpaceGroup::HallNumberFromSymmetryOperations(pbnm, orthorhombic, Centring::primitive, 62), 294uz);
  EXPECT_EQ(SKSpaceGroup::HallNumberFromSymmetryOperations(operationsOfSetting(292), orthorhombic, Centring::primitive,
                                                           std::nullopt),
            292uz);
}

TEST(SYMMETRY_STRUCTURE_FIXES, monoclinic_change_of_basis_relabels_c_and_i_centrings)
{
  // space group 5: C 1 2 1 (Hall 9, C-centred) and I 1 2 1 (Hall 11, I-centred) are the same group in two of the
  // monoclinic bases; the centring of a transformed setting used to stay C or I, so neither matched the other
  const std::vector<SKSeitzMatrix> cCentred = operationsOfSetting(9);
  const std::vector<SKSeitzMatrix> iCentred = operationsOfSetting(11);
  ASSERT_EQ(cCentred.size(), 4uz);
  ASSERT_EQ(iCentred.size(), 4uz);
  const auto& bases = SKRotationalChangeOfBasis::changeOfMonoclinicCentering;

  EXPECT_FALSE(
      SKSpaceGroup::getOriginShift(9, Centring::body, SKRotationalChangeOfBasis::identity, iCentred, 1e-2).has_value());
  EXPECT_TRUE(
      std::any_of(bases.begin(), bases.end(), [&](const SKRotationalChangeOfBasis& basis)
                  { return SKSpaceGroup::getOriginShift(9, Centring::body, basis, iCentred, 1e-2).has_value(); }));
  EXPECT_TRUE(
      std::any_of(bases.begin(), bases.end(), [&](const SKRotationalChangeOfBasis& basis)
                  { return SKSpaceGroup::getOriginShift(11, Centring::c_face, basis, cCentred, 1e-2).has_value(); }));

  // each setting is still identified as listed
  const double3x3 monoclinic(double3(5.0, 0.0, 0.0), double3(0.0, 6.0, 0.0), double3(-1.5, 0.0, 7.0));
  EXPECT_EQ(SKSpaceGroup::HallNumberFromSymmetryOperations(iCentred, monoclinic, Centring::body, 5), 11uz);
  EXPECT_EQ(SKSpaceGroup::HallNumberFromSymmetryOperations(cCentred, monoclinic, Centring::c_face, 5), 9uz);
}

TEST(SYMMETRY_STRUCTURE_FIXES, seitz_matrices_of_a_centrosymmetric_setting_are_not_doubled)
{
  // P -1 (Hall 2), C 1 2/c 1 (Hall 90) and F d -3 m with the inversion centre off the origin (Hall 525) and on it
  // (Hall 526): the encoded operations are the whole group; appended inverted copies doubled it with operations
  // that were not in the group
  for (const std::size_t hallNumber : {2uz, 90uz, 525uz, 526uz})
  {
    const SKSpaceGroupSetting setting = SKSpaceGroup(hallNumber).spaceGroupSetting();
    const std::vector<SKSeitzIntegerMatrix> matrices = setting.SeitzMatricesWithoutTranslation();
    const std::unordered_set<SKSeitzIntegerMatrix, SKSeitzIntegerMatrix::hashFunction> full =
        setting.fullSeitzMatrices().operations;
    EXPECT_EQ(matrices.size(), setting.order()) << "Hall " << hallNumber;
    EXPECT_EQ(full.size(), setting.order()) << "Hall " << hallNumber;
    for (const SKSeitzIntegerMatrix& m : matrices)
    {
      EXPECT_TRUE(full.contains(m)) << "Hall " << hallNumber;
    }
  }
}

TEST(SYMMETRY_STRUCTURE_FIXES, lattice_translation_and_inversion_centre_strings_are_in_twenty_fourths)
{
  // the stored components are in units of 1/24; dividing by 12 printed (1/2, 1/2, 0) as (1/1, 1/1, 0)
  EXPECT_EQ(SKSpaceGroup::latticeTranslationStrings(526),
            (std::vector<std::string>{"(0,0,0)", "(0,1/2,1/2)", "(1/2,0,1/2)", "(1/2,1/2,0)"}));
  EXPECT_EQ(SKSpaceGroup::latticeTranslationStrings(436),
            (std::vector<std::string>{"(0,0,0)", "(2/3,1/3,1/3)", "(1/3,2/3,2/3)", ""}));
  EXPECT_EQ(SKSpaceGroup::inversionCenterString(525), "(1/8,1/8,1/8)");
  EXPECT_EQ(SKSpaceGroup::inversionCenterString(526), "(0,0,0)");
}

TEST(SYMMETRY_STRUCTURE_FIXES, asymmetric_unit_bounds_are_open_at_one_and_may_lie_below_minus_one_half)
{
  // bounds in 1/48 with an odd value for an open boundary: {24, 49} is [1/2, 1); 49 used to be read as 49/48
  SKAsymmetricUnit upperHalf({24, 49}, {0, 24}, {0, 24});
  EXPECT_TRUE(upperHalf.contains(double3(0.5, 0.1, 0.1)));
  EXPECT_TRUE(upperHalf.contains(double3(0.999, 0.1, 0.1)));
  EXPECT_FALSE(upperHalf.contains(double3(1.0, 0.1, 0.1)));
  EXPECT_FALSE(upperHalf.contains(double3(0.0, 0.1, 0.1)));
  EXPECT_FALSE(upperHalf.contains(double3(0.49, 0.1, 0.1)));

  // {-48, 0} is [-1, 0], the whole cell: points whose image lies below -1/2 used to be rejected
  SKAsymmetricUnit wholeCell({-48, 0}, {0, 24}, {0, 24});
  for (const double x : {0.1, 0.25, 0.5, 0.75, -0.75, 1.3})
  {
    EXPECT_TRUE(wholeCell.contains(double3(x, 0.1, 0.1))) << x;
  }

  // {-32, 8} is [-2/3, 1/6]
  SKAsymmetricUnit skewed({-32, 8}, {0, 24}, {0, 24});
  EXPECT_TRUE(skewed.contains(double3(0.5, 0.1, 0.1)));   // image -1/2
  EXPECT_TRUE(skewed.contains(double3(0.4, 0.1, 0.1)));   // image -0.6
  EXPECT_TRUE(skewed.contains(double3(0.1, 0.1, 0.1)));   // itself
  EXPECT_FALSE(skewed.contains(double3(0.3, 0.1, 0.1)));  // 0.3 > 1/6 and -0.7 < -2/3
  EXPECT_FALSE(skewed.contains(double3(0.2, 0.1, 0.1)));
}

TEST(SYMMETRY_STRUCTURE_FIXES, found_space_group_carries_its_change_of_basis_and_wrapped_atoms)
{
  // the result used to carry a default-constructed (zero) change of basis instead of the one found
  const double3x3 cubic(double3(10.0, 0.0, 0.0), double3(0.0, 10.0, 0.0), double3(0.0, 0.0, 10.0));
  const auto generate = [&](std::size_t hallNumber, double3 site)
  {
    std::vector<std::tuple<double3, std::size_t, double>> atoms;
    for (const SKSeitzIntegerMatrix& m : SKSpaceGroup(hallNumber).spaceGroupSetting().fullSeitzMatrices().operations)
    {
      const double3 position = double3::fract(m * site);
      if (std::none_of(atoms.begin(), atoms.end(),
                       [&](const auto& atom)
                       {
                         const double3 d = std::get<0>(atom) - position;
                         return (d - double3(std::rint(d.x), std::rint(d.y), std::rint(d.z))).length_squared() < 1e-12;
                       }))
      {
        atoms.emplace_back(position, 1, 1.0);
      }
    }
    return atoms;
  };

  // P m -3 m from one atom at the origin, and P a -3 (the special case of the search) from its 8c site (x, x, x)
  for (const std::size_t number : {221uz, 205uz})
  {
    const std::size_t hallNumber = SKSpaceGroupDataBase::spaceGroupHallData[number].front();
    const double3 site = number == 221 ? double3(0.0, 0.0, 0.0) : double3(0.1, 0.1, 0.1);
    const std::optional<SKSpaceGroup::FoundSpaceGroupInfo> found =
        SKSpaceGroup::findSpaceGroup(cubic, generate(hallNumber, site), false, 1e-2);
    ASSERT_TRUE(found.has_value()) << number;
    EXPECT_EQ(SKSpaceGroup(found->HallNumber).spaceGroupSetting().number(), number);
    // a change of basis is a rotation with an inverse, not a zero matrix
    EXPECT_TRUE(found->changeOfBasis.rotationMatrix * found->changeOfBasis.inverseRotationMatrix ==
                SKRotationMatrix::identity)
        << number;
    for (const auto& atom : found->atoms)
    {
      const double3 p = std::get<0>(atom);
      EXPECT_TRUE(p.x >= -1e-2 && p.x < 1.0 + 1e-2 && p.y >= -1e-2 && p.y < 1.0 + 1e-2 && p.z >= -1e-2 &&
                  p.z < 1.0 + 1e-2)
          << number << " " << p.x << " " << p.y << " " << p.z;
    }
  }
}

TEST(SYMMETRY_STRUCTURE_FIXES, lattice_symmetry_counts_follow_the_metric)
{
  // the holohedry of each lattice at a precision of 1e-2: cubic 48, tetragonal 16, orthorhombic 8, triclinic 2.
  // A precision of 3 or more would have passed every angle test, and lengths differing by less than 3
  const double3x3 cubic(double3(10.0, 0.0, 0.0), double3(0.0, 10.0, 0.0), double3(0.0, 0.0, 10.0));
  const double3x3 tetragonal(double3(10.0, 0.0, 0.0), double3(0.0, 10.0, 0.0), double3(0.0, 0.0, 12.0));
  const double3x3 orthorhombic(double3(10.0, 0.0, 0.0), double3(0.0, 11.0, 0.0), double3(0.0, 0.0, 12.0));
  const double3x3 triclinic(double3(10.0, 0.0, 0.0), double3(1.0, 11.0, 0.0), double3(0.5, 1.5, 12.0));
  EXPECT_EQ(SKSymmetryCell::findLatticeSymmetry(cubic, 1e-2).rotations().size(), 48uz);
  EXPECT_EQ(SKSymmetryCell::findLatticeSymmetry(tetragonal, 1e-2).rotations().size(), 16uz);
  EXPECT_EQ(SKSymmetryCell::findLatticeSymmetry(orthorhombic, 1e-2).rotations().size(), 8uz);
  EXPECT_EQ(SKSymmetryCell::findLatticeSymmetry(triclinic, 1e-2).rotations().size(), 2uz);
}
