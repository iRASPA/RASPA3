#include <gtest/gtest.h>

import std;

import int3;
import double3;
import double3x3;
import float4;
import skdefinitions;
import skcell;
import skatom;
import skboundingbox;
import sksymmetrycell;
import skpointsymmetryset;
import skrotationmatrix;
import skseitzintegermatrix;
import skseitzmatrix;
import sksymmetryoperationset;
import skintegersymmetryoperationset;
import skspacegroupsetting;
import sktransformationmatrix;
import skasymmetricunit;
import skrotationalchangeofbasis;
import skspacegroupdatabase;
import skapolloniusdiagram;
import stringutils;
import skelement;
import skpointgroup;
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

TEST(SYMMETRY_STRUCTURE_FIXES, left_handed_cells_have_positive_volume_widths_and_replicas)
{
  // the same box with c reversed: a negative determinant, but the volume, the distances between opposite faces
  // and the number of replicas needed for a cutoff are those of the right-handed box (they came out negative)
  const double3x3 right(double3(10.0, 0.0, 0.0), double3(2.0, 12.0, 0.0), double3(1.0, 3.0, 14.0));
  const double3x3 left(double3(10.0, 0.0, 0.0), double3(2.0, 12.0, 0.0), double3(-1.0, -3.0, -14.0));
  ASSERT_LT(left.determinant(), 0.0);
  SKCell rightCell(right);
  SKCell leftCell(left);

  EXPECT_NEAR(rightCell.volume(), 10.0 * 12.0 * 14.0, 1e-9);
  EXPECT_NEAR(leftCell.volume(), rightCell.volume(), 1e-9);

  const double3 rightWidths = rightCell.perpendicularWidths();
  const double3 leftWidths = leftCell.perpendicularWidths();
  EXPECT_GT(rightWidths.x, 0.0);
  EXPECT_NEAR(leftWidths.x, rightWidths.x, 1e-9);
  EXPECT_NEAR(leftWidths.y, rightWidths.y, 1e-9);
  EXPECT_NEAR(leftWidths.z, rightWidths.z, 1e-9);

  const int3 rightReplicas = rightCell.numberOfReplicas(12.0);
  const int3 leftReplicas = leftCell.numberOfReplicas(12.0);
  EXPECT_EQ(rightReplicas.x, 3);  // the width along a is 10, twice the cutoff is 24
  EXPECT_EQ(leftReplicas.x, rightReplicas.x);
  EXPECT_EQ(leftReplicas.y, rightReplicas.y);
  EXPECT_EQ(leftReplicas.z, rightReplicas.z);
}

TEST(SYMMETRY_STRUCTURE_FIXES, find_primitive_trims_at_the_precision_the_group_was_found_at)
{
  using Atom = std::tuple<double3, std::size_t, double>;
  const double3x3 cubic(double3(10.0, 0.0, 0.0), double3(0.0, 10.0, 0.0), double3(0.0, 0.0, 10.0));
  const double requested = 0.02;

  // a P m -3 m orbit with noise of a quarter of the requested precision, which makes the search shrink its precision
  // (to about 0.017 here) and settle on P 1, and a second orbit of the same type displaced from the first by 0.018
  // with the same noise, so that every pair is exactly 0.018 apart: distinct at the precision the group was found
  // at, coincident at the requested one. The trim into the primitive cell used the requested precision and merged
  // each pair, returning 48 atoms for a 96-atom cell. Seeds are tried until the shrink lands below the displacement.
  const double displacement = 0.9 * requested;
  std::vector<Atom> atoms;
  std::optional<SKSpaceGroup::FoundSpaceGroupInfo> found;
  for (std::uint32_t seed = 1; seed <= 50 && !found; ++seed)
  {
    std::mt19937 rng(seed);
    std::normal_distribution<double> noise(0.0, 0.25 * requested / 10.0);
    atoms.clear();
    for (const SKSeitzIntegerMatrix& m : SKSpaceGroup(517).spaceGroupSetting().fullSeitzMatrices().operations)
    {
      const double3 position =
          double3::fract(m * double3(0.12, 0.23, 0.31) + double3(noise(rng), noise(rng), noise(rng)));
      atoms.emplace_back(position, 1, 1.0);
      atoms.emplace_back(double3::fract(position + double3(displacement / 10.0, 0.0, 0.0)), 1, 1.0);
    }
    found = SKSpaceGroup::findSpaceGroup(cubic, atoms, false, requested);
    if (found && !(found->symmetryPrecision < displacement && found->atoms.size() == atoms.size()))
    {
      found = std::nullopt;
    }
  }
  ASSERT_TRUE(found.has_value()) << "no seed made the search shrink its precision below the displacement";
  EXPECT_EQ(found->HallNumber, 1uz);

  const std::optional<SKSpaceGroup::FoundPrimitiveCellInfo> primitive =
      SKSpaceGroup::SKFindPrimitive(cubic, atoms, false, requested);
  ASSERT_TRUE(primitive.has_value());
  // P 1: the primitive cell is the cell itself and holds all 96 atoms
  EXPECT_NEAR(std::fabs(primitive->cell.unitCell().determinant()), std::fabs(cubic.determinant()), 1e-6);
  EXPECT_EQ(primitive->atoms.size(), atoms.size()) << "found precision " << found->symmetryPrecision;
}

TEST(SYMMETRY_STRUCTURE_FIXES, cell_angles_of_collinear_axes_are_not_nan)
{
  // the build uses -ffast-math, under which std::isnan folds to false, a NaN passes EXPECT_NEAR, and even a test of
  // the bit pattern is folded away; the text the C library prints for the value is beyond the optimizer's reach
  const auto isNaN = [](double value) { return std::to_string(value).find("nan") != std::string::npos; };

  // a metric tensor whose off-diagonal elements exceed the product of the two lengths by one ulp: the cosine is
  // 1 + 2e-16 (and -1 - 2e-16 for the antiparallel pair), of which acos is NaN unless clamped
  const double overOne = std::nextafter(1.0, 2.0);
  const SKSymmetryCell fromTensor = SKSymmetryCell::createFromMetricTensor(
      double3x3(double3(1.0, -overOne, 0.0), double3(-overOne, 1.0, overOne), double3(0.0, overOne, 1.0)));
  EXPECT_FALSE(isNaN(fromTensor.alpha()));
  EXPECT_NEAR(fromTensor.alpha(), 0.0, 1e-7);
  EXPECT_FALSE(isNaN(fromTensor.gamma()));
  EXPECT_NEAR(fromTensor.gamma(), std::numbers::pi, 1e-7);

  // parallel and antiparallel axes from their vectors: the cosine rounds to within an ulp of +-1 either way
  const double3 b(0.1 * 3.0, 0.2 * 3.0, 0.3 * 3.0);
  const double3 c(0.1 * 3.0 * 10.0 * 0.37, 0.2 * 3.0 * 10.0 * 0.37, 0.3 * 3.0 * 10.0 * 0.37);
  for (const double3x3& axes : {double3x3(double3(10.0, 0.0, 0.0), b, c), double3x3(double3(10.0, 0.0, 0.0), b, -c),
                                double3x3(double3(10.0, 0.0, 0.0), -overOne * b, b)})
  {
    const SKCell cell(axes);
    EXPECT_FALSE(isNaN(cell.alpha()));
    EXPECT_FALSE(isNaN(cell.beta()));
    EXPECT_FALSE(isNaN(cell.gamma()));
    const SKSymmetryCell fromCell = SKSymmetryCell::createFromUnitCell(axes);
    EXPECT_FALSE(isNaN(fromCell.alpha()));
    EXPECT_FALSE(isNaN(fromCell.beta()));
    EXPECT_FALSE(isNaN(fromCell.gamma()));
    EXPECT_NEAR(fromCell.alpha(), cell.alpha(), 1e-12);
  }
}

TEST(SYMMETRY_STRUCTURE_FIXES, angle_setters_keep_a_rotated_cell_oriented_and_setters_update_the_bounding_box)
{
  const double degree = std::numbers::pi / 180.0;
  // a monoclinic cell rotated by 30 degrees about z and then 40 degrees about x
  const double3x3 canonical =
      double3x3::fromCellParameters(10.0, 11.0, 12.0, 80.0 * degree, 90.0 * degree, 90.0 * degree);
  const double cz = std::cos(30.0 * degree), sz = std::sin(30.0 * degree);
  const double cx = std::cos(40.0 * degree), sx = std::sin(40.0 * degree);
  const double3x3 aboutZ(double3(cz, sz, 0.0), double3(-sz, cz, 0.0), double3(0.0, 0.0, 1.0));
  const double3x3 aboutX(double3(1.0, 0.0, 0.0), double3(0.0, cx, sx), double3(0.0, -sx, cx));
  const double3x3 rotation = aboutX * aboutZ;
  const double3x3 rotated = rotation * canonical;
  SKCell cell(rotated);

  // setting alpha used to put the new c axis in the canonical frame (a along x) while a and b stayed rotated, which
  // changed beta and gamma as well; now only alpha changes and the a and b axes stay where they were
  cell.setAlphaAngle(70.0 * degree);
  EXPECT_NEAR(cell.a(), 10.0, 1e-10);
  EXPECT_NEAR(cell.b(), 11.0, 1e-10);
  EXPECT_NEAR(cell.c(), 12.0, 1e-10);
  EXPECT_NEAR(cell.alpha(), 70.0 * degree, 1e-10);
  EXPECT_NEAR(cell.beta(), 90.0 * degree, 1e-10);
  EXPECT_NEAR(cell.gamma(), 90.0 * degree, 1e-10);
  EXPECT_NEAR((cell.unitCell()[0] - rotated[0]).length(), 0.0, 1e-10);
  EXPECT_NEAR((cell.unitCell()[1] - rotated[1]).length(), 0.0, 1e-10);

  // setting gamma rebuilds b and c; the cell is the rotated image of the cell with the new parameters
  cell.gamma(100.0 * degree);
  cell.setBetaAngle(85.0 * degree);
  const double3x3 expected =
      rotation * double3x3::fromCellParameters(10.0, 11.0, 12.0, 70.0 * degree, 85.0 * degree, 100.0 * degree);
  for (std::size_t i = 0; i < 3; ++i)
    for (std::size_t j = 0; j < 3; ++j) EXPECT_NEAR(cell.unitCell()[i][j], expected[i][j], 1e-9) << i << " " << j;
  EXPECT_NEAR((cell.unitCell()[0] - rotated[0]).length(), 0.0, 1e-10);

  // a cell in the canonical frame stays in it, with exact zeros
  SKCell canonicalCell(10.0, 11.0, 12.0, 90.0 * degree, 90.0 * degree, 90.0 * degree);
  canonicalCell.setAlphaAngle(80.0 * degree);
  EXPECT_EQ(canonicalCell.unitCell().ay, 0.0);
  EXPECT_EQ(canonicalCell.unitCell().az, 0.0);
  EXPECT_EQ(canonicalCell.unitCell().bz, 0.0);
  EXPECT_NEAR(canonicalCell.alpha(), 80.0 * degree, 1e-12);

  // the bounding box followed no setter: it encloses the full cell after each of them
  SKCell box(10.0, 10.0, 10.0, 90.0 * degree, 90.0 * degree, 90.0 * degree);
  box.setMaximumReplicas(int3(1, 0, 0));
  EXPECT_NEAR(box.boundingBox().maximum().x, 20.0, 1e-12);
  box.setMinimumReplicaY(-1);
  EXPECT_NEAR(box.boundingBox().minimum().y, -10.0, 1e-12);
  box.setLengthA(5.0);
  EXPECT_NEAR(box.boundingBox().maximum().x, 10.0, 1e-12);
  box.c(7.0);
  EXPECT_NEAR(box.boundingBox().maximum().z, 7.0, 1e-12);
  box.setGammaAngle(60.0 * degree);
  EXPECT_EQ(box.boundingBox().maximum(), box.enclosingBoundingBox().maximum());
  EXPECT_EQ(box.boundingBox().minimum(), box.enclosingBoundingBox().minimum());
  box.setUnitCell(double3x3(double3(3.0, 0.0, 0.0), double3(0.0, 4.0, 0.0), double3(0.0, 0.0, 5.0)));
  EXPECT_NEAR(box.boundingBox().maximum().x, 6.0, 1e-12);
  EXPECT_NEAR(box.boundingBox().minimum().y, -4.0, 1e-12);
  box.setBox(double3x3(double3(8.0, 0.0, 0.0), double3(0.0, 8.0, 0.0), double3(0.0, 0.0, 8.0)));
  EXPECT_NEAR(box.boundingBox().maximum().x, 8.0, 1e-12);
  EXPECT_NEAR(box.boundingBox().minimum().y, -4.0, 1e-12);
}

TEST(SYMMETRY_STRUCTURE_FIXES, symmetrize_without_operations_leaves_the_atoms_alone)
{
  using Atom = std::tuple<double3, std::size_t, double>;
  const double3x3 cubic(double3(10.0, 0.0, 0.0), double3(0.0, 10.0, 0.0), double3(0.0, 0.0, 10.0));
  const std::vector<Atom> atoms{Atom{double3(0.1, 0.2, 0.3), 1, 1.0}, Atom{double3(0.6, 0.7, 0.8), 2, 0.5}};

  // an empty set has no operation that maps an atom onto itself; the average over none was 0/0 (NaN coordinates,
  // which -ffast-math hides from std::isnan, hence the text check)
  const std::vector<Atom> unchanged = SKIntegerSymmetryOperationSet().symmetrize(cubic, atoms, 1e-2);
  ASSERT_EQ(unchanged.size(), atoms.size());
  for (std::size_t i = 0; i < atoms.size(); ++i)
  {
    const double3 position = std::get<0>(unchanged[i]);
    EXPECT_EQ(std::to_string(position.x).find("nan"), std::string::npos);
    EXPECT_EQ(position, std::get<0>(atoms[i]));
    EXPECT_EQ(std::get<1>(unchanged[i]), std::get<1>(atoms[i]));
    EXPECT_EQ(std::get<2>(unchanged[i]), std::get<2>(atoms[i]));
  }

  // with a group the identity always matches, so a general position is returned as it is and a position just off a
  // mirror of P m -3 m is put onto it
  SKIntegerSymmetryOperationSet operations = SKSpaceGroup(517).spaceGroupSetting().fullSeitzMatrices();
  const std::vector<Atom> general = operations.symmetrize(cubic, {Atom{double3(0.1, 0.2, 0.3), 1, 1.0}}, 1e-2);
  EXPECT_NEAR((std::get<0>(general.front()) - double3(0.1, 0.2, 0.3)).length(), 0.0, 1e-12);
  const std::vector<Atom> onMirror = operations.symmetrize(cubic, {Atom{double3(0.5004, 0.2, 0.3), 1, 1.0}}, 1e-2);
  EXPECT_NEAR(std::get<0>(onMirror.front()).x, 0.5, 1e-12);
  EXPECT_NEAR(std::get<0>(onMirror.front()).y, 0.2, 1e-12);
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

  // every setting's centring agrees with the lattice letter of its Hermann-Mauguin symbol; 38 B settings were
  // stored as primitive. R is left out: in rhombohedral axes an R lattice is primitive.
  const std::map<char, Centring> letters{{'P', Centring::primitive}, {'A', Centring::a_face}, {'B', Centring::b_face},
                                         {'C', Centring::c_face},    {'I', Centring::body},   {'F', Centring::face}};
  std::size_t bSettings = 0;
  for (std::size_t hallNumber = 1; hallNumber <= 530; ++hallNumber)
  {
    const SKSpaceGroupSetting& setting = SKSpaceGroupDataBase::spaceGroupData[hallNumber];
    const std::string symbol = setting.HMString();
    if (symbol.empty() || !letters.contains(symbol.front())) continue;
    EXPECT_EQ(setting.centring(), letters.at(symbol.front())) << "Hall " << hallNumber << " " << symbol;
    bSettings += symbol.front() == 'B';
  }
  EXPECT_EQ(bSettings, 40uz);
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

TEST(SYMMETRY_STRUCTURE_FIXES, poscar_skips_comment_lines_strips_bang_comments_and_rejects_cells_without_volume)
{
  // comment-only lines between the sections (they used to read as the end of the file), '!' comments (only '#'
  // was recognised), and the title line itself starting with '#'
  const std::string content =
      "# title\n"
      "1.0  ! the scaling factor\n"
      "# the lattice\n"
      "4.0 0.0 0.0\n"
      "0.0 5.0 0.0 # b\n"
      "0.0 0.0 6.0 ! c\n"
      "! two species\n"
      "Si O\n"
      "1 1\n"
      "Direct\n"
      "# the atoms\n"
      "0.1 0.2 0.3 ! Si\n"
      "   # indented comment\n"
      "0.4 0.5 0.6 # O\n";
  SKPOSCARParser parser(content);
  parser.startParsing();
  EXPECT_NEAR(parser.movies().front().front()->cell->unitCell().by, 5.0, 1e-12);
  const std::vector<std::tuple<double3, std::size_t, double>> atoms = parser.firstTestFrame();
  ASSERT_EQ(atoms.size(), 2uz);
  EXPECT_EQ(std::get<1>(atoms[0]), 14uz);
  EXPECT_EQ(std::get<1>(atoms[1]), 8uz);
  EXPECT_NEAR(std::get<0>(atoms[1]).z, 0.6, 1e-12);

  // a zero scaling factor and linearly dependent lattice vectors gave infinite or NaN coordinates
  const std::string zeroScale = "cell\n0.0\n4.0 0.0 0.0\n0.0 5.0 0.0\n0.0 0.0 6.0\nSi\n1\nDirect\n0.1 0.2 0.3\n";
  EXPECT_THROW(SKPOSCARParser(zeroScale).startParsing(), std::runtime_error);
  const std::string flatCell = "cell\n1.0\n4.0 0.0 0.0\n8.0 0.0 0.0\n0.0 0.0 6.0\nSi\n1\nCartesian\n1.0 2.0 3.0\n";
  EXPECT_THROW(SKPOSCARParser(flatCell).startParsing(), std::runtime_error);
  const std::string flatCellByVolume =
      "cell\n-100.0\n4.0 0.0 0.0\n8.0 0.0 0.0\n0.0 0.0 6.0\nSi\n1\nDirect\n0.1 0.2 0.3\n";
  EXPECT_THROW(SKPOSCARParser(flatCellByVolume).startParsing(), std::runtime_error);
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

TEST(SYMMETRY_STRUCTURE_FIXES, seitz_matrices_without_translation_are_the_distinct_rotations)
{
  // P -1 (Hall 2), C 1 2/c 1 (Hall 90) and F d -3 m in both origin choices (Hall 525, 526): the rotation parts,
  // each once and without translation, as many as the order divided by the centring translations. The function
  // used to return the full operations, translations and centring copies included (and before that twice them).
  for (const std::size_t hallNumber : {2uz, 90uz, 525uz, 526uz})
  {
    const SKSpaceGroupSetting setting = SKSpaceGroup(hallNumber).spaceGroupSetting();
    const std::vector<SKSeitzIntegerMatrix> rotations = setting.SeitzMatricesWithoutTranslation();
    const auto full = setting.fullSeitzMatrices().operations;
    EXPECT_EQ(full.size(), setting.order()) << "Hall " << hallNumber;
    EXPECT_EQ(rotations.size(), setting.order() / setting.latticeTranslations().size()) << "Hall " << hallNumber;
    for (std::size_t i = 0; i < rotations.size(); ++i)
    {
      EXPECT_EQ(rotations[i].translation, int3(0, 0, 0)) << "Hall " << hallNumber;
      for (std::size_t j = i + 1; j < rotations.size(); ++j)
        EXPECT_FALSE(rotations[i].rotation == rotations[j].rotation) << "Hall " << hallNumber;
      EXPECT_TRUE(std::any_of(full.begin(), full.end(),
                              [&](const SKSeitzIntegerMatrix& m) { return m.rotation == rotations[i].rotation; }))
          << "Hall " << hallNumber;
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

TEST(SYMMETRY_STRUCTURE_FIXES, smallest_primitive_cell_of_left_handed_and_small_cells)
{
  // a 2 x 1 x 1 supercell of a simple cubic lattice: two atoms, the primitive cell has half the volume
  using Atom = std::tuple<double3, std::size_t, double>;
  const std::vector<Atom> atoms{Atom{double3(0.0, 0.0, 0.0), 1, 1.0}, Atom{double3(0.5, 0.0, 0.0), 1, 1.0}};

  // left-handed basis (negative determinant): the signed volume made every candidate "larger" than the minimum
  const double3x3 leftHanded(double3(20.0, 0.0, 0.0), double3(0.0, 10.0, 0.0), double3(0.0, 0.0, -10.0));
  const double3x3 primitiveOfLeftHanded =
      SKSymmetryCell::findSmallestPrimitiveCell(atoms, atoms, leftHanded, false, 1e-2);
  EXPECT_NEAR(std::fabs(primitiveOfLeftHanded.determinant()), 1000.0, 1e-6);

  // a cell of two cubic Angstrom: the primitive cell of one cubic Angstrom fell under the absolute threshold of one
  const double3x3 small(double3(2.0, 0.0, 0.0), double3(0.0, 1.0, 0.0), double3(0.0, 0.0, 1.0));
  const double3x3 primitiveOfSmall = SKSymmetryCell::findSmallestPrimitiveCell(atoms, atoms, small, false, 1e-2);
  EXPECT_NEAR(std::fabs(primitiveOfSmall.determinant()), 1.0, 1e-9);

  // the ordinary case is unchanged
  const double3x3 cell(double3(20.0, 0.0, 0.0), double3(0.0, 10.0, 0.0), double3(0.0, 0.0, 10.0));
  EXPECT_NEAR(std::fabs(SKSymmetryCell::findSmallestPrimitiveCell(atoms, atoms, cell, false, 1e-2).determinant()),
              1000.0, 1e-6);
}

TEST(SYMMETRY_STRUCTURE_FIXES, stored_inversion_centres_match_the_inversion_operation)
{
  // inversion through c is the operation (-1, 2c): every centrosymmetric setting must list an operation with the
  // inversion rotation whose translation is twice the stored centre, modulo lattice translations (units of 1/24)
  const SKRotationMatrix inversion(int3(-1, 0, 0), int3(0, -1, 0), int3(0, 0, -1));
  for (std::size_t hallNumber = 1; hallNumber <= 530; ++hallNumber)
  {
    const SKSpaceGroupSetting& setting = SKSpaceGroupDataBase::spaceGroupData[hallNumber];
    const int3 centre = setting.inversionCenter();
    const SKSeitzIntegerMatrix expected(inversion, int3(2 * centre.x, 2 * centre.y, 2 * centre.z));
    const auto operations = setting.fullSeitzMatrices().operations;
    const bool centrosymmetric = std::any_of(operations.begin(), operations.end(),
                                             [&](const SKSeitzIntegerMatrix& m) { return m.rotation == inversion; });
    if (!centrosymmetric)
    {
      EXPECT_EQ(centre, int3(0, 0, 0)) << "Hall " << hallNumber;
      continue;
    }
    // lattice translations: the stored list, all in 1/24 units
    bool found = false;
    for (const int3& translation : setting.latticeTranslations())
    {
      const SKSeitzIntegerMatrix shifted(inversion, expected.translation + translation);
      found = found || operations.contains(shifted);
    }
    EXPECT_TRUE(found) << "Hall " << hallNumber << " " << setting.HMString() << ": centre (" << centre.x << ", "
                       << centre.y << ", " << centre.z << ")";
  }
}

TEST(SYMMETRY_STRUCTURE_FIXES, database_flags_agree_with_the_operations_and_settings_are_consistent)
{
  const SKRotationMatrix inversion(int3(-1, 0, 0), int3(0, -1, 0), int3(0, 0, -1));
  std::map<std::size_t, std::string> schoenfliesOfGroup;
  for (std::size_t hallNumber = 1; hallNumber <= 530; ++hallNumber)
  {
    const SKSpaceGroupSetting& setting = SKSpaceGroupDataBase::spaceGroupData[hallNumber];
    const auto operations = setting.fullSeitzMatrices().operations;
    const std::vector<int3> lattice = setting.latticeTranslations();
    const auto isLatticeTranslation = [&](const SKSeitzIntegerMatrix& m)
    {
      return std::any_of(lattice.begin(), lattice.end(),
                         [&](const int3& t) { return SKSeitzIntegerMatrix(m.rotation, t) == m; });
    };

    // inversion at the origin: the operation (-1, 0) is listed (P -1 and P m m m said false, three
    // non-centrosymmetric groups said true)
    const bool inversionAtOrigin = std::any_of(operations.begin(), operations.end(), [&](const SKSeitzIntegerMatrix& m)
                                               { return m.rotation == inversion && isLatticeTranslation(m); });
    EXPECT_EQ(setting.inversionAtOrigin(), inversionAtOrigin) << "Hall " << hallNumber << " " << setting.HMString();

    // symmorphic: every operation's translation is a lattice translation (P -3 1 c and P -3 m 1 were swapped)
    const bool symmorphic = std::all_of(operations.begin(), operations.end(), isLatticeTranslation);
    EXPECT_EQ(setting.symmorphicity() == Symmorphicity::symmorphic, symmorphic)
        << "Hall " << hallNumber << " " << setting.HMString();

    // the settings of one group share the Schoenflies symbol, which is never empty
    EXPECT_FALSE(setting.schoenflies().empty()) << "Hall " << hallNumber;
    const auto [it, inserted] = schoenfliesOfGroup.emplace(setting.number(), setting.schoenflies());
    EXPECT_EQ(it->second, setting.schoenflies()) << "Hall " << hallNumber << " " << setting.HMString();

    // a change of basis has a non-zero determinant (seven rhombohedral-axes settings had a zero column)
    EXPECT_NE(setting.transformationMatrix().determinant(), 0) << "Hall " << hallNumber;

    // the asymmetric unit lies within the cell and its images: no bound beyond 1 + 1/4 (three said 2)
    const SKAsymmetricUnit unit = setting.asymmetricUnit();
    for (const std::pair<int, int>& range : {unit.a, unit.b, unit.c})
    {
      EXPECT_LE(range.second, 60) << "Hall " << hallNumber;
      EXPECT_GE(range.first, -48) << "Hall " << hallNumber;
    }
  }
  EXPECT_EQ(SKSpaceGroupDataBase::spaceGroupData[1].schoenflies(), "C1^1");
  EXPECT_EQ(SKSpaceGroupDataBase::spaceGroupData[2].schoenflies(), "Ci^1");
  EXPECT_EQ(SKSpaceGroupDataBase::spaceGroupData[3].schoenflies(), "C2^1");

  // the first Hall number listed for a group is its standard setting (group 228 listed origin choice 1 first)
  for (std::size_t number = 1; number <= 230; ++number)
  {
    const std::vector<std::size_t>& hallNumbers = SKSpaceGroupDataBase::spaceGroupHallData[number];
    ASSERT_FALSE(hallNumbers.empty()) << number;
    EXPECT_TRUE(SKSpaceGroupDataBase::spaceGroupData[hallNumbers.front()].standardSetting()) << "group " << number;
  }
}

TEST(SYMMETRY_STRUCTURE_FIXES, point_group_symbol_element_periods_and_residue_tables)
{
  // D3h is -6m2, not -6m
  EXPECT_EQ(SKPointGroup::pointGroupData[26].symbol(), "-6m2");

  // beryllium is in period 2 and selenium in period 4 (both said 3)
  EXPECT_EQ(PredefinedElements::predefinedElements[4]._chemicalSymbol, "Be");
  EXPECT_EQ(PredefinedElements::predefinedElements[4]._period, 2uz);
  EXPECT_EQ(PredefinedElements::predefinedElements[34]._chemicalSymbol, "Se");
  EXPECT_EQ(PredefinedElements::predefinedElements[34]._period, 4uz);

  // the oxygen of water is oxygen
  EXPECT_EQ(PredefinedElements::residueDefinitionsElement.at("HOH+O"), "O");

  // PDB atom names in the bonded-atom table carry no padding, so exact lookups find them
  for (const auto& [key, bonded] : PredefinedElements::residueDefinitionsBondedAtoms)
  {
    EXPECT_EQ(key, trim(key));
    for (const std::string& name : bonded) EXPECT_EQ(name, trim(name)) << key;
  }
  EXPECT_EQ(PredefinedElements::residueDefinitionsBondedAtoms.at("ALA+CA"), (std::vector<std::string>{"N", "CB", "C"}));
}

TEST(SYMMETRY_STRUCTURE_FIXES, coincident_groups_merge_transitively)
{
  const double3x3 cell(double3(10.0, 0.0, 0.0), double3(0.0, 10.0, 0.0), double3(0.0, 0.0, 10.0));
  const double tolerance = 1.0e-3;
  // two early copies 1.05 tolerances apart, and a third between them within 0.55 tolerances of both: the third
  // used to be linked to the first only, which left the first two in separate groups
  const std::vector<double3> positions{double3(5.0, 5.0, 5.0), double3(5.0 + 1.05 * tolerance, 5.0, 5.0),
                                       double3(5.0 + 0.525 * tolerance, 5.0, 5.0), double3(2.0, 2.0, 2.0)};
  const std::vector<std::size_t> group =
      coincidentGroups(cell, cell.inverse(), double3(10.0, 10.0, 10.0), positions, tolerance);
  ASSERT_EQ(group.size(), 4uz);
  EXPECT_EQ(group[0], group[1]);
  EXPECT_EQ(group[1], group[2]);
  EXPECT_NE(group[3], group[0]);
  EXPECT_EQ(group[3], 3uz);

  // coincidence across the periodic boundary still counts
  const std::vector<double3> wrapped{double3(0.0001, 5.0, 5.0), double3(9.9995, 5.0, 5.0)};
  const std::vector<std::size_t> wrappedGroup =
      coincidentGroups(cell, cell.inverse(), double3(10.0, 10.0, 10.0), wrapped, tolerance);
  EXPECT_EQ(wrappedGroup[0], wrappedGroup[1]);
}

TEST(SYMMETRY_STRUCTURE_FIXES, largest_empty_sphere_radius_and_position_describe_the_same_vertex)
{
  // in the diagram of the entire space every vertex can lie inside the spheres, with a negative radius; the radius
  // used to start from 0 and the position from -1, so they could describe no vertex or different ones
  SKApolloniusDiagram diagram;
  diagram.vertices.push_back(SKApolloniusVertex{double3(1.0, 1.0, 1.0), -2.0, {}, {}, {}});
  diagram.vertices.push_back(SKApolloniusVertex{double3(3.0, 3.0, 3.0), -0.5, {}, {}, {}});
  diagram.vertices.push_back(SKApolloniusVertex{double3(5.0, 5.0, 5.0), -1.5, {}, {}, {}});
  EXPECT_DOUBLE_EQ(diagram.largestEmptySphereRadius(), -0.5);
  EXPECT_EQ(diagram.largestEmptySpherePosition(), double3(3.0, 3.0, 3.0));

  // without vertices there is no sphere: a radius of zero at the origin
  const SKApolloniusDiagram empty;
  EXPECT_DOUBLE_EQ(empty.largestEmptySphereRadius(), 0.0);
  EXPECT_EQ(empty.largestEmptySpherePosition(), double3(0.0, 0.0, 0.0));
}

TEST(SYMMETRY_STRUCTURE_FIXES, cell_setters_constructors_and_replica_indices_stay_consistent)
{
  // the a/b/c and angle setters left the full cell (unit cell times replicas) at its old value
  SKCell cell(10.0, 10.0, 10.0, std::numbers::pi / 2.0, std::numbers::pi / 2.0, std::numbers::pi / 2.0);
  cell.setMaximumReplicas(int3(1, 2, 0));  // 2 x 3 x 1 replicas
  cell.a(20.0);
  EXPECT_NEAR(cell.box().ax, 40.0, 1e-12);
  cell.alpha(std::numbers::pi / 3.0);
  // the b vector is scaled by 3 replicas, the c vector (one replica along z) is the unit cell's
  EXPECT_NEAR(cell.box().by, 3.0 * cell.unitCell().by, 1e-12);
  EXPECT_NEAR(cell.box().cy, cell.unitCell().cy, 1e-12);
  EXPECT_NEAR(cell.box().cz, cell.unitCell().cz, 1e-12);
  EXPECT_NEAR((cell.box() * cell.box().inverse()).ax, 1.0, 1e-12);

  // the bounding-box and super-cell constructors left the bounding box at its default
  const SKBoundingBox box(double3(0.0, 0.0, 0.0), double3(5.0, 6.0, 7.0));
  SKCell fromBox(box);
  EXPECT_EQ(fromBox.boundingBox().maximum(), fromBox.enclosingBoundingBox().maximum());
  EXPECT_EQ(fromBox.boundingBox().minimum(), fromBox.enclosingBoundingBox().minimum());
  SKCell superCell(cell, int3(0, 0, 0), int3(1, 1, 1));
  EXPECT_EQ(superCell.boundingBox().maximum(), superCell.enclosingBoundingBox().maximum());
  EXPECT_EQ(superCell.boundingBox().minimum(), superCell.enclosingBoundingBox().minimum());

  // the replica of an index is the one at that position of the render translation list (x fastest)
  SKCell replicas(10.0, 10.0, 10.0, std::numbers::pi / 2.0, std::numbers::pi / 2.0, std::numbers::pi / 2.0);
  replicas.setMinimumReplicas(int3(-1, 0, 0));
  replicas.setMaximumReplicas(int3(0, 2, 3));  // 2 x 3 x 4
  const std::vector<float4> translations = replicas.renderTranslationVectors();
  ASSERT_EQ(translations.size(), 24uz);
  for (int index = 0; index < 24; ++index)
  {
    const int3 replica = replicas.replicaFromIndex(index);
    EXPECT_EQ(replica, int3(static_cast<int>(translations[static_cast<std::size_t>(index)].x),
                            static_cast<int>(translations[static_cast<std::size_t>(index)].y),
                            static_cast<int>(translations[static_cast<std::size_t>(index)].z)))
        << index;
  }
}

TEST(SYMMETRY_STRUCTURE_FIXES, delaunay_reduction_sorts_near_equal_lengths_safely)
{
  // a cubic cell: all seven candidate vectors have the same length up to rounding, the case in which the former
  // "a + 1e-10 < b" comparator was not a strict weak ordering; the reduced cell is a right-handed cell of the same
  // volume with edges of the same length
  const double3x3 cubic(double3(10.0, 0.0, 0.0), double3(0.0, 10.0, 0.0), double3(0.0, 0.0, 10.0));
  const std::optional<double3x3> reducedCubic = SKSymmetryCell::computeDelaunayReducedCell(cubic, 1e-2);
  ASSERT_TRUE(reducedCubic.has_value());
  EXPECT_NEAR(reducedCubic->determinant(), 1000.0, 1e-9);
  for (std::size_t i = 0; i < 3; ++i) EXPECT_NEAR((*reducedCubic)[i].length(), 10.0, 1e-9) << i;

  // a skewed cell reduces to one of the same volume whose edges are no longer than the originals
  const double3x3 skewed(double3(10.0, 0.0, 0.0), double3(9.0, 1.0, 0.0), double3(0.0, 0.0, 10.0));
  const std::optional<double3x3> reducedSkewed = SKSymmetryCell::computeDelaunayReducedCell(skewed, 1e-2);
  ASSERT_TRUE(reducedSkewed.has_value());
  EXPECT_NEAR(std::fabs(reducedSkewed->determinant()), 100.0, 1e-9);
  for (std::size_t i = 0; i < 3; ++i) EXPECT_LE((*reducedSkewed)[i].length(), 10.0 + 1e-9) << i;

  // the two-dimensional reduction (the b axis kept) of the same cells
  const std::optional<double3x3> reduced2D = SKSymmetryCell::computeDelaunayReducedCell2D(cubic, 1e-2);
  ASSERT_TRUE(reduced2D.has_value());
  EXPECT_NEAR(std::fabs(reduced2D->determinant()), 1000.0, 1e-9);
}

TEST(SYMMETRY_STRUCTURE_FIXES, orthogonal_axis_directions_of_low_rotation_orders_terminate)
{
  // an order of zero wrapped the loop bound around to the largest std::size_t; orders zero and one sum no further
  // powers of the rotation, so the identity alone annihilates no axis
  EXPECT_TRUE(SKRotationMatrix::r_2_100.orthogonalToAxisDirection(0).empty());
  EXPECT_TRUE(SKRotationMatrix::r_2_100.orthogonalToAxisDirection(1).empty());

  // for the two-fold rotation about x the axes annihilated by 1 + R are the directions perpendicular to x
  const std::vector<int3> perpendicular = SKRotationMatrix::r_2_100.orthogonalToAxisDirection(2);
  ASSERT_FALSE(perpendicular.empty());
  for (const int3& axis : perpendicular) EXPECT_EQ(axis.x, 0) << axis.x << " " << axis.y << " " << axis.z;
}

TEST(SYMMETRY_STRUCTURE_FIXES, default_constructed_symmetry_objects_have_defined_members)
{
  // these members used to be left indeterminate by the default constructors
  const SKSymmetryCell cell;
  EXPECT_EQ(cell.a(), 0.0);
  EXPECT_EQ(cell.b(), 0.0);
  EXPECT_EQ(cell.c(), 0.0);
  EXPECT_DOUBLE_EQ(cell.alpha(), std::numbers::pi / 2.0);
  EXPECT_DOUBLE_EQ(cell.beta(), std::numbers::pi / 2.0);
  EXPECT_DOUBLE_EQ(cell.gamma(), std::numbers::pi / 2.0);

  SKAsymmetricAtom atom;
  EXPECT_EQ(atom.asymmetricIndex(), 0uz);
  SKAsymmetricAtom carbon("C", 6);
  EXPECT_EQ(carbon.asymmetricIndex(), 0uz);
  carbon.setAsymmetricIndex(7);
  SKAsymmetricAtom copy(carbon);
  EXPECT_EQ(copy.asymmetricIndex(), 7uz);

  SKAtomCopy atomCopy;
  EXPECT_EQ(atomCopy.asymmetricIndex(), -1);
}

TEST(SYMMETRY_STRUCTURE_FIXES, space_group_symbols_match_with_their_spaces_removed)
{
  // the comparison against the stored symbol without its spaces was lost, so a symbol written like "P121/c1" or
  // "Fm-3m" found nothing; several settings share a Hermann-Mauguin symbol, so the two spellings are compared with
  // each other rather than with a fixed Hall number
  for (const std::size_t hallNumber : {9uz, 81uz, 526uz})
  {
    const SKSpaceGroupSetting& setting = SKSpaceGroupDataBase::spaceGroupData[hallNumber];
    for (const auto& [symbol, lookup] : {std::pair{setting.HMString(), &SKSpaceGroup::HallNumberFromHMString},
                                         std::pair{setting.HallString(), &SKSpaceGroup::HallNumber}})
    {
      std::string stripped = symbol;
      std::erase(stripped, ' ');
      ASSERT_NE(stripped, symbol) << symbol;
      const std::optional<std::size_t> found = lookup(symbol);
      ASSERT_TRUE(found.has_value()) << symbol;
      EXPECT_EQ(lookup(stripped), found) << stripped;
    }
  }
  EXPECT_FALSE(SKSpaceGroup::HallNumberFromHMString("not a symbol").has_value());
}

TEST(SYMMETRY_STRUCTURE_FIXES, cif_symmetry_operations_accept_decimal_translations)
{
  // only integers and fractions were parsed; decimal translations made the whole operation unreadable
  const std::optional<SKSeitzMatrix> decimal = SKSpaceGroup::parseCifSymmetryOperation("0.5+x, 0.25+y, -z+0.75");
  ASSERT_TRUE(decimal.has_value());
  EXPECT_DOUBLE_EQ(decimal->translation.x, 0.5);
  EXPECT_DOUBLE_EQ(decimal->translation.y, 0.25);
  EXPECT_DOUBLE_EQ(decimal->translation.z, 0.75);
  EXPECT_EQ(decimal->rotation.int3x3_m.m33, -1);

  const std::optional<SKSeitzMatrix> third = SKSpaceGroup::parseCifSymmetryOperation("0.33333+x, .5+y, z");
  ASSERT_TRUE(third.has_value());
  EXPECT_NEAR(third->translation.x, 1.0 / 3.0, 1e-5);
  EXPECT_DOUBLE_EQ(third->translation.y, 0.5);

  // fractions and integers still parse, and a lone '.' does not
  const std::optional<SKSeitzMatrix> fraction = SKSpaceGroup::parseCifSymmetryOperation("x+1/2, y-1, z");
  ASSERT_TRUE(fraction.has_value());
  EXPECT_DOUBLE_EQ(fraction->translation.x, 0.5);
  EXPECT_DOUBLE_EQ(fraction->translation.y, -1.0);
  EXPECT_FALSE(SKSpaceGroup::parseCifSymmetryOperation("x+., y, z").has_value());
}

TEST(SYMMETRY_STRUCTURE_FIXES, operation_sets_return_themselves_and_asymmetric_atoms_lie_in_the_unit)
{
  // the full operations of a set are the set (the integer version returned an empty set, the other was undefined)
  const SKSpaceGroupSetting setting =
      SKSpaceGroup(SKSpaceGroupDataBase::spaceGroupHallData[221].front()).spaceGroupSetting();  // P m -3 m
  SKIntegerSymmetryOperationSet integerSet = setting.fullSeitzMatrices();
  EXPECT_EQ(integerSet.fullSeitzMatrices().operations.size(), setting.order());
  SKSymmetryOperationSet set({SKSeitzMatrix(SKRotationMatrix::identity, double3(0.0, 0.0, 0.0))});
  EXPECT_EQ(set.fullSeitzMatrices().operations.size(), 1uz);

  // the orbit of (0.9, 0.9, 0.9) in P m -3 m is one asymmetric atom; its representative is moved into the
  // asymmetric unit of the setting, which the disabled containment check never did
  const double3x3 cubic(double3(10.0, 0.0, 0.0), double3(0.0, 10.0, 0.0), double3(0.0, 0.0, 10.0));
  std::vector<std::tuple<double3, std::size_t, double>> atoms;
  for (const SKSeitzIntegerMatrix& m : integerSet.operations)
  {
    const double3 position = double3::fract(m * double3(0.9, 0.9, 0.9));
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
  ASSERT_EQ(atoms.size(), 8uz);
  const std::vector<std::tuple<double3, std::size_t, double>> asymmetric =
      integerSet.asymmetricAtoms(setting.asymmetricUnit(), atoms, cubic, false, 1e-2);
  ASSERT_EQ(asymmetric.size(), 1uz);
  SKAsymmetricUnit unit = setting.asymmetricUnit();
  EXPECT_TRUE(unit.contains(std::get<0>(asymmetric.front())));
  EXPECT_NEAR(std::get<0>(asymmetric.front()).x, 0.1, 1e-12);
  EXPECT_NEAR(std::get<0>(asymmetric.front()).y, 0.1, 1e-12);
  EXPECT_NEAR(std::get<0>(asymmetric.front()).z, 0.1, 1e-12);
}

TEST(SYMMETRY_STRUCTURE_FIXES, hermann_mauguin_lookup_honours_setting_qualifiers)
{
  // F d -3 m: Hall 525 is origin choice 1, 526 origin choice 2 (standard); R -3 m: Hall 458 hexagonal axes
  // (standard), 459 rhombohedral axes. A qualifier after a colon used to make the symbol unreadable, and without
  // one the first setting in the table was returned, so the second origin and the rhombohedral axes were unreachable
  EXPECT_EQ(SKSpaceGroup::HallNumberFromHMString("F d -3 m :2"), 526uz);
  EXPECT_EQ(SKSpaceGroup::HallNumberFromHMString("F d -3 m :1"), 525uz);
  EXPECT_EQ(SKSpaceGroup::HallNumberFromHMString("Fd-3m:1"), 525uz);
  EXPECT_EQ(SKSpaceGroup::HallNumberFromHMString("'F d -3 m :1'"), 525uz);
  EXPECT_EQ(SKSpaceGroup::HallNumberFromHMString("F d -3 m"), 526uz);
  EXPECT_EQ(SKSpaceGroup::HallNumberFromHMString("R -3 m :R"), 459uz);
  EXPECT_EQ(SKSpaceGroup::HallNumberFromHMString("R -3 m :H"), 458uz);
  EXPECT_EQ(SKSpaceGroup::HallNumberFromHMString("R -3 m"), 458uz);
  EXPECT_FALSE(SKSpaceGroup::HallNumberFromHMString("F d -3 m :3").has_value());
  // an axis-setting qualifier selects that setting: P b n m is P n m a in the "cab" setting (Hall 294)
  EXPECT_EQ(SKSpaceGroup::HallNumberFromHMString("P b n m :cab"), 294uz);
  EXPECT_EQ(SKSpaceGroup::HallNumberFromHMString("P b n m"), 294uz);
}

TEST(SYMMETRY_STRUCTURE_FIXES, centring_from_symbol_follows_the_database_for_rhombohedral_settings)
{
  // on rhombohedral axes the lattice is primitive in the database, and an H-prefixed symbol is the hexagonal-axes
  // setting stored as an R lattice; "R ..." used to give R always and "H ..." an H centring no setting has
  EXPECT_EQ(SKSpaceGroup::centringFromHMString("R -3 m"), Centring::r);
  EXPECT_EQ(SKSpaceGroup::centringFromHMString("R -3 m :H"), Centring::r);
  EXPECT_EQ(SKSpaceGroup::centringFromHMString("'R -3 m :R'"), Centring::primitive);
  EXPECT_EQ(SKSpaceGroup::centringFromHMString("H -3 m"), Centring::r);
  EXPECT_EQ(SKSpaceGroup::centringFromHMString("F d -3 m :2"), Centring::face);
  for (const std::size_t hallNumber : {434uz, 437uz, 445uz, 451uz, 453uz, 459uz, 461uz})
  {
    const SKSpaceGroupSetting& setting = SKSpaceGroupDataBase::spaceGroupData[hallNumber];
    EXPECT_EQ(SKSpaceGroup::centringFromHMString(setting.HMString() + " :" + setting.qualifier()), setting.centring())
        << hallNumber;
  }

  // so that the operations of R -3 m on rhombohedral axes are identified with the centring the symbol gives
  const double3x3 rhombohedral(double3(10.0, 0.0, 0.0), double3(-5.0, 8.66, 0.0), double3(0.0, 0.0, 10.0));
  const std::vector<SKSeitzMatrix> operations = operationsOfSetting(459);
  EXPECT_EQ(SKSpaceGroup::HallNumberFromSymmetryOperations(operations, rhombohedral,
                                                           SKSpaceGroup::centringFromHMString("R -3 m :R"), 166),
            459uz);
}

TEST(SYMMETRY_STRUCTURE_FIXES, hexagonal_centring_transformations_keep_the_c_axis)
{
  // the H cell is three primitive cells side by side in the ab-plane: determinant 3 (it was 9) and the inverse
  // leaves c alone (it divided it by 3)
  EXPECT_EQ(SKTransformationMatrix::primitiveToHexagonal.determinant(), 3);
  const double3x3 product = SKTransformationMatrix::hexagonalToPrimitive *
                            double3x3(SKTransformationMatrix::primitiveToHexagonal.transformation);
  for (std::size_t i = 0; i < 3; ++i)
    for (std::size_t j = 0; j < 3; ++j) EXPECT_NEAR(product[i][j], i == j ? 1.0 : 0.0, 1e-12) << i << " " << j;
  EXPECT_NEAR(SKTransformationMatrix::hexagonalToPrimitive.determinant(), 1.0 / 3.0, 1e-12);
}

TEST(SYMMETRY_STRUCTURE_FIXES, reverse_hexagonal_cell_is_three_rhombohedral_cells)
{
  // the matrix had determinant 1 (its columns were not a hexagonal cell on rhombohedral axes at all)
  EXPECT_EQ(SKTransformationMatrix::rhombohedralHexagonalToReverse.determinant(), 3);
  EXPECT_NEAR(SKTransformationMatrix::rhombohedralReverseToPrimitive.determinant(), 1.0 / 3.0, 1e-12);
  const double3x3 product = SKTransformationMatrix::rhombohedralReverseToPrimitive *
                            double3x3(SKTransformationMatrix::rhombohedralHexagonalToReverse.transformation);
  for (std::size_t i = 0; i < 3; ++i)
    for (std::size_t j = 0; j < 3; ++j) EXPECT_NEAR(product[i][j], i == j ? 1.0 : 0.0, 1e-12) << i << " " << j;

  // a rhombohedral cell (a = 5, alpha = 70 degrees) taken to the hexagonal cell: a_H = b_H at 120 degrees, c_H
  // perpendicular to both, three times the volume
  const double alpha = 70.0 * std::numbers::pi / 180.0;
  const double a = 5.0;
  const double3 aR(a, 0.0, 0.0);
  const double3 bR(a * std::cos(alpha), a * std::sin(alpha), 0.0);
  const double cx = a * std::cos(alpha);
  const double cy = a * (std::cos(alpha) - std::cos(alpha) * std::cos(alpha)) / std::sin(alpha);
  const double3 cR(cx, cy, std::sqrt(a * a - cx * cx - cy * cy));
  const double3x3 rhombohedral(aR, bR, cR);
  const double3x3 hexagonal = rhombohedral * SKTransformationMatrix::rhombohedralHexagonalToReverse;
  EXPECT_NEAR(hexagonal[0].length(), hexagonal[1].length(), 1e-12);
  EXPECT_NEAR(double3::dot(hexagonal[0], hexagonal[1]), -0.5 * hexagonal[0].length_squared(), 1e-10);
  EXPECT_NEAR(double3::dot(hexagonal[0], hexagonal[2]), 0.0, 1e-10);
  EXPECT_NEAR(double3::dot(hexagonal[1], hexagonal[2]), 0.0, 1e-10);
  EXPECT_NEAR(hexagonal.determinant(), 3.0 * rhombohedral.determinant(), 1e-10);

  // and it is the reverse setting: a_R sits at (1/3, 2/3, 1/3) of the hexagonal cell, not the obverse (2/3, 1/3, 1/3)
  const double3 aRInHexagonal = hexagonal.inverse() * aR;
  EXPECT_NEAR(aRInHexagonal.x, 1.0 / 3.0, 1e-12);
  EXPECT_NEAR(aRInHexagonal.y, 2.0 / 3.0, 1e-12);
  EXPECT_NEAR(aRInHexagonal.z, 1.0 / 3.0, 1e-12);
}

TEST(SYMMETRY_STRUCTURE_FIXES, niggli_cell_search_identifies_centred_groups_and_a_short_principal_axis)
{
  // F d -3 m from its 8a site (origin choice 2): a centred group, which the former direct match of the primitive
  // Niggli cell against primitive settings could not reach
  const double3x3 cubic(double3(10.0, 0.0, 0.0), double3(0.0, 10.0, 0.0), double3(0.0, 0.0, 10.0));
  std::vector<std::tuple<double3, std::size_t, double>> diamond;
  for (const SKSeitzIntegerMatrix& m : SKSpaceGroup(526).spaceGroupSetting().fullSeitzMatrices().operations)
  {
    const double3 position = double3::fract(m * double3(0.125, 0.125, 0.125));
    if (std::none_of(diamond.begin(), diamond.end(),
                     [&](const auto& atom)
                     {
                       const double3 d = std::get<0>(atom) - position;
                       return (d - double3(std::rint(d.x), std::rint(d.y), std::rint(d.z))).length_squared() < 1e-12;
                     }))
    {
      diamond.emplace_back(position, 1, 1.0);
    }
  }
  ASSERT_EQ(diamond.size(), 8uz);
  const std::optional<SKSpaceGroup::FoundNiggliCellInfo> fcc =
      SKSpaceGroup::findNiggliCell(cubic, diamond, false, 1e-2);
  ASSERT_TRUE(fcc.has_value());
  EXPECT_EQ(SKSpaceGroup(fcc->HallNumber).spaceGroupSetting().number(), 227uz);
  EXPECT_EQ(fcc->atoms.size(), 1uz);

  // one atom in a tetragonal cell whose c axis is the shortest: P 4/m m m, with the Niggli order putting the
  // four-fold axis first rather than third
  const double3x3 tetragonal(double3(10.0, 0.0, 0.0), double3(0.0, 10.0, 0.0), double3(0.0, 0.0, 6.0));
  const std::optional<SKSpaceGroup::FoundNiggliCellInfo> primitive = SKSpaceGroup::findNiggliCell(
      tetragonal, {std::tuple<double3, std::size_t, double>{double3(0.0, 0.0, 0.0), 1, 1.0}}, false, 1e-2);
  ASSERT_TRUE(primitive.has_value());
  EXPECT_EQ(SKSpaceGroup(primitive->HallNumber).spaceGroupSetting().number(), 123uz);
  EXPECT_NEAR(primitive->cell.c(), 6.0, 1e-9);
}

TEST(SYMMETRY_STRUCTURE_FIXES, symmetry_operations_in_a_reverse_or_permuted_setting_are_identified)
{
  // R -3 m on hexagonal axes (Hall 458) with its operations rewritten in the reverse setting, a rotation by 180
  // degrees about c; the matcher tried the identity basis only for trigonal groups and found nothing
  const std::vector<SKSeitzMatrix> obverse = operationsOfSetting(458);
  const std::vector<SKSeitzMatrix> reverse =
      SKSymmetryOperationSet(obverse)
          .changedBasis(SKTransformationMatrix(int3(-1, 0, 0), int3(0, -1, 0), int3(0, 0, 1)))
          .operations;
  const double3x3 hexagonal(double3(10.0, 0.0, 0.0), double3(-5.0, 8.6602540378, 0.0), double3(0.0, 0.0, 12.0));
  EXPECT_EQ(SKSpaceGroup::HallNumberFromSymmetryOperations(obverse, hexagonal, Centring::r, 166), 458uz);
  EXPECT_EQ(SKSpaceGroup::HallNumberFromSymmetryOperations(reverse, hexagonal, Centring::r, 166), 458uz);

  // P 4/m m m with its four-fold axis along a instead of c
  const std::size_t p4mmm = SKSpaceGroupDataBase::spaceGroupHallData[123].front();
  const std::vector<SKSeitzMatrix> permuted =
      SKSymmetryOperationSet(operationsOfSetting(p4mmm))
          .changedBasis(SKTransformationMatrix(int3(0, 1, 0), int3(0, 0, 1), int3(1, 0, 0)))
          .operations;
  const double3x3 tetragonal(double3(6.0, 0.0, 0.0), double3(0.0, 10.0, 0.0), double3(0.0, 0.0, 10.0));
  EXPECT_EQ(SKSpaceGroup::HallNumberFromSymmetryOperations(permuted, tetragonal, Centring::primitive, 123), p4mmm);
  EXPECT_EQ(SKSpaceGroup::HallNumberFromSymmetryOperations(permuted, tetragonal, Centring::primitive, std::nullopt),
            p4mmm);
}

TEST(SYMMETRY_STRUCTURE_FIXES, partial_occupancies_keep_every_element_on_a_shared_site)
{
  using Atom = std::tuple<double3, std::size_t, double>;
  const double3x3 cubic(double3(10.0, 0.0, 0.0), double3(0.0, 10.0, 0.0), double3(0.0, 0.0, 10.0));

  // two elements sharing one site (substitutional disorder): the primitive cell keeps both, with or without partial
  // occupancies (the flag used to merge them into the first element, dropping the second)
  const std::vector<Atom> shared{Atom{double3(0.25, 0.25, 0.25), 1, 0.5}, Atom{double3(0.25, 0.25, 0.25), 2, 0.5}};
  EXPECT_EQ(SKSymmetryCell::trim(shared, cubic, cubic, true, 1e-2).size(), 2uz);
  EXPECT_EQ(SKSymmetryCell::trim(shared, cubic, cubic, false, 1e-2).size(), 2uz);
  // a lattice translation copy of the same type is still trimmed away
  const std::vector<Atom> copies{Atom{double3(0.25, 0.25, 0.25), 1, 1.0}, Atom{double3(0.75, 0.25, 0.25), 1, 1.0}};
  const double3x3 halved(double3(5.0, 0.0, 0.0), double3(0.0, 10.0, 0.0), double3(0.0, 0.0, 10.0));
  EXPECT_EQ(SKSymmetryCell::trim(copies, cubic, halved, true, 1e-2).size(), 1uz);

  // two elements related by the inversion of P m -3 m: one asymmetric atom per element, so that expanding the
  // asymmetric unit puts both elements back on every site of the orbit
  SKIntegerSymmetryOperationSet operations =
      SKSpaceGroup(SKSpaceGroupDataBase::spaceGroupHallData[221].front()).spaceGroupSetting().fullSeitzMatrices();
  const SKAsymmetricUnit unit =
      SKSpaceGroup(SKSpaceGroupDataBase::spaceGroupHallData[221].front()).spaceGroupSetting().asymmetricUnit();
  std::vector<Atom> related{Atom{double3(0.1, 0.1, 0.1), 1, 0.5}, Atom{double3(0.9, 0.9, 0.9), 2, 0.5}};
  for (bool allowPartialOccupancies : {true, false})
  {
    const auto asymmetric = operations.asymmetricAtoms(unit, related, cubic, allowPartialOccupancies, 1e-2);
    ASSERT_EQ(asymmetric.size(), 2uz) << allowPartialOccupancies;
    EXPECT_EQ(std::get<1>(asymmetric[0]), 1uz);
    EXPECT_EQ(std::get<1>(asymmetric[1]), 2uz);
  }
  // the same element on both sites is one orbit
  std::vector<Atom> sameType{Atom{double3(0.1, 0.1, 0.1), 1, 1.0}, Atom{double3(0.9, 0.9, 0.9), 1, 1.0}};
  EXPECT_EQ(operations.asymmetricAtoms(unit, sameType, cubic, true, 1e-2).size(), 1uz);

  // the symmetry search itself still merges the elements of a shared site when partial occupancies are allowed
  EXPECT_TRUE(SKSymmetryCell::testTranslationalSymmetry(
      double3(0.5, 0.0, 0.0), {Atom{double3(0.25, 0.25, 0.25), 1, 0.5}, Atom{double3(0.75, 0.25, 0.25), 2, 0.5}}, cubic,
      true, 1e-2));
}

TEST(SYMMETRY_STRUCTURE_FIXES, asymmetric_unit_images_keep_their_shift_and_tolerate_round_off)
{
  // {12, 60} is [1/4, 5/4]: the image of 0.1 inside it is 1.1, not 0.1 wrapped into the cell
  const SKAsymmetricUnit shifted({12, 60}, {0, 24}, {0, 24});
  const std::optional<double3> image = shifted.imageInside(double3(0.1, 0.1, 0.1));
  ASSERT_TRUE(image.has_value());
  EXPECT_NEAR(image->x, 1.1, 1e-12);
  EXPECT_NEAR(image->y, 0.1, 1e-12);

  // round-off below a closed bound at 0 used to expel the point from [0, 1/2]
  const SKAsymmetricUnit half({0, 24}, {0, 24}, {0, 24});
  EXPECT_TRUE(half.contains(double3(-1e-17, 0.1, 0.1)));
  EXPECT_TRUE(half.contains(double3(0.5 + 1e-12, 0.1, 0.1)));
  EXPECT_FALSE(half.contains(double3(0.51, 0.1, 0.1)));
  // an open bound stays exact: [1/2, 1) excludes 1.0 (whose image 0.0 is below 1/2) but holds 1.0 - 1e-9
  const SKAsymmetricUnit open({24, 49}, {0, 24}, {0, 24});
  EXPECT_FALSE(open.contains(double3(1.0, 0.1, 0.1)));
  EXPECT_TRUE(open.contains(double3(1.0 - 1e-9, 0.1, 0.1)));

  // A 1 2 1 (Hall 10) has an asymmetric unit with z in [-1/2, 0]: the representative keeps a z in that range
  const SKSpaceGroupSetting setting = SKSpaceGroup(10).spaceGroupSetting();
  ASSERT_EQ(setting.asymmetricUnit().c, (std::pair<int, int>{-24, 0}));
  const double3x3 cell(double3(10.0, 0.0, 0.0), double3(0.0, 11.0, 0.0), double3(0.0, 0.0, 12.0));
  std::vector<std::tuple<double3, std::size_t, double>> atoms;
  for (const SKSeitzIntegerMatrix& m : setting.fullSeitzMatrices().operations)
  {
    const double3 position = double3::fract(m * double3(0.3, 0.2, 0.1));
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
  SKIntegerSymmetryOperationSet operations = setting.fullSeitzMatrices();
  const auto asymmetric = operations.asymmetricAtoms(setting.asymmetricUnit(), atoms, cell, false, 1e-2);
  ASSERT_EQ(asymmetric.size(), 1uz);
  const double3 representative = std::get<0>(asymmetric.front());
  EXPECT_TRUE(setting.asymmetricUnit().contains(representative));
  EXPECT_LE(representative.z, 0.0 + 1e-8);
  EXPECT_GE(representative.z, -0.5 - 1e-8);
}

TEST(SYMMETRY_STRUCTURE_FIXES, asymmetric_unit_keeps_a_point_that_is_already_inside)
{
  // a unit one cell wide holds both images of a point on its boundary; the image below used to win over the point
  const SKAsymmetricUnit symmetric({-24, 24}, {0, 24}, {0, 24});
  std::optional<double3> image = symmetric.imageInside(double3(0.5, 0.1, 0.1));
  ASSERT_TRUE(image.has_value());
  EXPECT_NEAR(image->x, 0.5, 1e-12);
  image = symmetric.imageInside(double3(-0.5, 0.1, 0.1));
  ASSERT_TRUE(image.has_value());
  EXPECT_NEAR(image->x, -0.5, 1e-12);

  const SKAsymmetricUnit wholeCell({-48, 0}, {0, 24}, {0, 24});
  image = wholeCell.imageInside(double3(0.0, 0.1, 0.1));
  ASSERT_TRUE(image.has_value());
  EXPECT_NEAR(image->x, 0.0, 1e-12);

  // [-1/2, 1) (Hall 14, 17, ...) holds 0.75 itself; it used to come back as its centred image -0.25
  const SKAsymmetricUnit wide({-24, 49}, {0, 24}, {0, 24});
  image = wide.imageInside(double3(0.75, 0.1, 0.1));
  ASSERT_TRUE(image.has_value());
  EXPECT_NEAR(image->x, 0.75, 1e-12);
  image = wide.imageInside(double3(0.3, 0.1, 0.1));
  ASSERT_TRUE(image.has_value());
  EXPECT_NEAR(image->x, 0.3, 1e-12);
  // [-1, 1) (Hall 22, 25, ...): points in [0, 1) stay where they are
  const SKAsymmetricUnit twoCells({-48, 49}, {0, 24}, {0, 24});
  image = twoCells.imageInside(double3(0.9, 0.1, 0.1));
  ASSERT_TRUE(image.has_value());
  EXPECT_NEAR(image->x, 0.9, 1e-12);

  // a point outside the unit still comes back as the image that is inside
  image = wide.imageInside(double3(1.25, 0.1, 0.1));
  ASSERT_TRUE(image.has_value());
  EXPECT_NEAR(image->x, 0.25, 1e-12);
  image = wide.imageInside(double3(-0.75, 0.1, 0.1));
  ASSERT_TRUE(image.has_value());
  EXPECT_NEAR(image->x, 0.25, 1e-12);

  // the representative of an atom already in a one-cell-wide unit is not moved by asymmetricAtoms (I 1 2 1, Hall 14)
  const SKSpaceGroupSetting setting = SKSpaceGroup(14).spaceGroupSetting();
  ASSERT_EQ(setting.asymmetricUnit().a, (std::pair<int, int>{-24, 49}));
  const double3x3 cell(double3(10.0, 0.0, 0.0), double3(0.0, 11.0, 0.0), double3(0.0, 0.0, 12.0));
  const double3 start(0.75, 0.1, 0.2);
  ASSERT_TRUE(setting.asymmetricUnit().contains(start));
  std::vector<std::tuple<double3, std::size_t, double>> atoms{{start, 1, 1.0}};
  SKIntegerSymmetryOperationSet operations = setting.fullSeitzMatrices();
  const auto asymmetric = operations.asymmetricAtoms(setting.asymmetricUnit(), atoms, cell, false, 1e-2);
  ASSERT_EQ(asymmetric.size(), 1uz);
  EXPECT_NEAR(std::get<0>(asymmetric.front()).x, 0.75, 1e-12);
}

TEST(SYMMETRY_STRUCTURE_FIXES, exact_symbol_matches_take_precedence_over_space_stripped_ones)
{
  // " P 3 2" (Hall 438, P 3 1 2) and " P 32" (Hall 432, P 3_2) strip to the same "p32"; the exact match used to
  // lose to the stripped match of the earlier entry. Likewise " P 6 2" (Hall 471) and " P 62" (Hall 465).
  EXPECT_EQ(SKSpaceGroup::HallNumber(" P 3 2"), 438uz);
  EXPECT_EQ(SKSpaceGroup::HallNumber(" P 32"), 432uz);
  EXPECT_EQ(SKSpaceGroup::HallNumber(" P 6 2"), 471uz);
  EXPECT_EQ(SKSpaceGroup::HallNumber(" P 62"), 465uz);
  // Hermann-Mauguin: "P 3 1 2" (group 149) and "P 31 1 2" (group 151) strip to the same "p3112"
  const std::optional<std::size_t> p312 = SKSpaceGroup::HallNumberFromHMString("P 3 1 2");
  const std::optional<std::size_t> p3112 = SKSpaceGroup::HallNumberFromHMString("P 31 1 2");
  ASSERT_TRUE(p312.has_value() && p3112.has_value());
  EXPECT_EQ(SKSpaceGroup(*p312).spaceGroupSetting().number(), 149uz);
  EXPECT_EQ(SKSpaceGroup(*p3112).spaceGroupSetting().number(), 151uz);
}

TEST(SYMMETRY_STRUCTURE_FIXES, shared_hall_symbols_of_group_68_are_reached_through_the_setting_qualifier)
{
  // the two origin-choice-1 settings of each pair have the same Hall symbol and the same operations (the Hall
  // notation does not distinguish them); the bare symbol gives the first, the axis qualifier selects the second
  for (const auto& [first, second, symbol, qualifier] :
       {std::tuple{322uz, 324uz, " C 2 2 -1ac", "ba-c"}, std::tuple{326uz, 328uz, " A 2 2 -1ab", "-cba"},
        std::tuple{330uz, 332uz, " B 2 2 -1ab", "a-cb"}})
  {
    const SKSpaceGroupSetting& one = SKSpaceGroupDataBase::spaceGroupData[first];
    const SKSpaceGroupSetting& other = SKSpaceGroupDataBase::spaceGroupData[second];
    EXPECT_EQ(one.HallString(), other.HallString());
    EXPECT_EQ(one.encodedSeitz(), other.encodedSeitz());
    EXPECT_EQ(SKSpaceGroup::HallNumber(symbol), first);
    EXPECT_EQ(SKSpaceGroup::HallNumber(std::string(symbol) + " :" + qualifier), second);
    EXPECT_EQ(SKSpaceGroup::HallNumber(std::string(symbol) + " :" + one.qualifier()), first);
  }
  // an origin choice as qualifier works for Hall symbols too
  EXPECT_EQ(SKSpaceGroup::HallNumber("-F 4vw 2vw 3 :2"), 526uz);
  EXPECT_FALSE(SKSpaceGroup::HallNumber("-F 4vw 2vw 3 :1").has_value());
}

TEST(SYMMETRY_STRUCTURE_FIXES, rotations_matching_no_point_group_give_no_point_group)
{
  // the identity with a four-fold rotation but without its square is no group: nothing, not point group number 0
  const SKPointSymmetrySet broken({SKRotationMatrix::identity, SKRotationMatrix::r_4_001});
  EXPECT_FALSE(SKPointGroup::fromPointSymmetry(broken).has_value());
  EXPECT_EQ(SKPointGroup(broken).number(), 0uz);

  // the identity and the inversion form point group -1 (number 2)
  const SKPointSymmetrySet ci({SKRotationMatrix::identity, SKRotationMatrix::inversionIdentity});
  const std::optional<SKPointGroup> pointGroup = SKPointGroup::fromPointSymmetry(ci);
  ASSERT_TRUE(pointGroup.has_value());
  EXPECT_EQ(pointGroup->number(), 2uz);
}

TEST(SYMMETRY_STRUCTURE_FIXES, hall_numbers_beyond_the_table_are_rejected)
{
  // the table has 531 entries; index 531 used to be read past its end
  EXPECT_NO_THROW(SKSpaceGroup(530));
  EXPECT_THROW(SKSpaceGroup(531), std::out_of_range);
  EXPECT_THROW(SKSpaceGroup::latticeTranslationStrings(531), std::out_of_range);
  EXPECT_THROW(SKSpaceGroup::inversionCenterString(10000), std::out_of_range);
  EXPECT_EQ(SKSpaceGroup::latticeTranslationStrings(530).size(), 4uz);
}
