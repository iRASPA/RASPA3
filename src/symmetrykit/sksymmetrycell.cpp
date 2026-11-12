module;

module sksymmetrycell;

import std;

import int3;
import int3x3;
import double3;
import double3x3;
import double4x3;

import skdefinitions;
import skrotationmatrix;
import skpointgroup;
import skpointsymmetryset;
import sktransformationmatrix;

namespace
{
using Atom = std::tuple<double3, std::size_t, double>;

// The atoms binned by their wrapped fractional position. The tests of a candidate symmetry operation ask, for every
// atom, whether an atom of its type lies within the precision of the atom's image, and read that from the bin of
// the image and the bins around it instead of from a pass over every atom, which made each test cost the square of
// the atom count. A bin is at least the grid's tolerance wide measured across the cell (the perpendicular width of
// the cell along an axis over the number of bins along it), so that two points within that tolerance of each other,
// by the component-wise wrapped difference the tests take, are in the same or in adjacent bins; and about an atom's
// share of the volume wide, so that few atoms are looked at. A query may ask for any tolerance up to the grid's.
struct AtomGrid
{
  const std::vector<Atom>& atoms;
  double3x3 unitCell;
  std::array<std::int32_t, 3> bins{1, 1, 1};
  // how near a bin's face, in bins, a point has to be for a match of it to lie in the bin beyond that face
  std::array<double, 3> margin{0.0, 0.0, 0.0};
  std::vector<std::size_t> first;    // where each bin's atoms start in 'ordered', with one entry past the last bin
  std::vector<std::size_t> ordered;  // the atom indices, bin by bin

  AtomGrid(const std::vector<Atom>& atoms, const double3x3& unitCell, double tolerance)
      : atoms(atoms), unitCell(unitCell)
  {
    const double volume = std::fabs(unitCell.determinant());
    if (volume > 0.0 && !atoms.empty())
    {
      const double3x3 inverse = unitCell.inverse();
      const double width = std::max(tolerance, std::cbrt(volume / double(atoms.size())));
      for (std::size_t k = 0; k < 3; ++k)
      {
        // the perpendicular width of the cell along axis k is one over the length of row k of the inverse
        const double rowLength =
            std::sqrt(inverse[0][k] * inverse[0][k] + inverse[1][k] * inverse[1][k] + inverse[2][k] * inverse[2][k]);
        // at most 64 along an axis, a quarter of a million bins in all; beyond that the bins hold more than an atom
        bins[k] = static_cast<std::int32_t>(std::clamp(std::floor(1.0 / (rowLength * width)), 1.0, 64.0));
        // the tolerance across the cell in bins, with a little to spare for the rounding of the comparisons
        margin[k] = (1.01 * tolerance * rowLength + 1.0e-12) * bins[k];
      }
    }

    first.assign(
        static_cast<std::size_t>(bins[0]) * static_cast<std::size_t>(bins[1]) * static_cast<std::size_t>(bins[2]) + 1,
        0);
    for (const Atom& atom : atoms) ++first[binOf(cellOf(std::get<0>(atom))) + 1];
    for (std::size_t b = 1; b < first.size(); ++b) first[b] += first[b - 1];
    std::vector<std::size_t> next(first.begin(), first.end() - 1);
    ordered.resize(atoms.size());
    for (std::size_t i = 0; i < atoms.size(); ++i) ordered[next[binOf(cellOf(std::get<0>(atoms[i])))]++] = i;
  }

  std::array<std::int32_t, 3> cellOf(const double3& position) const
  {
    const double3 wrapped = double3::fract(position);
    return {std::min(bins[0] - 1, static_cast<std::int32_t>(wrapped.x * bins[0])),
            std::min(bins[1] - 1, static_cast<std::int32_t>(wrapped.y * bins[1])),
            std::min(bins[2] - 1, static_cast<std::int32_t>(wrapped.z * bins[2]))};
  }

  std::size_t binOf(const std::array<std::int32_t, 3>& cell) const
  {
    return static_cast<std::size_t>((cell[2] * bins[1] + cell[1]) * bins[0] + cell[0]);
  }

  // The offset, atom minus position wrapped component-wise, to the nearest atom within the tolerance, of the type
  // unless any type will do.
  std::optional<double3> nearest(const double3& position, std::size_t type, bool anyType, double tolerance) const
  {
    // the bins to look in along each axis: the point's own, and the one beyond each face of it that the point is
    // within the margin of. A match is within the tolerance of the point on every axis, so a point further than
    // that from a face has its matches on its own side of it; the bins are wider than the tolerance, so for most
    // points that is one bin, where visiting all twenty-seven around the point cost more than the atoms in them.
    const double3 wrapped = double3::fract(position);
    std::array<std::array<std::int32_t, 3>, 3> around;
    std::array<std::size_t, 3> counts;
    for (std::size_t k = 0; k < 3; ++k)
    {
      const std::int32_t n = bins[k];
      const double scaled = wrapped[k] * n;
      const std::int32_t centre = std::min(n - 1, static_cast<std::int32_t>(scaled));
      const std::int32_t lower = (centre + n - 1) % n;
      const std::int32_t upper = (centre + 1) % n;
      around[k][0] = centre;
      counts[k] = 1;
      if (lower != centre && scaled - centre < margin[k]) around[k][counts[k]++] = lower;
      if (upper != centre && (upper != lower || counts[k] == 1) && 1.0 - (scaled - centre) < margin[k])
      {
        around[k][counts[k]++] = upper;
      }
    }

    double closestDistanceSquared = tolerance * tolerance;
    std::optional<double3> closest;
    for (std::size_t iz = 0; iz < counts[2]; ++iz)
    {
      for (std::size_t iy = 0; iy < counts[1]; ++iy)
      {
        for (std::size_t ix = 0; ix < counts[0]; ++ix)
        {
          const std::size_t bin = binOf({around[0][ix], around[1][iy], around[2][iz]});
          for (std::size_t slot = first[bin]; slot < first[bin + 1]; ++slot)
          {
            const Atom& atom = atoms[ordered[slot]];
            if (!anyType && std::get<1>(atom) != type) continue;
            double3 dr = std::get<0>(atom) - position;
            dr.x -= std::rint(dr.x);
            dr.y -= std::rint(dr.y);
            dr.z -= std::rint(dr.z);
            const double distanceSquared = (unitCell * dr).length_squared();
            if (distanceSquared < closestDistanceSquared)
            {
              closestDistanceSquared = distanceSquared;
              closest = dr;
            }
          }
        }
      }
    }
    return closest;
  }
};

// Whether the image of every atom under an operation lies within the precision of an atom of its type (of any atom,
// when partial occupancies are allowed).
template <typename Image>
bool everyImageOnAnAtom(const AtomGrid& grid, Image&& image, bool allowPartialOccupancies, double precision)
{
  for (const Atom& atom : grid.atoms)
  {
    if (!grid.nearest(image(std::get<0>(atom)), std::get<1>(atom), allowPartialOccupancies, precision)) return false;
  }
  return true;
}

// The translation of a symmetry operation averaged over all atoms: each atom's image under the rotation and the
// given translation is matched to the nearest atom of its type within the tolerance, and the mean of the remaining
// offsets is added to the translation. Without a match for every atom there is no refinement.
std::optional<double3> refinedTranslationOn(const AtomGrid& grid, const double3& translationVector,
                                            const SKRotationMatrix& rotationMatrix, bool allowPartialOccupancies,
                                            double tolerance)
{
  double3 offsetSum(0.0, 0.0, 0.0);
  for (const Atom& atom : grid.atoms)
  {
    const std::optional<double3> closest = grid.nearest(rotationMatrix * std::get<0>(atom) + translationVector,
                                                        std::get<1>(atom), allowPartialOccupancies, tolerance);
    if (!closest) return std::nullopt;
    offsetSum += *closest;
  }
  return translationVector + offsetSum / double(grid.atoms.size());
}

// Whether a transformed metric describes the cell within the precision: the lengths first, since most candidates
// fail on them, and the angles, which take arc-cosines, only for the ones that pass.
bool metricsSimilar(const SKSymmetryCell& metricCell, const double3x3& transformedMetricTensor,
                    double symmetryPrecision)
{
  const double3 lengthDifference(std::fabs(metricCell.a() - std::sqrt(transformedMetricTensor[0][0])),
                                 std::fabs(metricCell.b() - std::sqrt(transformedMetricTensor[1][1])),
                                 std::fabs(metricCell.c() - std::sqrt(transformedMetricTensor[2][2])));
  if (lengthDifference.length() > symmetryPrecision)
  {
    return false;
  }

  const SKSymmetryCell transformedMetricCell = SKSymmetryCell::createFromMetricTensor(transformedMetricTensor);
  const double3 angleDifference(std::sin(std::fabs(metricCell.alpha() - transformedMetricCell.alpha())),
                                std::sin(std::fabs(metricCell.beta() - transformedMetricCell.beta())),
                                std::sin(std::fabs(metricCell.gamma() - transformedMetricCell.gamma())));
  if (std::fabs(std::max({angleDifference.x, angleDifference.y, angleDifference.z})) > symmetryPrecision)
  {
    return false;
  }
  return true;
}
}  // namespace

SKSymmetryCell::SKSymmetryCell() {}

SKSymmetryCell::SKSymmetryCell(double a, double b, double c, double alpha, double beta, double gamma)
{
  _a = a;
  _b = b;
  _c = c;
  _alpha = alpha * std::numbers::pi / 180.0;
  _beta = beta * std::numbers::pi / 180.0;
  _gamma = gamma * std::numbers::pi / 180.0;
}

SKSymmetryCell SKSymmetryCell::createFromMetricTensor(double3x3 metricTensor)
{
  double A = metricTensor[0][0];
  double B = metricTensor[1][1];
  double C = metricTensor[2][2];

  SKSymmetryCell cell;
  cell._a = std::sqrt(A);
  cell._b = std::sqrt(B);
  cell._c = std::sqrt(C);
  // round-off can put the cosine a little beyond +-1 for (nearly) collinear axes, where acos gives NaN
  cell._alpha = std::acos(std::clamp(metricTensor[1][2] / (cell._b * cell._c), -1.0, 1.0));
  cell._beta = std::acos(std::clamp(metricTensor[0][2] / (cell._a * cell._c), -1.0, 1.0));
  cell._gamma = std::acos(std::clamp(metricTensor[0][1] / (cell._a * cell._b), -1.0, 1.0));
  return cell;
}

SKSymmetryCell SKSymmetryCell::createFromUnitCell(double3x3 unitCell)
{
  double3 column1 = unitCell[0];
  double3 column2 = unitCell[1];
  double3 column3 = unitCell[2];
  double length1 = column1.length();
  double length2 = column2.length();
  double length3 = column3.length();

  SKSymmetryCell cell;
  cell._a = length1;
  cell._b = length2;
  cell._c = length3;
  cell._alpha = std::acos(std::clamp(double3::dot(column2, column3) / (length2 * length3), -1.0, 1.0));
  cell._beta = std::acos(std::clamp(double3::dot(column1, column3) / (length1 * length3), -1.0, 1.0));
  cell._gamma = std::acos(std::clamp(double3::dot(column1, column2) / (length1 * length2), -1.0, 1.0));
  return cell;
}

SKSymmetryCell SKSymmetryCell::idealized(std::size_t pointGroupNumber, std::string qualifier)
{
  // int pointGroupNumber = spaceGroup.spaceGroupSetting().pointGroupNumber();
  Holohedry holohedry = SKPointGroup::pointGroupData[pointGroupNumber].holohedry();
  switch (holohedry)
  {
    case Holohedry::none:
      return *this;
    case Holohedry::triclinic:
    {
      double cg = std::cos(_gamma);
      double cb = std::cos(_beta);
      double ca = std::cos(_alpha);
      double sg = std::sin(_gamma);
      double temp = _c * std::sqrt(1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg) / sg;
      return SKSymmetryCell::createFromUnitCell(double3x3(double3(_a, 0.0, 0.0), double3(_b * cg, _b * sg, 0.0),
                                                          double3(_c * cb, _c * (ca - cb * cg) / sg, temp)));
    }
    case Holohedry::monoclinic:
      return SKSymmetryCell::createFromUnitCell(double3x3(double3(_a, 0.0, 0.0), double3(0.0, _b, 0.0),
                                                          double3(_c * std::cos(_beta), 0.0, _c * std::sin(_beta))));
    case Holohedry::orthorhombic:
      return SKSymmetryCell::createFromUnitCell(
          double3x3(double3(_a, 0.0, 0.0), double3(0.0, _b, 0.0), double3(0.0, 0.0, _c)));
    case Holohedry::tetragonal:
      return SKSymmetryCell::createFromUnitCell(
          double3x3(double3(0.5 * (_a + _b), 0.0, 0.0), double3(0.0, 0.5 * (_a + _b), 0.0), double3(0.0, 0.0, _c)));
    case Holohedry::trigonal:
      if (qualifier == "R")
      // if (spaceGroup.spaceGroupSetting().qualifier() == "R")
      {
        double avg = (_a + _b + _c) / 3.0;
        double angle = std::acos(std::clamp((std::cos(_gamma) + std::cos(_beta) + std::cos(_alpha)) / 3.0, -1.0, 1.0));
        // Reference, https://homepage.univie.ac.at/michael.leitner/lattice/struk/rgr.html
        double ahex = 2.0 * avg * std::sin(0.5 * angle);
        double chex = (_a + _b + _c) / 3.0 * std::sqrt(3.0 * (1.0 + 2.0 * std::cos(angle)));
        return SKSymmetryCell::createFromUnitCell(
            double3x3(double3(ahex / 2.0, -ahex / (2.0 * std::sqrt(3.0)), chex / 3.0),
                      double3(0.0, ahex / std::sqrt(3.0), chex / 3.0),
                      double3(-ahex / 2.0, -ahex / (2.0 * std::sqrt(3.0)), chex / 3.0)));
      }
      return SKSymmetryCell::createFromUnitCell(
          double3x3(double3(0.5 * (_a + _b), 0.0, 0.0),
                    double3(-(_a + _b) / 4.0, (_a + _b) / 4.0 * std::sqrt(3.0), 0.0), double3(0.0, 0.0, _c)));
    case Holohedry::hexagonal:
      return SKSymmetryCell::createFromUnitCell(
          double3x3(double3(0.5 * (_a + _b), 0.0, 0.0),
                    double3(-(_a + _b) / 4.0, (_a + _b) / 4.0 * std::sqrt(3.0), 0.0), double3(0.0, 0.0, _c)));
    case Holohedry::cubic:
    {
      double edge = (_a + _b + _c) / 3.0;
      return SKSymmetryCell::createFromUnitCell(
          double3x3(double3(edge, 0.0, 0.0), double3(0.0, edge, 0.0), double3(0.0, 0.0, edge)));
    }
    default:
      return SKSymmetryCell::createFromUnitCell(double3x3());
  }
}

double3x3 SKSymmetryCell::unitCell() const
{
  // exact zeros for right angles (see 'double3x3::fromCellParameters')
  return double3x3::fromCellParameters(_a, _b, _c, _alpha, _beta, _gamma);
}

double3x3 SKSymmetryCell::metricTensor()
{
  double half_xi = _b * _c * std::cos(_alpha);
  double half_eta = _a * _c * std::cos(_beta);
  double half_zeta = _a * _b * std::cos(_gamma);

  double3 v1 = double3(_a * _a, half_zeta, half_eta);
  double3 v2 = double3(half_zeta, _b * _b, half_xi);
  double3 v3 = double3(half_eta, half_xi, _c * _c);
  return double3x3(v1, v2, v3);
}

double SKSymmetryCell::volume()
{
  double cosAlpha = std::cos(_alpha);
  double cosBeta = std::cos(_beta);
  double cosGamma = std::cos(_gamma);
  double temp =
      1.0 - cosAlpha * cosAlpha - cosBeta * cosBeta - cosGamma * cosGamma + 2.0 * cosAlpha * cosBeta * cosGamma;
  return _a * _b * _c * std::sqrt(temp);
}

bool SKSymmetryCell::isOverlap(double3 a, double3 b, double3x3 lattice, double symmetryPrecision = 1e-2)
{
  // the signed fractional difference: a component-wise absolute value before the wrap flips the relative signs of
  // the components (0.01, -0.99 wraps to 0.01, 0.01 but |.| gives 0.01, 0.99 and then 0.01, -0.01), which in a
  // non-orthogonal cell changes the Cartesian distance
  double3 dr = a - b;
  dr.x -= std::rint(dr.x);
  dr.y -= std::rint(dr.y);
  dr.z -= std::rint(dr.z);
  if ((lattice * dr).length_squared() < symmetryPrecision * symmetryPrecision)
  {
    return true;
  }
  return false;
}

double3x3 SKSymmetryCell::findSmallestPrimitiveCell(
    const std::vector<std::tuple<double3, std::size_t, double>>& reducedAtoms,
    const std::vector<std::tuple<double3, std::size_t, double>>& atoms, double3x3 unitCell,
    bool allowPartialOccupancies, double symmetryPrecision = 1e-2)
{
  std::vector<double3> translationVectors{};

  if (reducedAtoms.size() > 0)
  {
    const AtomGrid grid(atoms, unitCell, symmetryPrecision);
    double3 origin = std::get<0>(reducedAtoms[0]);

    for (std::size_t i = 1; i < reducedAtoms.size(); i++)
    {
      double3 vec = std::get<0>(reducedAtoms[i]) - origin;

      if (everyImageOnAnAtom(
              grid, [&vec](const double3& position) { return position + vec; }, allowPartialOccupancies,
              symmetryPrecision))
      {
        double3 a = double3(vec.x - std::rint(vec.x), vec.y - std::rint(vec.y), vec.z - std::rint(vec.z));
        if (a.x < 0.0 - 1e-10)
        {
          a.x += 1.0;
        }
        if (a.y < 0.0 - 1e-10)
        {
          a.y += 1.0;
        }
        if (a.z < 0.0 - 1e-10)
        {
          a.z += 1.0;
        }
        translationVectors.push_back(a);
      }
    }

    translationVectors.push_back(double3(1.0, 0.0, 0.0));
    translationVectors.push_back(double3(0.0, 1.0, 0.0));
    translationVectors.push_back(double3(0.0, 0.0, 1.0));
  }

  std::size_t size = translationVectors.size();

  double3x3 smallestCell = unitCell;
  // the volume, not the signed determinant: a left-handed cell gave a negative minimum that no candidate undercut
  double initialVolume = std::fabs(unitCell.determinant());
  double minimumVolume = initialVolume;
  // a cell of lattice translations holds at least one atom, so anything smaller than half of that volume is a
  // (nearly) coplanar triple rather than a cell; the former absolute threshold of 1 Angstrom^3 rejected the
  // primitive cell of small or reduced-unit cells
  const double smallestPossibleVolume = initialVolume / static_cast<double>(std::max(std::size_t{1}, atoms.size()));

  for (std::size_t i = 0; i < size; i++)
  {
    const double3 tmpv1 = unitCell * translationVectors[i];
    for (std::size_t j = i + 1; j < size; j++)
    {
      const double3 tmpv2 = unitCell * translationVectors[j];
      const double faceArea12 = double3::cross(tmpv1, tmpv2).length();
      for (std::size_t k = j + 1; k < size; k++)
      {
        double3 tmpv3 = unitCell * translationVectors[k];
        double3x3 cell = double3x3(tmpv1, tmpv2, tmpv3);
        double volume = std::fabs(cell.determinant());

        // a translation, another and (nearly) their sum, read off noisy coordinates, span a sliver whose volume
        // passed the threshold below when the given basis was long and skewed: its smallest height, the volume over
        // the largest face, is within the precision, where no cell of a structure resolved at that precision is
        const double largestFaceArea =
            std::max({faceArea12, double3::cross(tmpv2, tmpv3).length(), double3::cross(tmpv1, tmpv3).length()});
        if (volume < 2.0 * symmetryPrecision * largestFaceArea) continue;

        if ((volume > 0.5 * smallestPossibleVolume) && (volume < minimumVolume))
        {
          minimumVolume = volume;
          smallestCell = double3x3(tmpv1, tmpv2, tmpv3);

          // initialVolume/volume is nearly integer
          if (std::size_t(std::rint(initialVolume / volume)) == (size - 2))
          {
            double3x3 relativeLattice = double3x3(translationVectors[i], translationVectors[j], translationVectors[k]);

            // inverse is nearly integer
            double3x3 temp = relativeLattice.inverse();
            const int3x3 temp2 = temp.toInt3x3();
            return unitCell * double3x3(temp2).inverse();
          }
        }
      }
    }
  }
  return smallestCell;
}

bool SKSymmetryCell::testTranslationalSymmetry(double3 translationVector,
                                               const std::vector<std::tuple<double3, std::size_t, double>>& atoms,
                                               double3x3 unitCell, bool allowPartialOccupancies,
                                               double precision = 1e-2)
{
  const AtomGrid grid(atoms, unitCell, precision);
  return everyImageOnAnAtom(
      grid, [&translationVector](const double3& position) { return position + translationVector; },
      allowPartialOccupancies, precision);
}

bool SKSymmetryCell::testSymmetry(double3 translationVector, SKRotationMatrix rotationMatrix,
                                  const std::vector<std::tuple<double3, std::size_t, double>>& atoms,
                                  double3x3 unitCell, bool allowPartialOccupancies, double precision = 1e-2)
{
  const AtomGrid grid(atoms, unitCell, precision);
  return everyImageOnAnAtom(
      grid, [&](const double3& position) { return rotationMatrix * position + translationVector; },
      allowPartialOccupancies, precision);
}

/// Computes translation vectors for symmetry operations
///
/// - parameter unitCell:            the unit cell
/// - parameter fractionalPositions: the fractional positions of the atomic configuration
/// - parameter rotationMatrix:      the symmetry elements
/// - parameter symmetryPrecision:   the precision of the search (default: 1e-2)
///
/// - returns: the list of translation vectors, including (0,0,0)
std::vector<double3> SKSymmetryCell::primitiveTranslationVectors(
    double3x3 unitCell, const std::vector<std::tuple<double3, std::size_t, double>>& reducedAtoms,
    const std::vector<std::tuple<double3, std::size_t, double>>& atoms, SKRotationMatrix rotationMatrix,
    bool allowPartialOccupancies, double symmetryPrecision = 1e-2)
{
  std::vector<double3> translationVectors{};

  if (reducedAtoms.size() > 0)
  {
    // wide enough for the refinement at twice the precision below
    const AtomGrid grid(atoms, unitCell, 2.0 * symmetryPrecision);
    double3 origin = rotationMatrix * std::get<0>(reducedAtoms[0]);

    for (std::size_t i = 0; i < reducedAtoms.size(); i++)
    {
      // the translation read off one pair of atoms carries the positional errors of both; averaged over every atom
      // it carries a fraction of them, which keeps the found operations, and the origin derived from them,
      // within the precision for noisy coordinates. A pair whose translation misses the precision by its own error
      // is refined at twice the precision first and then held to the precision.
      //
      // The refinement at the precision has a value exactly when every atom's image has an atom of its type within
      // the precision, which is the test of the operation; a separate test used to run the same search once more.
      double3 vec = std::get<0>(reducedAtoms[i]) - origin;
      std::optional<double3> refined =
          refinedTranslationOn(grid, vec, rotationMatrix, allowPartialOccupancies, symmetryPrecision);
      if (!refined)
      {
        refined = refinedTranslationOn(grid, vec, rotationMatrix, allowPartialOccupancies, 2.0 * symmetryPrecision);
        if (refined && !everyImageOnAnAtom(
                           grid, [&](const double3& position) { return rotationMatrix * position + *refined; },
                           allowPartialOccupancies, symmetryPrecision))
        {
          refined = std::nullopt;
        }
      }
      if (refined)
      {
        translationVectors.push_back(*refined);
      }
    }
  }
  return translationVectors;
}

std::optional<double3> SKSymmetryCell::refinedTranslation(
    double3 translationVector, SKRotationMatrix rotationMatrix,
    const std::vector<std::tuple<double3, std::size_t, double>>& atoms, double3x3 unitCell,
    bool allowPartialOccupancies, double tolerance)
{
  return refinedTranslationOn(AtomGrid(atoms, unitCell, tolerance), translationVector, rotationMatrix,
                              allowPartialOccupancies, tolerance);
}

/// Compute the Delaunay reduced cell
///
/// - parameter unitCell:          the original unit cell
/// - parameter symmetryPrecision: the precision of the cell
///
/// - returns: the Delaunay cell
///
///   We start with a lattice basis (b_i) 1≤ i ≤ n (n=2,3). This basis is extended by a factor b_n+1 = -(b_1 + ... + b_n
///   ) All scalar products  b_i . b_k (1 ≤ i < k ≤ n+1) are considered. The reduction is performed by mnimizing the
///   sum: sum = b_1^2 + ... + b_n+1^2. It can be shown that this sum can be reduced as long as one of the scalar
///   products is still positive. If e.g. the scalar product b_1 . b_2 is still positive, a transformation can be
///   performed such that the sum sum' of the transformed b_i^2 is smaller than sum: b_1' = -b_1 b_2' = b_2 b_3' = b_1 +
///   b_3 b_4' = b+1 + b_4
///
///   If all the scalar products are less than or equal to zero, the three shortest vectors forming the reduced basis
///   are contained in the set V = {b_1, b_2, b_3, b_4, b_1 + b_2, b_2 + b_3, b_3 + b_1} which corresponds to the
///   maximal set of faces of the Dirichlet domain (14 faces).
///
///   Reference: International Tables for Crystallography, Vol.A: Space Group Symmetry, page 747
namespace
{
// Sorts vectors by length, with lengths within 1e-10 (squared) of each other counting as equal. The vectors are
// ranked by exact length, consecutive ranks closer than the tolerance form one group, and std::sort then orders
// the vectors by group: a strict weak ordering that decides every pair as the former comparator "a + 1e-10 < b"
// did whenever that one was consistent, so the arrangement of equal lengths is unchanged. The former comparator
// was not a strict weak ordering (intransitive for lengths in steps smaller than the tolerance), which std::sort
// may answer with garbage.
void sortByLengthWithinTolerance(std::vector<double3>& vectors)
{
  std::vector<std::size_t> ranked(vectors.size());
  std::iota(ranked.begin(), ranked.end(), 0uz);
  std::sort(ranked.begin(), ranked.end(), [&vectors](std::size_t a, std::size_t b)
            { return vectors[a].length_squared() < vectors[b].length_squared(); });
  std::vector<std::size_t> group(vectors.size());
  for (std::size_t i = 1; i < ranked.size(); ++i)
  {
    const bool separate = vectors[ranked[i]].length_squared() - vectors[ranked[i - 1]].length_squared() >= 1e-10;
    group[ranked[i]] = group[ranked[i - 1]] + (separate ? 1 : 0);
  }
  std::vector<std::pair<std::size_t, double3>> keyed;
  keyed.reserve(vectors.size());
  for (std::size_t i = 0; i < vectors.size(); ++i) keyed.emplace_back(group[i], vectors[i]);
  std::sort(keyed.begin(), keyed.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
  for (std::size_t i = 0; i < vectors.size(); ++i) vectors[i] = keyed[i].second;
}
}  // namespace

std::optional<double3x3> SKSymmetryCell::computeDelaunayReducedCell(double3x3 unitCell, double symmetryPrecision = 1e-2)
{
  double3 additionalBasisVector = -(unitCell[0] + unitCell[1] + unitCell[2]);
  double4x3 extendedBasis = double4x3(unitCell[0], unitCell[1], unitCell[2], additionalBasisVector);

  bool somePositive = false;
  do
  {
    somePositive = false;
    // (i,j) in (0,1), (0,2), (0,3), (1,2), (1,3), (2,3); k,l denote the other two vectors
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>> basisIndices = {
        std::make_tuple(0, 1, 2, 3), std::make_tuple(0, 2, 1, 3), std::make_tuple(0, 3, 1, 2),
        std::make_tuple(1, 2, 0, 3), std::make_tuple(1, 3, 0, 2), std::make_tuple(2, 3, 0, 1)};
    for (const auto& [i, j, k, l] : basisIndices)
    {
      if (double3::dot(extendedBasis[i], extendedBasis[j]) > symmetryPrecision)
      {
        extendedBasis[k] += extendedBasis[i];
        extendedBasis[l] += extendedBasis[i];

        extendedBasis[i] = -extendedBasis[i];

        // start over (until all dotproducts are negative or zero)
        somePositive = true;
        break;
      }
    }
  } while (somePositive);

  // Search in the array {b1, b2, b3, b4, b1+b2, b2+b3, b3+b1}, sorted by length
  std::vector<double3> b = {extendedBasis[0],
                            extendedBasis[1],
                            extendedBasis[2],
                            extendedBasis[3],
                            extendedBasis[0] + extendedBasis[1],
                            extendedBasis[1] + extendedBasis[2],
                            extendedBasis[2] + extendedBasis[0]};

  sortByLengthWithinTolerance(b);

  // take the first two vectors, combined with a vector that has a non-zero, positive volume
  for (std::size_t i = 2; i < 7; i++)
  {
    double3x3 trialUnitCell = double3x3(b[0], b[1], b[i]);
    double volume = trialUnitCell.determinant();

    if (std::fabs(volume) > symmetryPrecision)
    {
      return (volume > 0) ? trialUnitCell : -trialUnitCell;
    }
  }

  return std::nullopt;
}

std::optional<double3x3> SKSymmetryCell::computeDelaunayReducedCell2D(double3x3 unitCell,
                                                                      double symmetryPrecision = 1e-2)
{
  double3x3 lattice2D = double3x3();
  lattice2D[0] = unitCell[0];
  lattice2D[1] = unitCell[2];

  double3x3 extendedBasis = double3x3(lattice2D[0], lattice2D[1], -(lattice2D[0] + lattice2D[1]));

  bool somePositive = false;
  do
  {
    somePositive = false;
    // (i,j) in (0,1), (0,2), (1,2); k denote the other two vectors
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> basisIndices = {
        std::make_tuple(0, 1, 2), std::make_tuple(0, 2, 1), std::make_tuple(1, 2, 0)};
    for (auto& [i, j, k] : basisIndices)
    {
      if (double3::dot(extendedBasis[i], extendedBasis[j]) > symmetryPrecision)
      {
        extendedBasis[k] += 2.0 * extendedBasis[i];
        extendedBasis[i] = -extendedBasis[i];

        // start over (until all dotproducts are negative or zero)
        somePositive = true;
        break;
      }
    }
  } while (somePositive);

  // Search in the set {b1, b2, b3, b1+b2}
  std::vector<double3> b = {extendedBasis[0], extendedBasis[1], extendedBasis[2], extendedBasis[0] + extendedBasis[1]};
  sortByLengthWithinTolerance(b);

  for (std::size_t i = 1; i < 4; i++)
  {
    double3x3 tmpmat = double3x3(b[0], unitCell[1], b[i]);

    if (std::fabs(tmpmat.determinant()) > symmetryPrecision)
    {
      extendedBasis[0] = b[0];
      extendedBasis[1] = b[i];
      break;
    }
  }

  double3x3 basis = double3x3();
  basis[0] = extendedBasis[0];
  basis[1] = unitCell[1];
  basis[2] = extendedBasis[1];

  double volume = basis.determinant();

  if (std::fabs(volume) < symmetryPrecision)
  {
    return std::nullopt;
  }

  if (volume < 0.0)
  {
    basis[1] = -basis[1];
  }

  return basis;
}

// Determining the lattice symmetry is equivalent to determining the Bravais type.
// Note: no default precision here. The definition carried '= 3.0', which no caller could see (the declaration has
// none), and which would have made checkMetricSimilarity's angle test, sin(|angle difference|) > precision, impossible
// to fail.
SKPointSymmetrySet SKSymmetryCell::findLatticeSymmetry(double3x3 reducedLattice, double symmetryPrecision)
{
  const std::vector<int3> latticeAxes = {
      int3(1, 0, 0),   int3(0, 1, 0),    int3(0, 0, 1),  int3(-1, 0, 0),  int3(0, -1, 0),  int3(0, 0, -1),
      int3(0, 1, 1),   int3(1, 0, 1),    int3(1, 1, 0),  int3(0, -1, -1), int3(-1, 0, -1), int3(-1, -1, 0),
      int3(0, 1, -1),  int3(-1, 0, 1),   int3(1, -1, 0), int3(0, -1, 1),  int3(1, 0, -1),  int3(-1, 1, 0),
      int3(1, 1, 1),   int3(-1, -1, -1), int3(-1, 1, 1), int3(1, -1, 1),  int3(1, 1, -1),  int3(1, -1, -1),
      int3(-1, 1, -1), int3(-1, -1, 1)};

  std::vector<SKRotationMatrix> pointSymmetries{};

  double3x3 latticeMetricMatrix = reducedLattice.transpose() * reducedLattice;
  // the cell of the lattice's metric, which every candidate is compared with, derived once: it used to be derived
  // again, with its square roots and arc-cosines, for each of the 6960 candidates
  const SKSymmetryCell metricCell = SKSymmetryCell::createFromMetricTensor(latticeMetricMatrix);

  // uses a stored list of all possible lattice vectors and loop over all possible permutations
  for (const int3& firstAxis : latticeAxes)
  {
    for (const int3& secondAxis : latticeAxes)
    {
      for (const int3& thirdAxis : latticeAxes)
      {
        SKRotationMatrix axes = SKRotationMatrix(firstAxis, secondAxis, thirdAxis);
        int determinant = axes.determinant();

        // if the determinant is 1 or -1 we have a (proper) rotation  (6960 proper rotations)
        if (determinant == 1 || determinant == -1)
        {
          double3x3 transformationMatrix = double3x3(axes.int3x3_m);

          // the inverse of a rotation matrix is its transpose, so we use the transpose here
          double3x3 newLattice = reducedLattice * transformationMatrix;
          double3x3 transformedLatticeMetricMatrix = newLattice.transpose() * newLattice;

          if (metricsSimilar(metricCell, transformedLatticeMetricMatrix, symmetryPrecision))
          {
            pointSymmetries.push_back(axes);
          }
        }
      }
    }
  }

  // distinct axis triples are distinct matrices; the former conjugation by lattice^-1 * lattice was the identity
  // up to inversion noise, rounded back to the same integers
  return SKPointSymmetrySet(pointSymmetries);
}

bool SKSymmetryCell::checkMetricSimilarity(double3x3 transformedMetricTensor, double3x3 metricTensor,
                                           double symmetryPrecision = 1e-2)
{
  return metricsSimilar(SKSymmetryCell::createFromMetricTensor(metricTensor), transformedMetricTensor,
                        symmetryPrecision);
}

std::vector<std::tuple<double3, std::size_t, double>> SKSymmetryCell::trim(
    const std::vector<std::tuple<double3, std::size_t, double>>& atoms, double3x3 from, double3x3 to,
    [[maybe_unused]] bool allowPartialOccupancies, double symmetryPrecision = 1e-2)
{
  double3x3 changeOfBasis = to.inverse() * from;

  std::vector<std::tuple<double3, std::size_t, double>> trimmedAtoms{};
  std::transform(
      atoms.begin(), atoms.end(), std::back_inserter(trimmedAtoms),
      [changeOfBasis](const std::tuple<double3, std::size_t, double>& atom) -> std::tuple<double3, std::size_t, double>
      {
        return std::make_tuple(double3::fract(changeOfBasis * std::get<0>(atom)), std::get<1>(atom), std::get<2>(atom));
      });

  std::vector<std::size_t> overlapTable = std::vector<std::size_t>(trimmedAtoms.size());
  std::fill(overlapTable.begin(), overlapTable.end(), -1);

  std::vector<std::tuple<double3, std::size_t, double>> result{};

  for (std::size_t i = 0; i < trimmedAtoms.size(); i++)
  {
    overlapTable[i] = i;
    for (std::size_t j = 0; j < trimmedAtoms.size(); j++)
    {
      // only atoms of one type are translation copies of each other. Partial occupancies let the symmetry search
      // (testTranslationalSymmetry, testSymmetry) treat a site shared by several elements as one site, but the
      // primitive cell must keep every element on it
      if (std::get<1>(trimmedAtoms[i]) == std::get<1>(trimmedAtoms[j]))
      {
        // signed difference, see isOverlap: an absolute value here changes the distance in a non-orthogonal cell
        double3 dr = std::get<0>(trimmedAtoms[i]) - std::get<0>(trimmedAtoms[j]);
        dr.x -= std::rint(dr.x);
        dr.y -= std::rint(dr.y);
        dr.z -= std::rint(dr.z);
        if ((to * dr).length_squared() < symmetryPrecision * symmetryPrecision)
        {
          if (overlapTable[j] == j)
          {
            overlapTable[i] = j;
            break;
          }
        }
      }
    }
  }

  for (std::size_t i = 0; i < trimmedAtoms.size(); i++)
  {
    if (overlapTable[i] == i)
    {
      result.push_back(trimmedAtoms[i]);
    }
  }

  return result;
}

const double epsilon = 1.0e-5;

bool SKSymmetryCell::isSmallerThen(double x, double y) { return x < (y - epsilon); }

bool SKSymmetryCell::isLargerThen(double x, double y) { return isSmallerThen(y, x); }

bool SKSymmetryCell::isSmallerEqualThen(double x, double y) { return !(y < (x - epsilon)); }

bool SKSymmetryCell::isLargerEqualThen(double x, double y) { return !(x < (y - epsilon)); }

bool SKSymmetryCell::isEqualTo(double x, double y) { return !((x < (y - epsilon)) || (y < (x - epsilon))); }

bool SKSymmetryCell::isLargerThenZeroXiEtaZeta(double xi, double eta, double zeta)
{
  std::size_t n_positive = 0;
  std::size_t n_zero = 0;

  if (isSmallerThen(0, xi))
  {
    n_positive += 1;
  }
  else if (!isSmallerThen(xi, 0.0))
  {
    n_zero += 1;
  }

  if (isSmallerThen(0.0, eta))
  {
    n_positive += 1;
  }
  else if (!isSmallerThen(eta, 0.0))
  {
    n_zero += 1;
  }

  if (isSmallerThen(0.0, zeta))
  {
    n_positive += 1;
  }
  else if (!isSmallerThen(zeta, 0.0))
  {
    n_zero += 1;
  }

  return ((n_positive == 3) || (n_zero == 0 && n_positive == 1));
}

template <typename T>
int sign(T val)
{
  return (T(0) < val) - (val < T(0));
}

/// Ref: I. Krivy, B. Gruber,  "A Unified Algorithm for Determining the Reduced (Niggli) Cell",  Acta Cryst. (1976).
/// A32, 297-298
///    R. W. Grosse-Kunstleve, N. K. Sauter and P. D. Adams, "Numerically stable algorithms for the computation of
///    reduced unit cells", Acta Cryst. (2004). A60, 1-6
std::optional<std::pair<SKSymmetryCell, SKTransformationMatrix>>
SKSymmetryCell::computeReducedNiggliCellAndChangeOfBasisMatrix()
{
  int counter = 0;

  // step 0:
  double A = (_a * _a);
  double B = (_b * _b);
  double C = (_c * _c);
  double xi = (2.0 * _b * _c * std::cos(_alpha));
  double eta = (2.0 * _a * _c * std::cos(_beta));
  double zeta = (2.0 * _a * _b * std::cos(_gamma));

  // Every branch condition below is homogeneous of degree one in (A, B, C, xi, eta, zeta), so scaling all six
  // quantities by a common positive factor leaves the chosen integer change-of-basis matrices unchanged. We
  // normalize by a characteristic squared length (the geometric mean of A, B, C) so that the fixed absolute
  // tolerance `epsilon` behaves as a relative tolerance. This follows the recommendation of Grosse-Kunstleve
  // et al. (2004) to scale the tolerance by (a b c)^(1/3), which in squared-length units equals (A B C)^(1/3).
  // The scale is only used for the comparisons; the metric quantities are restored to absolute units before
  // constructing the reduced cell, and the returned change-of-basis matrix is unaffected.
  double scale = std::cbrt(A * B * C);
  if (!(scale > 0.0))
  {
    scale = 1.0;
  }
  A /= scale;
  B /= scale;
  C /= scale;
  xi /= scale;
  eta /= scale;
  zeta /= scale;

  SKTransformationMatrix changeOfBasisMatrix = SKTransformationMatrix::identity;

algorithmStart:
{
  counter += 1;
  if (counter > 10000)
  {
    return std::nullopt;
  }

  // step 1
  if (SKSymmetryCell::isLargerThen(A, B) ||
      (SKSymmetryCell::isEqualTo(A, B) && (SKSymmetryCell::isLargerThen(std::fabs(xi), std::fabs(eta)))))
  {
    SKTransformationMatrix matrixC = SKTransformationMatrix(int3(0, -1, 0), int3(-1, 0, 0), int3(0, 0, -1));
    // assert(matrixC.determinant() == 1);
    changeOfBasisMatrix = changeOfBasisMatrix * matrixC;

    // Swap x, y and ensures proper sign of determinant
    std::swap(A, B);
    std::swap(xi, eta);
  }

  // step 2
  if (SKSymmetryCell::isLargerThen(B, C) ||
      (SKSymmetryCell::isEqualTo(B, C) && (SKSymmetryCell::isLargerThen(std::fabs(eta), std::fabs(zeta)))))
  {
    SKTransformationMatrix matrixC = SKTransformationMatrix(int3(-1, 0, 0), int3(0, 0, -1), int3(0, -1, 0));
    // assert(matrixC.determinant() == 1);
    changeOfBasisMatrix = changeOfBasisMatrix * matrixC;

    // Swap y, z and ensures proper sign of determinant
    std::swap(B, C);
    std::swap(eta, zeta);

    goto algorithmStart;
  }

  // step 3
  if (SKSymmetryCell::isLargerThenZeroXiEtaZeta(xi, eta, zeta))
  {
    std::array<int, 3> f = {1, 1, 1};
    if (SKSymmetryCell::isSmallerThen(xi, 0.0))
    {
      f[0] = -1;
    }
    if (SKSymmetryCell::isSmallerThen(eta, 0.0))
    {
      f[1] = -1;
    }
    if (SKSymmetryCell::isSmallerThen(zeta, 0.0))
    {
      f[2] = -1;
    }
    SKTransformationMatrix matrixC = SKTransformationMatrix(int3(f[0], 0, 0), int3(0, f[1], 0), int3(0, 0, f[2]));
    // assert(matrixC.determinant() == 1);
    changeOfBasisMatrix = changeOfBasisMatrix * matrixC;

    xi = std::fabs(xi);
    eta = std::fabs(eta);
    zeta = std::fabs(zeta);
  }
  else  // step 4:
  {
    std::size_t p = 0;
    std::array<int, 3> f = {1, 1, 1};
    if (SKSymmetryCell::isLargerThen(xi, 0.0))
    {
      f[0] = -1;
    }
    else if (!SKSymmetryCell::isSmallerThen(xi, 0.0))
    {
      p = 0;
    }
    if (SKSymmetryCell::isLargerThen(eta, 0.0))
    {
      f[1] = -1;
    }
    else if (!SKSymmetryCell::isSmallerThen(eta, 0.0))
    {
      p = 1;
    }
    if (SKSymmetryCell::isLargerThen(zeta, 0.0))
    {
      f[2] = -1;
    }
    else if (!SKSymmetryCell::isSmallerThen(zeta, 0.0))
    {
      p = 2;
    }
    if (f[0] * f[1] * f[2] < 0)
    {
      f[p] = -1;
    }
    SKTransformationMatrix matrixC = SKTransformationMatrix(int3(f[0], 0, 0), int3(0, f[1], 0), int3(0, 0, f[2]));
    // assert(matrixC.determinant() == 1);
    changeOfBasisMatrix = changeOfBasisMatrix * matrixC;

    xi = -std::fabs(xi);
    eta = -std::fabs(eta);
    zeta = -std::fabs(zeta);
  }

  // step 5
  if ((SKSymmetryCell::isLargerThen(std::fabs(xi), B)) ||
      (SKSymmetryCell::isEqualTo(xi, B) && SKSymmetryCell::isSmallerThen(2.0 * eta, zeta)) ||
      (SKSymmetryCell::isEqualTo(xi, -B) && SKSymmetryCell::isSmallerThen(zeta, 0.0)))
  {
    SKTransformationMatrix matrixC = SKTransformationMatrix(int3(1, 0, 0), int3(0, 1, 0), int3(0, -sign(xi), 1));
    // assert(matrixC.determinant() == 1);
    changeOfBasisMatrix = changeOfBasisMatrix * matrixC;

    C = B + C - xi * sign(xi);
    eta = eta - zeta * sign(xi);
    xi = xi - 2.0 * B * sign(xi);

    goto algorithmStart;
  }

  // step 6
  if ((SKSymmetryCell::isLargerThen(std::fabs(eta), A)) ||
      (SKSymmetryCell::isEqualTo(eta, A) && SKSymmetryCell::isSmallerThen(2.0 * xi, zeta)) ||
      (SKSymmetryCell::isEqualTo(eta, -A) && SKSymmetryCell::isSmallerThen(zeta, 0.0)))
  {
    SKTransformationMatrix matrixC = SKTransformationMatrix(int3(1, 0, 0), int3(0, 1, 0), int3(-sign(eta), 0, 1));
    // assert(matrixC.determinant() == 1);
    changeOfBasisMatrix = changeOfBasisMatrix * matrixC;

    C = A + C - eta * sign(eta);
    xi = xi - zeta * sign(eta);
    eta = eta - 2.0 * A * sign(eta);

    goto algorithmStart;
  }

  // step7
  if ((SKSymmetryCell::isLargerThen(std::fabs(zeta), A)) ||
      (SKSymmetryCell::isEqualTo(zeta, A) && SKSymmetryCell::isSmallerThen(2.0 * xi, eta)) ||
      (SKSymmetryCell::isEqualTo(zeta, -A) && SKSymmetryCell::isSmallerThen(eta, 0.0)))
  {
    SKTransformationMatrix matrixC = SKTransformationMatrix(int3(1, 0, 0), int3(-sign(zeta), 1, 0), int3(0, 0, 1));
    // assert(matrixC.determinant() == 1);
    changeOfBasisMatrix = changeOfBasisMatrix * matrixC;

    B = A + B - zeta * sign(zeta);
    xi = xi - eta * sign(zeta);
    zeta = zeta - 2.0 * A * sign(zeta);

    goto algorithmStart;
  }

  // step 8
  if (SKSymmetryCell::isSmallerThen(xi + eta + zeta + A + B, 0.0) ||
      (SKSymmetryCell::isEqualTo(xi + eta + zeta + A + B, 0.0) &&
       SKSymmetryCell::isLargerThen(2.0 * (A + eta) + zeta, 0.0)))
  {
    SKTransformationMatrix matrixC = SKTransformationMatrix(int3(1, 0, 0), int3(0, 1, 0), int3(1, 1, 1));
    // assert(matrixC.determinant() == 1);
    changeOfBasisMatrix = changeOfBasisMatrix * matrixC;

    C = A + B + C + xi + eta + zeta;
    xi = 2.0 * B + xi + zeta;
    eta = 2.0 * A + eta + zeta;

    goto algorithmStart;
  }
}

  // Restore absolute (squared-length) units. Angles are scale-invariant, but lengths must be rescaled.
  A *= scale;
  B *= scale;
  C *= scale;
  xi *= scale;
  eta *= scale;
  zeta *= scale;

  auto safeAcosDegrees = [](double cosine)
  { return std::acos(std::clamp(cosine, -1.0, 1.0)) * 180.0 / std::numbers::pi; };

  SKSymmetryCell cell = SKSymmetryCell(std::sqrt(A), std::sqrt(B), std::sqrt(C),
                                       safeAcosDegrees(xi / (2.0 * std::sqrt(B) * std::sqrt(C))),
                                       safeAcosDegrees(eta / (2.0 * std::sqrt(A) * std::sqrt(C))),
                                       safeAcosDegrees(zeta / (2.0 * std::sqrt(A) * std::sqrt(B))));

  return std::make_pair(cell, changeOfBasisMatrix);
}
