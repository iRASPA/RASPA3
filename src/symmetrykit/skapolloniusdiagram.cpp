module;

module skapolloniusdiagram;

import std;

import int3;
import double3;
import double3x3;
import skvoronoi;

// Implementation notes
// --------------------
// Two tangency problems arise, and both collapse onto the same small linear system. Requiring a
// sphere of radius t centred at c to touch site i from outside gives |c - x_i| = r_i + t. Squaring
// yields |c|² - 2 c·x_i + |x_i|² = t² + 2 r_i t + r_i², and subtracting the equation for one chosen
// site cancels both the |c|² and the t² term, leaving conditions that are *linear* in the four
// unknowns (c, t):
//
//     2 (x_i - x_0)·c + 2 (r_i - r_0) t = (|x_i|² - r_i²) - (|x_0|² - r_0²).
//
// Four sites give three such rows, which is exactly enough to write c as an affine function of t;
// substituting back into any one tangency condition leaves a single quadratic in t. That is the
// vertex case. For a point on the trisector curve of three sites, two rows come from tangency and
// the third from the sampling plane, and the same solve applies. So one routine,
// solveTangencyRows, does the algebra for both, and neither case needs iteration or a high-degree
// predicate.
//
// The expensive part is not the algebra but deciding which quadruples of sites to try. Candidates
// are drawn from radical-diagram adjacency: the radical diagram is cheap, coincides with the
// Apollonius diagram when the radii are equal, and stays combinatorially close as they spread.
// Every candidate is then certified against an exhaustive neighbourhood search, so adjacency only
// influences which vertices are *found*, never whether a found vertex is real.

namespace
{
// Three rows in four unknowns leave a line, and a quadratic meets a line at most twice, so a tangency has
// at most two solutions and they go in a fixed place. Millions of quadruples are solved and nearly all of
// them are thrown away again, so what a solution is returned in is worth as much as the algebra: a vector
// here would spend more time in the allocator than in the arithmetic.
struct TangentSpheres
{
  std::array<SKApolloniusTangentSphere, 2> spheres{};
  std::size_t count{0};

  void add(const SKApolloniusTangentSphere& sphere) { spheres[count++] = sphere; }
  const SKApolloniusTangentSphere* begin() const { return spheres.data(); }
  const SKApolloniusTangentSphere* end() const { return spheres.data() + count; }
  bool empty() const { return count == 0; }
};

// Solves three rows that are linear in the four unknowns (c, t) together with the tangency condition for
// the reference site, returning the at most two (c, t), restricted to t >= 0 unless negative clearance is
// admitted. Rows are given as `rows[i]·c + timeCoefficients[i] * t = constants[i]`.
//
// Three rows in four unknowns leave a line, and the tangency condition, being quadratic, meets that line at
// most twice. The line is taken from the null space of the rows rather than by writing c as a function of t
// and substituting, because that function does not always exist. Three collinear site centres contribute
// two rows with parallel coefficients of c, leaving those coefficients of rank two; what the rows then fix
// is t, with c free to move along a direction, which is the circle of tangent spheres that three collinear
// spheres share. It is the same line either way, so reading it off the null space treats the arrangements
// alike, and the arrangement that has no c as a function of t is not special: it is simply the one whose
// null vector has no component along t.
TangentSpheres solveTangencyRows(const std::array<double3, 3>& rows, const double3& timeCoefficients,
                                 const double3& constants, const double3& referenceCentre, double referenceRadius,
                                 bool allowNegativeRadius = false)
{
  // The rows as three of four columns, the last column holding the coefficients of t.
  std::array<std::array<double, 4>, 3> system;
  double scale = 0.0;
  for (std::size_t i = 0; i < 3; ++i)
  {
    system[i] = {rows[i].x, rows[i].y, rows[i].z, timeCoefficients[i]};
    scale = std::max(scale, std::sqrt(rows[i].length_squared() + timeCoefficients[i] * timeCoefficients[i]));
  }

  // The determinant of the rows with one column struck out. Alternating in sign these are the components of
  // the null vector, the generalised cross product of the three rows.
  auto minorWithout = [&](std::size_t skipped)
  {
    std::array<std::size_t, 3> kept{};
    for (std::size_t column = 0, written = 0; column < 4; ++column)
      if (column != skipped) kept[written++] = column;
    auto entry = [&](std::size_t row, std::size_t column) { return system[row][kept[column]]; };
    return entry(0, 0) * (entry(1, 1) * entry(2, 2) - entry(1, 2) * entry(2, 1)) -
           entry(0, 1) * (entry(1, 0) * entry(2, 2) - entry(1, 2) * entry(2, 0)) +
           entry(0, 2) * (entry(1, 0) * entry(2, 1) - entry(1, 1) * entry(2, 0));
  };

  std::array<double, 4> along{minorWithout(0), -minorWithout(1), minorWithout(2), -minorWithout(3)};
  double alongNorm = std::sqrt(along[0] * along[0] + along[1] * along[1] + along[2] * along[2] + along[3] * along[3]);

  // Every minor vanishing means the rows have rank below three and leave a surface rather than a line, so
  // there is no isolated solution to find: four concyclic centres of equal radii, say, or for the trisector
  // case a sampling plane that holds the curve instead of cutting it.
  if (scale <= 0.0 || alongNorm < 1.0e-12 * scale * scale * scale) return {};
  for (double& component : along) component /= alongNorm;

  // A point on the line, found by striking out the column whose minor is largest, which is the best
  // conditioned of the four ways to do it, and holding that unknown at zero.
  std::size_t held = 0;
  for (std::size_t column = 1; column < 4; ++column)
    if (std::abs(along[column]) > std::abs(along[held])) held = column;
  std::array<std::size_t, 3> kept{};
  for (std::size_t column = 0, written = 0; column < 4; ++column)
    if (column != held) kept[written++] = column;

  // double3x3 is built from columns, so each column of the reduced rows becomes one of its columns.
  double3x3 reduced(double3(system[0][kept[0]], system[1][kept[0]], system[2][kept[0]]),
                    double3(system[0][kept[1]], system[1][kept[1]], system[2][kept[1]]),
                    double3(system[0][kept[2]], system[1][kept[2]], system[2][kept[2]]));
  double3 solved = reduced.inverse() * constants;

  std::array<double, 4> base{};
  base[held] = 0.0;
  for (std::size_t column = 0; column < 3; ++column) base[kept[column]] = solved[column];

  double3 baseCentre(base[0], base[1], base[2]);
  double baseRadius = base[3];
  double3 alongCentre(along[0], along[1], along[2]);
  double alongRadius = along[3];

  // Substituting c = baseCentre + s * alongCentre and t = baseRadius + s * alongRadius into the unsquared
  // tangency |c - x_ref| = r_ref + t, squared.
  double3 offset = baseCentre - referenceCentre;
  double reach = referenceRadius + baseRadius;
  double a = double3::dot(alongCentre, alongCentre) - alongRadius * alongRadius;
  double b = 2.0 * (double3::dot(offset, alongCentre) - reach * alongRadius);
  double c = double3::dot(offset, offset) - reach * reach;

  std::array<double, 2> roots{};
  std::size_t rootCount = 0;
  if (std::abs(a) < 1.0e-14)
  {
    if (std::abs(b) > 1.0e-14 * scale) roots[rootCount++] = -c / b;
  }
  else
  {
    double discriminant = b * b - 4.0 * a * c;
    if (discriminant < 0.0) return {};
    double squareRoot = std::sqrt(discriminant);
    roots[rootCount++] = (-b + squareRoot) / (2.0 * a);
    roots[rootCount++] = (-b - squareRoot) / (2.0 * a);
  }

  TangentSpheres spheres;
  for (std::size_t root = 0; root < rootCount; ++root)
  {
    double s = roots[root];
    double t = baseRadius + s * alongRadius;
    if (t < 0.0 && !allowNegativeRadius) continue;

    // The rows hold the tangency conditions squared, so a root only satisfies the unsquared condition
    // |c - x| = r + t when r + t is not negative. With t confined to be non-negative this is automatic;
    // once negative t is admitted it has to be imposed, or the extraneous roots of the squaring are
    // taken for vertices.
    if (referenceRadius + t < 0.0) continue;

    spheres.add(SKApolloniusTangentSphere{baseCentre + s * alongCentre, t});
  }
  return spheres;
}

// The row contributed by requiring tangency to `centre`/`radius` relative to a reference site.
void tangencyRow(const double3& centre, double radius, const double3& referenceCentre, double referenceRadius,
                 double3& row, double& timeCoefficient, double& constant)
{
  row = 2.0 * (centre - referenceCentre);
  timeCoefficient = 2.0 * (radius - referenceRadius);
  constant = (double3::dot(centre, centre) - radius * radius) -
             (double3::dot(referenceCentre, referenceCentre) - referenceRadius * referenceRadius);
}

// The points of the trisector of three sites at one value of the clearance.
//
// Imposing all three tangencies at a fixed clearance leaves two linear conditions on the centre, whose
// intersection is a line, and one quadratic, so a trisector meets each clearance level in at most two
// points, one on each of the two branches that make up the curve. Sweeping the clearance therefore
// traverses the whole curve, and, usefully, the sweep runs out at a finite clearance exactly when the
// curve is closed: an open trisector reaches every clearance above its minimum.
std::vector<double3> trisectorPointsAtClearance(const std::array<double3, 3>& centres,
                                                const std::array<double, 3>& radii, double clearance)
{
  double referenceDistance = radii[0] + clearance;
  if (referenceDistance < 0.0) return {};

  std::array<double3, 2> rows;
  std::array<double, 2> constants;
  for (std::size_t i = 0; i < 2; ++i)
  {
    double timeCoefficient = 0.0;
    tangencyRow(centres[i + 1], radii[i + 1], centres[0], radii[0], rows[i], timeCoefficient, constants[i]);
    constants[i] -= timeCoefficient * clearance;
  }

  double3 direction = double3::cross(rows[0], rows[1]);
  double directionLength = direction.length();
  double scale = std::max(rows[0].length(), rows[1].length());
  if (scale <= 0.0 || directionLength < 1.0e-12 * scale * scale) return {};

  // Any point on the line where the two planes meet; the third row pins it down along the line.
  double3x3 matrix(double3(rows[0].x, rows[1].x, direction.x), double3(rows[0].y, rows[1].y, direction.y),
                   double3(rows[0].z, rows[1].z, direction.z));
  if (std::abs(matrix.determinant()) < 1.0e-30) return {};
  double3 base = matrix.inverse() * double3(constants[0], constants[1], 0.0);

  double3 offset = base - centres[0];
  double a = double3::dot(direction, direction);
  double b = 2.0 * double3::dot(offset, direction);
  double c = double3::dot(offset, offset) - referenceDistance * referenceDistance;
  double discriminant = b * b - 4.0 * a * c;
  if (discriminant < 0.0) return {};

  double root = std::sqrt(discriminant);
  std::vector<double3> points;
  points.push_back(base + ((-b + root) / (2.0 * a)) * direction);
  if (root > 1.0e-9 * (std::abs(b) + a)) points.push_back(base + ((-b - root) / (2.0 * a)) * direction);
  return points;
}

// A point of the bisector sheet of two sites, at one clearance and one angle about the axis joining their
// centres. The sheet is one branch of a hyperboloid of revolution about that axis, so at each clearance
// the two tangent spheres meet in a circle, and the circle shrinks to a point at the apex below which the
// sheet does not reach.
std::optional<double3> bisectorSheetPoint(const double3& firstCentre, double firstRadius,
                                         const double3& secondCentre, double secondRadius, double clearance,
                                         double angle)
{
  double firstDistance = firstRadius + clearance;
  double secondDistance = secondRadius + clearance;
  if (firstDistance < 0.0 || secondDistance < 0.0) return std::nullopt;

  double3 axis = secondCentre - firstCentre;
  double separation = axis.length();
  if (separation < 1.0e-12) return std::nullopt;
  axis = (1.0 / separation) * axis;

  double along = (separation * separation + firstDistance * firstDistance - secondDistance * secondDistance) /
                 (2.0 * separation);
  double radialSquared = firstDistance * firstDistance - along * along;
  if (radialSquared < 0.0) return std::nullopt;
  double radial = std::sqrt(radialSquared);

  double3 reference = std::abs(axis.x) < 0.9 ? double3(1.0, 0.0, 0.0) : double3(0.0, 1.0, 0.0);
  double3 inPlane = double3::cross(axis, reference);
  inPlane = (1.0 / inPlane.length()) * inPlane;
  double3 across = double3::cross(axis, inPlane);

  return firstCentre + along * axis + (radial * std::cos(angle)) * inPlane + (radial * std::sin(angle)) * across;
}

// Identifies a vertex by its four sites in a translation-invariant way, so that every periodic
// copy, and every one of the four sites that generates it, yields the same key.
//
// The four sites alone are not an identity. A quadruple admits up to two distinct tangent spheres,
// the two roots of the quadratic, and both can be empty and so both be genuine vertices; Kamarianakis
// distinguishes them as v_ijkl and v_ikjl by the orientation of the tetrahedron of tangency points.
// Keying on the sites alone therefore collides two different vertices and silently drops one, which
// leaves the arcs between them unpaired. `orientation` separates them.
struct SiteTuple
{
  std::array<std::int64_t, 16> data{};
  std::int64_t orientation{0};
  bool operator<(const SiteTuple& other) const
  {
    if (data != other.data) return data < other.data;
    return orientation < other.orientation;
  }
  bool operator==(const SiteTuple& other) const { return data == other.data && orientation == other.orientation; }
};

// For each position, one of the positions it coincides with to within `tolerance`, the same one for all of
// a coinciding set, so that what shares an answer is what is the same point. Coincidence is periodic: a
// point near a face of the cell coincides with one near the opposite face.
//
// A grid of the tolerance's own size is laid over the cell and the points are dropped into it, so a point
// can only coincide with the points in its own cell of the grid or in one of the twenty-six around it, and
// what has to be searched does not grow with how many points coincide. Nothing is gained by looking further
// than the first match, the rest of the set being reached through it.
//
// The grid may be coarser than the tolerance and still be right, since what it has to do is put every pair
// closer than the tolerance within a cell of each other, which a coarser grid also does; that is what lets
// the division count be capped, and it is capped so that a very large cell cannot ask for more divisions
// than there are numbers to count them with. Only the occupied cells are held, so the count costs nothing.
std::vector<std::size_t> coincidentGroups(const double3x3& unitCell, const double3x3& inverseCell,
                                          const double3& perpendicularWidths, const std::vector<double3>& positions,
                                          double tolerance)
{
  constexpr std::int32_t mostDivisions = 1 << 20;
  auto divisionsAlong = [&](double width)
  { return std::clamp(static_cast<std::int32_t>(std::min(width / tolerance, 1.0e9)), 1, mostDivisions); };
  int3 divisions(divisionsAlong(perpendicularWidths.x), divisionsAlong(perpendicularWidths.y),
                 divisionsAlong(perpendicularWidths.z));

  auto cellOf = [&](const double3& position)
  {
    double3 fractional = double3::fract(inverseCell * position);
    return int3(std::min(divisions.x - 1, static_cast<std::int32_t>(fractional.x * divisions.x)),
                std::min(divisions.y - 1, static_cast<std::int32_t>(fractional.y * divisions.y)),
                std::min(divisions.z - 1, static_cast<std::int32_t>(fractional.z * divisions.z)));
  };

  struct CellHash
  {
    std::size_t operator()(const std::array<std::int32_t, 3>& cell) const
    {
      std::size_t hash = 1469598103934665603ull;
      for (std::int32_t coordinate : cell) hash = (hash ^ static_cast<std::size_t>(coordinate)) * 1099511628211ull;
      return hash;
    }
  };
  std::unordered_map<std::array<std::int32_t, 3>, std::vector<std::size_t>, CellHash> occupants;

  std::vector<std::size_t> group(positions.size());
  for (std::size_t p = 0; p < positions.size(); ++p)
  {
    group[p] = p;
    int3 own = cellOf(positions[p]);

    bool found = false;
    for (std::int32_t ox = -1; ox <= 1 && !found; ++ox)
      for (std::int32_t oy = -1; oy <= 1 && !found; ++oy)
        for (std::int32_t oz = -1; oz <= 1 && !found; ++oz)
        {
          std::array<std::int32_t, 3> neighbour{(own.x + ox + divisions.x) % divisions.x,
                                                (own.y + oy + divisions.y) % divisions.y,
                                                (own.z + oz + divisions.z) % divisions.z};
          auto entry = occupants.find(neighbour);
          if (entry == occupants.end()) continue;
          for (std::size_t other : entry->second)
          {
            double3 delta = inverseCell * (positions[p] - positions[other]);
            delta = double3(delta.x - std::round(delta.x), delta.y - std::round(delta.y),
                            delta.z - std::round(delta.z));
            if ((unitCell * delta).length() < tolerance)
            {
              group[p] = other;
              found = true;
              break;
            }
          }
        }

    occupants[{own.x, own.y, own.z}].push_back(p);
  }
  return group;
}

// For the sets that are asked about a key far more often than they are read in order. Ordinary FNV-1a
// over the words of the key.
struct SiteWordsHash
{
  std::size_t operator()(const std::vector<std::int64_t>& words) const
  {
    std::size_t hash = 1469598103934665603ull;
    for (std::int64_t word : words) hash = (hash ^ static_cast<std::size_t>(word)) * 1099511628211ull;
    return hash;
  }
};

struct SiteTupleHash
{
  std::size_t operator()(const SiteTuple& key) const
  {
    std::size_t hash = 1469598103934665603ull;
    for (std::int64_t word : key.data) hash = (hash ^ static_cast<std::size_t>(word)) * 1099511628211ull;
    return (hash ^ static_cast<std::size_t>(key.orientation)) * 1099511628211ull;
  }
};

// Sorts (site, image) pairs and rebases the images on the first pair, which removes the lattice
// translation. `count` pairs are used; the remainder of the key stays zero. When `permutation` is
// given it receives the sorted order, so that a caller can evaluate the orientation of the tangency
// tetrahedron in the same canonical order and get a sign that does not depend on how the quadruple
// happened to be enumerated.
SiteTuple canonicalTuple(const std::size_t* indices, const int3* images, std::size_t count, int3& appliedOffset,
                         std::array<std::size_t, 4>* permutation = nullptr)
{
  std::array<std::array<std::int64_t, 5>, 4> pairs{};
  for (std::size_t i = 0; i < count; ++i)
    pairs[i] = {static_cast<std::int64_t>(indices[i]), images[i].x, images[i].y, images[i].z,
                static_cast<std::int64_t>(i)};
  std::sort(pairs.begin(), pairs.begin() + static_cast<std::ptrdiff_t>(count),
            [](const auto& lhs, const auto& rhs)
            {
              return std::tie(lhs[0], lhs[1], lhs[2], lhs[3]) < std::tie(rhs[0], rhs[1], rhs[2], rhs[3]);
            });

  appliedOffset = int3(static_cast<std::int32_t>(pairs[0][1]), static_cast<std::int32_t>(pairs[0][2]),
                       static_cast<std::int32_t>(pairs[0][3]));

  SiteTuple key;
  for (std::size_t i = 0; i < count; ++i)
  {
    key.data[4 * i + 0] = pairs[i][0];
    key.data[4 * i + 1] = pairs[i][1] - appliedOffset.x;
    key.data[4 * i + 2] = pairs[i][2] - appliedOffset.y;
    key.data[4 * i + 3] = pairs[i][3] - appliedOffset.z;
    if (permutation != nullptr) (*permutation)[i] = static_cast<std::size_t>(pairs[i][4]);
  }
  return key;
}

// The edges leaving a vertex: which triples of its sites carry one, and in which direction.
//
// Wang et al. give the tangent to the trisector of three sites at a point v as the cross product of
// the differences of the unit vectors from v to the sites, each such difference being the normal of
// one bisector sheet at v. That fixes the line but not which of its two halves is in the diagram.
// Moving off v along the trisector keeps the three sites equidistant while every other site becomes
// nearer or farther; only where all of them become farther does the tangent sphere stay empty and the
// arc belong to the diagram. The clearance to site m has gradient u_m, the unit vector from m to v, so
// the gap between site m and the triple grows along a direction d at the rate (u_m - u_i)·d, and its
// sign is the test. A half on which no site encroaches is an edge.
//
// With four sites each triple omits one site, whose gap grows one way and shrinks the other, so exactly
// one half of each of the four trisectors survives and there are four branches. With more, a triple can
// lose both halves, to different sites, and only those triples whose tangency points span a facet of
// the hull of all of them keep one.
//
// A rate too small to have a reliable sign leaves the half undecided: the gap is then governed by
// curvature rather than slope, which this test cannot see. Rather than guess, such a branch is admitted
// and reported through `ambiguous`, so that a diagram resting on a guess cannot claim to be consistent.
std::vector<SKApolloniusBranch> vertexBranches(const double3& centre, const std::vector<double3>& siteCentres,
                                              std::size_t& ambiguous)
{
  std::size_t count = siteCentres.size();
  std::vector<double3> unitToSite(count);
  for (std::size_t s = 0; s < count; ++s)
  {
    double3 delta = centre - siteCentres[s];
    double length = delta.length();
    unitToSite[s] = (length > 1.0e-12) ? (1.0 / length) * delta : double3(0.0, 0.0, 0.0);
  }

  constexpr double decidableRate = 1.0e-9;
  std::vector<SKApolloniusBranch> branches;
  for (std::size_t i = 0; i < count; ++i)
    for (std::size_t j = i + 1; j < count; ++j)
      for (std::size_t k = j + 1; k < count; ++k)
      {
        double3 tangent = double3::cross(unitToSite[i] - unitToSite[j], unitToSite[j] - unitToSite[k]);
        double length = tangent.length();
        if (length < 1.0e-14) continue;  // the three bisector normals are dependent: no trisector here
        tangent = (1.0 / length) * tangent;

        for (double sign : {1.0, -1.0})
        {
          double3 direction = sign * tangent;
          bool outgoing = true;
          bool undecided = false;
          for (std::size_t m = 0; m < count && outgoing; ++m)
          {
            if (m == i || m == j || m == k) continue;
            double rate = double3::dot(unitToSite[m] - unitToSite[i], direction);
            if (rate < -decidableRate)
              outgoing = false;
            else if (rate <= decidableRate)
              undecided = true;
          }
          if (!outgoing) continue;
          if (undecided) ++ambiguous;
          branches.push_back(SKApolloniusBranch{{i, j, k}, direction});
        }
      }
  return branches;
}

// Orientation of the tetrahedron formed by the points at which a tangent sphere touches its four
// sites, taken in canonical site order. This is the sign that distinguishes the two vertices a
// quadruple can support.
std::int64_t tangencyOrientation(const double3& centre, double radius, const std::array<double3, 4>& siteCentres,
                                 const std::array<std::size_t, 4>& permutation)
{
  std::array<double3, 4> tangency;
  for (std::size_t i = 0; i < 4; ++i)
  {
    double3 delta = siteCentres[permutation[i]] - centre;
    double length = delta.length();
    tangency[i] = (length > 1.0e-12) ? centre + (radius / length) * delta : centre;
  }
  double determinant = double3::dot(tangency[1] - tangency[0],
                                    double3::cross(tangency[2] - tangency[0], tangency[3] - tangency[0]));
  return (determinant >= 0.0) ? 1 : -1;
}


// A cell list over site images, used for the exhaustive emptiness certification. The certification
// must not depend on the candidate adjacency, or a spurious vertex could slip through.
struct SiteGrid
{
  double3x3 cell;
  double3x3 inverseCell;
  int3 binCount;
  double3 binWidth;
  std::vector<double3> positions;
  std::vector<std::vector<std::size_t>> bins;

  void build(const double3x3& unitCell, const double3x3& unitCellInverse, const double3& perpendicularWidths,
             const std::vector<double3>& fractionalPositions, double volume)
  {
    cell = unitCell;
    inverseCell = unitCellInverse;
    std::size_t count = fractionalPositions.size();
    double target = std::cbrt(volume / std::max(1.0, static_cast<double>(count) / 4.0));
    binCount = int3(std::max(1, static_cast<std::int32_t>(perpendicularWidths.x / target)),
                    std::max(1, static_cast<std::int32_t>(perpendicularWidths.y / target)),
                    std::max(1, static_cast<std::int32_t>(perpendicularWidths.z / target)));
    binWidth = double3(perpendicularWidths.x / binCount.x, perpendicularWidths.y / binCount.y,
                       perpendicularWidths.z / binCount.z);

    positions.resize(count);
    bins.assign(static_cast<std::size_t>(binCount.x) * static_cast<std::size_t>(binCount.y) *
                    static_cast<std::size_t>(binCount.z),
                {});
    for (std::size_t i = 0; i < count; ++i)
    {
      double3 fractional = double3::fract(fractionalPositions[i]);
      positions[i] = cell * fractional;
      std::int32_t bx = std::min(binCount.x - 1, static_cast<std::int32_t>(fractional.x * binCount.x));
      std::int32_t by = std::min(binCount.y - 1, static_cast<std::int32_t>(fractional.y * binCount.y));
      std::int32_t bz = std::min(binCount.z - 1, static_cast<std::int32_t>(fractional.z * binCount.z));
      bins[static_cast<std::size_t>((bz * binCount.y + by) * binCount.x + bx)].push_back(i);
    }
  }

  // Visits every site image whose centre lies within `searchRadius` of `wrappedCentre`, reporting
  // the site index, its image position, and the lattice image it came from.
  template <typename Visitor>
  void forEachNear(const double3& wrappedCentre, double searchRadius, Visitor&& visit) const
  {
    auto split = [](std::int32_t coordinate, std::int32_t extent)
    {
      std::int32_t image = (coordinate >= 0) ? coordinate / extent : -((-coordinate + extent - 1) / extent);
      return std::pair<std::int32_t, std::int32_t>{coordinate - image * extent, image};
    };

    double3 fractional = inverseCell * wrappedCentre;
    std::int32_t cx = std::min(binCount.x - 1, static_cast<std::int32_t>(fractional.x * binCount.x));
    std::int32_t cy = std::min(binCount.y - 1, static_cast<std::int32_t>(fractional.y * binCount.y));
    std::int32_t cz = std::min(binCount.z - 1, static_cast<std::int32_t>(fractional.z * binCount.z));

    int3 span(static_cast<std::int32_t>(std::ceil(searchRadius / binWidth.x)) + 1,
              static_cast<std::int32_t>(std::ceil(searchRadius / binWidth.y)) + 1,
              static_cast<std::int32_t>(std::ceil(searchRadius / binWidth.z)) + 1);

    for (std::int32_t oz = -span.z; oz <= span.z; ++oz)
      for (std::int32_t oy = -span.y; oy <= span.y; ++oy)
        for (std::int32_t ox = -span.x; ox <= span.x; ++ox)
        {
          auto [bx, lx] = split(cx + ox, binCount.x);
          auto [by, ly] = split(cy + oy, binCount.y);
          auto [bz, lz] = split(cz + oz, binCount.z);
          int3 image(lx, ly, lz);
          double3 shift = cell * double3(lx, ly, lz);
          for (std::size_t j : bins[static_cast<std::size_t>((bz * binCount.y + by) * binCount.x + bx)])
            visit(j, positions[j] + shift, image);
        }
  }

  // Whether any site image overlaps the sphere of radius `radius` about `wrappedCentre`, which is to say
  // whether any site comes nearer to it than r_j + radius. This is the certification a candidate vertex
  // has to pass, and nearly every candidate proposed fails it, so what the test costs is decided by how
  // quickly a failure is found. The bins are therefore taken in shells outwards from the one the centre
  // is in, and the search stops at the first site that overlaps: an intruder is usually close by, being a
  // site that was nearly a site of the vertex, and is met in the first shell or two. A sphere that is
  // genuinely empty costs the full search, as it must.
  bool intrudes(const double3& wrappedCentre, double radius, const std::vector<double>& radii, double maximumRadius,
                double tolerance) const
  {
    auto split = [](std::int32_t coordinate, std::int32_t extent)
    {
      std::int32_t image = (coordinate >= 0) ? coordinate / extent : -((-coordinate + extent - 1) / extent);
      return std::pair<std::int32_t, std::int32_t>{coordinate - image * extent, image};
    };

    double searchRadius = std::max(0.0, radius + maximumRadius);
    double3 fractional = inverseCell * wrappedCentre;
    std::int32_t cx = std::min(binCount.x - 1, static_cast<std::int32_t>(fractional.x * binCount.x));
    std::int32_t cy = std::min(binCount.y - 1, static_cast<std::int32_t>(fractional.y * binCount.y));
    std::int32_t cz = std::min(binCount.z - 1, static_cast<std::int32_t>(fractional.z * binCount.z));

    int3 span(static_cast<std::int32_t>(std::ceil(searchRadius / binWidth.x)) + 1,
              static_cast<std::int32_t>(std::ceil(searchRadius / binWidth.y)) + 1,
              static_cast<std::int32_t>(std::ceil(searchRadius / binWidth.z)) + 1);
    std::int32_t widest = std::max({span.x, span.y, span.z});

    for (std::int32_t shell = 0; shell <= widest; ++shell)
      for (std::int32_t oz = -std::min(shell, span.z); oz <= std::min(shell, span.z); ++oz)
        for (std::int32_t oy = -std::min(shell, span.y); oy <= std::min(shell, span.y); ++oy)
          for (std::int32_t ox = -std::min(shell, span.x); ox <= std::min(shell, span.x); ++ox)
          {
            if (std::max({std::abs(ox), std::abs(oy), std::abs(oz)}) != shell) continue;

            auto [bx, lx] = split(cx + ox, binCount.x);
            auto [by, ly] = split(cy + oy, binCount.y);
            auto [bz, lz] = split(cz + oz, binCount.z);
            double3 shift = cell * double3(lx, ly, lz);
            for (std::size_t j : bins[static_cast<std::size_t>((bz * binCount.y + by) * binCount.x + bx)])
              if ((wrappedCentre - (positions[j] + shift)).length() < radii[j] + radius - tolerance) return true;
          }
    return false;
  }

  // Clearance min_j(|x - x_j| - r_j) at a point, together with the site image that attains it.
  // `upperBound` must be an upper bound on the answer so the search radius can be sized exactly.
  double clearance(const double3& wrappedPoint, const std::vector<double>& radii, double upperBound,
                   double maximumRadius, std::size_t& nearestIndex, int3& nearestImage) const
  {
    double best = upperBound;
    nearestIndex = std::numeric_limits<std::size_t>::max();
    forEachNear(wrappedPoint, upperBound + maximumRadius,
                [&](std::size_t j, const double3& image, const int3& latticeImage)
                {
                  double value = (wrappedPoint - image).length() - radii[j];
                  if (value < best)
                  {
                    best = value;
                    nearestIndex = j;
                    nearestImage = latticeImage;
                  }
                });
    return best;
  }
};

}  // namespace

namespace
{
// The spheres tangent to four sites, as `skApolloniusTangentSpheres` gives them but without allocating,
// which is what the searches over quadruples use.
TangentSpheres tangentSpheresOf(const std::array<double3, 4>& centres, const std::array<double, 4>& radii,
                                bool allowNegativeRadius)
{
  std::array<double3, 3> rows;
  double3 timeCoefficients(0.0, 0.0, 0.0);
  double3 constants(0.0, 0.0, 0.0);
  double timeArray[3];
  double constantArray[3];
  for (std::size_t i = 0; i < 3; ++i)
    tangencyRow(centres[i + 1], radii[i + 1], centres[0], radii[0], rows[i], timeArray[i], constantArray[i]);

  timeCoefficients = double3(timeArray[0], timeArray[1], timeArray[2]);
  constants = double3(constantArray[0], constantArray[1], constantArray[2]);
  TangentSpheres solutions =
      solveTangencyRows(rows, timeCoefficients, constants, centres[0], radii[0], allowNegativeRadius);

  // solveTangencyRows can only impose the unsquared tangency for the site it was given as reference.
  // The other three were folded into differences of squares, so each needs the same test, and each
  // needs its distance to come out right rather than merely right in magnitude.
  TangentSpheres kept;
  for (const SKApolloniusTangentSphere& sphere : solutions)
  {
    bool tangent = true;
    for (std::size_t i = 0; i < 4 && tangent; ++i)
    {
      double tangentDistance = radii[i] + sphere.radius;
      if (tangentDistance < 0.0 || std::abs((sphere.centre - centres[i]).length() - tangentDistance) >
                                       1.0e-6 * std::max(1.0, tangentDistance))
        tangent = false;
    }
    if (tangent) kept.add(sphere);
  }
  return kept;
}

}  // namespace

std::vector<SKApolloniusTangentSphere> skApolloniusTangentSpheres(const std::array<double3, 4>& centres,
                                                                 const std::array<double, 4>& radii,
                                                                 bool allowNegativeRadius)
{
  TangentSpheres spheres = tangentSpheresOf(centres, radii, allowNegativeRadius);
  return std::vector<SKApolloniusTangentSphere>(spheres.begin(), spheres.end());
}

std::optional<SKApolloniusTangentSphere> skApolloniusTrisectorPoint(const std::array<double3, 3>& centres,
                                                                   const std::array<double, 3>& radii,
                                                                   const double3& from, const double3& to, double t,
                                                                   bool allowNegativeRadius)
{
  // Two rows from tangency to sites 1 and 2 relative to site 0, and a third from the plane cutting
  // the chord at parameter t, which selects one point on the curve.
  double3 chord = to - from;
  double chordLength = chord.length();
  if (chordLength < 1.0e-12) return std::nullopt;
  double3 normal = chord / chordLength;
  double3 planePoint = from + t * chord;

  std::array<double3, 3> rows;
  double timeArray[3];
  double constantArray[3];
  for (std::size_t i = 0; i < 2; ++i)
    tangencyRow(centres[i + 1], radii[i + 1], centres[0], radii[0], rows[i], timeArray[i], constantArray[i]);
  rows[2] = normal;
  timeArray[2] = 0.0;
  constantArray[2] = double3::dot(normal, planePoint);

  TangentSpheres solutions =
      solveTangencyRows(rows, double3(timeArray[0], timeArray[1], timeArray[2]),
                        double3(constantArray[0], constantArray[1], constantArray[2]), centres[0], radii[0],
                        allowNegativeRadius);

  // The rows carry the tangency to sites 1 and 2 squared, so a root of them need not satisfy the
  // unsquared condition. With the clearance confined to be non-negative it always does; once the curve
  // is followed inside the sites the wrong branch of each square root becomes a solution too, and has
  // to be thrown out.
  //
  // The plane meets the curve in up to two points, on opposite branches; the one nearest the chord is
  // kept, which is the branch the arc between the two vertices runs along.
  std::optional<SKApolloniusTangentSphere> best;
  double bestDistance = std::numeric_limits<double>::max();
  for (const SKApolloniusTangentSphere& sphere : solutions)
  {
    bool tangent = true;
    for (std::size_t i = 0; i < 3 && tangent; ++i)
    {
      double tangentDistance = radii[i] + sphere.radius;
      if (tangentDistance < 0.0 || std::abs((sphere.centre - centres[i]).length() - tangentDistance) >
                                       1.0e-6 * std::max(1.0, tangentDistance))
        tangent = false;
    }
    if (!tangent) continue;

    double distance = (sphere.centre - planePoint).length();
    if (distance < bestDistance)
    {
      bestDistance = distance;
      best = sphere;
    }
  }
  return best;
}

double SKApolloniusDiagram::largestEmptySphereRadius() const
{
  double maximum = 0.0;
  for (const SKApolloniusVertex& vertex : vertices) maximum = std::max(maximum, vertex.radius);
  return maximum;
}

double3 SKApolloniusDiagram::largestEmptySpherePosition() const
{
  double maximum = -1.0;
  double3 position(0.0, 0.0, 0.0);
  for (const SKApolloniusVertex& vertex : vertices)
    if (vertex.radius > maximum)
    {
      maximum = vertex.radius;
      position = vertex.position;
    }
  return position;
}

SKApolloniusDiagram SKApolloniusDiagram::create(const double3x3& unitCell,
                                                const std::vector<double3>& fractionalPositions,
                                                const std::vector<double>& radii, std::size_t neighbourRings,
                                                SKApolloniusRegion region)
{
  // Over all of space a vertex of negative clearance is admitted, and a site intrudes on a candidate
  // only if its own clearance there is smaller, which is the same test read with a signed radius.
  bool allowNegativeRadius = (region == SKApolloniusRegion::EntireSpace);

  SKApolloniusDiagram diagram;
  std::size_t siteCount = fractionalPositions.size();
  // One site in a periodic cell is already an infinite set of spheres, and the four a vertex needs may
  // be four images of it, so what has to be non-empty is the cell and not the count of sites in it.
  if (siteCount == 0 || radii.size() != siteCount) return diagram;

  double3x3 inverseCell = unitCell.inverse();
  double volume = std::abs(unitCell.determinant());

  // Perpendicular widths of the cell, needed to size the bins. The columns of the cell matrix are the
  // lattice vectors, so they are what the widths are built from; taking the rows instead would give the
  // widths of a different cell whenever the cell is not orthogonal.
  double3 a = unitCell[0];
  double3 b = unitCell[1];
  double3 c = unitCell[2];
  double3 perpendicularWidths(volume / double3::cross(b, c).length(), volume / double3::cross(c, a).length(),
                              volume / double3::cross(a, b).length());

  double maximumRadius = 0.0;
  for (double radius : radii) maximumRadius = std::max(maximumRadius, radius);

  // An upper bound on the radius of any empty sphere. Four sites in near-degenerate position have a
  // common tangent sphere that is astronomically large, and such a sphere is never a vertex here
  // because periodicity guarantees a site image within the lattice covering radius of every point.
  // Half the sum of the lattice vector lengths bounds all four half-diagonals and so bounds the
  // covering radius. Discarding candidates above it keeps the emptiness search radius finite, which
  // it otherwise would not be.
  double emptyRadiusBound = 0.5 * (a.length() + b.length() + c.length());

  SiteGrid grid;
  grid.build(unitCell, inverseCell, perpendicularWidths, fractionalPositions, volume);

  // Candidate adjacency from the radical diagram. It coincides with Apollonius adjacency for equal
  // radii and stays close as they spread, which makes it a strong seed; correctness of the result
  // does not rest on it, only the chance of finding every vertex does.
  SKVoronoi radical(unitCell, fractionalPositions, radii);
  std::vector<std::vector<std::pair<std::size_t, int3>>> adjacency(siteCount);
  {
    std::vector<SKVoronoiCell> radicalCells = radical.computeAllCells();
    for (std::size_t i = 0; i < siteCount; ++i)
      for (const SKVoronoiFace& face : radicalCells[i].faces)
        adjacency[i].push_back({face.neighborIndex, face.neighborImage});
  }

  // Widening by extra rings composes the image shifts, so a neighbour of a neighbour is recorded at
  // the image where it actually sits relative to the original site.
  for (std::size_t ring = 1; ring < std::max<std::size_t>(1, neighbourRings); ++ring)
  {
    std::vector<std::vector<std::pair<std::size_t, int3>>> widened = adjacency;
    for (std::size_t i = 0; i < siteCount; ++i)
    {
      for (const auto& [j, imageJ] : adjacency[i])
        for (const auto& [k, imageK] : adjacency[j])
        {
          int3 image(imageJ.x + imageK.x, imageJ.y + imageK.y, imageJ.z + imageK.z);
          if (k == i && image.x == 0 && image.y == 0 && image.z == 0) continue;
          widened[i].push_back({k, image});
        }
      std::sort(widened[i].begin(), widened[i].end(),
                [](const auto& lhs, const auto& rhs)
                {
                  return std::tie(lhs.first, lhs.second.x, lhs.second.y, lhs.second.z) <
                         std::tie(rhs.first, rhs.second.x, rhs.second.y, rhs.second.z);
                });
      widened[i].erase(std::unique(widened[i].begin(), widened[i].end(),
                                   [](const auto& lhs, const auto& rhs)
                                   {
                                     return lhs.first == rhs.first && lhs.second.x == rhs.second.x &&
                                            lhs.second.y == rhs.second.y && lhs.second.z == rhs.second.z;
                                   }),
                       widened[i].end());
    }
    adjacency = std::move(widened);
  }

  // Enumerate quadruples: each site together with three of its candidate neighbours. Every vertex
  // is reachable from all four of its sites, so it will be proposed several times and deduplicated
  // by its canonical site tuple.
  std::map<SiteTuple, std::size_t> vertexLookup;
  constexpr double emptinessTolerance = 1.0e-7;

  // Quadruples already solved, keyed by their canonical form, from which the lattice translation has been
  // removed. A quadruple is reached from each of its four sites, and in the sweep below from every box that
  // its vertex is near, so most of the ones proposed have been solved already. What a quadruple gives is
  // fixed by the shape it makes: translating it translates its tangent spheres, which are wrapped back into
  // the cell and canonicalised to the same keys, and leaves the emptiness of those spheres alone. So a
  // quadruple seen a second time can add nothing, and both the tangency and the certification of it are
  // spared.
  std::unordered_set<SiteTuple, SiteTupleHash> solvedQuadruples;
  auto quadrupleIsNew = [&solvedQuadruples](const std::size_t* indices, const int3* images)
  {
    int3 unusedOffset;
    return solvedQuadruples.insert(canonicalTuple(indices, images, 4, unusedOffset)).second;
  };

  for (std::size_t i = 0; i < siteCount; ++i)
  {
    const std::vector<std::pair<std::size_t, int3>>& neighbours = adjacency[i];
    std::size_t count = neighbours.size();
    double3 centreI = grid.positions[i];

    for (std::size_t p = 0; p < count; ++p)
      for (std::size_t q = p + 1; q < count; ++q)
        for (std::size_t r = q + 1; r < count; ++r)
        {
          std::array<std::size_t, 4> siteIndices{i, neighbours[p].first, neighbours[q].first, neighbours[r].first};
          std::array<int3, 4> siteImages{int3(0, 0, 0), neighbours[p].second, neighbours[q].second,
                                         neighbours[r].second};

          if (!quadrupleIsNew(siteIndices.data(), siteImages.data())) continue;

          std::array<double3, 4> centres;
          std::array<double, 4> tangentRadii;
          centres[0] = centreI;
          tangentRadii[0] = radii[i];
          for (std::size_t s = 1; s < 4; ++s)
          {
            centres[s] = grid.positions[siteIndices[s]] +
                         unitCell * double3(siteImages[s].x, siteImages[s].y, siteImages[s].z);
            tangentRadii[s] = radii[siteIndices[s]];
          }

          for (const SKApolloniusTangentSphere& sphere : tangentSpheresOf(centres, tangentRadii, allowNegativeRadius))
          {
            if (sphere.radius > emptyRadiusBound) continue;

            // Certify emptiness exhaustively. Any site closer than r_j + t would overlap the
            // tangent sphere, which would mean this is not a vertex of the diagram. The search
            // radius is exact: a site can only intrude if its centre is within t + r_max.
            double3 wrappedCentre = unitCell * double3::fract(inverseCell * sphere.centre);
            double3 latticeShift = wrappedCentre - sphere.centre;
            if (grid.intrudes(wrappedCentre, sphere.radius, radii, maximumRadius, emptinessTolerance)) continue;

            // Re-express the four sites relative to the wrapped centre, so that stored images are
            // consistent with the stored position and every periodic copy canonicalises alike.
            double3 shiftFractional = inverseCell * latticeShift;
            int3 shift(static_cast<std::int32_t>(std::llround(shiftFractional.x)),
                       static_cast<std::int32_t>(std::llround(shiftFractional.y)),
                       static_cast<std::int32_t>(std::llround(shiftFractional.z)));
            std::array<int3, 4> shiftedImages;
            for (std::size_t s = 0; s < 4; ++s)
              shiftedImages[s] = int3(siteImages[s].x + shift.x, siteImages[s].y + shift.y, siteImages[s].z + shift.z);

            int3 unusedOffset;
            std::array<std::size_t, 4> permutation{};
            SiteTuple key = canonicalTuple(siteIndices.data(), shiftedImages.data(), 4, unusedOffset, &permutation);
            std::array<double3, 4> shiftedCentres;
            for (std::size_t s = 0; s < 4; ++s) shiftedCentres[s] = centres[s] + latticeShift;
            key.orientation = tangencyOrientation(wrappedCentre, sphere.radius, shiftedCentres, permutation);
            if (vertexLookup.contains(key)) continue;

            // Branches are left until the cotangent sets are known, since a vertex found here as a
            // quadruple may turn out to be part of a larger set and its edges depend on the whole of it.
            vertexLookup.emplace(key, diagram.vertices.size());
            diagram.vertices.push_back(
                SKApolloniusVertex{wrappedCentre, sphere.radius,
                                   std::vector<std::size_t>(siteIndices.begin(), siteIndices.end()),
                                   std::vector<int3>(shiftedImages.begin(), shiftedImages.end()),
                                   {}});
          }
        }
  }

  // Each branch of a vertex runs along the trisector of three of its sites, a one-dimensional locus
  // whose two ends are vertices. So each triple should be met by exactly two branches, and that pairing
  // is the edge. The branch is carried alongside the vertex, since it holds the direction the arc leaves
  // in, and since a vertex may meet one triple with two branches when the trisector passes through it.
  struct IncidentBranch
  {
    std::size_t vertex;
    std::size_t branch;
    int3 offset;  // lattice shift canonicalisation applied to reach the key's frame
  };
  auto buildTripleLookup = [&]()
  {
    std::map<SiteTuple, std::vector<IncidentBranch>> lookup;
    for (std::size_t v = 0; v < diagram.vertices.size(); ++v)
    {
      const SKApolloniusVertex& vertex = diagram.vertices[v];
      for (std::size_t b = 0; b < vertex.branches.size(); ++b)
      {
        std::array<std::size_t, 3> indices;
        std::array<int3, 3> images;
        for (std::size_t s = 0; s < 3; ++s)
        {
          indices[s] = vertex.siteIndices[vertex.branches[b].sites[s]];
          images[s] = vertex.siteImages[vertex.branches[b].sites[s]];
        }
        int3 offset;
        SiteTuple key = canonicalTuple(indices.data(), images.data(), 3, offset);
        lookup[key].push_back(IncidentBranch{v, b, offset});
      }
    }
    return lookup;
  };

  // Completeness pass, following the vertex location of Wang et al., Section 5.2. Enumerating
  // candidate quadruples from radical adjacency is a heuristic and does miss vertices. Rather than
  // widen it, which converges slowly and costs cubically, the cell is swept: a box is subdivided
  // until few enough sites can be nearest anywhere inside it to take those sites in quadruples, and
  // every quadruple of them is solved. Since every vertex has four nearest sites, and the sites that
  // can be nearest in a box containing it include them, the quadruple that gives it is among those
  // solved, so no vertex can escape the sweep whenever the box holding it is finished. Boxes near a
  // face or an edge of the diagram come down to two or three sites and cost nothing, which is what
  // keeps the refinement local.
  //
  // The sweep also collects the triples it passes, which is how the trisectors that carry no vertex are
  // later found: they cannot be reached from the vertices, so they need a source of candidates that
  // does not depend on vertices existing.
  std::unordered_set<SiteTuple, SiteTupleHash> candidateTriples;
  {
    // Wang et al. describe a cube carrying the sites that may dominate part of it, and reach that set by
    // a Dijkstra-like propagation over cubes through a queue of distance events. What matters about the
    // propagation is that a cube is told which sites to reckon with by its parent rather than by a fresh
    // search: a box lies inside its parent, so the sites that concern it are among those that concerned
    // the parent, and the whole sweep then costs one distance per site per box instead of a range query
    // per box. The queue itself is unnecessary, since a box knows its own parent.
    //
    // For the parent's list to answer for the child it has to be drawn a little wider than the test the
    // child applies to it, and the margin that closes is four times the circumradius. Writing d_i for the
    // clearance to site i and S(c, r, a) for {i : d_i(c) <= min_j d_j(c) + a r}, the clearance is
    // 1-Lipschitz, so for a child of centre c' and circumradius r/2 inside a parent of centre c and
    // circumradius r,
    //
    //     i in S(c', r/2, 4)  =>  d_i(c) <= d_i(c') + r <= min_j d_j(c') + 2r + r <= min_j d_j(c) + 4r,
    //
    // that is S(c', r/2, 4) is contained in S(c, r, 4). The margin four is the fixed point of that
    // recursion and so is the one that lets the lists nest all the way down. Theorem 6 needs only
    // S(c, r, 2) to decide the box, which is read off the wider list at no cost.
    struct CandidateSite
    {
      std::size_t index;
      int3 image;
      double3 position;  // of that image, so descending the tree costs one distance and no lattice work
    };

    // The list belongs to the level rather than to the box. All eight children of a box reckon with the
    // same sites, and the boxes are taken depth first, so while a box of some level is being processed no
    // other box of that level has been reached and the one list per level is the only one any of them
    // needs. That is what keeps the sweep free of allocation: a box costs one distance per site on its
    // list and nothing else.
    // A box is its corner and its level, since every box of a level has the same size: the size and the
    // circumradius that go with it are worked out once per level rather than once per box. It also carries
    // how many sites could be nearest in its parent, which is what tells it whether splitting is achieving
    // anything.
    struct Cube
    {
      double3 originFractional;
      std::size_t level;
      std::size_t parentNearestCount;
    };

    // Half the longest diagonal, the radius of the sphere about the centre that encloses the box.
    auto circumradius = [&](const double3& sizeFractional)
    {
      double largest = 0.0;
      for (double sx : {-0.5, 0.5})
        for (double sy : {-0.5, 0.5})
          for (double sz : {-0.5, 0.5})
            largest = std::max(largest, (unitCell * double3(sx * sizeFractional.x, sy * sizeFractional.y,
                                                           sz * sizeFractional.z))
                                            .length());
      return largest;
    };

    // The list a box of the coarsest level starts from, which is the only one that has to come out of the
    // grid. The query has to reach every site within 4r of the nearest one, and how far that is is not
    // known until the nearest is, so it widens until it does.
    auto gatherFromGrid = [&](const double3& centre, double radius, std::vector<CandidateSite>& out)
    {
      out.clear();
      double searchRadius = maximumRadius + 4.0 * radius + 1.0;
      double smallest = 0.0;
      bool found = false;
      for (std::size_t attempt = 0; attempt < 8; ++attempt)
      {
        found = false;
        smallest = std::numeric_limits<double>::max();
        grid.forEachNear(centre, searchRadius,
                         [&](std::size_t j, const double3& image, const int3&)
                         {
                           found = true;
                           smallest = std::min(smallest, (centre - image).length() - radii[j]);
                         });
        if (found && searchRadius >= smallest + 4.0 * radius + maximumRadius) break;
        searchRadius = found ? (smallest + 4.0 * radius + maximumRadius + 1.0e-6) : (2.0 * searchRadius);
      }
      if (!found) return;

      double threshold = smallest + 4.0 * radius;
      grid.forEachNear(centre, searchRadius,
                       [&](std::size_t j, const double3& image, const int3& imageOffset)
                       {
                         if ((centre - image).length() - radii[j] <= threshold)
                           out.push_back(CandidateSite{j, imageOffset, image});
                       });
    };

    // `nearby` is the list of the box the quadruple came from, which is at hand and in cache. Nearly every
    // quadruple proposed meets at no vertex, and the site that spoils it is usually one of these, so trying
    // them before the grid is what most candidates cost in total. It is only a shortcut: a candidate that
    // survives it is still certified against everything.
    auto certifyAndStore = [&](const std::array<std::size_t, 4>& quadIndices, const std::array<int3, 4>& quadImages,
                               const std::vector<CandidateSite>& nearby)
    {
      if (!quadrupleIsNew(quadIndices.data(), quadImages.data())) return;

      std::array<double3, 4> quadCentres;
      std::array<double, 4> quadRadii;
      for (std::size_t s = 0; s < 4; ++s)
      {
        quadCentres[s] =
            grid.positions[quadIndices[s]] + unitCell * double3(quadImages[s].x, quadImages[s].y, quadImages[s].z);
        quadRadii[s] = radii[quadIndices[s]];
      }

      for (const SKApolloniusTangentSphere& sphere : tangentSpheresOf(quadCentres, quadRadii, allowNegativeRadius))
      {
        if (sphere.radius > emptyRadiusBound) continue;

        bool spoiled = false;
        for (std::size_t s = 0; s < nearby.size() && !spoiled; ++s)
          spoiled = (sphere.centre - nearby[s].position).length() <
                    radii[nearby[s].index] + sphere.radius - emptinessTolerance;
        if (spoiled) continue;

        double3 wrappedCentre = unitCell * double3::fract(inverseCell * sphere.centre);
        double3 latticeShift = wrappedCentre - sphere.centre;
        if (grid.intrudes(wrappedCentre, sphere.radius, radii, maximumRadius, emptinessTolerance)) continue;

        double3 shiftFractional = inverseCell * latticeShift;
        int3 shift(static_cast<std::int32_t>(std::llround(shiftFractional.x)),
                   static_cast<std::int32_t>(std::llround(shiftFractional.y)),
                   static_cast<std::int32_t>(std::llround(shiftFractional.z)));
        std::array<int3, 4> shiftedImages;
        for (std::size_t s = 0; s < 4; ++s)
          shiftedImages[s] = int3(quadImages[s].x + shift.x, quadImages[s].y + shift.y, quadImages[s].z + shift.z);

        int3 unusedOffset;
        std::array<std::size_t, 4> permutation{};
        SiteTuple newKey = canonicalTuple(quadIndices.data(), shiftedImages.data(), 4, unusedOffset, &permutation);
        std::array<double3, 4> shiftedCentres;
        for (std::size_t s = 0; s < 4; ++s) shiftedCentres[s] = quadCentres[s] + latticeShift;
        newKey.orientation = tangencyOrientation(wrappedCentre, sphere.radius, shiftedCentres, permutation);
        if (vertexLookup.contains(newKey)) continue;

        vertexLookup.emplace(newKey, diagram.vertices.size());
        diagram.vertices.push_back(
            SKApolloniusVertex{wrappedCentre, sphere.radius,
                               std::vector<std::size_t>(quadIndices.begin(), quadIndices.end()),
                               std::vector<int3>(shiftedImages.begin(), shiftedImages.end()),
                               {}});
      }
    };

    // A floor on refinement, at a box far smaller than any site separation, so that the recursion has a
    // bound that does not depend on the sites being in any particular position. The rules below finish a
    // box long before it, and on the structures this has been run on it is never reached; it is here for
    // the configuration that would defeat them.
    constexpr double smallestBoxCircumradius = 1.0e-3;

    // The coarsest boxes are sized by the geometry and not as a fraction of the cell, so that the list a
    // box has to reckon with is set by how the sites are packed rather than by how large the cell is.
    // Dividing the cell a fixed number of times would grow that list with the cell, which is the one way
    // this sweep could fail to be linear in the number of sites.
    constexpr double initialBoxEdge = 3.0;
    int3 initialDivisions(std::max(1, static_cast<std::int32_t>(std::llround(perpendicularWidths.x / initialBoxEdge))),
                          std::max(1, static_cast<std::int32_t>(std::llround(perpendicularWidths.y / initialBoxEdge))),
                          std::max(1, static_cast<std::int32_t>(std::llround(perpendicularWidths.z / initialBoxEdge))));
    double3 initialSize(1.0 / initialDivisions.x, 1.0 / initialDivisions.y, 1.0 / initialDivisions.z);

    std::vector<Cube> stack;
    std::vector<std::vector<CandidateSite>> levelReachable(1);
    std::vector<double3> levelSize{initialSize};
    std::vector<double> levelRadius{circumradius(initialSize)};
    auto reachLevel = [&](std::size_t level)
    {
      while (levelSize.size() <= level)
      {
        levelSize.push_back(0.5 * levelSize.back());
        levelRadius.push_back(0.5 * levelRadius.back());
      }
      if (levelReachable.size() <= level) levelReachable.resize(level + 1);
    };

    for (std::int32_t ix = 0; ix < initialDivisions.x; ++ix)
      for (std::int32_t iy = 0; iy < initialDivisions.y; ++iy)
        for (std::int32_t iz = 0; iz < initialDivisions.z; ++iz)
          stack.push_back(Cube{double3(static_cast<double>(ix) * initialSize.x, static_cast<double>(iy) * initialSize.y,
                                       static_cast<double>(iz) * initialSize.z),
                               0, std::numeric_limits<std::size_t>::max()});

    std::vector<double> clearances;
    std::vector<std::pair<std::size_t, int3>> nearestSites;  // S(c, r, 2): can be nearest inside the box
    std::vector<std::int64_t> finishedKey;
    std::unordered_set<std::vector<std::int64_t>, SiteWordsHash> finishedSets;

    while (!stack.empty())
    {
      Cube cube = stack.back();
      stack.pop_back();

      // Room for this box's own list and for the one it may write for its children, taken before either is
      // referred to so that neither reference is left dangling by the growth.
      reachLevel(cube.level + 1);

      double3 centre = unitCell * (cube.originFractional + 0.5 * levelSize[cube.level]);
      double radius = levelRadius[cube.level];

      // A box of the coarsest level is the one that has no parent to inherit from.
      if (cube.level == 0) gatherFromGrid(centre, radius, levelReachable[0]);
      const std::vector<CandidateSite>& reachable = levelReachable[cube.level];

      clearances.resize(reachable.size());
      double smallest = std::numeric_limits<double>::max();
      for (std::size_t s = 0; s < reachable.size(); ++s)
      {
        clearances[s] = (centre - reachable[s].position).length() - radii[reachable[s].index];
        smallest = std::min(smallest, clearances[s]);
      }
      if (reachable.empty()) continue;

      // Over the free region alone, a box whose clearance is negative throughout holds no vertex and no
      // point of any trisector the diagram keeps, since the clearance cannot rise by more than the
      // circumradius anywhere in the box. Such a box lies inside the sites and is not explored.
      if (!allowNegativeRadius && smallest + radius < 0.0) continue;

      nearestSites.clear();
      for (std::size_t s = 0; s < reachable.size(); ++s)
        if (clearances[s] < smallest + 2.0 * radius)
          nearestSites.push_back({reachable[s].index, reachable[s].image});

      // Where a box is down to a few possible nearest sites there is nothing left to narrow: whatever the
      // diagram does inside it, it does with those sites and no others. So the box is finished here rather
      // than refined, by taking its sites in threes and in fours.
      //
      // Every quadruple is solved, and each vertex inside the box is tangent to four sites that are nearest
      // there, so it is tangent to four of these and is among the solutions; the certification throws out
      // the quadruples that meet at no vertex of the diagram, as it does everywhere else. Every triple is
      // recorded, and a trisector inside the box likewise runs on three of these, which is what the search
      // for the closed ones later needs. And since a child box lies inside this one, the same holds of
      // anything in any of them: refining could turn up nothing that solving these does not.
      //
      // Since that argument holds of a box of any size, a box may be finished whenever finishing it is
      // cheaper than splitting it, and when to do so is a question of cost alone and not of correctness.
      // Two things answer it. One is that the box is down to few sites, six being where the fifteen
      // quadruples of six cost about what the level of refinement they replace costs. The other is that
      // splitting the parent did not reduce the count, so splitting this box is unlikely to either, and the
      // quadruples are better paid for once here than in eight children and again in theirs. This second is
      // what a degeneracy looks like from inside the sweep, and having it is what lets the sweep leave
      // degeneracies alone rather than refine at them: sites exactly cotangent stay cotangent however far a
      // box around them is split, and Wang et al.'s condition of refining until four sites remain is
      // precisely the one that would never be met. It is held to sets whose quadruples are few enough to be
      // worth it, a coarse box holding hundreds of sites with no degeneracy anywhere near it.
      constexpr std::size_t fewSites = 6;
      constexpr std::size_t affordableSites = 16;
      bool splittingIsHelping = nearestSites.size() < cube.parentNearestCount;
      if (nearestSites.size() <= fewSites || radius < smallestBoxCircumradius ||
          (!splittingIsHelping && nearestSites.size() <= affordableSites))
      {
        if (nearestSites.size() > affordableSites) ++diagram.verification.degenerateSweepBoxes;

        const std::vector<std::pair<std::size_t, int3>>& sites = nearestSites;
        if (sites.size() >= 3)
        {
          // Boxes far outnumber the sets of sites they finish with: the boxes along one edge of the diagram
          // are finished by the same three sites, those around one vertex by the same four, and a box in one
          // corner of the cell by the same set as its periodic twin in another. So the set is asked about
          // before it is worked with, and only a set not seen before is taken in threes and fours. What the
          // set is asked in is its canonical form, sorted and with the lattice translation removed, which is
          // what makes the twins the same set.
          std::sort(nearestSites.begin(), nearestSites.end(),
                    [](const auto& lhs, const auto& rhs)
                    {
                      return std::tie(lhs.first, lhs.second.x, lhs.second.y, lhs.second.z) <
                             std::tie(rhs.first, rhs.second.x, rhs.second.y, rhs.second.z);
                    });
          int3 base = nearestSites.front().second;
          finishedKey.clear();
          for (const auto& [index, image] : nearestSites)
            finishedKey.insert(finishedKey.end(), {static_cast<std::int64_t>(index), image.x - base.x,
                                                  image.y - base.y, image.z - base.z});
          if (!finishedSets.insert(finishedKey).second) continue;

          for (std::size_t i = 0; i < sites.size(); ++i)
            for (std::size_t j = i + 1; j < sites.size(); ++j)
              for (std::size_t k = j + 1; k < sites.size(); ++k)
              {
                std::array<std::size_t, 3> tripleIndices{sites[i].first, sites[j].first, sites[k].first};
                std::array<int3, 3> tripleImages{sites[i].second, sites[j].second, sites[k].second};
                int3 tripleOffset;
                candidateTriples.insert(canonicalTuple(tripleIndices.data(), tripleImages.data(), 3, tripleOffset));

                for (std::size_t l = k + 1; l < sites.size(); ++l)
                  certifyAndStore({sites[i].first, sites[j].first, sites[k].first, sites[l].first},
                                  {sites[i].second, sites[j].second, sites[k].second, sites[l].second}, reachable);
              }
        }
        continue;
      }

      // The sites the children reckon with, which are those on this box's list that they can still be
      // concerned by, written to the level below for all eight of them to read.
      std::vector<CandidateSite>& childReachable = levelReachable[cube.level + 1];
      childReachable.clear();
      for (std::size_t s = 0; s < reachable.size(); ++s)
        if (clearances[s] <= smallest + 4.0 * radius) childReachable.push_back(reachable[s]);

      const double3& half = levelSize[cube.level + 1];
      for (std::size_t ix = 0; ix < 2; ++ix)
        for (std::size_t iy = 0; iy < 2; ++iy)
          for (std::size_t iz = 0; iz < 2; ++iz)
            stack.push_back(Cube{cube.originFractional + double3(static_cast<double>(ix) * half.x,
                                                                 static_cast<double>(iy) * half.y,
                                                                 static_cast<double>(iz) * half.z),
                                 cube.level + 1, nearestSites.size()});
    }
  }

  // Cotangent sets, and the branches that follow from them.
  //
  // A sphere touching more than four sites is one vertex of the diagram belonging to all of them, but its
  // tangency is satisfied by each quadruple drawn from the set, so the two passes above have constructed
  // it once per quadruple. Those copies sit at the same point and are gathered here into the single vertex
  // they are, carrying the union of their sites. Nothing else can coincide: the radius of a certified
  // tangent sphere is the clearance at its centre, so its centre determines it, and two vertices at one
  // point can only be one vertex found twice.
  //
  // Coincidence is judged on the minimum image, since a vertex on a cell face may be wrapped to either
  // side of it. Which copies to compare comes from a grid of the tolerance's own size laid over the cell:
  // a copy can only coincide with one in its own cell of that grid or in a neighbouring one, and the cells
  // wrap, so a vertex on a face of the unit cell meets the copies of it that came out on the far side. Which
  // matters because a degenerate configuration can produce a hundred copies of one vertex and there is no
  // bound on how many: comparing the copies with one another would be quadratic in a quantity that is not
  // the size of the answer.
  {
    constexpr double coincidenceTolerance = 1.0e-6;
    std::size_t rawCount = diagram.vertices.size();
    std::vector<std::size_t> parent(rawCount);
    for (std::size_t v = 0; v < rawCount; ++v) parent[v] = v;
    auto find = [&parent](std::size_t node)
    {
      while (parent[node] != node) node = parent[node] = parent[parent[node]];
      return node;
    };

    std::vector<double3> rawPositions;
    rawPositions.reserve(rawCount);
    for (const SKApolloniusVertex& vertex : diagram.vertices) rawPositions.push_back(vertex.position);
    std::vector<std::size_t> groupOf =
        coincidentGroups(unitCell, inverseCell, perpendicularWidths, rawPositions, coincidenceTolerance);
    for (std::size_t v = 0; v < rawCount; ++v)
    {
      std::size_t left = find(v);
      std::size_t right = find(groupOf[v]);
      if (left != right) parent[left] = right;
    }

    std::map<std::size_t, std::vector<std::size_t>> groups;
    for (std::size_t v = 0; v < rawCount; ++v) groups[find(v)].push_back(v);

    std::vector<SKApolloniusVertex> merged;
    merged.reserve(groups.size());
    for (const auto& [root, members] : groups)
    {
      // The first member's frame is the one the set is expressed in; the others are shifted onto it, so
      // that a copy wrapped to the far face of the cell contributes its sites at the images they occupy
      // as seen from here.
      const SKApolloniusVertex& leader = diagram.vertices[members.front()];
      std::vector<std::pair<std::size_t, int3>> cotangent;
      for (std::size_t member : members)
      {
        const SKApolloniusVertex& copy = diagram.vertices[member];
        double3 delta = inverseCell * (leader.position - copy.position);
        int3 shift(static_cast<std::int32_t>(std::llround(delta.x)), static_cast<std::int32_t>(std::llround(delta.y)),
                   static_cast<std::int32_t>(std::llround(delta.z)));
        for (std::size_t s = 0; s < copy.siteIndices.size(); ++s)
          cotangent.push_back({copy.siteIndices[s], int3(copy.siteImages[s].x + shift.x, copy.siteImages[s].y + shift.y,
                                                        copy.siteImages[s].z + shift.z)});
      }
      std::sort(cotangent.begin(), cotangent.end(),
                [](const auto& lhs, const auto& rhs)
                {
                  return std::tie(lhs.first, lhs.second.x, lhs.second.y, lhs.second.z) <
                         std::tie(rhs.first, rhs.second.x, rhs.second.y, rhs.second.z);
                });
      cotangent.erase(std::unique(cotangent.begin(), cotangent.end(),
                                  [](const auto& lhs, const auto& rhs)
                                  {
                                    return lhs.first == rhs.first && lhs.second.x == rhs.second.x &&
                                           lhs.second.y == rhs.second.y && lhs.second.z == rhs.second.z;
                                  }),
                      cotangent.end());

      SKApolloniusVertex vertex{leader.position, leader.radius, {}, {}, {}};
      std::vector<double3> siteCentres;
      for (const auto& [index, image] : cotangent)
      {
        vertex.siteIndices.push_back(index);
        vertex.siteImages.push_back(image);
        siteCentres.push_back(grid.positions[index] + unitCell * double3(image.x, image.y, image.z));
      }
      vertex.branches = vertexBranches(vertex.position, siteCentres, diagram.verification.ambiguousBranches);
      if (vertex.siteIndices.size() > 4) ++diagram.verification.degenerateVertices;
      merged.push_back(std::move(vertex));
    }
    diagram.vertices = std::move(merged);
  }

  std::map<SiteTuple, std::vector<IncidentBranch>> tripleLookup = buildTripleLookup();

  for (const auto& [key, incident] : tripleLookup)
  {
    if (incident.size() < 2)
    {
      // One vertex on a triple means its arc does not end at a second vertex inside the free region.
      // The subdivision sweep above is exhaustive over vertices of non-negative clearance, so the
      // missing end cannot be one this construction was entitled to find: it is a vertex of negative
      // clearance, lying inside the spheres, at which the arc leaves the free region. Such an arc is
      // clipped rather than absent, and is counted apart from a genuine pairing failure.
      ++diagram.verification.truncatedTriples;
      continue;
    }

    // The triple as seen in the frame of the key itself, i.e. with zero applied offset.
    std::array<std::size_t, 3> baseIndices;
    std::array<int3, 3> baseImages;
    std::array<double3, 3> baseCentres;
    std::array<double, 3> tangentRadii;
    for (std::size_t s = 0; s < 3; ++s)
    {
      baseIndices[s] = static_cast<std::size_t>(key.data[4 * s + 0]);
      baseImages[s] = int3(static_cast<std::int32_t>(key.data[4 * s + 1]),
                           static_cast<std::int32_t>(key.data[4 * s + 2]),
                           static_cast<std::int32_t>(key.data[4 * s + 3]));
      baseCentres[s] = grid.positions[baseIndices[s]] +
                       unitCell * double3(baseImages[s].x, baseImages[s].y, baseImages[s].z);
      tangentRadii[s] = radii[baseIndices[s]];
    }

    // A trisector carries an ordering of its points, and a vertex bounds an edge together with its
    // neighbour in that ordering, not with an arbitrary partner. What orders it is the curve's own
    // parameter, which is come by in closed form, because a trisector is always a plane conic.
    //
    // Subtracting the tangency conditions pairwise leaves, for i = 1 and 2,
    //
    //     n_i . x + d_i t = K_i,   n_i = 2(c_i - c_0),  d_i = 2(r_i - r_0),
    //                              K_i = |c_i|^2 - |c_0|^2 - r_i^2 + r_0^2,
    //
    // one equation for each pair, linear in the position x and in the clearance t alike. Unless all three
    // radii are alike, t can be eliminated between the two, and what is left is a single linear equation in
    // x: the plane the curve lies in. Within that plane t is an affine function of position, so what remains
    // of the tangency to site 0, |x - c_0| = r_0 + t, is a quadratic. Laying the plane's first axis along the
    // gradient of t leaves that quadratic already in principal position, with coefficients 1 - g^2 and 1 for
    // g the size of the gradient, so there is no eigenproblem to solve: the conic is an ellipse for g below
    // one and the branch of a hyperbola for g above it, and either way its own parameter is at hand.
    //
    // Nothing in this is a special arrangement to be recognised, which is the point of doing it this way.
    // Collinear centres are g = 0: the clearance does not vary within the plane, the ellipse is a circle,
    // and the angle round it orders the vertices. That is the arrangement no ordering by clearance can
    // handle, since the clearance is the same at every point of it. Nearly collinear centres are a small g
    // and a nearly round ellipse, ordered by the very same formula, so nothing turns on how nearly collinear
    // they are, which is what a threshold on collinearity would have turned on. Radii nearly alike are a
    // large g and a nearly straight branch, and again the same formula, approaching the straight line it
    // ought to. Only radii exactly alike stand apart, and only because there is then no t to eliminate: the
    // two equations are planes in their own right and the curve is the line they cut in.
    enum class Shape
    {
      Line,     // three equal radii: the ordinary Voronoi edge, ordered along its direction
      Ellipse,  // ordered by the angle round it, and closed, so its last vertex meets its first
      Branch    // one branch of a hyperbola, or a parabola, ordered across the axis of the conic
    };

    // The two pairwise equations, the one with the larger difference of radii first so that the clearance is
    // eliminated using the better conditioned of the two.
    std::array<double3, 2> difference;
    std::array<double, 2> radiusDifference;
    std::array<double, 2> constants;
    for (std::size_t i = 0; i < 2; ++i)
    {
      difference[i] = 2.0 * (baseCentres[i + 1] - baseCentres[0]);
      radiusDifference[i] = 2.0 * (tangentRadii[i + 1] - tangentRadii[0]);
      constants[i] = double3::dot(baseCentres[i + 1], baseCentres[i + 1]) -
                     double3::dot(baseCentres[0], baseCentres[0]) - tangentRadii[i + 1] * tangentRadii[i + 1] +
                     tangentRadii[0] * tangentRadii[0];
    }
    if (std::abs(radiusDifference[1]) > std::abs(radiusDifference[0]))
    {
      std::swap(difference[0], difference[1]);
      std::swap(radiusDifference[0], radiusDifference[1]);
      std::swap(constants[0], constants[1]);
    }

    Shape shape = Shape::Line;
    double3 lineDirection(0.0, 0.0, 0.0);
    double3 planeOrigin(0.0, 0.0, 0.0);
    double3 firstAxis(0.0, 0.0, 0.0);   // along the gradient of the clearance within the plane
    double3 secondAxis(0.0, 0.0, 0.0);  // across it
    double centreOffset = 0.0;          // where the conic is centred along the first axis
    double firstSemiAxis = 0.0;
    double secondSemiAxis = 0.0;
    // The clearance is affine along the first axis, t = originClearance + rate * u, which is what makes
    // the bottleneck of an arc a matter of the least u it reaches and no search at all.
    double rate = 0.0;
    double originClearance = 0.0;
    double reach = 0.0;
    double flattening = 0.0;
    double extent = 0.0;
    double heightSquared = 0.0;  // |planeOrigin - c_0|^2, the offset of the plane from site 0

    double radiusScale = std::max({tangentRadii[0], tangentRadii[1], tangentRadii[2], 1.0});
    if (std::abs(radiusDifference[0]) < 1.0e-14 * radiusScale)
    {
      lineDirection = double3::cross(difference[0], difference[1]);
      double straightness = difference[0].length() * difference[1].length();
      if (lineDirection.length() < 1.0e-12 * straightness)
      {
        // Equal radii on collinear centres: the two planes are parallel and there is no curve to order.
        ++diagram.verification.unpairedTriples;
        continue;
      }
      lineDirection = lineDirection / lineDirection.length();
    }
    else
    {
      double3 planeNormal = radiusDifference[0] * difference[1] - radiusDifference[1] * difference[0];
      double planeScale = std::abs(radiusDifference[0]) * difference[1].length() +
                          std::abs(radiusDifference[1]) * difference[0].length();
      if (planeNormal.length() < 1.0e-12 * planeScale)
      {
        // The two equations say the same thing, so the points of equal clearance are a surface rather than
        // a curve and there is no ordering to make. No vertex can arise on such a triple to begin with.
        ++diagram.verification.unpairedTriples;
        continue;
      }
      double planeConstant =
          (radiusDifference[0] * constants[1] - radiusDifference[1] * constants[0]) / planeNormal.length();
      planeNormal = planeNormal / planeNormal.length();

      // The plane, taken with its origin at the point of it nearest site 0, so that the offset from that
      // site is along the normal and drops out of the conic's linear terms.
      planeOrigin = baseCentres[0] + (planeConstant - double3::dot(planeNormal, baseCentres[0])) * planeNormal;
      double3 fromSite = planeOrigin - baseCentres[0];

      double3 gradient =
          -(difference[0] - double3::dot(difference[0], planeNormal) * planeNormal) / radiusDifference[0];
      rate = gradient.length();
      firstAxis = (rate > 0.0) ? gradient / rate
                               : double3::cross(planeNormal, (std::abs(planeNormal.x) < 0.9) ? double3(1.0, 0.0, 0.0)
                                                                                             : double3(0.0, 1.0, 0.0));
      firstAxis = firstAxis / firstAxis.length();
      secondAxis = double3::cross(planeNormal, firstAxis);

      originClearance = (constants[0] - double3::dot(difference[0], planeOrigin)) / radiusDifference[0];
      reach = tangentRadii[0] + originClearance;
      heightSquared = double3::dot(fromSite, fromSite);
      flattening = 1.0 - rate * rate;

      if (rate < 1.0)
      {
        // (1 - g^2)(u - u_c)^2 + v^2 = C, an ellipse, and a circle when the clearance does not vary at all.
        shape = Shape::Ellipse;
        centreOffset = reach * rate / flattening;
        extent = reach * reach / flattening - heightSquared;
        if (extent <= 0.0)
        {
          // No real curve at all, so the vertices reported on it cannot be there.
          ++diagram.verification.unpairedTriples;
          continue;
        }
        secondSemiAxis = std::sqrt(extent);
        firstSemiAxis = std::sqrt(extent / flattening);
      }
      else
      {
        // Above one the conic opens: each branch of it, and a parabola too, gives the first coordinate as a
        // function of the second, so the second orders the curve outright. Only one of the two branches can
        // be a trisector, the other being where the tangency holds with its sign reversed, so there is one
        // curve to order here as much as in the other cases.
        shape = Shape::Branch;
        // The same constants as the ellipse's, and used the same way to place a point on the curve. At
        // exactly one they are not there to be had, the conic being a parabola, which is carried by its
        // own form below rather than by a case of its own.
        if (std::abs(flattening) > 1.0e-12)
        {
          centreOffset = reach * rate / flattening;
          extent = reach * reach / flattening - heightSquared;
        }
      }
    }

    struct OrderedVertex
    {
      std::size_t vertex;
      int3 offset;
      double along;  // position along the curve, increasing the way the ordering runs
      bool forward;  // whether the edge leaves along the direction the ordering runs
    };

    std::vector<OrderedVertex> ordered;
    ordered.reserve(incident.size());
    for (const IncidentBranch& entry : incident)
    {
      // Undo the offset that canonicalisation applied, which puts every incident vertex in the
      // single frame in which they all share this copy of the triple.
      const SKApolloniusVertex& incidentVertex = diagram.vertices[entry.vertex];
      double3 position = incidentVertex.position - unitCell * double3(entry.offset.x, entry.offset.y, entry.offset.z);
      double3 outgoing = incidentVertex.branches[entry.branch].direction;

      double along = 0.0;
      bool forward = false;
      if (shape == Shape::Line)
      {
        along = double3::dot(position, lineDirection);
        forward = double3::dot(outgoing, lineDirection) > 0.0;
      }
      else
      {
        double first = double3::dot(position - planeOrigin, firstAxis) - centreOffset;
        double second = double3::dot(position - planeOrigin, secondAxis);
        double firstStep = double3::dot(outgoing, firstAxis);
        double secondStep = double3::dot(outgoing, secondAxis);
        if (shape == Shape::Ellipse)
        {
          // The angle of the point on the ellipse, taken on the ellipse's own scale in each direction so
          // that a flattened one is measured as soundly as a round one, and the tangent there.
          along = std::atan2(second / secondSemiAxis, first / firstSemiAxis);
          forward = secondSemiAxis * secondSemiAxis * first * secondStep -
                        firstSemiAxis * firstSemiAxis * second * firstStep >
                    0.0;
        }
        else
        {
          along = second;
          forward = secondStep > 0.0;
        }
      }

      ordered.push_back(OrderedVertex{entry.vertex, entry.offset, along, forward});
    }
    if (ordered.size() < 2)
    {
      ++diagram.verification.unpairedTriples;
      continue;
    }
    std::sort(ordered.begin(), ordered.end(),
              [](const OrderedVertex& lhs, const OrderedVertex& rhs) { return lhs.along < rhs.along; });

    // A point of the curve at a given value of that same parameter, and the clearance there. Both are in
    // closed form: the position because the conic is in principal position in the plane's own axes, and
    // the clearance because it is affine along the first of them, t = originClearance + g u. So a sample
    // costs a cosine at worst and needs no tangency system solved for it.
    //
    // The straight line of three equal radii is the one case with no such t to be had, the clearance
    // there being the distance to a site rather than an affine function, and it is taken as it comes.
    double3 lineBase(0.0, 0.0, 0.0);
    double branchSign = 1.0;
    {
      const OrderedVertex& reference = ordered.front();
      double3 referencePosition =
          diagram.vertices[reference.vertex].position -
          unitCell * double3(reference.offset.x, reference.offset.y, reference.offset.z);
      if (shape == Shape::Line)
      {
        lineBase = referencePosition - reference.along * lineDirection;
      }
      else if (shape == Shape::Branch && std::abs(flattening) > 1.0e-12)
      {
        // Which of the hyperbola's two branches is the trisector, the other holding the tangency with
        // its sign reversed. The whole curve is on one of them, so one vertex settles it.
        double first = double3::dot(referencePosition - planeOrigin, firstAxis) - centreOffset;
        branchSign = (first < 0.0) ? -1.0 : 1.0;
      }
    }

    auto curvePoint = [&](double along) -> std::optional<std::pair<double3, double>>
    {
      if (shape == Shape::Line)
      {
        double3 position = lineBase + along * lineDirection;
        return std::pair{position, (position - baseCentres[0]).length() - tangentRadii[0]};
      }

      double first = 0.0;
      double second = 0.0;
      if (shape == Shape::Ellipse)
      {
        first = firstSemiAxis * std::cos(along);
        second = secondSemiAxis * std::sin(along);
      }
      else
      {
        second = along;
        if (std::abs(flattening) > 1.0e-12)
        {
          double squared = (extent - second * second) / flattening;
          if (squared < 0.0) return std::nullopt;
          first = branchSign * std::sqrt(squared);
        }
        else
        {
          // Equal to one the conic is a parabola, which has no centre to measure from: 2 h u = v^2 +
          // |c_0 - o|^2 - h^2 for h the reach, straight from the tangency to site 0.
          if (std::abs(reach) < 1.0e-12) return std::nullopt;
          first = (second * second + heightSquared - reach * reach) / (2.0 * reach);
        }
      }

      double alongFirstAxis = centreOffset + first;
      return std::pair{planeOrigin + alongFirstAxis * firstAxis + second * secondAxis,
                       originClearance + rate * alongFirstAxis};
    };

    // The clearance along the curve is smooth and turns at one point at most: the far end of the ellipse,
    // the waist of the open branch, or the foot of the perpendicular from site 0 on the straight line. The
    // least clearance on an arc is therefore at one of its two ends or at that point, so putting it among
    // the samples makes the bottleneck exact instead of as fine as the sampling is.
    std::array<double, 3> turningPoints{};
    std::size_t turningCount = 0;
    if (shape == Shape::Line)
    {
      turningPoints[turningCount++] = double3::dot(baseCentres[0] - lineBase, lineDirection);
    }
    else if (shape == Shape::Ellipse)
    {
      // Where the first coordinate is least, which is where an affine clearance is least. A circle, the
      // mark of collinear centres, has the same clearance everywhere and no turn to look for.
      if (rate > 0.0)
      {
        turningPoints[turningCount++] = -std::numbers::pi;
        turningPoints[turningCount++] = std::numbers::pi;
        turningPoints[turningCount++] = 3.0 * std::numbers::pi;
      }
    }
    else
    {
      turningPoints[turningCount++] = 0.0;
    }

    // Consecutive vertices along the curve delimit its arcs, including the arc that wraps from the
    // last back to the first, which is a real edge on a closed trisector and is the pair of infinite
    // ends on an open one. An arc is in the diagram only if the edges leaving both of its endpoints
    // run into it, which is what tells the two apart and needs no classification of the curve: on an
    // open trisector the extreme vertices send their edges outwards, so the wrapping arc is refused.
    for (std::size_t p = 0; p < ordered.size(); ++p)
    {
      std::size_t q = (p + 1) % ordered.size();
      if (q == p) continue;

      // Going from p to q runs with the ordering, so p must leave that way and q must leave against it
      // for the arc between them to be an edge.
      if (!ordered[p].forward) continue;
      if (ordered[q].forward) continue;

      std::size_t fromVertex = ordered[p].vertex;
      std::size_t toVertex = ordered[q].vertex;
      const int3& fromOffset = ordered[p].offset;
      const int3& toOffset = ordered[q].offset;

      // The two halves of a trisector that runs through a vertex are consecutive in this ordering, and
      // the arc between them is the point itself rather than an edge.
      if (fromVertex == toVertex && fromOffset.x == toOffset.x && fromOffset.y == toOffset.y &&
          fromOffset.z == toOffset.z)
        continue;

      int3 toImage(fromOffset.x - toOffset.x, fromOffset.y - toOffset.y, fromOffset.z - toOffset.z);

      std::array<std::size_t, 3> indices = baseIndices;
      std::array<int3, 3> images;
      std::array<double3, 3> centres;
      for (std::size_t s = 0; s < 3; ++s)
      {
        images[s] = int3(baseImages[s].x + fromOffset.x, baseImages[s].y + fromOffset.y,
                         baseImages[s].z + fromOffset.z);
        centres[s] = grid.positions[indices[s]] + unitCell * double3(images[s].x, images[s].y, images[s].z);
      }

      double3 from = diagram.vertices[fromVertex].position;
      double3 to = diagram.vertices[toVertex].position + unitCell * double3(toImage.x, toImage.y, toImage.z);

      // Sample the arc to get its bottleneck and length, walking the curve in its own parameter: the very
      // one the vertices were ordered by, so that the samples land on the stretch of curve between them
      // and nowhere else. Cutting the chord into equal parts instead reaches only the part of a strongly
      // curved arc that projects onto the chord, and where the cutting plane meets the curve twice it
      // takes whichever point is nearer the chord, which on such an arc is the wrong one. Either way it
      // reads a bottleneck wider than the arc really has, and a probe is then let through a wall.
      //
      // The stretch of the curve inside the sites is followed too. An arc of the free space may leave
      // it in the middle and return: that is what the window between two cages is, wide at both ends
      // and shut where it passes the ring of atoms. Its bottleneck is the clearance at that ring, and
      // it is negative, which is how a passage a probe cannot pass is told apart from one it can. Stop
      // at the boundary of the free region instead and both look alike.
      // The samples are kept, since the narrowest of them is where the passage is tightest and the
      // window across it is cut there; a sample that could not be placed is left out of that choice,
      // as it is left out of the bottleneck.
      constexpr std::size_t sampleCount = 16;
      std::array<double3, sampleCount + 1 + turningPoints.size()> curve;
      std::array<double, sampleCount + 1 + turningPoints.size()> clearance;

      // The arc runs from p to q the way the ordering runs. On the closed ellipse the last vertex meets
      // the first the long way round, which is a turn of the parameter past its wrap.
      double alongFrom = ordered[p].along;
      double alongTo = ordered[q].along;
      bool spanned = true;
      if (alongTo < alongFrom)
      {
        if (shape == Shape::Ellipse)
          alongTo += 2.0 * std::numbers::pi;
        else
          spanned = false;  // an open curve is ordered one way only, so there is no such arc to walk
      }

      // The curve is placed in the frame of the key, and this arc in the frame where its own `from`
      // vertex sits at the position it is stored at.
      double3 toArcFrame = unitCell * double3(fromOffset.x, fromOffset.y, fromOffset.z);

      curve[0] = from;
      clearance[0] = diagram.vertices[fromVertex].radius;
      curve[sampleCount] = to;
      clearance[sampleCount] = diagram.vertices[toVertex].radius;

      double length = 0.0;
      bool sampled = spanned;
      for (std::size_t s = 1; s <= sampleCount; ++s)
      {
        double parameter = static_cast<double>(s) / static_cast<double>(sampleCount);
        std::optional<std::pair<double3, double>> point =
            spanned ? curvePoint(alongFrom + parameter * (alongTo - alongFrom)) : std::nullopt;
        if (point.has_value() && s < sampleCount)
        {
          curve[s] = point->first + toArcFrame;
          clearance[s] = point->second;
        }
        else if (s < sampleCount)
        {
          curve[s] = from + parameter * (to - from);
          clearance[s] = std::numeric_limits<double>::max();
          sampled = false;
        }
        length += (curve[s] - curve[s - 1]).length();
      }

      // The turn of the clearance, where the arc reaches it, which is the bottleneck exactly. Its own
      // neighbours are kept alongside it, since the samples no longer run in order once it is added and
      // the direction across the passage is read off them below.
      std::size_t sampleTotal = sampleCount + 1;
      std::array<double3, turningPoints.size()> turningChord{};
      for (std::size_t k = 0; k < turningCount && spanned; ++k)
      {
        double turn = turningPoints[k];
        if (turn <= alongFrom || turn >= alongTo) continue;
        std::optional<std::pair<double3, double>> point = curvePoint(turn);
        if (!point.has_value())
        {
          sampled = false;
          continue;
        }
        double step = 1.0e-4 * std::max(std::abs(alongTo - alongFrom), 1.0e-6);
        std::optional<std::pair<double3, double>> before = curvePoint(std::max(turn - step, alongFrom));
        std::optional<std::pair<double3, double>> after = curvePoint(std::min(turn + step, alongTo));
        turningChord[sampleTotal - (sampleCount + 1)] =
            (before.has_value() && after.has_value()) ? after->first - before->first : to - from;
        curve[sampleTotal] = point->first + toArcFrame;
        clearance[sampleTotal] = point->second;
        ++sampleTotal;
      }
      if (!sampled) ++diagram.verification.unsampledArcs;

      std::size_t narrowest = 0;
      for (std::size_t s = 1; s < sampleTotal; ++s)
        if (clearance[s] < clearance[narrowest]) narrowest = s;

      // The tangent of the trisector at the narrowest sample. The curve is the intersection of two
      // clearance bisectors, whose normals at a point are the differences of the unit vectors from the
      // sites to it, so the tangent is the cross product of those two differences: exact, and no more
      // work than the chord between neighbouring samples would be. It is oriented by that chord, since
      // the cross product fixes the line and not which way along it the arc runs. Where the two normals
      // are parallel, which is where the curve is not cut out by them transversally, the chord is all
      // there is.
      double3 chord = (narrowest > sampleCount)
                          ? turningChord[narrowest - (sampleCount + 1)]
                          : curve[std::min(narrowest + 1, sampleCount)] - curve[narrowest > 0 ? narrowest - 1 : 0];
      std::array<double3, 3> fromSite;
      for (std::size_t s = 0; s < 3; ++s)
      {
        double3 offset = curve[narrowest] - centres[s];
        double offsetLength = offset.length();
        fromSite[s] = offsetLength > 0.0 ? offset / offsetLength : double3(0.0, 0.0, 0.0);
      }
      double3 tangent = double3::cross(fromSite[0] - fromSite[1], fromSite[0] - fromSite[2]);
      if (tangent.length() < 1.0e-10) tangent = chord;
      if (double3::dot(tangent, chord) < 0.0) tangent = -tangent;

      double tangentLength = tangent.length();
      double3 direction = tangentLength > 0.0 ? tangent / tangentLength : double3(0.0, 0.0, 0.0);

      diagram.edges.push_back(SKApolloniusEdge{fromVertex, toVertex, toImage, indices, images, clearance[narrowest],
                                               curve[narrowest], direction, length, false});
    }
  }

  // Trisectors that carry no vertex.
  //
  // Everything above reaches the diagram through its vertices, so a trisector on which no fourth site
  // ever intrudes is invisible to it: there is no vertex to start from. Such a curve is closed, and it
  // is a real Apollonius edge, a ring-shaped channel whose narrowest point is a bottleneck a probe must
  // pass. Wang et al. describe them in Sections 4.2 and 5.5 and close the gap by replicating the
  // offending site into four perturbed copies, which reintroduces vertices artificially. They are found
  // directly here instead, by sweeping the clearance along each candidate trisector: the sweep decides
  // whether the curve is bounded, and where it is bounded it also traverses it, so the same pass proves
  // the curve closed, proves no site intrudes, and measures the bottleneck and length.
  //
  // The widest clearance each ring reaches is kept alongside it: together with the bottleneck it gives the
  // band of clearances the ring occupies, which is what decides later whether the ring bounds a face of
  // its own or is a hole in a larger one.
  std::map<std::size_t, double> ringTopClearance;
  {
    double smallestRadius = std::numeric_limits<double>::max();
    for (double radius : radii) smallestRadius = std::min(smallestRadius, radius);
    double lowestClearance = allowNegativeRadius ? -smallestRadius : 0.0;

    constexpr std::size_t bracketSamples = 192;
    constexpr std::size_t traverseSamples = 96;
    constexpr std::size_t bisectionSteps = 40;
    double intrusionTolerance = 1.0e-7;

    // In the order the triples canonicalise in, so that the rings come out in an order that is the
    // structure's and not the hash table's.
    std::vector<SiteTuple> ringCandidates(candidateTriples.begin(), candidateTriples.end());
    std::sort(ringCandidates.begin(), ringCandidates.end());

    for (const SiteTuple& key : ringCandidates)
    {
      if (tripleLookup.contains(key)) continue;  // the curve already carries vertices

      std::array<std::size_t, 3> indices;
      std::array<int3, 3> images;
      std::array<double3, 3> centres;
      std::array<double, 3> tripleRadii;
      for (std::size_t s = 0; s < 3; ++s)
      {
        indices[s] = static_cast<std::size_t>(key.data[4 * s + 0]);
        images[s] = int3(static_cast<std::int32_t>(key.data[4 * s + 1]),
                         static_cast<std::int32_t>(key.data[4 * s + 2]),
                         static_cast<std::int32_t>(key.data[4 * s + 3]));
        centres[s] = grid.positions[indices[s]] + unitCell * double3(images[s].x, images[s].y, images[s].z);
        tripleRadii[s] = radii[indices[s]];
      }

      // An open trisector reaches every clearance above its minimum, so if the curve still exists at the
      // largest clearance any empty sphere can have, it is not closed and is not what is sought here.
      if (!trisectorPointsAtClearance(centres, tripleRadii, emptyRadiusBound).empty()) continue;

      // Bracket the clearances the curve occupies.
      double firstPresent = 0.0;
      double lastPresent = 0.0;
      bool anyPresent = false;
      for (std::size_t sample = 0; sample <= bracketSamples; ++sample)
      {
        double clearance = lowestClearance + (emptyRadiusBound - lowestClearance) *
                                                 static_cast<double>(sample) / static_cast<double>(bracketSamples);
        if (trisectorPointsAtClearance(centres, tripleRadii, clearance).empty()) continue;
        if (!anyPresent) firstPresent = clearance;
        lastPresent = clearance;
        anyPresent = true;
      }
      if (!anyPresent) continue;

      // Refine the two ends of that range. They are where the branches meet, so the lower end is the
      // bottleneck of the ring and the upper end its widest point.
      double clearanceBelow = std::max(lowestClearance, firstPresent - (emptyRadiusBound - lowestClearance) /
                                                                          static_cast<double>(bracketSamples));
      double minimumClearance = firstPresent;
      for (std::size_t step = 0; step < bisectionSteps; ++step)
      {
        double middle = 0.5 * (clearanceBelow + minimumClearance);
        if (trisectorPointsAtClearance(centres, tripleRadii, middle).empty())
          clearanceBelow = middle;
        else
          minimumClearance = middle;
      }
      double clearanceAbove =
          lastPresent + (emptyRadiusBound - lowestClearance) / static_cast<double>(bracketSamples);
      double maximumClearance = lastPresent;
      for (std::size_t step = 0; step < bisectionSteps; ++step)
      {
        double middle = 0.5 * (maximumClearance + clearanceAbove);
        if (trisectorPointsAtClearance(centres, tripleRadii, middle).empty())
          clearanceAbove = middle;
        else
          maximumClearance = middle;
      }
      if (maximumClearance <= minimumClearance) continue;

      // Traverse the ring, on both branches, and reject it the moment any site is nearer than the
      // clearance the three defining sites share. Those three sit exactly at that clearance, so the
      // tolerance excludes them without having to identify them.
      bool intruded = false;
      double length = 0.0;
      std::array<std::optional<double3>, 2> previousOnBranch;
      std::array<double3, 2> firstOnBranch;
      std::array<std::optional<double3>, 2> secondOnBranch;
      for (std::size_t sample = 0; sample <= traverseSamples && !intruded; ++sample)
      {
        double clearance = minimumClearance + (maximumClearance - minimumClearance) *
                                                  static_cast<double>(sample) / static_cast<double>(traverseSamples);
        std::vector<double3> points = trisectorPointsAtClearance(centres, tripleRadii, clearance);
        for (std::size_t branch = 0; branch < points.size() && !intruded; ++branch)
        {
          double3 wrapped = unitCell * double3::fract(inverseCell * points[branch]);
          grid.forEachNear(wrapped, std::max(0.0, clearance + maximumRadius),
                           [&](std::size_t j, const double3& image, const int3&)
                           {
                             if ((wrapped - image).length() - radii[j] < clearance - intrusionTolerance)
                               intruded = true;
                           });

          if (previousOnBranch[branch].has_value())
          {
            length += (points[branch] - previousOnBranch[branch].value()).length();
            if (!secondOnBranch[branch].has_value()) secondOnBranch[branch] = points[branch];
          }
          else
          {
            firstOnBranch[branch] = points[branch];
          }
          previousOnBranch[branch] = points[branch];
        }
      }
      if (intruded) continue;

      // The two branches meet at both ends, so closing the ring joins their far ends to each other.
      if (previousOnBranch[0].has_value() && previousOnBranch[1].has_value())
      {
        length += (previousOnBranch[0].value() - previousOnBranch[1].value()).length();
        length += (firstOnBranch[0] - firstOnBranch[1]).length();
      }

      // The narrowest point of a ring is the pinch where its two branches meet, at the lowest clearance
      // it reaches. The ring runs from one branch to the other through that point, so the chord between
      // the branches just above the pinch is the direction it runs in there.
      double3 pinch = previousOnBranch[0].has_value() ? firstOnBranch[0] : firstOnBranch[1];
      double3 chord(0.0, 0.0, 0.0);
      if (previousOnBranch[0].has_value() && previousOnBranch[1].has_value())
      {
        pinch = 0.5 * (firstOnBranch[0] + firstOnBranch[1]);
        if (secondOnBranch[0].has_value() && secondOnBranch[1].has_value())
          chord = secondOnBranch[0].value() - secondOnBranch[1].value();
      }
      double chordLength = chord.length();

      ringTopClearance[diagram.edges.size()] = maximumClearance;
      diagram.edges.push_back(SKApolloniusEdge{std::numeric_limits<std::size_t>::max(),
                                               std::numeric_limits<std::size_t>::max(), int3(0, 0, 0), indices, images,
                                               minimumClearance, pinch,
                                               chordLength > 0.0 ? chord / chordLength : double3(0.0, 0.0, 0.0), length,
                                               true});
      ++diagram.verification.vertexlessLoops;
    }
  }

  // Faces, one per site of each pair so that every cell owns its own, and cells.
  std::map<SiteTuple, std::size_t> faceLookup;
  auto faceKey = [](std::size_t owner, std::size_t other, const int3& relativeImage)
  {
    SiteTuple key;
    key.data[0] = static_cast<std::int64_t>(owner);
    key.data[4] = static_cast<std::int64_t>(other);
    key.data[5] = relativeImage.x;
    key.data[6] = relativeImage.y;
    key.data[7] = relativeImage.z;
    return key;
  };

  // Group the edges by the bisector surface they lie on. That surface is not yet a face: the edges on
  // one surface can enclose several regions of it that are separate faces of the diagram (Wang et al.,
  // Section 5.4), so each group is split below into its connected parts. Which two of an edge's three
  // sites the surface belongs to is kept with it, since a ring needs to know its third site to work out
  // which region of the surface it bounds.
  struct PatchEdge
  {
    std::size_t edgeIndex;
    std::size_t ownerLocal;  // position of the owning site among the edge's three
    std::size_t otherLocal;
  };
  struct PrimitivePatch
  {
    std::size_t owner;
    std::size_t other;
    int3 relativeImage;
    std::vector<PatchEdge> edges;
  };
  std::vector<PrimitivePatch> patches;
  std::map<SiteTuple, std::size_t> patchLookup;

  for (std::size_t e = 0; e < diagram.edges.size(); ++e)
  {
    const SKApolloniusEdge& edge = diagram.edges[e];
    for (std::size_t s = 0; s < 3; ++s)
      for (std::size_t o = 0; o < 3; ++o)
      {
        if (s == o) continue;
        int3 relative(edge.siteImages[o].x - edge.siteImages[s].x, edge.siteImages[o].y - edge.siteImages[s].y,
                      edge.siteImages[o].z - edge.siteImages[s].z);
        SiteTuple key = faceKey(edge.siteIndices[s], edge.siteIndices[o], relative);
        auto [iterator, inserted] = patchLookup.try_emplace(key, patches.size());
        if (inserted) patches.push_back(PrimitivePatch{edge.siteIndices[s], edge.siteIndices[o], relative, {}});
        patches[iterator->second].edges.push_back(PatchEdge{e, s, o});
      }
  }

  // Which region of a bisector surface a ring bounds.
  //
  // A ring is a closed curve confined to the band of clearances it spans, and its own trisector is the
  // only place on the surface where the third site is exactly as near, so the parts of the surface below
  // and above that band each lie wholly on one side of the ring. The diagram keeps the side where the
  // third site is the farther one, so testing that site just outside the band settles which side that is.
  // If either part beyond the band is kept, the region next to the ring reaches past it and the ring is
  // one boundary component of a larger region, a hole punched in it. If neither is, the region is the disc
  // the ring encloses and the ring bounds a face on its own, which is the face hump of Wang et al.,
  // Section 5.5.
  auto ringBoundsOwnFace = [&](const PatchEdge& patchEdge)
  {
    const SKApolloniusEdge& ring = diagram.edges[patchEdge.edgeIndex];
    auto found = ringTopClearance.find(patchEdge.edgeIndex);
    if (found == ringTopClearance.end()) return true;

    std::size_t thirdLocal = 3 - patchEdge.ownerLocal - patchEdge.otherLocal;
    auto centreOf = [&](std::size_t local)
    {
      return grid.positions[ring.siteIndices[local]] +
             unitCell * double3(ring.siteImages[local].x, ring.siteImages[local].y, ring.siteImages[local].z);
    };
    double3 thirdCentre = centreOf(thirdLocal);
    double thirdRadius = radii[ring.siteIndices[thirdLocal]];

    double margin = std::max(1.0e-6, 1.0e-3 * (found->second - ring.bottleneckRadius));
    for (double clearance : {ring.bottleneckRadius - margin, found->second + margin})
    {
      std::optional<double3> point =
          bisectorSheetPoint(centreOf(patchEdge.ownerLocal), radii[ring.siteIndices[patchEdge.ownerLocal]],
                             centreOf(patchEdge.otherLocal), radii[ring.siteIndices[patchEdge.otherLocal]], clearance,
                             0.0);
      if (!point.has_value()) continue;
      if ((point.value() - thirdCentre).length() - thirdRadius > clearance) return false;
    }
    return true;
  };

  // Two edges of a surface belong to the same face when they meet at a vertex, so the faces are the
  // connected components of the surface's boundary. Rings meet nothing and so cannot be placed that way;
  // each is instead assigned to the region it bounds, which is either a region already found, as one more
  // boundary component of it, or a region of its own.
  for (const PrimitivePatch& patch : patches)
  {
    std::map<std::size_t, std::size_t> parent;
    auto find = [&parent](std::size_t node)
    {
      while (parent[node] != node) node = parent[node] = parent[parent[node]];
      return node;
    };
    for (const PatchEdge& patchEdge : patch.edges)
    {
      if (diagram.edges[patchEdge.edgeIndex].isLoop) continue;
      for (std::size_t endpoint : {diagram.edges[patchEdge.edgeIndex].from, diagram.edges[patchEdge.edgeIndex].to})
        parent.try_emplace(endpoint, endpoint);
    }
    for (const PatchEdge& patchEdge : patch.edges)
    {
      if (diagram.edges[patchEdge.edgeIndex].isLoop) continue;
      std::size_t left = find(diagram.edges[patchEdge.edgeIndex].from);
      std::size_t right = find(diagram.edges[patchEdge.edgeIndex].to);
      if (left != right) parent[left] = right;
    }

    std::map<std::size_t, std::size_t> componentFace;
    for (const PatchEdge& patchEdge : patch.edges)
    {
      if (diagram.edges[patchEdge.edgeIndex].isLoop) continue;
      std::size_t root = find(diagram.edges[patchEdge.edgeIndex].from);
      auto [iterator, inserted] = componentFace.try_emplace(root, diagram.faces.size());
      if (inserted)
        diagram.faces.push_back(SKApolloniusFace{patch.owner, patch.other, patch.relativeImage, {}, false});
      diagram.faces[iterator->second].edgeIndices.push_back(patchEdge.edgeIndex);
    }

    for (const PatchEdge& patchEdge : patch.edges)
    {
      if (!diagram.edges[patchEdge.edgeIndex].isLoop) continue;

      if (!componentFace.empty() && !ringBoundsOwnFace(patchEdge))
      {
        // A hole belongs to the region around it. Deciding which region that is only takes work when the
        // surface carries several, since then the hole is inside one of them and the others are elsewhere
        // on the same surface; that has not been seen, and rather than guess it is counted so that it
        // cannot pass unnoticed.
        if (componentFace.size() > 1) ++diagram.verification.ringsOfUncertainFace;
        diagram.faces[componentFace.begin()->second].edgeIndices.push_back(patchEdge.edgeIndex);
        continue;
      }

      diagram.faces.push_back(
          SKApolloniusFace{patch.owner, patch.other, patch.relativeImage, {patchEdge.edgeIndex}, true});
    }
  }

  // A patch closes when every vertex on its boundary is met by exactly two of its edges. A ring carries no
  // vertex, so it neither satisfies nor violates that condition, and a face bounded only by rings is closed
  // by construction.
  for (SKApolloniusFace& face : diagram.faces)
  {
    std::map<std::size_t, std::size_t> incidence;
    for (std::size_t e : face.edgeIndices)
    {
      if (diagram.edges[e].isLoop) continue;
      ++incidence[diagram.edges[e].from];
      ++incidence[diagram.edges[e].to];
    }
    face.isClosed = !face.edgeIndices.empty();
    for (const auto& [vertex, degree] : incidence)
      if (degree != 2) face.isClosed = false;
    if (!face.isClosed) ++diagram.verification.unclosedFaces;
  }

  diagram.cells.resize(siteCount);
  for (std::size_t i = 0; i < siteCount; ++i) diagram.cells[i] = SKApolloniusCell{i, {}, {}, true};
  for (std::size_t f = 0; f < diagram.faces.size(); ++f)
  {
    diagram.cells[diagram.faces[f].site1].faceIndices.push_back(f);
    diagram.cells[diagram.faces[f].site1].isEmpty = false;
  }
  for (std::size_t v = 0; v < diagram.vertices.size(); ++v)
    for (std::size_t site : diagram.vertices[v].siteIndices) diagram.cells[site].vertexIndices.push_back(v);
  for (SKApolloniusCell& cell : diagram.cells)
  {
    std::sort(cell.vertexIndices.begin(), cell.vertexIndices.end());
    cell.vertexIndices.erase(std::unique(cell.vertexIndices.begin(), cell.vertexIndices.end()),
                             cell.vertexIndices.end());
  }

  // Verification: a vertex carries one edge along each of its branches, four of them in the general
  // case and more where the configuration is degenerate.
  std::vector<std::size_t> valence(diagram.vertices.size(), 0);
  for (const SKApolloniusEdge& edge : diagram.edges)
  {
    if (edge.isLoop) continue;  // a closed trisector has no endpoints to contribute valence to
    ++valence[edge.from];
    ++valence[edge.to];
  }
  diagram.verification.vertexCount = diagram.vertices.size();
  for (std::size_t v = 0; v < diagram.vertices.size(); ++v)
    if (valence[v] == diagram.vertices[v].branches.size()) ++diagram.verification.verticesOfFullValence;

  // Vertices left sharing a position would be a cotangent set that the gathering above failed to bring
  // together, which would leave each copy pairing along triples the others also claim. The gathering is
  // meant to make this impossible, so any survivor is a defect and is reported as one.
  std::vector<double3> mergedPositions;
  mergedPositions.reserve(diagram.vertices.size());
  for (const SKApolloniusVertex& vertex : diagram.vertices) mergedPositions.push_back(vertex.position);
  std::vector<std::size_t> coincidentWith =
      coincidentGroups(unitCell, inverseCell, perpendicularWidths, mergedPositions, 1.0e-6);
  for (std::size_t v = 0; v < coincidentWith.size(); ++v)
    if (coincidentWith[v] != v) ++diagram.verification.coincidentVertices;

  return diagram;
}
