module;

module skvoronoi;

import std;

import int3;
import double3;
import double3x3;

// Implementation notes
// --------------------
// The interface is fractional-native, but each Voronoi cell is constructed as a convex
// polyhedron in Cartesian displacement coordinates x = r - r_i (relative to the site):
// candidate separations are mapped through the cell matrix h once, after which every
// bisector is the plain Euclidean plane x·Δ = ½|Δ|². This is algebraically identical to
// cutting in fractional space with the metric tensor M = hᵀh (the plane distances are the
// same dot products), but avoids one matrix-vector product per candidate and per vertex
// bound. Fractional vertex positions are recovered through h⁻¹ on output only.
//
// Candidates are found with a cell list over the fractional unit cube: bins are walked
// outward in Chebyshev shells (wrapping periodically, so lattice images are enumerated on
// the fly), and the walk terminates once the Cartesian lower bound of the next shell
// exceeds twice the largest vertex distance of the current cell — beyond that no plane
// can intersect it.

namespace
{
struct PolyFace
{
  std::vector<std::size_t> vertexIndices;
  bool hasSource{false};
  std::size_t neighborIndex{};
  int3 neighborImage{};
  double rsq{};
  double3 delta{};  // Cartesian displacement towards the neighbor [Å]
};

struct Polyhedron
{
  std::vector<double3> vertices;  // Cartesian, relative to the site [Å]
  std::vector<PolyFace> faces;
};

struct Candidate
{
  double rsq;
  double3 delta;  // Cartesian displacement towards the neighbor [Å]
  std::size_t index;
  int3 image;
};

// Reusable scratch buffers; one instance is shared across all plane cuts of a cell (and
// across cells in computeAllCells), so steady-state cutting performs no heap allocation.
struct Workspace
{
  std::vector<double> distance;
  std::vector<std::int8_t> side;
  std::vector<std::size_t> clipped;  // clipped polygon under construction
  struct EdgeCut
  {
    std::size_t a;
    std::size_t b;
    std::size_t index;
  };
  std::vector<EdgeCut> edgeCuts;
  std::vector<std::pair<double, std::size_t>> loop;  // (angle, vertex) for the closing face
  std::vector<std::size_t> remap;
  std::vector<double3> compactedVertices;
  std::vector<Candidate> shellCandidates;
};

Polyhedron makeInitialCube(double halfWidth)
{
  Polyhedron poly;
  poly.vertices.reserve(64);
  poly.faces.reserve(32);
  poly.vertices = {double3(-halfWidth, -halfWidth, -halfWidth), double3(halfWidth, -halfWidth, -halfWidth),
                   double3(halfWidth, halfWidth, -halfWidth),   double3(-halfWidth, halfWidth, -halfWidth),
                   double3(-halfWidth, -halfWidth, halfWidth),  double3(halfWidth, -halfWidth, halfWidth),
                   double3(halfWidth, halfWidth, halfWidth),    double3(-halfWidth, halfWidth, halfWidth)};
  poly.faces.resize(6);
  poly.faces[0].vertexIndices = {0, 3, 2, 1};  // z = -h
  poly.faces[1].vertexIndices = {4, 5, 6, 7};  // z = +h
  poly.faces[2].vertexIndices = {0, 1, 5, 4};  // y = -h
  poly.faces[3].vertexIndices = {2, 3, 7, 6};  // y = +h
  poly.faces[4].vertexIndices = {0, 4, 7, 3};  // x = -h
  poly.faces[5].vertexIndices = {1, 2, 6, 5};  // x = +h
  return poly;
}

void compactVertices(Polyhedron &poly, Workspace &ws)
{
  ws.remap.assign(poly.vertices.size(), std::numeric_limits<std::size_t>::max());
  ws.compactedVertices.clear();
  for (PolyFace &face : poly.faces)
  {
    for (std::size_t &index : face.vertexIndices)
    {
      if (ws.remap[index] == std::numeric_limits<std::size_t>::max())
      {
        ws.remap[index] = ws.compactedVertices.size();
        ws.compactedVertices.push_back(poly.vertices[index]);
      }
      index = ws.remap[index];
    }
  }
  std::swap(poly.vertices, ws.compactedVertices);
}

// Cuts the polyhedron by the halfspace { x : x·n <= b }. The closing face inherits the
// neighbor bookkeeping of `source`. Returns true when the plane actually cut something.
bool cutByPlane(Polyhedron &poly, const double3 &n, double b, const PolyFace &source, Workspace &ws)
{
  const double tolerance = 1.0e-10 * std::max(1.0, std::fabs(b));
  const std::size_t oldVertexCount = poly.vertices.size();

  ws.distance.resize(oldVertexCount);
  ws.side.resize(oldVertexCount);
  bool anyOutside = false;
  bool anyInside = false;
  for (std::size_t i = 0; i < oldVertexCount; ++i)
  {
    ws.distance[i] = double3::dot(poly.vertices[i], n) - b;
    ws.side[i] = ws.distance[i] > tolerance ? 1 : (ws.distance[i] < -tolerance ? -1 : 0);
    anyOutside = anyOutside || (ws.side[i] > 0);
    anyInside = anyInside || (ws.side[i] < 0);
  }
  if (!anyOutside) return false;
  if (!anyInside)
  {
    // Every remaining vertex is on the far side of the plane: the radical plane of a
    // heavier neighbour has clipped the cell away entirely. This is a hidden site (its
    // power cell is empty), exactly as voro++ reports for a fully engulfed atom. Signal
    // the empty cell by clearing the polyhedron.
    poly.vertices.clear();
    poly.faces.clear();
    return true;
  }

  // Clip every face polygon in place; interpolated vertices are cached per polyhedron
  // edge (linear scan — a cut crosses only a handful of edges) so adjacent faces share
  // them.
  ws.edgeCuts.clear();
  for (PolyFace &face : poly.faces)
  {
    std::size_t count = face.vertexIndices.size();
    bool touched = false;
    for (std::size_t k = 0; k < count && !touched; ++k) touched = (ws.side[face.vertexIndices[k]] >= 0);
    if (!touched) continue;  // face strictly inside; unchanged

    ws.clipped.clear();
    for (std::size_t k = 0; k < count; ++k)
    {
      std::size_t a = face.vertexIndices[k];
      std::size_t b2 = face.vertexIndices[(k + 1) % count];
      if (ws.side[a] <= 0) ws.clipped.push_back(a);
      if (static_cast<int>(ws.side[a]) * static_cast<int>(ws.side[b2]) < 0)
      {
        std::size_t lo = std::min(a, b2), hi = std::max(a, b2);
        std::size_t newIndex = std::numeric_limits<std::size_t>::max();
        for (const Workspace::EdgeCut &cut : ws.edgeCuts)
        {
          if (cut.a == lo && cut.b == hi)
          {
            newIndex = cut.index;
            break;
          }
        }
        if (newIndex == std::numeric_limits<std::size_t>::max())
        {
          double t = ws.distance[a] / (ws.distance[a] - ws.distance[b2]);
          poly.vertices.push_back(poly.vertices[a] + t * (poly.vertices[b2] - poly.vertices[a]));
          newIndex = poly.vertices.size() - 1;
          ws.edgeCuts.push_back({lo, hi, newIndex});
        }
        ws.clipped.push_back(newIndex);
      }
    }
    face.vertexIndices.assign(ws.clipped.begin(), ws.clipped.end());
  }
  std::erase_if(poly.faces, [](const PolyFace &face) { return face.vertexIndices.size() < 3; });

  // The cross-section consists of the interpolated vertices plus any original vertices
  // lying exactly on the plane.
  ws.loop.clear();
  for (const Workspace::EdgeCut &cut : ws.edgeCuts) ws.loop.push_back({0.0, cut.index});
  for (std::size_t i = 0; i < oldVertexCount; ++i)
  {
    if (ws.side[i] == 0) ws.loop.push_back({0.0, i});
  }

  if (ws.loop.size() >= 3)
  {
    // Order the cross-section right-handed around the outward normal n. The cross-section
    // of a convex polyhedron is a convex polygon, so an angular sort around its centroid
    // is sufficient. A monotone pseudo-angle (quadrant + normalized slope) replaces
    // atan2; it induces the same ordering at a fraction of the cost.
    double3 axis = double3::normalize(n);
    double3 helper = std::fabs(axis.x) < 0.9 ? double3(1.0, 0.0, 0.0) : double3(0.0, 1.0, 0.0);
    double3 u = double3::normalize(double3::cross(helper, axis));
    double3 v = double3::cross(axis, u);  // u × v = axis

    double3 center(0.0, 0.0, 0.0);
    for (const auto &[angle, index] : ws.loop) center += poly.vertices[index];
    center /= static_cast<double>(ws.loop.size());

    for (auto &[angle, index] : ws.loop)
    {
      double3 p = poly.vertices[index] - center;
      double x = double3::dot(p, u);
      double y = double3::dot(p, v);
      double r = std::fabs(x) + std::fabs(y);
      double q = (r > 0.0) ? y / r : 0.0;               // diamond angle in [-1,1]
      angle = (x >= 0.0) ? q : ((q >= 0.0) ? 2.0 - q : -2.0 - q);  // monotone in [-2,2)
    }
    std::sort(ws.loop.begin(), ws.loop.end());

    PolyFace closing = source;
    closing.vertexIndices.reserve(ws.loop.size());
    for (const auto &[angle, index] : ws.loop) closing.vertexIndices.push_back(index);
    poly.faces.push_back(std::move(closing));
  }

  compactVertices(poly, ws);
  return true;
}

double maximumRadiusSquared(const Polyhedron &poly)
{
  double maximum = 0.0;
  for (const double3 &x : poly.vertices)
  {
    maximum = std::max(maximum, x.length_squared());
  }
  return maximum;
}
}  // namespace

SKVoronoi::SKVoronoi(const double3x3 &unitCell, const std::vector<double3> &fractionalPositions,
                     const std::vector<double> &radii)
    : _unitCell(unitCell)
{
  double determinant = unitCell.determinant();
  if (std::fabs(determinant) < 1.0e-12)
  {
    throw std::runtime_error("SKVoronoi: singular unit cell");
  }
  if (fractionalPositions.empty())
  {
    throw std::runtime_error("SKVoronoi: no positions given");
  }
  if (!radii.empty() && radii.size() != fractionalPositions.size())
  {
    throw std::runtime_error("SKVoronoi: radii size does not match positions size");
  }

  // Squared radii; the unweighted diagram is the special case of all radii zero.
  _radiiSquared.assign(fractionalPositions.size(), 0.0);
  if (!radii.empty())
  {
    for (std::size_t i = 0; i < radii.size(); ++i) _radiiSquared[i] = radii[i] * radii[i];
    _maximumRadiusSquared = *std::max_element(_radiiSquared.begin(), _radiiSquared.end());
  }

  _inverseUnitCell = unitCell.inverse();
  _metricTensor = unitCell.transpose() * unitCell;
  _cellVolume = std::fabs(determinant);

  double3 a = double3(unitCell[0].x, unitCell[0].y, unitCell[0].z);
  double3 b = double3(unitCell[1].x, unitCell[1].y, unitCell[1].z);
  double3 c = double3(unitCell[2].x, unitCell[2].y, unitCell[2].z);
  _perpendicularWidths = double3(_cellVolume / double3::cross(b, c).length(),
                                 _cellVolume / double3::cross(c, a).length(),
                                 _cellVolume / double3::cross(a, b).length());

  _positions.reserve(fractionalPositions.size());
  _cartesianPositions.reserve(fractionalPositions.size());
  for (const double3 &position : fractionalPositions)
  {
    _positions.push_back(double3::fract(position));
    _cartesianPositions.push_back(_unitCell * _positions.back());
  }

  // Size the cell-list grid to roughly four sites per bin, with bin counts proportional
  // to the perpendicular widths so that the bins are approximately metrically cubic.
  double targetBinSize = std::cbrt(_cellVolume / std::max(1.0, static_cast<double>(_positions.size()) / 4.0));
  _gridSize = int3(std::max(1, static_cast<int>(_perpendicularWidths.x / targetBinSize)),
                   std::max(1, static_cast<int>(_perpendicularWidths.y / targetBinSize)),
                   std::max(1, static_cast<int>(_perpendicularWidths.z / targetBinSize)));
  _minimumBinWidth = std::min({_perpendicularWidths.x / static_cast<double>(_gridSize.x),
                               _perpendicularWidths.y / static_cast<double>(_gridSize.y),
                               _perpendicularWidths.z / static_cast<double>(_gridSize.z)});

  _bins.assign(static_cast<std::size_t>(_gridSize.x) * static_cast<std::size_t>(_gridSize.y) *
                   static_cast<std::size_t>(_gridSize.z),
               {});
  for (std::size_t i = 0; i < _positions.size(); ++i)
  {
    int bx = std::min(_gridSize.x - 1, static_cast<int>(_positions[i].x * static_cast<double>(_gridSize.x)));
    int by = std::min(_gridSize.y - 1, static_cast<int>(_positions[i].y * static_cast<double>(_gridSize.y)));
    int bz = std::min(_gridSize.z - 1, static_cast<int>(_positions[i].z * static_cast<double>(_gridSize.z)));
    _bins[static_cast<std::size_t>((bz * _gridSize.y + by) * _gridSize.x + bx)].push_back(i);
  }

  // The Voronoi cell of any site is contained in the Voronoi cell of the lattice itself,
  // which lies within a ball of radius half the longest body diagonal (the covering
  // radius bound); use that as the half-width of the initial bounding cube.
  double maximumBodyDiagonal = 0.0;
  for (double sx : {-1.0, 1.0})
  {
    for (double sy : {-1.0, 1.0})
    {
      maximumBodyDiagonal = std::max(maximumBodyDiagonal, (_unitCell * double3(sx, sy, 1.0)).length());
    }
  }
  // A larger radius pushes this site's faces outward by up to ½(rᵢ²−rⱼ²)/|Δ|, so the power
  // cell of the biggest atom can extend beyond the geometric bound; add that margin.
  _initialHalfWidth = 0.5 * maximumBodyDiagonal + 1.0 + std::sqrt(_maximumRadiusSquared);
}

SKVoronoiCell SKVoronoi::computeCell(std::size_t siteIndex) const
{
  const double3 site = _positions[siteIndex];
  const double3 siteCartesian = _cartesianPositions[siteIndex];

  const int3 siteBin(std::min(_gridSize.x - 1, static_cast<int>(site.x * static_cast<double>(_gridSize.x))),
                     std::min(_gridSize.y - 1, static_cast<int>(site.y * static_cast<double>(_gridSize.y))),
                     std::min(_gridSize.z - 1, static_cast<int>(site.z * static_cast<double>(_gridSize.z))));

  // Reused across calls (per thread) so that steady-state construction does not allocate.
  static thread_local Workspace ws;

  const double siteRadiusSquared = _radiiSquared[siteIndex];

  Polyhedron poly = makeInitialCube(_initialHalfWidth);
  double maxRadiusSquared = maximumRadiusSquared(poly);

  // Farthest a cutting site can be and still trim the current cell. For the unweighted
  // diagram a bisector at Cartesian distance D sits at D/2, so D must be below twice the
  // circumradius R. With radical planes the plane between this site and a neighbour of
  // radius rⱼ sits at ½(D + (rᵢ²−rⱼ²)/D); a heavier neighbour pulls it inward, so the
  // search must extend to the D where even the heaviest possible neighbour's plane just
  // reaches R: D = R + sqrt(R² + max(0, rmax² − rᵢ²)). This reduces to 2R when unweighted.
  auto searchRadiusSquaredOf = [&](double circumRadiusSquared)
  {
    double c = std::max(0.0, _maximumRadiusSquared - siteRadiusSquared);
    double s = std::sqrt(circumRadiusSquared) + std::sqrt(circumRadiusSquared + c);
    return s * s;
  };
  double searchRadiusSquared = searchRadiusSquaredOf(maxRadiusSquared);

  // Walk bin shells outward: shell k holds all bins with Chebyshev offset k from the
  // site's bin. Offsets beyond the grid size wrap around and carry a periodic image, so
  // arbitrary images are enumerated without a preset range. Every bin in shell k has all
  // of its sites at Cartesian distance at least (k-1)·(minimum bin width) from the site,
  // so the walk can stop once that bound exceeds the search radius above.
  auto binAndImage = [this](int coordinate, int gridExtent) -> std::pair<int, int>
  {
    int image = (coordinate >= 0) ? coordinate / gridExtent : -((-coordinate + gridExtent - 1) / gridExtent);
    return {coordinate - image * gridExtent, image};
  };

  std::vector<Candidate> &shellCandidates = ws.shellCandidates;
  for (int k = 0;; ++k)
  {
    double lowerBound = static_cast<double>(k - 1) * _minimumBinWidth;
    if (k > 0 && lowerBound * lowerBound > searchRadiusSquared) break;

    shellCandidates.clear();
    for (int ox = -k; ox <= k; ++ox)
    {
      for (int oy = -k; oy <= k; ++oy)
      {
        for (int oz = -k; oz <= k; ++oz)
        {
          if (std::max({std::abs(ox), std::abs(oy), std::abs(oz)}) != k) continue;

          auto [bx, lx] = binAndImage(siteBin.x + ox, _gridSize.x);
          auto [by, ly] = binAndImage(siteBin.y + oy, _gridSize.y);
          auto [bz, lz] = binAndImage(siteBin.z + oz, _gridSize.z);

          // One image-offset transform per bin; candidates then need plain dot products.
          double3 imageShift =
              _unitCell * double3(static_cast<double>(lx), static_cast<double>(ly), static_cast<double>(lz)) -
              siteCartesian;

          for (std::size_t j : _bins[static_cast<std::size_t>((bz * _gridSize.y + by) * _gridSize.x + bx)])
          {
            if (j == siteIndex && lx == 0 && ly == 0 && lz == 0) continue;
            double3 delta = _cartesianPositions[j] + imageShift;
            double rsq = delta.length_squared();
            if (rsq < 1.0e-16)
            {
              throw std::runtime_error("SKVoronoi: coinciding sites");
            }
            shellCandidates.push_back({rsq, delta, j, int3(lx, ly, lz)});
          }
        }
      }
    }
    std::sort(shellCandidates.begin(), shellCandidates.end(),
              [](const Candidate &lhs, const Candidate &rhs) { return lhs.rsq < rhs.rsq; });

    for (const Candidate &candidate : shellCandidates)
    {
      // Candidates in this shell are sorted by distance; once one is beyond the search
      // radius so are the rest (the circumradius only shrinks as cuts proceed).
      if (candidate.rsq > searchRadiusSquared) break;

      PolyFace source;
      source.hasSource = true;
      source.neighborIndex = candidate.index;
      source.neighborImage = candidate.image;
      source.rsq = candidate.rsq;
      source.delta = candidate.delta;

      // Radical (power) plane offset; ½|Δ|² for the unweighted case (all radii zero).
      double offset = 0.5 * (candidate.rsq + siteRadiusSquared - _radiiSquared[candidate.index]);

      if (cutByPlane(poly, candidate.delta, offset, source, ws))
      {
        if (poly.faces.empty()) break;  // hidden site: radical cell is empty
        maxRadiusSquared = maximumRadiusSquared(poly);
        searchRadiusSquared = searchRadiusSquaredOf(maxRadiusSquared);
      }
    }
    if (poly.faces.empty()) break;
  }

  // A hidden site (engulfed by a heavier neighbour in the power metric) has an empty cell,
  // exactly as voro++ reports; return a degenerate zero-volume cell for it.
  if (poly.vertices.empty())
  {
    SKVoronoiCell cell;
    cell.siteIndex = siteIndex;
    cell.sitePositionFractional = site;
    cell.volume = 0.0;
    cell.centroidCartesian = double3(0.0, 0.0, 0.0);
    return cell;
  }

  // The initial cube provably contains the cell, so all of its faces must have been cut.
  if (std::any_of(poly.faces.begin(), poly.faces.end(), [](const PolyFace &face) { return !face.hasSource; }))
  {
    throw std::runtime_error("SKVoronoi: cell construction failed (initial bounding cube face survived)");
  }

  // Assemble the result; fractional vertices are recovered through the inverse cell matrix.
  SKVoronoiCell cell;
  cell.siteIndex = siteIndex;
  cell.sitePositionFractional = site;
  cell.verticesCartesian = std::move(poly.vertices);
  cell.verticesFractional.reserve(cell.verticesCartesian.size());
  for (const double3 &x : cell.verticesCartesian)
  {
    cell.verticesFractional.push_back(_inverseUnitCell * x);
  }

  cell.faces.reserve(poly.faces.size());
  for (PolyFace &face : poly.faces)
  {
    SKVoronoiFace outputFace;
    outputFace.neighborIndex = face.neighborIndex;
    outputFace.neighborImage = face.neighborImage;
    outputFace.neighborDistance = std::sqrt(face.rsq);
    outputFace.vertexIndices = std::move(face.vertexIndices);
    outputFace.normalCartesian = double3::normalize(face.delta);

    double3 areaVector(0.0, 0.0, 0.0);
    const double3 &p0 = cell.verticesCartesian[outputFace.vertexIndices[0]];
    for (std::size_t k = 1; k + 1 < outputFace.vertexIndices.size(); ++k)
    {
      const double3 &p1 = cell.verticesCartesian[outputFace.vertexIndices[k]];
      const double3 &p2 = cell.verticesCartesian[outputFace.vertexIndices[k + 1]];
      areaVector += double3::cross(p1 - p0, p2 - p0);
    }
    outputFace.area = 0.5 * areaVector.length();

    cell.faces.push_back(std::move(outputFace));
  }

  // Volume and centroid from the fan of tetrahedra spanned by the site (which lies strictly
  // inside the cell) and the triangulated faces.
  double volume = 0.0;
  double3 centroid(0.0, 0.0, 0.0);
  for (const SKVoronoiFace &face : cell.faces)
  {
    double signedVolume = 0.0;
    double3 signedCentroid(0.0, 0.0, 0.0);
    const double3 &p0 = cell.verticesCartesian[face.vertexIndices[0]];
    for (std::size_t k = 1; k + 1 < face.vertexIndices.size(); ++k)
    {
      const double3 &p1 = cell.verticesCartesian[face.vertexIndices[k]];
      const double3 &p2 = cell.verticesCartesian[face.vertexIndices[k + 1]];
      double tetrahedronVolume = double3::dot(p0, double3::cross(p1, p2)) / 6.0;
      signedVolume += tetrahedronVolume;
      signedCentroid += tetrahedronVolume * 0.25 * (p0 + p1 + p2);
    }
    if (signedVolume < 0.0)
    {
      signedVolume = -signedVolume;
      signedCentroid = -signedCentroid;
    }
    volume += signedVolume;
    centroid += signedCentroid;
  }
  cell.volume = volume;
  cell.centroidCartesian = centroid / volume;

  return cell;
}

std::vector<SKVoronoiCell> SKVoronoi::computeAllCells() const
{
  std::vector<SKVoronoiCell> cells;
  cells.reserve(_positions.size());
  for (std::size_t i = 0; i < _positions.size(); ++i)
  {
    cells.push_back(computeCell(i));
  }
  return cells;
}
