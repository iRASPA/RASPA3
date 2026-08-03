module;

module voronoi_network;

import std;

import int3;
import double3;
import double3x3;
import unit_cell;
import skvoronoi;

constexpr double mergeTolerance = 0.02;  // Å; zeo++ VOR_NODE_MERGE_THRESHOLD

// Squared distance from the origin to the segment [a, b].
double distanceSquaredOriginToSegment(const double3& a, const double3& b)
{
  double3 ab = b - a;
  double denominator = double3::dot(ab, ab);
  double t = (denominator > 0.0) ? std::clamp(-double3::dot(a, ab) / denominator, 0.0, 1.0) : 0.0;
  double3 closest = a + t * ab;
  return double3::dot(closest, closest);
}

// At most two spheres are tangent to four, so the answer fits in a pair of them. Which matters
// because the walk below solves this millions of times and a vector would allocate at every one.
struct TangentSpheres
{
  std::array<ApolloniusSphere, 2> spheres{};
  std::size_t count{0};

  void add(const ApolloniusSphere& sphere) { spheres[count++] = sphere; }
  const ApolloniusSphere* begin() const { return spheres.data(); }
  const ApolloniusSphere* end() const { return spheres.data() + count; }
};

// `weights` are the |x_i|^2 - r_i^2 of the four, which depend on one sphere each and not on the
// quadruple, so the walk below works them out once per atom instead of once per quadruple it appears in.
TangentSpheres tangentSpheresOf(const std::array<double3, 4>& centres, const std::array<double, 4>& radii,
                                const std::array<double, 4>& weights)
{
  // A sphere of radius t centred at c that touches all four from outside satisfies
  // |c - x_i| = r_i + t. Squaring gives |c|^2 - 2 c.x_i + |x_i|^2 = t^2 + 2 r_i t + r_i^2, and
  // subtracting the i = 0 equation cancels both |c|^2 and t^2, leaving three equations that are
  // linear in c and t:
  //     2 (x_i - x_0) . c + 2 (r_i - r_0) t = (|x_i|^2 - r_i^2) - (|x_0|^2 - r_0^2)
  double3 rows[3]{2.0 * (centres[1] - centres[0]), 2.0 * (centres[2] - centres[0]),
                  2.0 * (centres[3] - centres[0])};
  double3 rightHandSide(weights[1] - weights[0], weights[2] - weights[0], weights[3] - weights[0]);
  double3 timeCoefficient(2.0 * (radii[1] - radii[0]), 2.0 * (radii[2] - radii[0]), 2.0 * (radii[3] - radii[0]));

  // Both unknowns are solved from the same three rows, by Cramer's rule written with cross products.
  // The two right-hand sides then share the three cofactor vectors and the one determinant, and the
  // whole solve costs a single division; forming the inverse matrix instead would compute the
  // determinant a second time and divide nine times, which the walk below cannot afford to do a
  // million times over.
  double3 cofactor0 = double3::cross(rows[1], rows[2]);
  double3 cofactor1 = double3::cross(rows[2], rows[0]);
  double3 cofactor2 = double3::cross(rows[0], rows[1]);
  double determinant = double3::dot(rows[0], cofactor0);

  // Coplanar or coincident centres leave the centre undetermined along a direction. The test is
  // |det| < 1e-12 s^3 for s the longest row, squared throughout so that neither the row lengths nor
  // the scale need a root taken of them.
  double scaleSquared =
      std::max({rows[0].length_squared(), rows[1].length_squared(), rows[2].length_squared()});
  if (scaleSquared <= 0.0 ||
      determinant * determinant < 1.0e-24 * scaleSquared * scaleSquared * scaleSquared)
    return {};

  double inverseDeterminant = 1.0 / determinant;
  // centre at t = 0
  double3 base = inverseDeterminant * (rightHandSide.x * cofactor0 + rightHandSide.y * cofactor1 +
                                       rightHandSide.z * cofactor2);
  // centre moves by -t * direction
  double3 direction = inverseDeterminant * (timeCoefficient.x * cofactor0 + timeCoefficient.y * cofactor1 +
                                            timeCoefficient.z * cofactor2);

  // Substituting c = base - t * direction into |c - x_0| = r_0 + t gives a quadratic in t.
  double3 offset = base - centres[0];
  double a = double3::dot(direction, direction) - 1.0;
  double b = -2.0 * (double3::dot(offset, direction) + radii[0]);
  double c = double3::dot(offset, offset) - radii[0] * radii[0];

  std::array<double, 2> roots{};
  std::size_t rootCount = 0;
  if (std::abs(a) < 1.0e-14)
  {
    if (std::abs(b) > 1.0e-14) roots[rootCount++] = -c / b;
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
  for (std::size_t r = 0; r < rootCount; ++r)
  {
    double t = roots[r];
    if (t < 0.0) continue;
    spheres.add(ApolloniusSphere{base - t * direction, t});
  }
  return spheres;
}

std::vector<ApolloniusSphere> apolloniusTangentSpheres(const std::array<double3, 4>& centres,
                                                       const std::array<double, 4>& radii)
{
  std::array<double, 4> weights{};
  for (std::size_t i = 0; i < 4; ++i) weights[i] = double3::dot(centres[i], centres[i]) - radii[i] * radii[i];

  TangentSpheres spheres = tangentSpheresOf(centres, radii, weights);
  return std::vector<ApolloniusSphere>(spheres.begin(), spheres.end());
}

double VoronoiNetwork::largestIncludedSphereDiameter() const
{
  double maximum = 0.0;
  for (const VoronoiNode& node : nodes) maximum = std::max(maximum, node.maximalRadius);
  return 2.0 * maximum;
}

VoronoiNetwork VoronoiNetwork::create(const UnitCell& unitCell,
                                      const std::vector<double3>& fractionalPositions,
                                      const std::vector<double>& radii)
{
  VoronoiNetwork network;
  network.unitCell = unitCell;
  network.atomRadii = radii;

  const double3x3 cell = unitCell.cell;
  const double3x3 inverseCell = unitCell.inverseCell;

  // Radical (power) Voronoi diagram: the tessellation depends on the atom radii, matching
  // zeo++/voro++ semantics for radii-dependent pore analysis.
  SKVoronoi voronoi(cell, fractionalPositions, radii);
  network.atomPositionsFractional = voronoi.positions();
  const std::size_t numberOfAtoms = network.atomPositionsFractional.size();
  network.atomNodeVectors.assign(numberOfAtoms, {});

  // Tolerance-based vertex merging, matching zeo++/voro++ voronoi_network semantics: a
  // vertex is mapped onto the first stored node whose position differs by less than
  // mergeTolerance in every Cartesian component (minimum image). Grid quantisation is not
  // good enough here: two vertices within tolerance can land in different grid cells,
  // splitting one physical node in two. A split node drops the constricting atom from some
  // of its edges, which overestimates bottleneck radii (Df) and node radii (Di).
  double volume = unitCell.volume;
  double3 perpendicularWidths(volume / double3::cross(cell[1], cell[2]).length(),
                              volume / double3::cross(cell[2], cell[0]).length(),
                              volume / double3::cross(cell[0], cell[1]).length());
  // Bin width of ~1 Å (always > mergeTolerance) keeps the ±1-bin neighbour search exact
  // while the bin count stays proportional to the cell volume.
  int3 binCount(std::max(1, static_cast<std::int32_t>(perpendicularWidths.x)),
                std::max(1, static_cast<std::int32_t>(perpendicularWidths.y)),
                std::max(1, static_cast<std::int32_t>(perpendicularWidths.z)));
  std::vector<std::vector<std::size_t>> nodeBins(
      static_cast<std::size_t>(binCount.x) * static_cast<std::size_t>(binCount.y) *
      static_cast<std::size_t>(binCount.z));

  auto binIndexOf = [&](const double3& wrappedFractional) -> std::size_t
  {
    std::int32_t bx = std::min(binCount.x - 1, static_cast<std::int32_t>(wrappedFractional.x * binCount.x));
    std::int32_t by = std::min(binCount.y - 1, static_cast<std::int32_t>(wrappedFractional.y * binCount.y));
    std::int32_t bz = std::min(binCount.z - 1, static_cast<std::int32_t>(wrappedFractional.z * binCount.z));
    return static_cast<std::size_t>(bx) +
           static_cast<std::size_t>(binCount.x) *
               (static_cast<std::size_t>(by) + static_cast<std::size_t>(binCount.y) * static_cast<std::size_t>(bz));
  };

  // With one or two bins along an axis the -1 and +1 offsets wrap onto bins that are already
  // covered, so the scan range is narrowed to visit each neighbouring bin exactly once.
  auto offsetLow = [](std::int32_t count) { return count >= 3 ? -1 : 0; };
  auto offsetHigh = [](std::int32_t count) { return count >= 2 ? 1 : 0; };

  auto findOrAddNode = [&](const double3& unwrappedFractional, double contributionRadius,
                           std::size_t atomIndex) -> std::size_t
  {
    double3 wrapped = double3::fract(unwrappedFractional);
    std::int32_t bx = std::min(binCount.x - 1, static_cast<std::int32_t>(wrapped.x * binCount.x));
    std::int32_t by = std::min(binCount.y - 1, static_cast<std::int32_t>(wrapped.y * binCount.y));
    std::int32_t bz = std::min(binCount.z - 1, static_cast<std::int32_t>(wrapped.z * binCount.z));

    for (std::int32_t dz = offsetLow(binCount.z); dz <= offsetHigh(binCount.z); ++dz)
      for (std::int32_t dy = offsetLow(binCount.y); dy <= offsetHigh(binCount.y); ++dy)
        for (std::int32_t dx = offsetLow(binCount.x); dx <= offsetHigh(binCount.x); ++dx)
        {
          std::int32_t nx = (bx + dx + binCount.x) % binCount.x;
          std::int32_t ny = (by + dy + binCount.y) % binCount.y;
          std::int32_t nz = (bz + dz + binCount.z) % binCount.z;
          std::size_t bin = static_cast<std::size_t>(nx) +
                            static_cast<std::size_t>(binCount.x) *
                                (static_cast<std::size_t>(ny) +
                                 static_cast<std::size_t>(binCount.y) * static_cast<std::size_t>(nz));
          for (std::size_t candidate : nodeBins[bin])
          {
            double3 difference = network.nodes[candidate].fractional - wrapped;
            difference.x -= std::round(difference.x);
            difference.y -= std::round(difference.y);
            difference.z -= std::round(difference.z);
            double3 cartesian = cell * difference;
            if (std::abs(cartesian.x) < mergeTolerance && std::abs(cartesian.y) < mergeTolerance &&
                std::abs(cartesian.z) < mergeTolerance)
            {
              VoronoiNode& node = network.nodes[candidate];
              node.radius = std::min(node.radius, contributionRadius);
              if (std::find(node.atomIndices.begin(), node.atomIndices.end(), atomIndex) == node.atomIndices.end())
              {
                node.atomIndices.push_back(atomIndex);
              }
              return candidate;
            }
          }
        }

    std::size_t index = network.nodes.size();
    VoronoiNode node;
    node.fractional = wrapped;
    node.position = cell * wrapped;
    node.radius = contributionRadius;
    node.atomIndices = {atomIndex};
    network.nodes.push_back(std::move(node));
    nodeBins[binIndexOf(wrapped)].push_back(index);
    return index;
  };

  // Undirected edges keyed by (min node, max node, oriented lattice shift). Every cell that
  // shares an edge registers its own bottleneck contribution, so entries are collected into a
  // flat vector and reduced afterwards. Keeping the shift as three separate key components
  // avoids packing them into one integer, and sorting once is cheaper than maintaining a tree
  // with an allocation per edge.
  struct EdgeEntry
  {
    std::array<std::int64_t, 5> key;  // from node, to node, and the three lattice shift components
    double radius;
    double length;
  };
  std::vector<EdgeEntry> edgeEntries;

  std::vector<SKVoronoiCell> cells = voronoi.computeAllCells();

  std::size_t numberOfFaceEdges = 0;
  for (const SKVoronoiCell& polyhedron : cells)
    for (const SKVoronoiFace& face : polyhedron.faces) numberOfFaceEdges += face.vertexIndices.size();
  edgeEntries.reserve(numberOfFaceEdges);

  for (std::size_t i = 0; i < cells.size(); ++i)
  {
    const SKVoronoiCell& polyhedron = cells[i];
    double3 siteFractional = network.atomPositionsFractional[i];
    double radiusI = radii[i];

    // Fractional coordinate of each vertex (unwrapped) and its node index.
    std::vector<double3> vertexFractional(polyhedron.verticesCartesian.size());
    std::vector<std::size_t> vertexNode(polyhedron.verticesCartesian.size());
    for (std::size_t v = 0; v < polyhedron.verticesCartesian.size(); ++v)
    {
      const double3& relative = polyhedron.verticesCartesian[v];  // relative to the site
      double distanceToAtom = relative.length();
      vertexFractional[v] = siteFractional + inverseCell * relative;
      vertexNode[v] = findOrAddNode(vertexFractional[v], distanceToAtom - radiusI, i);

      network.atomNodeVectors[i].push_back({vertexNode[v], relative});
    }

    for (const SKVoronoiFace& face : polyhedron.faces)
    {
      std::size_t count = face.vertexIndices.size();
      for (std::size_t k = 0; k < count; ++k)
      {
        std::size_t p = face.vertexIndices[k];
        std::size_t q = face.vertexIndices[(k + 1) % count];

        double3 fractionalP = vertexFractional[p];
        double3 fractionalQ = vertexFractional[q];

        // Lattice shift between the two endpoint nodes along this edge. Using the merged
        // node positions (rather than floor() of each cell's unwrapped vertex coordinate)
        // makes the shift image-invariant, so all cells that share a Voronoi edge register
        // it under the same (node,node,shift) key. That is essential: the edge's bottleneck
        // radius is the minimum contribution over those cells, and a split key would drop
        // the constricting atom and overestimate the throat.
        double3 nodeA = network.nodes[vertexNode[p]].fractional;
        double3 nodeB = network.nodes[vertexNode[q]].fractional;
        double3 shift = (fractionalQ - fractionalP) - (nodeB - nodeA);
        int3 delta(static_cast<int>(std::lround(shift.x)), static_cast<int>(std::lround(shift.y)),
                   static_cast<int>(std::lround(shift.z)));

        double3 relativeP = polyhedron.verticesCartesian[p];
        double3 relativeQ = polyhedron.verticesCartesian[q];
        double length = (relativeQ - relativeP).length();

        // Drop degenerate self-edges. A self-loop with a non-zero lattice shift is a
        // legitimate periodic edge (e.g. every edge of a single-atom lattice), but only if
        // its two endpoints are geometrically distinct: when they merge to the same node
        // yet straddle a cell boundary the segment is (near) zero-length, which would
        // otherwise fake a wide percolating channel with bottleneck equal to the node
        // radius. Skip both the zero-shift case and the zero-length case.
        if (vertexNode[p] == vertexNode[q] &&
            ((delta.x == 0 && delta.y == 0 && delta.z == 0) || length < mergeTolerance))
          continue;

        // Bottleneck contribution of atom i: closest approach of the edge segment to the
        // atom centre (segment expressed relative to the atom).
        double bottleneck = std::sqrt(distanceSquaredOriginToSegment(relativeP, relativeQ)) - radiusI;

        std::size_t a = vertexNode[p];
        std::size_t b = vertexNode[q];
        int3 orientedDelta = delta;
        if (a > b)
        {
          std::swap(a, b);
          orientedDelta = int3(-delta.x, -delta.y, -delta.z);
        }
        std::array<std::int64_t, 5> edgeKey{
            static_cast<std::int64_t>(a),               static_cast<std::int64_t>(b),
            static_cast<std::int64_t>(orientedDelta.x), static_cast<std::int64_t>(orientedDelta.y),
            static_cast<std::int64_t>(orientedDelta.z)};
        edgeEntries.push_back(EdgeEntry{edgeKey, bottleneck, length});
      }
    }
  }

  // A stable sort groups the contributions to each edge while leaving the first-registered one
  // at the front of its run, so the retained length is the one insertion order would have kept.
  std::stable_sort(edgeEntries.begin(), edgeEntries.end(),
                   [](const EdgeEntry& lhs, const EdgeEntry& rhs) { return lhs.key < rhs.key; });

  std::size_t numberOfUniqueEdges = 0;
  for (std::size_t i = 0; i < edgeEntries.size(); ++i)
    if (i == 0 || edgeEntries[i].key != edgeEntries[i - 1].key) ++numberOfUniqueEdges;
  network.edges.reserve(2 * numberOfUniqueEdges);

  for (std::size_t start = 0; start < edgeEntries.size();)
  {
    // The bottleneck is the tightest contribution over the cells that share the edge.
    std::size_t end = start;
    double radius = edgeEntries[start].radius;
    while (end < edgeEntries.size() && edgeEntries[end].key == edgeEntries[start].key)
    {
      radius = std::min(radius, edgeEntries[end].radius);
      ++end;
    }

    const std::array<std::int64_t, 5>& key = edgeEntries[start].key;
    std::size_t a = static_cast<std::size_t>(key[0]);
    std::size_t b = static_cast<std::size_t>(key[1]);
    int3 delta(static_cast<std::int32_t>(key[2]), static_cast<std::int32_t>(key[3]),
               static_cast<std::int32_t>(key[4]));
    double length = edgeEntries[start].length;

    network.edges.push_back(VoronoiEdge{a, b, delta, radius, length});
    network.edges.push_back(VoronoiEdge{b, a, int3(-delta.x, -delta.y, -delta.z), radius, length});

    start = end;
  }

  // Every radius above is a clearance to the atoms whose cells bound that node or edge. The
  // tessellation is radical, so adjacency is decided by the power distance |x-x_j|^2 - r_j^2
  // while the room for a probe is the clearance |x-x_j| - r_j; the two disagree once the radii
  // differ, and an atom that bounds nothing here can still be the nearest in clearance. Left
  // uncorrected the radii are therefore too large, and an edge can be reported as wider than a
  // sphere can actually pass. Recomputing each radius as a minimum over all nearby atoms makes
  // the bottlenecks a rigorous lower bound: a sphere of that radius really does fit everywhere
  // along the edge. For equal radii the radical diagram is the ordinary Voronoi diagram and the
  // bounding atoms are already the nearest ones, so this pass leaves the network untouched.
  double maximumAtomRadius = 0.0;
  for (double radius : radii) maximumAtomRadius = std::max(maximumAtomRadius, radius);

  // Cell list over the atoms, sized for roughly four atoms per bin.
  double atomTargetBinSize = std::cbrt(volume / std::max(1.0, static_cast<double>(numberOfAtoms) / 4.0));
  int3 atomBinCount(std::max(1, static_cast<std::int32_t>(perpendicularWidths.x / atomTargetBinSize)),
                    std::max(1, static_cast<std::int32_t>(perpendicularWidths.y / atomTargetBinSize)),
                    std::max(1, static_cast<std::int32_t>(perpendicularWidths.z / atomTargetBinSize)));
  double3 atomBinWidth(perpendicularWidths.x / static_cast<double>(atomBinCount.x),
                       perpendicularWidths.y / static_cast<double>(atomBinCount.y),
                       perpendicularWidths.z / static_cast<double>(atomBinCount.z));

  std::vector<double3> atomPositionsCartesian(numberOfAtoms);
  std::vector<std::vector<std::size_t>> atomBins(static_cast<std::size_t>(atomBinCount.x) *
                                                 static_cast<std::size_t>(atomBinCount.y) *
                                                 static_cast<std::size_t>(atomBinCount.z));
  for (std::size_t j = 0; j < numberOfAtoms; ++j)
  {
    const double3& fractional = network.atomPositionsFractional[j];
    atomPositionsCartesian[j] = cell * fractional;
    std::int32_t bx = std::min(atomBinCount.x - 1, static_cast<std::int32_t>(fractional.x * atomBinCount.x));
    std::int32_t by = std::min(atomBinCount.y - 1, static_cast<std::int32_t>(fractional.y * atomBinCount.y));
    std::int32_t bz = std::min(atomBinCount.z - 1, static_cast<std::int32_t>(fractional.z * atomBinCount.z));
    atomBins[static_cast<std::size_t>((bz * atomBinCount.y + by) * atomBinCount.x + bx)].push_back(j);
  }

  // Splits an out-of-range bin coordinate into a wrapped index and the image it came from.
  auto binAndImage = [](std::int32_t coordinate, std::int32_t extent)
  {
    std::int32_t image = (coordinate >= 0) ? coordinate / extent : -((-coordinate + extent - 1) / extent);
    return std::pair<std::int32_t, std::int32_t>{coordinate - image * extent, image};
  };

  // Half the longest diagonal of a bin, which is the radius of the sphere about a bin's centre that
  // contains it. Together with the centre it says how near a bin can come to a point, and so which
  // bins are worth looking inside. Stepping one bin along an axis moves by a fixed vector, so a bin's
  // centre is reached by adding whole steps rather than by a coordinate transform of its own.
  double binCircumradius = 0.0;
  double3 binStepX = cell * double3(1.0 / static_cast<double>(atomBinCount.x), 0.0, 0.0);
  double3 binStepY = cell * double3(0.0, 1.0 / static_cast<double>(atomBinCount.y), 0.0);
  double3 binStepZ = cell * double3(0.0, 0.0, 1.0 / static_cast<double>(atomBinCount.z));
  for (double sx : {-0.5, 0.5})
    for (double sy : {-0.5, 0.5})
      for (double sz : {-0.5, 0.5})
        binCircumradius = std::max(binCircumradius, (sx * binStepX + sy * binStepY + sz * binStepZ).length());

  // Visits every atom image whose centre lies within searchRadius of the given wrapped point.
  //
  // The bins enumerated are a box around the point, but the region asked for is a ball, and in a box
  // wide enough to hold a ball most of the bins are outside it. So each bin is measured before it is
  // opened, by the distance from the point to its centre less its circumradius, which is how near
  // anything in it can be. This is measured on the centre and circumradius rather than per axis
  // because the axes of a triclinic cell are not perpendicular, and a bound taken along them
  // separately would not bound the distance.
  //
  // That measurement has to be cheaper than the atoms it saves, or the box is merely paid for twice
  // over. So it is a squared length against a squared reach, with no root, and the centre it needs is
  // built by accumulating the axis steps down the loop nest, with no transform and nothing else per
  // bin. The wrapped bin index and image are only worked out for the bins that survive it.
  auto forEachAtomNear = [&](const double3& wrappedCentre, double searchRadius, auto&& visit)
  {
    double3 centreFractional = inverseCell * wrappedCentre;
    std::int32_t cx = std::min(atomBinCount.x - 1, static_cast<std::int32_t>(centreFractional.x * atomBinCount.x));
    std::int32_t cy = std::min(atomBinCount.y - 1, static_cast<std::int32_t>(centreFractional.y * atomBinCount.y));
    std::int32_t cz = std::min(atomBinCount.z - 1, static_cast<std::int32_t>(centreFractional.z * atomBinCount.z));

    // A displacement of searchRadius spans at most ceil(searchRadius / bin width) bins; the
    // extra bin absorbs the position of the centre within its own bin.
    int3 span(static_cast<std::int32_t>(std::ceil(searchRadius / atomBinWidth.x)) + 1,
              static_cast<std::int32_t>(std::ceil(searchRadius / atomBinWidth.y)) + 1,
              static_cast<std::int32_t>(std::ceil(searchRadius / atomBinWidth.z)) + 1);
    double binReach = searchRadius + binCircumradius;
    double binReachSquared = binReach * binReach;

    // The centre of the point's own bin, as seen from the point.
    double3 ownBinCentre = (static_cast<double>(cx) + 0.5) * binStepX + (static_cast<double>(cy) + 0.5) * binStepY +
                           (static_cast<double>(cz) + 0.5) * binStepZ - wrappedCentre;

    for (std::int32_t oz = -span.z; oz <= span.z; ++oz)
    {
      double3 alongZ = ownBinCentre + static_cast<double>(oz) * binStepZ;
      for (std::int32_t oy = -span.y; oy <= span.y; ++oy)
      {
        double3 alongZY = alongZ + static_cast<double>(oy) * binStepY;
        for (std::int32_t ox = -span.x; ox <= span.x; ++ox)
        {
          double3 binCentre = alongZY + static_cast<double>(ox) * binStepX;
          if (double3::dot(binCentre, binCentre) > binReachSquared) continue;

          auto [bx, lx] = binAndImage(cx + ox, atomBinCount.x);
          auto [by, ly] = binAndImage(cy + oy, atomBinCount.y);
          auto [bz, lz] = binAndImage(cz + oz, atomBinCount.z);
          std::size_t bin = static_cast<std::size_t>((bz * atomBinCount.y + by) * atomBinCount.x + bx);
          if (atomBins[bin].empty()) continue;

          double3 imageShift =
              cell * double3(static_cast<double>(lx), static_cast<double>(ly), static_cast<double>(lz));
          for (std::size_t j : atomBins[bin])
          {
            visit(j, atomPositionsCartesian[j] + imageShift);
          }
        }
      }
    }
  };

  // The stored radius is an upper bound on the true clearance, so only atoms within
  // (radius + largest atom radius) of the feature can lower it.
  for (VoronoiNode& node : network.nodes)
  {
    double clearance = node.radius;
    forEachAtomNear(node.position, node.radius + maximumAtomRadius,
                    [&](std::size_t j, const double3& atomImage)
                    { clearance = std::min(clearance, (node.position - atomImage).length() - radii[j]); });
    node.radius = clearance;
  }

  // Edges were emitted as consecutive forward/backward pairs sharing one radius.
  for (std::size_t i = 0; i + 1 < network.edges.size(); i += 2)
  {
    VoronoiEdge& forward = network.edges[i];
    double3 p = network.nodes[forward.from].position;
    double3 q = network.nodes[forward.to].position + cell * double3(static_cast<double>(forward.delta.x),
                                                                    static_cast<double>(forward.delta.y),
                                                                    static_cast<double>(forward.delta.z));

    // The segment can leave the cell, so shift it to sit around its own wrapped midpoint.
    double3 midpoint = 0.5 * (p + q);
    double3 wrappedMidpoint = cell * double3::fract(inverseCell * midpoint);
    double3 shift = wrappedMidpoint - midpoint;
    double halfLength = 0.5 * (q - p).length();

    double clearance = forward.radius;
    forEachAtomNear(wrappedMidpoint, halfLength + forward.radius + maximumAtomRadius,
                    [&](std::size_t j, const double3& atomImage)
                    {
                      double distance = std::sqrt(
                          distanceSquaredOriginToSegment(p + shift - atomImage, q + shift - atomImage));
                      clearance = std::min(clearance, distance - radii[j]);
                    });
    forward.radius = clearance;
    network.edges[i + 1].radius = clearance;
  }

  // The radii above are now honest, but they are measured at the radical vertices, and that is
  // not where the free space is widest. The clearance min_j(|x - x_j| - r_j) peaks at an
  // Apollonius vertex, the centre of a sphere tangent to four atoms; the two coincide only when
  // the surrounding radii are equal, which is why the correction above leaves Di short. Each node
  // therefore walks from its vertex to the peak of its own pocket, using the closed form above to
  // supply the steps.
  struct Candidate
  {
    double3 centre;
    double radius;
    double clearance;
    double weight;  // |x|^2 - r^2, which the tangency solve needs and which is the atom's alone
  };
  std::vector<Candidate> candidates;
  constexpr std::size_t maximumApolloniusIterations = 16;
  constexpr std::size_t apolloniusCandidateCount = 8;
  constexpr double apolloniusConvergenceTolerance = 1.0e-12;

  for (VoronoiNode& node : network.nodes)
  {
    node.maximalRadius = node.radius;
    node.maximalPosition = node.position;
    if (node.radius <= 0.0) continue;

    // One solve is not enough: the four atoms tightest at the radical vertex are usually not the
    // four the maximal sphere ends up touching. So walk. Each step jumps to the Apollonius vertex
    // of the currently tightest four and keeps it only if the true clearance there is larger.
    //
    // Recording the clearance rather than the tangent radius is what makes this safe: a clearance
    // is by definition the radius of a sphere that touches nothing, so no step can ever claim
    // room that is not there, and no separate emptiness test is needed. Monotonicity does the
    // rest. The walk cannot cycle because the clearance strictly increases, and it cannot wander
    // into a neighbouring pocket because reaching one means crossing a bottleneck where the
    // clearance is lower. It halts at a point no tangent solve improves, the local peak.
    double3 point = node.position;
    for (std::size_t iteration = 0; iteration < maximumApolloniusIterations; ++iteration)
    {
      candidates.clear();
      forEachAtomNear(point, node.maximalRadius + maximumAtomRadius,
                      [&](std::size_t j, const double3& atomImage)
                      {
                        candidates.push_back(Candidate{atomImage, radii[j],
                                                       (point - atomImage).length() - radii[j],
                                                       double3::dot(atomImage, atomImage) - radii[j] * radii[j]});
                      });
      if (candidates.size() < 4) break;

      // Which four atoms the peak touches is not known in advance, and simply taking the four
      // tightest stalls short of it. Try every quadruple drawn from the nearest few instead; the
      // radius test below discards almost all of them before any real work is done.
      std::size_t considered = std::min(apolloniusCandidateCount, candidates.size());
      std::partial_sort(candidates.begin(), candidates.begin() + static_cast<std::ptrdiff_t>(considered),
                        candidates.end(),
                        [](const Candidate& lhs, const Candidate& rhs) { return lhs.clearance < rhs.clearance; });

      double bestClearance = node.maximalRadius;
      double3 bestPoint = point;
      for (std::size_t i0 = 0; i0 < considered; ++i0)
        for (std::size_t i1 = i0 + 1; i1 < considered; ++i1)
          for (std::size_t i2 = i1 + 1; i2 < considered; ++i2)
            for (std::size_t i3 = i2 + 1; i3 < considered; ++i3)
            {
              std::array<double3, 4> tangentCentres{candidates[i0].centre, candidates[i1].centre,
                                                    candidates[i2].centre, candidates[i3].centre};
              std::array<double, 4> tangentRadii{candidates[i0].radius, candidates[i1].radius, candidates[i2].radius,
                                                 candidates[i3].radius};
              std::array<double, 4> tangentWeights{candidates[i0].weight, candidates[i1].weight,
                                                   candidates[i2].weight, candidates[i3].weight};

              for (const ApolloniusSphere& sphere : tangentSpheresOf(tangentCentres, tangentRadii, tangentWeights))
              {
                // The sphere touches all four, so the clearance at its centre cannot exceed its
                // radius; a shorter radius cannot improve on what we already have.
                if (sphere.radius <= bestClearance) continue;
                // Keep the jump inside the pocket rather than trusting a distant solution.
                if ((sphere.centre - point).length() > node.radius) continue;

                // Almost every sphere that gets this far is spoiled by an atom that is already in
                // hand, since the atoms tightest at the point the sphere was solved from are the
                // ones most likely to reach into it. A clearance over any set of atoms is an upper
                // bound on the clearance over all of them, so a bound that already fails to improve
                // settles the sphere without going to the grid for the rest.
                double reachable = sphere.radius;
                for (std::size_t s = 0; s < considered; ++s)
                  reachable =
                      std::min(reachable, (sphere.centre - candidates[s].centre).length() - candidates[s].radius);
                if (reachable <= bestClearance) continue;

                double3 wrappedCentre = cell * double3::fract(inverseCell * sphere.centre);
                double clearance = sphere.radius;
                forEachAtomNear(wrappedCentre, sphere.radius + maximumAtomRadius,
                                [&](std::size_t j, const double3& atomImage)
                                { clearance = std::min(clearance, (wrappedCentre - atomImage).length() - radii[j]); });

                if (clearance > bestClearance)
                {
                  bestClearance = clearance;
                  bestPoint = wrappedCentre;
                }
              }
            }

      if (bestClearance <= node.maximalRadius + apolloniusConvergenceTolerance) break;
      node.maximalRadius = bestClearance;
      node.maximalPosition = bestPoint;
      point = bestPoint;
    }
  }

  return network;
}
