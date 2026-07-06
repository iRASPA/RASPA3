module;

module voronoi_network;

import std;

import int3;
import double3;
import double3x3;
import simulationbox;
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

double VoronoiNetwork::largestIncludedSphereDiameter() const
{
  double maximum = 0.0;
  for (const VoronoiNode& node : nodes) maximum = std::max(maximum, node.radius);
  return 2.0 * maximum;
}

VoronoiNetwork VoronoiNetwork::create(const SimulationBox& simulationBox,
                                      const std::vector<double3>& fractionalPositions,
                                      const std::vector<double>& radii)
{
  VoronoiNetwork network;
  network.simulationBox = simulationBox;
  network.atomRadii = radii;

  const double3x3 cell = simulationBox.cell;
  const double3x3 inverseCell = simulationBox.inverseCell;

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
  double volume = simulationBox.volume;
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

  auto findOrAddNode = [&](const double3& unwrappedFractional, double contributionRadius,
                           std::size_t atomIndex) -> std::size_t
  {
    double3 wrapped = double3::fract(unwrappedFractional);
    std::int32_t bx = std::min(binCount.x - 1, static_cast<std::int32_t>(wrapped.x * binCount.x));
    std::int32_t by = std::min(binCount.y - 1, static_cast<std::int32_t>(wrapped.y * binCount.y));
    std::int32_t bz = std::min(binCount.z - 1, static_cast<std::int32_t>(wrapped.z * binCount.z));

    for (std::int32_t dz = -1; dz <= 1; ++dz)
      for (std::int32_t dy = -1; dy <= 1; ++dy)
        for (std::int32_t dx = -1; dx <= 1; ++dx)
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

  // Undirected edges keyed by (min node, max node, oriented lattice shift); the value is
  // the running minimum bottleneck radius and the edge length.
  struct EdgeAccumulator
  {
    int3 delta;
    double radius;
    double length;
  };
  std::map<std::array<std::int64_t, 4>, EdgeAccumulator> edgeLookup;

  std::vector<SKVoronoiCell> cells = voronoi.computeAllCells();

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
        std::array<std::int64_t, 4> edgeKey{static_cast<std::int64_t>(a), static_cast<std::int64_t>(b),
                                            0, 0};
        // Fold the oriented lattice shift into the remaining key slots.
        edgeKey[2] = (static_cast<std::int64_t>(orientedDelta.x) << 42) ^
                     (static_cast<std::int64_t>(orientedDelta.y) << 21) ^
                     static_cast<std::int64_t>(orientedDelta.z);

        auto it = edgeLookup.find(edgeKey);
        if (it == edgeLookup.end())
        {
          edgeLookup[edgeKey] = EdgeAccumulator{orientedDelta, bottleneck, length};
        }
        else
        {
          it->second.radius = std::min(it->second.radius, bottleneck);
        }
      }
    }
  }

  network.edges.reserve(2 * edgeLookup.size());
  for (const auto& [key, accumulator] : edgeLookup)
  {
    std::size_t a = static_cast<std::size_t>(key[0]);
    std::size_t b = static_cast<std::size_t>(key[1]);
    network.edges.push_back(VoronoiEdge{a, b, accumulator.delta, accumulator.radius, accumulator.length});
    network.edges.push_back(VoronoiEdge{
        b, a, int3(-accumulator.delta.x, -accumulator.delta.y, -accumulator.delta.z), accumulator.radius,
        accumulator.length});
  }

  return network;
}
