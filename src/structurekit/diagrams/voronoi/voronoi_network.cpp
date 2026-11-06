module;

module voronoi_network;

import std;

import int3;
import double3;
import double3x3;
import unit_cell;
import skvoronoi;
import skapolloniusdiagram;

constexpr double mergeTolerance = 0.02;  // Å; zeo++ VOR_NODE_MERGE_THRESHOLD

// How far along the segment [a, b] its point nearest the origin is, from 0 at a to 1 at b.
double parameterNearestOriginOnSegment(const double3& a, const double3& b)
{
  double3 ab = b - a;
  double denominator = double3::dot(ab, ab);
  return (denominator > 0.0) ? std::clamp(-double3::dot(a, ab) / denominator, 0.0, 1.0) : 0.0;
}

// Squared distance from the origin to the segment [a, b].
double distanceSquaredOriginToSegment(const double3& a, const double3& b)
{
  double3 closest = a + parameterNearestOriginOnSegment(a, b) * (b - a);
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

  // Coplanar or coincident centres leave the centre undetermined along a direction as a function of t,
  // and the solve below, which writes it as one, has nothing to divide by. The test is |det| < 1e-12 s^3
  // for s the longest row, squared throughout so that neither the row lengths nor the scale need a root
  // taken of them.
  //
  // That is not to say there is no tangent sphere. Four centres in a plane with radii that differ have a
  // pair of them, one either side of the plane and of one radius: the rows then fix t, and it is the centre
  // that is free, along the normal. Returning nothing for those took a stepping stone away from the walk
  // wherever four atoms of two kinds lie in a plane, which in a crystal is everywhere. They are handed to
  // the diagram's own solve, which takes the line of solutions from the null space of the rows and so has
  // no such case; it is slower, and this is the one arrangement that needs it.
  double scaleSquared =
      std::max({rows[0].length_squared(), rows[1].length_squared(), rows[2].length_squared()});
  if (scaleSquared <= 0.0) return {};
  if (determinant * determinant < 1.0e-24 * scaleSquared * scaleSquared * scaleSquared)
  {
    TangentSpheres spheres;
    for (const SKApolloniusTangentSphere& sphere : skApolloniusTangentSpheres(centres, radii))
    {
      spheres.add(ApolloniusSphere{sphere.centre, sphere.radius});
    }
    return spheres;
  }

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

    // The two roots, (-b + root) / 2a first, each from the sum that does not cancel. Written out as
    // that formula, the root that stays finite as a goes to zero is a difference of two numbers nearly
    // alike, and it is the tangent sphere there is: the other is the one going off to infinity.
    double sum = -0.5 * (b + std::copysign(std::sqrt(discriminant), b));
    if (sum == 0.0)
    {
      roots[rootCount++] = 0.0;
    }
    else
    {
      roots[rootCount++] = (b >= 0.0) ? c / sum : sum / a;
      roots[rootCount++] = (b >= 0.0) ? sum / a : c / sum;
    }
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
                                      const std::vector<double>& givenRadii)
{
  VoronoiNetwork network;
  network.unitCell = unitCell;
  // Empty radii ask the diagram for the ordinary, unweighted Voronoi construction. From here on there is
  // one radius per atom whatever was asked, since everything below and everything downstream indexes the
  // radii by atom and takes the largest of them, and an unweighted atom is one of radius zero.
  const std::vector<double> radii =
      givenRadii.empty() ? std::vector<double>(fractionalPositions.size(), 0.0) : givenRadii;
  network.atomRadii = radii;

  const double3x3 cell = unitCell.cell;
  const double3x3 inverseCell = unitCell.inverseCell;

  // Radical (power) Voronoi diagram: the tessellation depends on the atom radii, matching
  // zeo++/voro++ semantics for radii-dependent pore analysis.
  SKVoronoi voronoi(cell, fractionalPositions, radii);
  network.atomPositionsFractional = voronoi.positions();
  const std::size_t numberOfAtoms = network.atomPositionsFractional.size();
  network.atomNodeVectors.assign(numberOfAtoms, {});

  // Gathering the vertices of the cells into nodes, after zeo++/voro++ voronoi_network. A vertex of the
  // diagram is a vertex of every cell that meets there, four of them or more, each of which has worked
  // it out for itself; and where the structure is nearly symmetric, a vertex that symmetry would make
  // one comes apart into several a few hundredths of an Ångström from each other. Left as they are, a
  // split node drops the constricting atom from some of its edges, which overestimates bottleneck
  // radii (Df) and node radii (Di). So vertices nearer each other than mergeTolerance are one node.
  //
  // Nearer is a distance. zeo++ asks instead that each Cartesian component of the difference be less
  // than the tolerance, which is a cube round the node and reaches 1.7 times as far along its diagonals
  // as along the axes: two vertices 0.03 Å apart are then one node or two according to which way the
  // crystal happens to be turned.
  //
  // And one node is everything joined by such pairs, however many steps it takes to get from one vertex
  // of it to another. zeo++ maps each vertex onto the first node already made that is near enough, and
  // makes a node of it otherwise. That depends on the order the vertices come in: of three in a row,
  // each near the next and the ends not near each other, the middle one goes with whichever end came
  // first, and the node stands where its first vertex stood. The order is that of the atoms in the
  // input and of the corners within each cell, neither of which is anything about the crystal, so the
  // same crystal with its atoms listed otherwise, or turned, came out with another network. Joined up
  // instead, the nodes are a property of the vertices alone, and a node stands at the mean of its
  // vertices, which turns with the crystal.
  //
  // What that gives up is a bound on how far across a node can be: a row of vertices each within the
  // tolerance of the next is one node however long the row. There is no such row in a diagram short of
  // a vertex of very many cells come apart, and gathering that into one node is what the tolerance is
  // for.
  double volume = unitCell.volume;
  double3 perpendicularWidths(volume / double3::cross(cell[1], cell[2]).length(),
                              volume / double3::cross(cell[2], cell[0]).length(),
                              volume / double3::cross(cell[0], cell[1]).length());
  // Bin width of ~1 Å (always > mergeTolerance) keeps the ±1-bin neighbour search exact
  // while the bin count stays proportional to the cell volume.
  int3 binCount(std::max(1, static_cast<std::int32_t>(perpendicularWidths.x)),
                std::max(1, static_cast<std::int32_t>(perpendicularWidths.y)),
                std::max(1, static_cast<std::int32_t>(perpendicularWidths.z)));

  // With one or two bins along an axis the -1 and +1 offsets wrap onto bins that are already
  // covered, so the scan range is narrowed to visit each neighbouring bin exactly once.
  auto offsetLow = [](std::int32_t count) { return count >= 3 ? -1 : 0; };
  auto offsetHigh = [](std::int32_t count) { return count >= 2 ? 1 : 0; };

  std::vector<SKVoronoiCell> cells = voronoi.computeAllCells();

  // Every vertex of every cell, numbered through: those of cell i from firstVertexOfCell[i] on.
  std::vector<std::size_t> firstVertexOfCell(cells.size() + 1, 0);
  for (std::size_t i = 0; i < cells.size(); ++i)
  {
    firstVertexOfCell[i + 1] = firstVertexOfCell[i] + cells[i].verticesCartesian.size();
  }
  const std::size_t numberOfCellVertices = firstVertexOfCell.back();

  std::vector<double3> vertexWrapped(numberOfCellVertices);   // fractional, in the home cell
  std::vector<double> vertexClearance(numberOfCellVertices);  // to the surface of the cell's own atom
  std::vector<int3> vertexBin(numberOfCellVertices);
  std::vector<std::vector<std::size_t>> vertexBins(static_cast<std::size_t>(binCount.x) *
                                                   static_cast<std::size_t>(binCount.y) *
                                                   static_cast<std::size_t>(binCount.z));
  auto binIndex = [&](std::int32_t bx, std::int32_t by, std::int32_t bz) -> std::size_t
  {
    return static_cast<std::size_t>(bx) +
           static_cast<std::size_t>(binCount.x) *
               (static_cast<std::size_t>(by) + static_cast<std::size_t>(binCount.y) * static_cast<std::size_t>(bz));
  };
  for (std::size_t i = 0; i < cells.size(); ++i)
  {
    for (std::size_t v = 0; v < cells[i].verticesCartesian.size(); ++v)
    {
      const std::size_t vertex = firstVertexOfCell[i] + v;
      const double3& relative = cells[i].verticesCartesian[v];  // relative to the site
      const double3 wrapped = double3::fract(network.atomPositionsFractional[i] + inverseCell * relative);
      vertexWrapped[vertex] = wrapped;
      vertexClearance[vertex] = relative.length() - radii[i];
      vertexBin[vertex] = int3(std::min(binCount.x - 1, static_cast<std::int32_t>(wrapped.x * binCount.x)),
                               std::min(binCount.y - 1, static_cast<std::int32_t>(wrapped.y * binCount.y)),
                               std::min(binCount.z - 1, static_cast<std::int32_t>(wrapped.z * binCount.z)));
      vertexBins[binIndex(vertexBin[vertex].x, vertexBin[vertex].y, vertexBin[vertex].z)].push_back(vertex);
    }
  }

  // From one vertex to another of the same home cell by the shortest way round, in fractions.
  auto shortestWay = [](const double3& from, const double3& to)
  {
    double3 difference = to - from;
    difference.x -= std::round(difference.x);
    difference.y -= std::round(difference.y);
    difference.z -= std::round(difference.z);
    return difference;
  };

  // The vertices joined up: each pointing at another of its node until one points at itself.
  std::vector<std::size_t> joinedTo(numberOfCellVertices);
  std::iota(joinedTo.begin(), joinedTo.end(), std::size_t{0});
  auto headOf = [&](std::size_t vertex)
  {
    while (joinedTo[vertex] != vertex)
    {
      joinedTo[vertex] = joinedTo[joinedTo[vertex]];
      vertex = joinedTo[vertex];
    }
    return vertex;
  };
  for (std::size_t vertex = 0; vertex < numberOfCellVertices; ++vertex)
  {
    const int3& bin = vertexBin[vertex];
    for (std::int32_t dz = offsetLow(binCount.z); dz <= offsetHigh(binCount.z); ++dz)
      for (std::int32_t dy = offsetLow(binCount.y); dy <= offsetHigh(binCount.y); ++dy)
        for (std::int32_t dx = offsetLow(binCount.x); dx <= offsetHigh(binCount.x); ++dx)
        {
          std::size_t neighbouringBin =
              binIndex((bin.x + dx + binCount.x) % binCount.x, (bin.y + dy + binCount.y) % binCount.y,
                       (bin.z + dz + binCount.z) % binCount.z);
          for (std::size_t other : vertexBins[neighbouringBin])
          {
            // Each pair once, and none that is joined already.
            if (other >= vertex) continue;
            std::size_t head = headOf(vertex);
            std::size_t otherHead = headOf(other);
            if (head == otherHead) continue;

            double3 cartesian = cell * shortestWay(vertexWrapped[vertex], vertexWrapped[other]);
            if (double3::dot(cartesian, cartesian) < mergeTolerance * mergeTolerance)
            {
              // The lower number heads, so that a node is headed by the first of its vertices.
              joinedTo[std::max(head, otherHead)] = std::min(head, otherHead);
            }
          }
        }
  }

  // A node to each head, in the order of the heads, standing at the mean of its vertices. The mean is taken
  // of the ways from the head, which are short, so that a node lying across a face of the cell is not
  // averaged from both sides of it. Its radius is the least clearance any of its cells gave it.
  constexpr std::size_t noNode = std::numeric_limits<std::size_t>::max();
  std::vector<std::size_t> nodeOfVertex(numberOfCellVertices, noNode);
  std::vector<double3> meanWay;
  std::vector<std::size_t> numberOfVertices;
  for (std::size_t i = 0; i < cells.size(); ++i)
  {
    for (std::size_t vertex = firstVertexOfCell[i]; vertex < firstVertexOfCell[i + 1]; ++vertex)
    {
      const std::size_t head = headOf(vertex);
      if (nodeOfVertex[head] == noNode)
      {
        nodeOfVertex[head] = network.nodes.size();
        VoronoiNode node;
        node.radius = vertexClearance[vertex];
        network.nodes.push_back(std::move(node));
        meanWay.push_back(double3(0.0, 0.0, 0.0));
        numberOfVertices.push_back(0);
      }
      const std::size_t nodeIndex = nodeOfVertex[head];
      nodeOfVertex[vertex] = nodeIndex;

      VoronoiNode& node = network.nodes[nodeIndex];
      node.radius = std::min(node.radius, vertexClearance[vertex]);
      if (std::find(node.atomIndices.begin(), node.atomIndices.end(), i) == node.atomIndices.end())
      {
        node.atomIndices.push_back(i);
      }
      meanWay[nodeIndex] += shortestWay(vertexWrapped[head], vertexWrapped[vertex]);
      ++numberOfVertices[nodeIndex];
    }
  }
  for (std::size_t vertex = 0; vertex < numberOfCellVertices; ++vertex)
  {
    if (headOf(vertex) != vertex) continue;
    const std::size_t nodeIndex = nodeOfVertex[vertex];
    VoronoiNode& node = network.nodes[nodeIndex];
    node.fractional = double3::fract(vertexWrapped[vertex] +
                                     meanWay[nodeIndex] / static_cast<double>(numberOfVertices[nodeIndex]));
    node.position = cell * node.fractional;
  }

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
      vertexFractional[v] = siteFractional + inverseCell * relative;
      vertexNode[v] = nodeOfVertex[firstVertexOfCell[i] + v];

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

  // No clearance in a periodic structure exceeds half the cell's longest diagonal: every point
  // lies in some image of the cell, within that distance of one of its corners, and the images
  // of any one atom sit at the same offset from every corner. The bound keeps the searches below
  // finite, as a tangent sphere solved from four nearly coplanar atoms can be arbitrarily wide.
  double maximumClearance = 0.5 * unitCell.longestDiagonal();

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

    // Where the edge is narrowest goes with how narrow it is. The radius is the least, over the atoms, of
    // the distance from an atom to the segment less its radius, and the least over the atoms of a least
    // over the segment is the least clearance anywhere on the segment: it is the clearance at the point of
    // the segment nearest the atom that sets it. So the edge has a narrowest point after all, and that
    // point is kept, with the direction of the segment, which is the way the passage runs there.
    //
    // The radius coming in was taken on the segment as each cell saw it, between its own vertices, and
    // this pass is on the segment between the nodes those vertices were gathered into. The two are the
    // same segment unless a node is several vertices, and then they differ by about the distance vertices
    // are gathered over and the radius kept is the smaller. So the clearance at the point kept is the
    // radius, or more than it by that little.
    double clearance = forward.radius;
    double narrowest = std::numeric_limits<double>::max();
    double along = 0.0;
    forEachAtomNear(wrappedMidpoint, halfLength + forward.radius + maximumAtomRadius,
                    [&](std::size_t j, const double3& atomImage)
                    {
                      double3 a = p + shift - atomImage;
                      double3 b = q + shift - atomImage;
                      double t = parameterNearestOriginOnSegment(a, b);
                      double distance = (a + t * (b - a)).length() - radii[j];
                      if (distance < narrowest)
                      {
                        narrowest = distance;
                        along = t;
                      }
                      clearance = std::min(clearance, distance);
                    });
    forward.radius = clearance;
    network.edges[i + 1].radius = clearance;

    // In the frame where `from` is in the home cell, which for the edge taken the other way round is one
    // lattice translation away and runs the other way.
    if (narrowest < std::numeric_limits<double>::max() && halfLength > 0.0)
    {
      VoronoiEdge& backward = network.edges[i + 1];
      forward.bottleneckPosition = p + along * (q - p);
      forward.bottleneckDirection = (q - p) / (2.0 * halfLength);
      forward.hasBottleneckGeometry = true;
      backward.bottleneckPosition = forward.bottleneckPosition - (q - network.nodes[forward.to].position);
      backward.bottleneckDirection = -forward.bottleneckDirection;
      backward.hasBottleneckGeometry = true;
    }
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
  constexpr std::size_t apolloniusCandidateLimit = 20;  // the most taken once ties are let in
  constexpr double apolloniusCandidateTie = 1.0e-9;     // Å; nearer than this in clearance is as near
  constexpr double apolloniusConvergenceTolerance = 1.0e-12;

  // Whether the segment between two points lies in the void: no atom reaches it. Measured on the
  // atoms around its midpoint, as the edges' bottlenecks are.
  auto segmentIsFree = [&](const double3& a, const double3& b)
  {
    double3 midpoint = 0.5 * (a + b);
    double3 wrappedMidpoint = cell * double3::fract(inverseCell * midpoint);
    double3 shift = wrappedMidpoint - midpoint;
    double halfLength = 0.5 * (b - a).length();
    bool free = true;
    forEachAtomNear(wrappedMidpoint, halfLength + maximumAtomRadius,
                    [&](std::size_t j, const double3& atomImage)
                    {
                      if (!free) return;
                      double squared = distanceSquaredOriginToSegment(a + shift - atomImage, b + shift - atomImage);
                      if (squared <= radii[j] * radii[j]) free = false;
                    });
    return free;
  };

  // A vertex inside an atom is no place to start from: no segment out of it is free, so the walk
  // could not leave it, and the pocket it would have belonged to has its own vertices to walk from.
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
    // rest: the walk cannot cycle because the clearance strictly increases. What it does not do
    // by itself is keep the walk in its own pore, since a jump is a jump and not a path, and a
    // tangent sphere solved from atoms on this side of a wall can sit on the other side of it. So
    // a jump is taken only along a segment that stays clear of every atom, which keeps the walk
    // inside the free space it started in and lets it go as far as that space does. A bound on the
    // jump's length cannot do the same job: tied to the vertex's own clearance it stopped the walk
    // dead wherever a vertex sat near a wall, and the peak of a pocket is routinely several such
    // clearances from the vertices around it.
    //
    // The candidates are the atoms whose surfaces lie within the current clearance of the point,
    // and one more atom radius beyond that: the four the next sphere touches are not all within
    // reach of the point it is solved from, and a walk that only sees the atoms already touching
    // it has nothing to climb with. It halts at a point no tangent solve improves, the local peak.
    double3 point = node.position;
    for (std::size_t iteration = 0; iteration < maximumApolloniusIterations; ++iteration)
    {
      candidates.clear();
      forEachAtomNear(point, node.maximalRadius + 2.0 * maximumAtomRadius,
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
      //
      // The nearest few are not cut off between two atoms that are equally near. Where the structure has
      // symmetry they come in shells, six at one clearance and six at the next, and a cut through a shell
      // takes whichever of its atoms rounding put first: the quadruple the peak needs is then among the
      // candidates or not by accident, and of six nodes that symmetry makes alike four reached their peak
      // and two stayed where they were. So an atom as near as the last one taken is taken too, up to a
      // limit that holds a shell of twelve beyond the first eight.
      std::size_t sorted = std::min(apolloniusCandidateLimit, candidates.size());
      std::partial_sort(candidates.begin(), candidates.begin() + static_cast<std::ptrdiff_t>(sorted),
                        candidates.end(),
                        [](const Candidate& lhs, const Candidate& rhs) { return lhs.clearance < rhs.clearance; });
      std::size_t considered = std::min(apolloniusCandidateCount, sorted);
      while (considered < sorted &&
             candidates[considered].clearance - candidates[considered - 1].clearance < apolloniusCandidateTie)
      {
        ++considered;
      }

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
                // A sphere wider than any clearance the cell admits is not empty, and its centre is
                // that far from the atoms it was solved from: measuring the clearance there, over a
                // search ball that wide, would be a walk over the whole of space.
                if (sphere.radius > maximumClearance) continue;

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

                // Only a sphere that improves is worth the walk to it, and only one that can be
                // walked to through the void is taken.
                if (clearance > bestClearance && segmentIsFree(point, sphere.centre))
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
