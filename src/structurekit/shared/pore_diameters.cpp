module;

module pore_diameters;

import std;

import int3;
import voronoi_network;

// Union-find over network nodes that also tracks the integer lattice offset of each node
// relative to its set representative. Merging two nodes that are already connected with a
// different net offset signals a path that percolates through the periodic boundary.
struct PeriodicUnionFind
{
  std::vector<std::size_t> parent;
  std::vector<int3> offset;  // lattice offset of node relative to parent

  explicit PeriodicUnionFind(std::size_t n) : parent(n), offset(n, int3(0, 0, 0))
  {
    std::iota(parent.begin(), parent.end(), std::size_t{0});
  }

  std::pair<std::size_t, int3> find(std::size_t x)
  {
    int3 accumulated(0, 0, 0);
    while (parent[x] != x)
    {
      accumulated += offset[x];
      x = parent[x];
    }
    return {x, accumulated};
  }
};

PercolatingPath widestPercolatingPath(const VoronoiNetwork& network, std::span<const std::size_t> nodes)
{
  PercolatingPath path;
  if (network.nodes.empty() || network.edges.empty() || nodes.empty()) return path;

  std::vector<char> inSubset(network.nodes.size(), 0);
  for (std::size_t node : nodes) inSubset[node] = 1;

  // Process edges from widest to narrowest; the first edge that closes a loop with a
  // non-zero net lattice offset marks the percolation threshold (the free sphere).
  std::vector<std::size_t> order;
  order.reserve(network.edges.size());
  for (std::size_t i = 0; i < network.edges.size(); ++i)
  {
    const VoronoiEdge& edge = network.edges[i];
    if (inSubset[edge.from] != 0 && inSubset[edge.to] != 0) order.push_back(i);
  }
  std::sort(order.begin(), order.end(),
            [&](std::size_t a, std::size_t b) { return network.edges[a].radius > network.edges[b].radius; });

  // The components are read off once every edge as wide as the bottleneck has been added, not when the
  // first of them closed the loop: a node behind a neck exactly that wide is reachable at the bottleneck's
  // radius, so it is on the path, and which of two tied edges the sort put first cannot decide that. The same
  // goes for a second channel closing a loop of its own at that radius, so every loop closed at it is kept.
  PeriodicUnionFind unionFind(network.nodes.size());
  std::vector<std::size_t> percolatingMembers;

  for (std::size_t index : order)
  {
    const VoronoiEdge& edge = network.edges[index];
    // an edge with no clearance is buried in the atoms (the network keeps them): no sphere passes it, so a
    // structure that percolates only through such edges does not percolate
    if (edge.radius <= 0.0) break;
    if (path.percolates && edge.radius < path.radius) break;

    auto [rootFrom, accFrom] = unionFind.find(edge.from);
    auto [rootTo, accTo] = unionFind.find(edge.to);

    if (rootFrom == rootTo)
    {
      int3 net = accFrom + edge.delta - accTo;
      if (net.x != 0 || net.y != 0 || net.z != 0)
      {
        if (!path.percolates)
        {
          path.percolates = true;
          path.limitingEdge = index;
          path.radius = edge.radius;
        }
        percolatingMembers.push_back(edge.from);
      }
    }
    else
    {
      unionFind.parent[rootTo] = rootFrom;
      unionFind.offset[rootTo] = accFrom + edge.delta - accTo;
    }
  }

  if (!path.percolates) return path;

  std::set<std::size_t> percolatingRoots;
  for (std::size_t member : percolatingMembers) percolatingRoots.insert(unionFind.find(member).first);
  for (std::size_t node : nodes)
  {
    if (percolatingRoots.contains(unionFind.find(node).first)) path.componentNodes.push_back(node);
  }

  return path;
}

PoreDiameters PoreDiameters::compute(const VoronoiNetwork& network)
{
  PoreDiameters diameters;
  diameters.includedSphereDiameter = network.largestIncludedSphereDiameter();

  std::vector<std::size_t> allNodes(network.nodes.size());
  std::iota(allNodes.begin(), allNodes.end(), std::size_t{0});

  PercolatingPath path = widestPercolatingPath(network, allNodes);
  if (!path.percolates) return diameters;

  diameters.freeSphereDiameter = 2.0 * path.radius;

  // Dif: the largest node radius reachable at the free sphere's width, in whichever channel it is.
  double maximumIncluded = 0.0;
  for (std::size_t node : path.componentNodes)
  {
    maximumIncluded = std::max(maximumIncluded, network.nodes[node].maximalRadius);
  }
  diameters.includedAlongFreePathDiameter = 2.0 * maximumIncluded;

  return diameters;
}
