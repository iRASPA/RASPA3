module;

module voronoi_pore_diameters;

import std;

import int3;
import double3;
import skspacegroupdatabase;
import atom;
import framework;
import forcefield;
import units;
import voronoi_network;

// Union-find over Voronoi nodes that also tracks the integer lattice offset of each node
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

PoreDiameters PoreDiameters::compute(const VoronoiNetwork& network)
{
  PoreDiameters diameters;
  diameters.includedSphereDiameter = network.largestIncludedSphereDiameter();

  if (network.nodes.empty() || network.edges.empty()) return diameters;

  // Process edges from widest to narrowest; the first edge that closes a loop with a
  // non-zero net lattice offset marks the percolation threshold (the free sphere).
  std::vector<std::size_t> order(network.edges.size());
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::sort(order.begin(), order.end(),
            [&](std::size_t a, std::size_t b) { return network.edges[a].radius > network.edges[b].radius; });

  PeriodicUnionFind unionFind(network.nodes.size());
  bool percolates = false;
  double freeRadius = 0.0;
  std::size_t percolatingRoot = 0;

  for (std::size_t index : order)
  {
    const VoronoiEdge& edge = network.edges[index];
    auto [rootFrom, accFrom] = unionFind.find(edge.from);
    auto [rootTo, accTo] = unionFind.find(edge.to);

    if (rootFrom == rootTo)
    {
      int3 net = accFrom + edge.delta - accTo;
      if (net.x != 0 || net.y != 0 || net.z != 0)
      {
        percolates = true;
        freeRadius = edge.radius;
        percolatingRoot = rootFrom;
        break;
      }
    }
    else
    {
      unionFind.parent[rootTo] = rootFrom;
      unionFind.offset[rootTo] = accFrom + edge.delta - accTo;
    }
  }

  if (!percolates) return diameters;

  diameters.freeSphereDiameter = 2.0 * freeRadius;

  // Dif: the largest node radius reachable within the component that first percolated.
  double maximumIncluded = 0.0;
  for (std::size_t i = 0; i < network.nodes.size(); ++i)
  {
    if (unionFind.find(i).first == percolatingRoot)
    {
      maximumIncluded = std::max(maximumIncluded, network.nodes[i].radius);
    }
  }
  diameters.includedAlongFreePathDiameter = 2.0 * maximumIncluded;

  return diameters;
}

void VoronoiPoreDiameters::run(const ForceField& forceField, const Framework& framework)
{
  std::chrono::system_clock::time_point time_begin = std::chrono::system_clock::now();

  std::vector<double3> fractionalPositions;
  std::vector<double> radii;
  fractionalPositions.reserve(framework.unitCellAtoms.size());
  radii.reserve(framework.unitCellAtoms.size());
  for (const Atom& atom : framework.unitCellAtoms)
  {
    fractionalPositions.push_back(framework.simulationBox.inverseCell * atom.position);
    std::size_t type = static_cast<std::size_t>(atom.type);
    radii.push_back(0.5 * forceField(type, type).sizeParameter());
  }

  VoronoiNetwork network = VoronoiNetwork::create(framework.simulationBox, fractionalPositions, radii);
  result = PoreDiameters::compute(network);

  std::chrono::duration<double> timing = std::chrono::system_clock::now() - time_begin;

  std::ofstream myfile;
  myfile.open(framework.name + ".voronoi.res.txt");
  std::print(myfile, "# Pore diameters from the Voronoi network (Di, Df, Dif)\n");
  std::print(myfile, "# Framework: {}\n", framework.name);
  std::print(myfile, "# Space-group Hall-number: {}\n", framework.spaceGroupHallNumber);
  std::print(myfile, "# Number of framework atoms: {}\n", framework.unitCellAtoms.size());
  std::print(myfile, "# Number of Voronoi nodes: {}\n", network.nodes.size());
  std::print(myfile, "# Number of Voronoi edges (directed): {}\n", network.edges.size());
  std::print(myfile, "# CPU Timing: {} [s]\n", timing.count());
  std::print(myfile, "Di (largest included sphere):            {} [Å]\n", result.includedSphereDiameter);
  std::print(myfile, "Df (largest free sphere):                {} [Å]\n", result.freeSphereDiameter);
  std::print(myfile, "Dif (included sphere along free path):   {} [Å]\n", result.includedAlongFreePathDiameter);
  myfile.close();
}
