module;

export module voronoi_channels;

import std;

import int3;
import framework;
import forcefield;
import voronoi_network;

// Identification of channels and pockets from the Voronoi network for a given probe
// radius, and their dimensionality.
//
// After pruning nodes/edges that are too narrow for the probe, each connected component
// is examined for periodic self-connection: if a depth-first walk reaches a node it has
// already visited but at a different periodic image, the component percolates and is a
// channel; otherwise it is an isolated pocket. The dimensionality (1/2/3) is the rank of
// the set of lattice-translation vectors along which the component connects to itself.

export struct VoronoiPore
{
  bool isChannel{false};       // true = channel (percolates), false = pocket
  int dimensionality{0};       // 0 pocket, 1/2/3 channel
  std::vector<std::size_t> nodeIndices;
};

export struct ChannelAnalysis
{
  std::vector<VoronoiPore> pores;
  std::vector<std::int32_t> nodePoreId;  // pore index per node, or -1 if pruned
  std::size_t numberOfChannels{0};
  std::size_t numberOfPockets{0};

  static ChannelAnalysis compute(const VoronoiNetwork& network, double probeRadius);
};

export struct VoronoiChannels
{
  ChannelAnalysis result;

  void run(const ForceField& forceField, const Framework& framework, std::string probePseudoAtom);
};

// Rank (0..3) of a set of integer lattice vectors, i.e. the dimensionality they span.
export int latticeVectorRank(const std::vector<int3>& vectors);
