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

  // Per node, the lattice translation carrying it into the frame of the first node of its pore, so that
  // the positions `nodes[i].position + cell * nodeLatticeOffset[i]` over a pore are one connected lift of
  // it rather than one representative of each periodic family. Zero for a pruned node.
  //
  // The walk below accumulates this already, to decide whether a component closes on a translate of
  // itself; it is kept because anything integrating over the boundary of a bounded pore has to assemble
  // that boundary in a single frame first, the pieces belonging to different lifts otherwise. Within a
  // pocket it is well defined: a second path to a node that disagreed with the first would be a loop with
  // a non-zero net translation, which is what makes a component a channel instead.
  std::vector<int3> nodeLatticeOffset;

  static ChannelAnalysis compute(const VoronoiNetwork& network, double probeRadius);
};

export struct VoronoiChannels
{
  ChannelAnalysis result;

  void run(const ForceField& forceField, const Framework& framework, std::string probePseudoAtom);

  // The same, over a network already built on the framework's own radii. Building the network costs
  // more than this analysis does, so a caller that has one already, as the pore-diameter pass does,
  // should not pay for a second identical one.
  void run(const ForceField& forceField, const Framework& framework, std::string probePseudoAtom,
           const VoronoiNetwork& network);
};

// Rank (0..3) of a set of integer lattice vectors, i.e. the dimensionality they span.
export int latticeVectorRank(const std::vector<int3>& vectors);
