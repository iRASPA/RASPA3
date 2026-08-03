module;

export module voronoi_channels;

import std;

import crystal;
import pair_interactions;
import voronoi_network;

// Which components of a network are channels and which are pockets is decided by the same walk whichever
// network it is, so the analysis itself lives apart from any of them. Re-exported because a caller asking
// this route for the channels is asking for that type back.
export import channel_analysis;

// The radical (Voronoi) route to them: the network built on the atoms of a framework, and the channels
// picked out of it for a probe.
export struct VoronoiChannels
{
  ChannelAnalysis result;

  void run(const PairInteractions& interactions, const Crystal& framework, std::string probePseudoAtom);

  // The same, over a network already built on the framework's own radii. Building the network costs
  // more than this analysis does, so a caller that has one already, as the pore-diameter pass does,
  // should not pay for a second identical one.
  void run(const PairInteractions& interactions, const Crystal& framework, std::string probePseudoAtom,
           const VoronoiNetwork& network);
};
