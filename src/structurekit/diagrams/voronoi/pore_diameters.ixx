module;

export module voronoi_pore_diameters;

import std;

import crystal;
import pair_interactions;
import voronoi_network;

// Di, Df and Dif are read off a network by arithmetic that belongs to no diagram, and they are read the same
// way here, from the Apollonius diagram and from a roadmap of samples. Re-exported because a caller asking
// this route for the diameters is asking for that type back.
export import pore_diameters;

// The radical (Voronoi) route to them: the network built on the atoms of a framework, and the diameters read
// off it. It is the cheap route and the approximate one. A radical vertex is not where the clearance peaks
// unless the atom radii are all equal, so Di comes from a local ascent away from the vertex, and an edge's
// radius is a bound over the whole segment rather than the narrowest point of the passage, so Df is a lower
// bound on the real bottleneck.
export struct VoronoiPoreDiameters
{
  PoreDiameters result;
  VoronoiNetwork network;  // the network the diameters were read from, kept for the analyses that follow

  void run(const PairInteractions& interactions, const Crystal& framework);
};
