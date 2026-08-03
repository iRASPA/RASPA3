module;

export module apollonius_network;

import std;

import double3;
import unit_cell;
import skapolloniusdiagram;
import voronoi_network;

// The pore network read off the Apollonius diagram.
//
// The pore analyses of zeo++ all run on the same object: a graph whose nodes are the local maxima of
// the clearance min_j(|x - x_j| - r_j) and whose edges are the passages between them, each annotated
// with the radius of the largest sphere that fits there. That graph is what `VoronoiNetwork` holds,
// and it is what is built here, only from the additively weighted diagram rather than the radical one.
// The Apollonius diagram is a Voronoi diagram, of the clearance rather than the distance, so the two
// networks are the same kind of object and every analysis downstream is unchanged.
//
// What differs is that the Apollonius network is exact where the radical one is not. A vertex of the
// radical diagram is equidistant in power, not in clearance, so unless the radii are all equal it is
// not where the clearance peaks: the largest included sphere sits at a nearby Apollonius vertex, and
// the radical network can only approach it by a local ascent. Likewise the narrowest point of a
// radical edge is not the narrowest point of the passage. Here the vertices are the peaks themselves,
// tangent to four sites and certified to overlap none, and each arc's bottleneck is the smallest
// tangent-sphere radius along the trisector curve. Di, Df and Dif then come out of the graph directly.
export struct ApolloniusPoreNetwork
{
  VoronoiNetwork network;

  std::size_t numberOfVertices{0};  // vertices of the diagram, one node each
  std::size_t numberOfArcs{0};      // arcs between two vertices, one undirected edge each

  // Closed trisectors carrying no vertex: ring-shaped passages that no fourth site ever intrudes on.
  // A ring has no endpoints, so a graph of nodes and edges has nowhere to put it and it is left out.
  // A ring that winds through the periodic boundary is a channel, and dropping it can only make Df
  // smaller, never larger, so the count is reported rather than hidden.
  std::size_t numberOfRings{0};

  // The diagram's own account of how it came out. A missing peak or a missing passage can only make a
  // pore diameter smaller, never larger, so the counters are carried along and written out rather than
  // reduced to a single flag.
  SKApolloniusDiagram::Verification verification;

  // Whether the part of the diagram this network is made of came out consistent: no two vertices left
  // at one point, every vertex carrying an edge along each of its branches, every triple of sites
  // paired at both ends of its arc save for the arcs the free region clipped, and every arc's
  // narrowest point measured.
  //
  // This is weaker than the diagram's own `verification.isComplete()`, which also asks that every
  // bisector patch close into a cycle and that no branch direction come out a tie. Both fail routinely
  // on a symmetric framework, where a sphere tangent to five sites at once is common: more than two
  // edges of a patch then meet at one vertex, and a fifth site can recede at zero rate along a branch.
  // Neither counts a node or an arc, and the patches are never read here. A branch sent the wrong way
  // by such a tie would be: it would leave its triple paired once or three times, which is checked.
  bool networkIsComplete() const;

  // The size of the diagram and what its verification came to, as comment lines, for the header of an
  // analysis' output file. A peak or a passage the diagram failed to find lowers what is read off it,
  // so that is reported with the numbers rather than left to be discovered elsewhere.
  void writeHeader(std::ostream& stream) const;

  static ApolloniusPoreNetwork create(const UnitCell& unitCell,
                                      const std::vector<double3>& fractionalPositions,
                                      const std::vector<double>& radii);
};
