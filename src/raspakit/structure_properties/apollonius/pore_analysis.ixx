module;

export module apollonius_pore_analysis;

import std;

import framework;
import forcefield;
import voronoi_pore_diameters;
import voronoi_channels;
import pore_window;

// The zeo++ pore properties computed from the Apollonius diagram: the pore diameters Di, Df and Dif,
// and the channel/pocket analysis for a probe.
//
// Both come out of one pore network by the very analyses that run on the radical Voronoi network, so
// what changes is not how the properties are read but what they are read from. Building the diagram
// costs far more than either analysis, which is why they are done together here rather than in one
// facade each.
//
// The diagram is the exact answer where the radical network is an approximation of it. Its vertices
// are the true maxima of the clearance, so Di is the deepest of them rather than the best a local
// ascent could find; its arcs carry the true narrowest point of each passage, so Df is the real
// bottleneck of the percolating path rather than the narrowest radical edge along it. The price is
// the construction, which is orders of magnitude slower than a radical tessellation.
//
// The diagram also says where each bottleneck is, not only how wide it is, which is what lets the window
// across it be measured and Df be given a shape rather than a single number; see module pore_window.
export struct ApolloniusPoreAnalysis
{
  PoreDiameters diameters;
  ChannelAnalysis channels;
  PoreWindow freeSphere;               // the window at the bottleneck that sets Df
  std::vector<ChannelWindow> windows;  // the narrowest window of each channel

  void run(const ForceField& forceField, const Framework& framework, std::string probePseudoAtom);
};
