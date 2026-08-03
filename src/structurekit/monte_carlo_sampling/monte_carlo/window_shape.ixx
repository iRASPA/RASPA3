module;

export module mc_window_shape;

import std;

import pore_window;
import sampled_structure;
import sampled_roadmap;

// The shape of the narrowest window of a pore system, by sampling.
//
// Df is the inscribed circle of the window at the bottleneck and it says nothing about how far the window
// reaches the other ways. A slot 3 Å across and 9 Å long and a round hole 3 Å across have the same Df and
// admit quite different molecules. So the plane across the passage at the bottleneck is measured: the
// shortest and longest free chord through it, and the largest-area ellipse that fits in it.
//
// What this needs of a route, and what has until now confined it to the Apollonius diagram, is not the
// width of the bottleneck but its *place*: a point on the narrowest cross-section and the direction the
// passage runs in there, without which there is no plane to measure in. A radical Voronoi edge carries a
// width and no such point, its arc being solved for rather than walked along. A sampled hop carries both,
// because the width of a hop was found by sliding a sphere along it and the narrowest point is where the
// sphere was smallest -- which is the one thing sampling gets for free that the exact construction does
// not. So the same `PoreWindow::measure` runs here unchanged, and reports the second route to these numbers
// that the two can be checked against each other on.
//
// The measurement in the plane is itself exact given the plane: rays are cast from the bottleneck out to
// the first atom each meets, over the whole turn. Both the widths and the ellipse are lower bounds, each
// direction being followed only as far as that first atom.
export struct MC_WindowShape
{
  // At the bottleneck of the widest path across the whole crystal, which is the Df the diameters report.
  PoreWindow freeSphere;
  double freeSphereDiameter{0.0};  // Å

  // And per channel, at that channel's own bottleneck. Numbered as in the sampled channel analysis.
  std::vector<ChannelWindow> channels;

  bool percolates{false};

  SampledRoadmap roadmap;
  double seconds{0.0};

  void run(const SampledStructure &structure, std::optional<std::size_t> numberOfIterations,
           std::optional<std::size_t> numberOfInnerSteps);

  void run(const SampledStructure &structure, const SampledRoadmap &roadmap);
};
