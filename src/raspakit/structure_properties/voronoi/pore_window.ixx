module;

export module pore_window;

import std;

import double3;
import voronoi_network;
import voronoi_channels;

// The shape of the narrowest window of a channel: what Df says with one number, said with two.
//
// Df is the diameter of the largest sphere that can travel through the structure, so it is the
// inscribed circle of the tightest window along the best percolating path and nothing more. A slit and
// a round hole of the same Df are the same number, though a flat molecule passes one and not the other,
// and the crystallographic tables of zeolite apertures have always quoted two numbers for this reason.
// The window is already located by the diagram, since the bottleneck of an arc is a point on it and the
// arc has a tangent there, so measuring the window costs a plane's worth of arithmetic on top of a
// structure that has been built anyway.
//
// What is measured. The plane through the bottleneck perpendicular to the passage cuts each atom that
// reaches it in a disc, of radius sqrt(r² - d²) for an atom of radius r whose centre stands d from the
// plane, and the free cross-section is what those discs leave. Two things are read from it, both
// centred on the bottleneck itself, which is the incentre of the window:
//
//   - the free width, the length of the free chord through the bottleneck, at its smallest and largest
//     over all directions in the plane. This is the aperture measured across, as the tables do it, and
//     the smallest of them is at least Df because the inscribed circle is one candidate chord.
//   - the largest-area ellipse that fits, over every orientation and elongation. Two numbers for the
//     shape of the window rather than for two of its chords, and the more useful of the pair when the
//     cross-section is not convex, where a long chord can pass through a place too tight to be of any
//     use.
//
// Both are lower bounds and cannot overstate the room. Each direction is followed only as far as the
// first atom it meets, so the region measured is the part of the cross-section visible from the
// bottleneck. That region is contained in the free cross-section whatever shape the latter has, so an
// ellipse inscribed in it is inscribed in the window: the ellipse reported is always free, and is the
// largest only when the window is star shaped about its incentre, which an aperture ringed by atoms is.
//
// What they are not. Neither is a criterion for whether a non-spherical molecule passes. That depends
// on the molecule's own shape and on the freedom it has to turn along the way, which is a question in
// its configuration space rather than in this plane. These are a description of the window.
//
// Nor is every bottleneck a window. The plane is the one the diagram gives, across the arc at its
// narrowest point, and where that plane is not an aperture ringed by atoms the measurement says so
// rather than inventing one; see `clipped`. Of the two zeolites whose windows are the textbook cases,
// this puts FAU's at 6.8 by 6.8 Å across its twelve-ring, perpendicular to <111> as it should be, and
// MFI's at 4.6 by 5.0 Å perpendicular to b, along the straight channel. Both agree with the
// crystallographic tables once their smaller oxygen radius is accounted for.
export struct PoreWindow
{
  bool measured{false};  // false when the passage is shut at the bottleneck, leaving no window to measure

  double3 position{};  // Cartesian, the narrowest point of the passage
  double3 normal{};    // unit, the direction the passage runs in there; the window lies across it

  double freeRadius{0.0};  // clearance at the bottleneck [Å]; twice this is the Df of the passage

  double smallestFreeWidth{0.0};  // free chord through the bottleneck, over all directions [Å]
  double largestFreeWidth{0.0};

  double minorAxis{0.0};          // diameters of the largest-area inscribed ellipse [Å]
  double majorAxis{0.0};
  double3 majorAxisDirection{};   // unit, in the plane of the window

  // A direction in which no atom was met within the reach the cell allows, so the window is open that
  // way and what is reported above is cut off there rather than measured. It means the plane across the
  // bottleneck is not an aperture ringed by atoms: either the passage has no ring around it, being a
  // saddle in a pore that is open in the other directions, or the arc crosses its ring obliquely and the
  // plane perpendicular to it runs off down the channel. The smallest free width is unaffected, being
  // bounded by the inscribed circle either way; the largest and the ellipse are not to be read.
  bool clipped{false};

  std::size_t boundingAtoms{0};  // atoms that bound the window, the ring around it

  // Measures the window in the plane through `position` perpendicular to `normal`, over the atoms of
  // the network and their periodic images.
  static PoreWindow measure(const VoronoiNetwork& network, const double3& position, const double3& normal);
};

// The narrowest window of one channel, at the bottleneck of its own widest percolating path. That
// bottleneck is where the channel's Df is set, so the two sit together.
export struct ChannelWindow
{
  std::size_t poreIndex{0};  // index into ChannelAnalysis::pores
  int dimensionality{0};
  double freeSphereDiameter{0.0};  // the Df of this channel alone [Å]
  std::size_t limitingEdge{0};     // index into network.edges of its bottleneck
  PoreWindow window;
};

// The window at the bottleneck that sets the structure's own Df, taken over every node of the network.
// This is the one that sits next to Df itself, and like Df it is a property of the structure rather than
// of a probe. `measured` comes back false where the structure does not percolate at all.
export PoreWindow freeSphereWindow(const VoronoiNetwork& network);

// The windows of every channel found by the channel analysis. Pockets are left out, having no
// percolating path to be limited by.
//
// A channel whose bottleneck the network does not locate is reported with `window.measured` false:
// which edge is the bottleneck is always known, but where along that edge it sits is known only for a
// network built from an Apollonius diagram. See VoronoiEdge::hasBottleneckGeometry.
export std::vector<ChannelWindow> channelWindows(const VoronoiNetwork& network, const ChannelAnalysis& channels);
