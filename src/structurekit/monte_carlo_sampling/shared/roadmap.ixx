module;

export module sampled_roadmap;

import std;

import double3;
import voronoi_network;
import channel_analysis;
import sampled_structure;
import sampling_backend;

// The void of a crystal as a graph, found by throwing points at it.
//
// Every geometric question about a pore system -- how much room there is, how much of it a probe can reach,
// how wide the widest sphere that gets through is, which way the channels run, which pockets are shut --
// is a question about the shape of one region of space, and the other routes here answer them by building
// that region: a Voronoi tessellation, an Apollonius diagram, a grid of clearances. Each construction is
// exact and each is a great deal of machinery whose failure modes have nothing to do with the question.
//
// This is the other way round it. A point is thrown at the cell; if it lands in the void it is kept, with
// the room around it. Two kept points close enough together are joined when a sphere can be slid in a
// straight line from one to the other, and how big that sphere may be is the width of the join. What comes
// out is a roadmap: a graph on which every path is a path a sphere can really travel, and whose bottleneck
// is a real one. It is the same shape of thing a diagram produces -- `VoronoiNetwork` is the type all of
// them are kept in -- so everything that reads a pore network reads this one too, and the answers can be
// set beside each other.
//
// Everything read off it is a lower bound, and this is the property worth understanding about the whole
// family. A sphere the roadmap says fits, fits, because a point of the void was found with that much room
// around it. A path the roadmap says is open, is open, because it is a chain of straight hops each of which
// was checked exactly. What sampling can miss is a deeper pocket that no point landed in and a wider way
// round that no chain of hops happened to take, and both of those can only make an answer smaller.
//
// Two stages carry the answers past what the points alone would give.
//
// A point that no neighbour of its in the roadmap beats is at the top of its own pocket as far as the
// sample can see, so it is walked further uphill until the room around it stops growing. That lands on the
// true deepest point of the pocket rather than the best point that happened to be drawn, and a few thousand
// points already put the largest included sphere within a ten-thousandth of an Ångström.
//
// The distinct places those walks end are the pockets of the crystal, and every pair of them is then tried
// as a way from one to the other. This is what a roadmap of samples cannot get by having more points in it:
// a hop between two points that merely fell near the middles of neighbouring pockets misses the middle of
// the window between them by about the spacing of the sample, however fine that spacing is made, whereas
// the way between the two pocket centres does not, being the way the passage was built around. Where a
// passage curves and the straight line between the centres clips a wall, the line is bent -- its most
// hemmed-in point is walked out sideways and the two halves are put through the same again.
export struct SampledRoadmap
{
  // The graph itself. Nodes `[0, numberOfVoidSamples)` are the points that landed in the void, in the order
  // they were drawn; nodes `[firstPocketNode, ...)` are the deepest point of each distinct pocket. Every
  // hop is stored both ways round, as the diagram builders store theirs, and carries the narrowest point of
  // the way it stands for in `bottleneckPosition` and `bottleneckDirection`.
  VoronoiNetwork network;

  // How the roadmap falls apart into pieces, and which of them run off across the crystal. Every property
  // read off a roadmap needs this and it is the same walk each time, so it is done once here rather than
  // four times over by the callers, who would otherwise be free to disagree about it.
  ChannelAnalysis components;

  // Per pore of `components`, whether it is a piece of the void the sample was fine enough to resolve. A
  // channel always is. A pocket is when a walk uphill ended inside it, that walk having found a genuine
  // local maximum of the room available; a piece with no such maximum in it is one or two points at the
  // very edge of the void that no hop happened to reach, and calling it a pocket would be counting a stray
  // as a cavity -- and, where the structure is being blocked, writing a sphere in a working pore.
  std::vector<char> poreIsResolved;

  std::size_t numberOfSamples{0};       // points thrown
  std::size_t numberOfVoidSamples{0};   // of those, the ones that landed in the void
  std::size_t numberOfLinks{0};         // hops between the nodes, counted once each
  std::size_t numberOfPeaks{0};         // samples no neighbour beats, each walked uphill
  std::size_t numberOfPocketCentres{0};  // distinct places those walks ended
  std::size_t firstPocketNode{0};       // where those enter `network.nodes`
  std::size_t numberOfJoins{0};         // ways found between components the rest of the roadmap left apart

  std::size_t numberOfChannels{0};         // components that run off across the crystal
  std::size_t numberOfPockets{0};          // and those that close on themselves and were resolved
  std::size_t numberOfUnresolvedPieces{0};  // pieces too thinly sampled to call either way

  double linkDistance{0.0};           // Å, how far apart two samples may be and still be joined
  std::size_t numberOfInnerSteps{0};  // uphill steps allowed a peak

  // Which backend measured the spheres, which is also the suffix of every file written from this roadmap,
  // so that running both leaves two sets of results side by side rather than one on top of the other.
  std::string backend{"cpu"};
  double seconds{0.0};

  // The share of the points that landed in the void. This is the void fraction of the structure at the
  // contact radii it was given, and it comes out of building the roadmap whether it was asked for or not.
  double voidFraction{0.0};

  // `numberOfIterations` points, each walked at most `numberOfInnerSteps` steps uphill if it turns out to
  // be a peak. The link distance is set from the sample size, there being one spacing at which the roadmap
  // is neither disconnected nor needlessly dense.
  //
  // The algorithm is here and the arithmetic is in the backend, which is what lets the processor and the
  // GPU run the same construction rather than two of them that have to be kept in step.
  static SampledRoadmap build(const SampledStructure &structure, const SamplingBackend &backend,
                              std::optional<std::size_t> numberOfIterations,
                              std::optional<std::size_t> numberOfInnerSteps);

  // The nodes that stand for the void the sample covers, which is the sample points and not the pocket
  // centres added on top of them. Anything counting volume has to count these and only these, the pocket
  // centres being extra nodes at places already counted.
  std::size_t numberOfVolumeNodes() const { return this->numberOfVoidSamples; }

  // Whether `node` belongs to a piece of the void that runs off across the crystal, which is where a probe
  // let in from outside can get to.
  bool isReachable(std::size_t node) const;

  // Whether it belongs to a pocket the sample resolved, which is where such a probe cannot get and a
  // molecule already inside would stay. A node in an unresolved piece is neither.
  bool isShutIn(std::size_t node) const;

  // The comment lines a report of a sampled route is headed with, below the structure's own.
  void writeHeader(std::ostream &stream) const;
};
